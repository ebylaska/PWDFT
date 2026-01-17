#ifndef _nwpw_cscf_mixing_HPP_
#define _nwpw_cscf_mixing_HPP_

#pragma once

/* nwpw_cscf_mixing.hpp
   Author - Eric Bylaska
        this class is used to keep perform kerker updates


*/

#include <iostream>

#include <cstring>
#include <cmath>
#include <algorithm>
#include "blas.h"
#include "Parallel.hpp"
#include "CGrid.hpp"
#include "nwpw_ckerker.hpp"

namespace pwdft {

const int ROWS = 40;
const int COLS = 40;


static void setElement(double A[], int row, int col, double value) {
    int index = (row - 1) * COLS + (col - 1);  // Convert 1-based (row, col) to 0-based 1D index
    A[index] = value;
}

static double getElement(double A[], int row, int col) {
    int index = (row - 1) * COLS + (col - 1);  // Convert 1-based (row, col) to 0-based 1D index
    return A[index];  // Convert to 0-based index by subtracting 1
}


/*
static void setElement(double A[], int row, int col, double value) {
    A[(col - 1) * COLS + (row - 1)] = value;   // column-major
}
static double getElement(double A[], int row, int col) {
    return A[(col - 1) * COLS + (row - 1)];
}
*/


class nwpw_cscf_mixing : public nwpw_ckerker {

   bool dokerker,isband;
   int algorithm,max_m,n2ft3d,nsize,ispin,ierr;
   double alpha_max,alpha,beta,convergence_threshold;
   double *rho_list = nullptr;
   Parallel *parall;

   int max_list=40;
   int ipiv[40];
   double  c[40],d[40],B[40*40];
   double  A[40*40],w[40],w0;
   double  Binv[40*40];
   double scf_error_prev = 0.0;
   int    scf_it = 0;
   int    scf_it_max = 0;
   //int indxf[max_m],indxu[max_m];
   int *indxf = nullptr;
   int *indxu = nullptr;;

    int m;
    int  one = 1;
    double rone=1.0;
    double mrone=-1.0;


public:
   /* constructors */
   /*************************************************
    *                                               *
    *        nwpw_cscf_mixing::nwpw_scf_mixing       *
    *                                               *
    *************************************************/
   nwpw_cscf_mixing(CGrid *mygrid0, const double g0, const int algorithm0, const double alpha0, const double beta0, const int max_m0, 
                   const int ispin0, const int nsize0, double *rho_in) : nwpw_ckerker(mygrid0, g0)
   {
      parall = mygrid0->c3db::parall;
      algorithm = algorithm0;
      alpha = alpha_max = alpha0;
      beta  = beta0;
      max_m = max_m0;
      n2ft3d = nsize0;
      nsize  = nsize0*ispin0;
      ispin  = ispin0;

      // std::cout << "nwpw_scf_mixing algorithm=" << algorithm << std::endl;
      /* simple mixing */
      if (algorithm==0)
      {
         rho_list = new (std::nothrow) double[nsize*2]();
         reset_mix(rho_in);
      }

      /* Broyden mixing */
      else if (algorithm==1)
      {
         rho_list = new (std::nothrow) double[nsize*8]();
         reset_mix(rho_in);
      }

      /* Johnson mixing */
      else if (algorithm==2)
      {
         rho_list = new (std::nothrow) double[nsize* (5+(2*max_m))]();
         indxf = new (std::nothrow) int[max_m]();
         indxu = new (std::nothrow) int[max_m]();

         reset_mix(rho_in);
         w0 = 0.01;
      }

      /* Anderson mixing */
      else if (algorithm==3)
      {
         rho_list = new (std::nothrow) double[nsize*4]();
         reset_mix(rho_in);
      }

      /* local Thomas-Fermi mixing */
      else if (algorithm==4)
      {
         rho_list = new (std::nothrow) double[nsize*3]();
         reset_mix(rho_in);
      }
      std::memcpy(rho_list, rho_in, nsize*sizeof(double));
     

   }



   /* destructor */
   /*******************************************
    *                                         *
    *    nwpw_cscf_mixing::~nwpw_cscf_mixing    *
    *                                         *
    *******************************************/
   ~nwpw_cscf_mixing() 
   {
      m = 0;
      if (rho_list) delete [] rho_list;
      if (indxf)    delete [] indxf;
      if (indxu)    delete [] indxu;
   }
   /*******************************************
    *                                         *
    *           nwpw_cscf_mixing::reset_mix    *
    *                                         *
    *******************************************/
   void reset_mix(double *rho_in)
   {
      m = 1;
      std::memcpy(rho_list,rho_in,nsize*sizeof(double));
      if (indxf && indxu) 
      {
         for (int i=0; i<max_m; ++i)
         {
            indxf[i] = (5+i)*nsize;
            indxu[i] = (5+max_m+i)*nsize;
         }
      }
   }

   /*******************************************
    *                                         *
    *        nwpw_cscf_mixing::set_alpha      *
    *                                         *
    *******************************************/
   void set_alpha(double a) { alpha = a; }

   /*******************************************
    *                                         *
    *        nwpw_cscf_mixing::get_alpha      *
    *                                         *
    *******************************************/
   double get_alpha() const { return alpha; }


   /*******************************************
    *                                         *
    *          nwpw_cscf_mixing::mix          *
    *                                         *
    *******************************************/
   void mix(double *vout, double *vnew, const double deltae, double *scf_error0) 
   {
      int nsize_global = parall->ISumAll(1, nsize);

      /* simple mixing */
      if (algorithm==0)
      {
         double *rr = rho_list;
         double *ff = rho_list+nsize;
         std::memcpy(ff,rr,nsize*sizeof(double));

         // vnew += alpha*ff;
         DAXPY_PWDFT(nsize,mrone,vout,one,ff,one);
         DSCAL_PWDFT(nsize,mrone,ff,one);

         // scf_error = sqrt(<F1|F1>)
         // Compute dot product on local rank
         double scf_error = DDOT_PWDFT(nsize,ff,one,ff,one);

         // Sum dot products across all ranks
         *scf_error0 = std::sqrt(parall->SumAll(1, scf_error) / (double)nsize_global);

         if (!std::isfinite(*scf_error0))
         {
            alpha *= 0.25;
            reset_mix(rr);
            return;
         }

         scf_it += 1;
         if (scf_it>=2)
         {
            if (*scf_error0 > scf_error_prev * 1.2) 
            {
               alpha *= 0.5;
               scf_it_max = 0;
            }
            else if (*scf_error0 < scf_error_prev * 0.8) 
            {
               scf_it_max += 1;
               if (scf_it_max > 5)
                  alpha = std::min(alpha_max, alpha * 1.1);
            }
            else
            {
               scf_it_max = 0;
            }
         }
         scf_error_prev = *scf_error0;
         alpha = std::clamp(alpha, 1e-4, alpha_max);
         
 
         // kerker stuff here
         for (auto ms=0; ms<ispin; ++ms)
            kerker_G(ff + ms*n2ft3d);
       
         std::memcpy(vnew,rr,nsize*sizeof(double));
         DAXPY_PWDFT(nsize,alpha,ff,one,vnew,one);  //vnew = vnew+alpha*f
         std::memcpy(rr,vnew,nsize*sizeof(double)); //vm=vnew
      }

      /* Broyden mixing */
      if (algorithm==1)
      {
         double *V0 = rho_list;
         double *V1 = rho_list + nsize;
         double *Vout0 = rho_list + 2*nsize;
         double *Vout1 = rho_list + 3*nsize;
         double *F0    = rho_list + 4*nsize;
         double *F1    = rho_list + 5*nsize;
         double *Vbar0 = rho_list + 6*nsize;
         double *Vbar1 = rho_list + 7*nsize;

         if (m==1)
         {
            // Vout0 = vout 
            std::memcpy(Vout0,vout,nsize*sizeof(double));

            // F0 = vout - V0
            std::memcpy(F0,V0,nsize*sizeof(double));
            DSCAL_PWDFT(nsize,mrone,F0,one);
            DAXPY_PWDFT(nsize,rone,vout,one,F0,one);

            // Apply Kerker filter to F0 to smooth residuals
            for (auto ms=0; ms<ispin; ++ms)
               kerker_G(F0 + ms*n2ft3d);

            // V1 = V0 + alpha*F0
            std::memcpy(V1,V0,nsize*sizeof(double));
            DAXPY_PWDFT(nsize,alpha,F0,one,V1,one);

            // scf_error = sqrt(<F0|F0>)
            double scf_error = DDOT_PWDFT(nsize,F0,one,F0,one);
            //*scf_error0 = std::sqrt(parall->SumAll(1,scf_error));
            //*scf_error0 = std::sqrt(parall->SumAll(1,scf_error))/((double) nsize_global);
            *scf_error0 = std::sqrt(parall->SumAll(1, scf_error) / (double)nsize_global);

            //scf_error = sqrt(scf_error);

            // vnew = V1 
            std::memcpy(vnew,V1,nsize*sizeof(double));
         }
         else //(m.gt.1) 
         {
            // vout1 = vout 
            std::memcpy(Vout1,vout,nsize*sizeof(double));

            // F1 = vout - V1 
            std::memcpy(F1,V1,nsize*sizeof(double));
            DSCAL_PWDFT(nsize,mrone,F1,one);
            DAXPY_PWDFT(nsize,rone,vout,one,F1,one);

            // Apply Kerker filter to F1 to smooth residuals
            for (auto ms=0; ms<ispin; ++ms)
               kerker_G(F1 + ms*n2ft3d);

            // scf_error = sqrt(<F1|F1>)
            double scf_error = DDOT_PWDFT(nsize,F1,one,F1,one);
            //*scf_error0 = std::sqrt(parall->SumAll(1,scf_error))/((double) nsize_global);
            *scf_error0 = std::sqrt(parall->SumAll(1, scf_error) / (double)nsize_global);

            //scf_error = sqrt(scf_error);
           
            // Beta = <F1|F1-F0>/<F1-F0/F1-F0> 
            std::memcpy(Vbar1,F1,nsize*sizeof(double));
            DAXPY_PWDFT(nsize,mrone,F0,one,Vbar1,one);
            double sum0 = DDOT_PWDFT(nsize,F1,one,Vbar1,one);
            double sum1 = DDOT_PWDFT(nsize,Vbar1,one,Vbar1,one);
            sum0 = parall->SumAll(1,sum0);
            sum1 = parall->SumAll(1,sum1);
            if (sum1 < 1e-10) 
                beta = 0.0;
            else
               beta = sum0/sum1;
            beta = std::clamp(beta, -1.0, 1.0);


            // Vbar1 = (1-Beta)*Vout1 + Beta*Vout0
            std::memcpy(Vbar1,Vout0,nsize*sizeof(double));
            DSCAL_PWDFT(nsize,beta,Vbar1,one);
            double ombeta = 1.0-beta;
            DAXPY_PWDFT(nsize,ombeta,Vout1,one,Vbar1,one);

            // Vbar0 = (1-Beta)*V1 + Beta*V0 
            std::memcpy(Vbar0,V0,nsize*sizeof(double));
            DSCAL_PWDFT(nsize,beta,Vbar0,one);
            DAXPY_PWDFT(nsize,ombeta,V1,one,Vbar0,one);

            // F0 = F1, Vout0 = Vout1, V0 = V1
            std::memcpy(F0,F1,nsize*sizeof(double));
            std::memcpy(Vout0,Vout1,nsize*sizeof(double));
            std::memcpy(V0,V1,nsize*sizeof(double));

            // V1 = (1-alpha)*Vbar0 + alpha*Vbar1
            std::memcpy(V1,Vbar1,nsize*sizeof(double));
            DSCAL_PWDFT(nsize,alpha,V1,one);
            double omalpha = 1.0-alpha;
            DAXPY_PWDFT(nsize,omalpha,Vbar0,one,V1,one);
 
            // vnew = V1
            std::memcpy(vnew,V1,nsize*sizeof(double));
         } 
         ++m;
      }


      /* Johnson mixing */
      /* Johnson mixing */
      if (algorithm==2)
      {
         double *V0 = rho_list;
         double *V1 = rho_list + nsize;
         double *F0 = rho_list + 2*nsize;
         double *F1 = rho_list + 3*nsize;
         double *dV = rho_list + 4*nsize;

         if (m==1)
         {
            // F0 = vout - V0
            std::memcpy(F0,V0,nsize*sizeof(double));
            DSCAL_PWDFT(nsize,mrone,F0,one);
            DAXPY_PWDFT(nsize,rone,vout,one,F0,one);

            // Apply Kerker filter to F0 to smooth residuals
            for (auto ms=0; ms<ispin; ++ms)
               kerker_G(F0 + ms*n2ft3d);

            // scf_error = sqrt(<F1|F1>)
            double scf_error = DDOT_PWDFT(nsize,F0,one,F0,one);
            //*scf_error0 = std::sqrt(parall->SumAll(1,scf_error))/((double) nsize_global);
            *scf_error0 = std::sqrt(parall->SumAll(1, scf_error) / (double)nsize_global);

            //scf_error = sqrt(scf_error);

            // V1 = V0 + alpha*F0
            std::memcpy(V1,V0,nsize*sizeof(double));
            DAXPY_PWDFT(nsize,alpha,F0,one,V1,one);

            // vnew = V1
            std::memcpy(vnew,V1,nsize*sizeof(double));

            std::memset(A,0,max_list*max_list*sizeof(double));

            //double sumer1 =  DDOT_PWDFT(nsize,V1,one,V1,one);
            //sumer1 = parall->SumAll(1,sumer1);
            //if (parall->is_master()) {
            //   std::cout << "m=" << m << " <F1|F1> = " << sumer1 << std::endl;
            //}

            ++m;
         }
         else //(m.gt.1)
         {
            // F1 = vout - V1
            std::memcpy(F1,V1,nsize*sizeof(double));
            //std::memcpy(F1,vnew,nsize*sizeof(double));
            DSCAL_PWDFT(nsize,mrone,F1,one);
            DAXPY_PWDFT(nsize,rone,vout,one,F1,one);

             // Apply Kerker filter to F0 to smooth residuals
             for (auto ms=0; ms<ispin; ++ms)
                kerker_G(F1 + ms*n2ft3d);
             //double sumer0 =  DDOT_PWDFT(nsize,F0,one,F0,one);
             //double sumer1 =  DDOT_PWDFT(nsize,F1,one,F1,one);
             //double sumer2 =  DDOT_PWDFT(nsize,vout,one,vout,one);
             //sumer0 = parall->SumAll(1,sumer0);
             //sumer1 = parall->SumAll(1,sumer1);
             //sumer2 = parall->SumAll(1,sumer2);
             //if (parall->is_master()) {
             //  std::cout << "<F0|F0> = " << sumer0 << std::endl;
             //  std::cout << "<F1|F1> = " << sumer1 << std::endl;
             //  std::cout << "<vout|vout> = " << sumer2 << std::endl;
             //}

            // scf_error = sqrt(<F1|F1>)
            double scf_error = DDOT_PWDFT(nsize,F1,one,F1,one);
            //*scf_error0 = std::sqrt(parall->SumAll(1,scf_error))/((double) nsize_global);
            *scf_error0 = std::sqrt(parall->SumAll(1, scf_error) / (double)nsize_global);
            //scf_error = std::sqrt(scf_error);

            //dF = dF(m-1), U = U(m-1)
            int ih = m - 2;
            double *dF = rho_list + indxf[ih];
            double *U  = rho_list + indxu[ih];
            //double *dF = rho_list + (5+m-1)*nsize;
            //double *U  = rho_list + (5+max_m+m-1)*nsize;

            // dF = (F1-F0)
            std::memcpy(dF,F1,nsize*sizeof(double));
            DAXPY_PWDFT(nsize,mrone,F0,one,dF,one);

            //sumer0 =  DDOT_PWDFT(nsize,dF,one,dF,one);
            //sumer0 = parall->SumAll(1,sumer0);
            //if (parall->is_master()) {
            //   std::cout << "<dF|dF> = " << sumer0 << std::endl;
            //}

            // dV  = (V1-V0)
            std::memcpy(dV,V1,nsize*sizeof(double));
            DAXPY_PWDFT(nsize,mrone,V0,one,dV,one);

            // U  = alpha*dF + dV
            std::memcpy(U,dV,nsize*sizeof(double));
            DAXPY_PWDFT(nsize,alpha,dF,one,U,one);

            // Define A,c and B
            int mm = m-1;
            for (auto i=1; i<m; ++i)
            {
               //define dFi here
               double *dFi = rho_list + indxf[i-1];
               //double *dFi = rho_list + indxf[i](5+i-1)*nsize;
               double sum0 = DDOT_PWDFT(nsize,dFi,one,dF,one);
               double sum1 = DDOT_PWDFT(nsize,dFi,one,F1,one);
               sum0 = parall->SumAll(1,sum0);
               sum1 = parall->SumAll(1,sum1);
               setElement(A,i,m-1,sum0);
               setElement(A,m-1,i,sum0);
               c[i-1] = sum1;
            }

            std::memset(B,0,max_list*max_list*sizeof(double));

           // double small = 0.0;
            for (auto i=1; i<=mm; ++i)
            {
               for (auto j=1; j<=mm; ++j)
               {
                  setElement(B,i,j,getElement(A,i,j));
               //   small += std::abs(getElement(A,i,j));
               }
            }
            //small /= (double)(mm*mm);
            //small = std::max(small, 1e-12);

            bool bad_diag = false;
            double small = 0.0;
            for (auto i=1; i<=mm; ++i)
            {
                double Aii = getElement(A,i,i);
                small += std::abs(Aii);

                double row_sum = 0.0;
                for (auto j=1; j<=mm; ++j)
                {
                   if (j != i) row_sum += std::abs(getElement(A,i,j));
                }

                
                //if (((Aii < 0.0) || (Aii < row_sum)) )  bad_diag = true;
                if (Aii < 0.0)  bad_diag = true;
            

            //if (parall->is_master()) {
            //    std::cout << "i=" << i << " Aii="<< Aii << " row_sum=" << row_sum << std::endl;
            //}

            }
            small /= mm;
            small = std::max(small, 1e-12);



            //if (parall->is_master()) {
            //   std::cout << "Johnson m=" << m
            //             << " mm=" << mm
            //             << " small=" << small << std::endl;
            //}
            //if (parall->is_master()) {
            //   std::cout << "A matrix:" << std::endl;
            //   for (int i=1; i<=mm; ++i) {
            //      for (int j=1; j<=mm; ++j)
            //         std::cout << std::setw(14) << getElement(A,i,j);
            //      std::cout << std::endl;
            //   }
            //}


            if (bad_diag)
            {
               // linear fallback
               //std::memcpy(V1, V0, nsize*sizeof(double));
               DAXPY_PWDFT(nsize, alpha, F1, one, V1, one);
               std::memcpy(vnew, V1, nsize*sizeof(double));

               reset_mix(V1);
               return;
            }


            for (auto i=1; i<m; ++i)
            for (auto j=1; j<m; ++j)
            {
               setElement(Binv,i,j,0.0);
            }
            for (auto i=1;i<m;++i)
            {
               setElement(Binv,i,i,1.0*small);
            }

            //int mm = m-1;
            ierr = 0;
            DGESV_PWDFT(mm,mm,B,max_list,ipiv,Binv,max_list,ierr);


            // Define d
            for (auto i=1;i<m;++i)
            {
               d[i-1] = 0.0;
               for (auto j=1; j<m; ++j)
                 d[i-1] -= (c[j-1]/small)*getElement(Binv,j,i);
            }

            // F0 = F1,  V0 = V1
            std::memcpy(F0,F1,nsize*sizeof(double));
            std::memcpy(V0,V1,nsize*sizeof(double));

            // V1 = V0 + alpha*F0 - Sum(i=1,m-1) d(i)*U(i)
            DAXPY_PWDFT(nsize,alpha,F0,one,V1,one);

            for (auto i=1; i<m; ++i)
            {
               //define U here
               //call nwpw_list_ptr(1,(5+max_m +i),U)
               double *U = rho_list + indxu[i-1];
               DAXPY_PWDFT(nsize,(d[i-1]),U,one,V1,one);
            }


            if (m<max_m)
               ++m;
            else
            {
              // Shift A matrix
              for (auto j=1; j<m-1; ++j)
              for (auto i=1; i<m-1; ++i)
                setElement(A,i,j, getElement(A,i+1,j+1));

               // Shift dF and U
               int itmpf = indxf[0];
               int itmpu = indxu[0];
               for (int i=0; i<(max_m-1); ++i)
               {
                  indxf[i] = indxf[i+1];
                  indxu[i] = indxu[i+1];
               }
               indxf[max_m-1] = itmpf;
               indxu[max_m-1] = itmpu;
               //call nwpw_list_shift_range(1,(5+1),(5+max_m))
               //call nwpw_list_shift_range(1,(5+max_m+1),(5+2*max_m))

            }

            std::memcpy(vnew,V1,nsize*sizeof(double));

         }
      }



      /* Anderson density mixing */
      if (algorithm==3)
      {
         if ((m==2) || (m==1) || (deltae>0.0))
         {
            double *rr = rho_list;
            double *ff = rho_list + nsize;
            double *tt = rho_list + 2*nsize;

            // ff=vout-vm
            std::memcpy(ff,vout,nsize*sizeof(double));
            DAXPY_PWDFT(nsize,mrone,rr,one,ff,one);

            double scf_error = DDOT_PWDFT(nsize,ff,one,ff,one);
            //*scf_error0 = std::sqrt(parall->SumAll(1,scf_error))/((double) nsize_global);
            *scf_error0 = std::sqrt(parall->SumAll(1, scf_error) / (double)nsize_global);

            //scf_error = std::sqrt(scf_error);
       
            for (auto ms=0; ms<ispin; ++ms)
               kerker_G(ff + ms*n2ft3d);
       
            std::memcpy(vnew,rr,nsize*sizeof(double));
            DAXPY_PWDFT(nsize,alpha,ff,one,vnew,one);
       
            std::memcpy(tt,vout,nsize*sizeof(double));
            std::memcpy(ff,rr,nsize*sizeof(double));
            std::memcpy(rr,vnew,nsize*sizeof(double));
         }
         else
         {
            double *rr = rho_list;           // vm 
            double *ss = rho_list + nsize;   // vm1
            double *tt = rho_list + 2*nsize; // vout1
            double *ff = rho_list + 3*nsize; //kk

            std::memcpy(ff,vout,nsize*sizeof(double));

            DAXPY_PWDFT(nsize,mrone,rr,one,ff,one);  // ff_ptr = vout - rr_ptr
            DAXPY_PWDFT(nsize,mrone,ss,one,tt,one);  // tt_ptr = vout1 - vm1

            double scf_error = DDOT_PWDFT(nsize,ff,one,ff,one);
            //*scf_error0 = std::sqrt(parall->SumAll(1,scf_error))/((double) nsize_global);
            *scf_error0 = std::sqrt(parall->SumAll(1, scf_error) / (double)nsize_global);

            //scf_error = std::sqrt(scf_error);

            for (auto ms=0; ms<ispin; ++ms)
            {
               int shift = ms*n2ft3d;

               kerker_G(ff + shift);
               kerker_G(tt + shift);

               // generate beta 
               double p00 = DDOT_PWDFT(n2ft3d,ff+shift,one,ff+shift,one); // p00 = <ff_ptr|ff_ptr>
               double p01 = DDOT_PWDFT(n2ft3d,ff+shift,one,tt+shift,one); // p01 = <ff_ptr|tt_ptr>
               double p11 = DDOT_PWDFT(n2ft3d,tt+shift,one,tt+shift,one); // p11 = <tt_ptr|tt_ptr>
               p00 = parall->SumAll(1,p00);
               p01 = parall->SumAll(1,p01);
               p11 = parall->SumAll(1,p11);
               double r00 = p00-2.00*p01+p11;
               beta = (p00-p01)/r00;

               if ((r00<0.0) || (beta<(-1.5))) beta = 1.0e-3;
               if (beta>1.5)  beta = 1.5;

               double ombeta = 1.0-beta;

               //*** vnew = (1-beta)*vm + beta*vm1 + alpha*((1-beta)*fm + beta*fm1) ***
               std::memcpy(vnew+shift,rr+shift,n2ft3d*sizeof(double));

               DSCAL_PWDFT(n2ft3d,ombeta,vnew+shift,one);
               DAXPY_PWDFT(n2ft3d,beta,ss+shift,one,vnew+shift,one);
               DSCAL_PWDFT(n2ft3d,ombeta,ff+shift,one);
               DAXPY_PWDFT(n2ft3d,beta,tt+shift,one,ff+shift,one);
               DAXPY_PWDFT(n2ft3d,alpha,ff+shift,one,vnew+shift,one);
            }
            std::memcpy(ss,rr,  nsize*sizeof(double));
            std::memcpy(rr,vnew,nsize*sizeof(double));
            std::memcpy(tt,vout,nsize*sizeof(double));
         }
         ++m; 
      }

      /* local Thomas Fermi mixing */
      if (algorithm==4)
      {
         const double twothirds = 2.0/3.0;
         double *rr = rho_list;
         double *ff = rho_list+nsize;
         double *tf = rho_list+2*nsize;
         std::memcpy(ff,rr,nsize*sizeof(double));

         // compute the  residual ff = vout - vold 
         DAXPY_PWDFT(nsize,mrone,vout,one,ff,one);
         DSCAL_PWDFT(nsize,mrone,ff,one);

         // scf_error = sqrt(<ff|ff>)
         double scf_error = DDOT_PWDFT(nsize,ff,one,ff,one);
         //*scf_error0 = std::sqrt(parall->SumAll(1,scf_error))/((double) nsize_global);
         *scf_error0 = std::sqrt(parall->SumAll(1, scf_error) / (double)nsize_global);

         //scf_error = std::sqrt(scf_error);

         // Apply Kerker filter to smooth residual ff
         for (auto ms=0; ms<ispin; ++ms)
            kerker_G(ff + ms*n2ft3d);

         // Apply TF mixing
         // tf = ff/(1+alpha*rho(n-1)**(2/3))= ff/(1+alpha*rr**
         // rho(n) = rho(n-1) + beta*tf
         //for (auto i=0; i<nsize; ++i)
         //   tf[i] = ff[i]/(1.0 + alpha*std::pow(rr[i],twothirds));
         for (auto i=0; i<nsize; ++i)
         {
            double rhoi = std::max(rr[i], 1e-12);
            tf[i] = ff[i] / (1.0 + alpha * std::pow(rhoi, 2.0/3.0));
         }


         // Apply Kerker filter to smooth residual ff
         //for (auto ms=0; ms<ispin; ++ms)
         //   kerker_G(tf + ms*n2ft3d);

         std::memcpy(vnew,rr,nsize*sizeof(double)); // 
         DAXPY_PWDFT(nsize,beta,tf,one,vnew,one);  // vnew(n) = rho(n-1) + beta*tf
         std::memcpy(rr,vnew,nsize*sizeof(double)); //vm=vnew
      }
   }
};

} // namespace pwdft

#endif

#ifndef _nwpw_scf_mixing_HPP_
#define _nwpw_scf_mixing_HPP_

#pragma once

/* nwpw_scf_mixing.hpp
   Author - Eric Bylaska
        this class is used to keep perform kerker updates


*/

#include <cmath>
#include "blas.h"
#include "Parallel.hpp"
#include "PGrid.hpp"
#include "nwpw_kerker.hpp"

namespace pwdft {

class nwpw_scf_mixing : public nwpw_kerker {

   bool dokerker,isband;
   int algorithm,max_m,n2ft3d,nsize,ispin,ierr;
   double alpha,beta,convergence_threshold;
   double *rho_list;
   Parallel *parall;

   int max_list=40;
   int ipiv[40];
   double  c[40],d[40],B[40*40];
   double  A[40*40],w[40],w0;
   double  Binv[40*40];

    int m;
    int  one = 1;
    double rone=1.0;
    double mrone=-1.0;

public:
   /* constructors */
   /*************************************************
    *                                               *
    *        nwpw_scf_mixing::nwpw_scf_mixing       *
    *                                               *
    *************************************************/
   nwpw_scf_mixing(PGrid *mygrid0, const double g0, const int algorithm0, const double alpha0, const int max_m0, 
                   const int ispin0, const int nsize0, double *rho_in) : nwpw_kerker(mygrid0, g0)
   {
      parall = mygrid0->d3db::parall;
      algorithm = algorithm0;
      alpha = alpha0;
      max_m = max_m0;
      n2ft3d = nsize0;
      nsize  = nsize0*ispin0;
      ispin  = ispin0;

      // std::cout << "nwpw_scf_mixing algorithm=" << algorithm << std::endl;
      /* simple mixing */
      if (algorithm==0)
      {
         rho_list = new (std::nothrow) double[nsize*2]();
         nwpw_scf_mixing_reset(rho_in);
      }

      /* Broyden mixing */
      else if (algorithm==1)
      {
         rho_list = new (std::nothrow) double[nsize*8]();
         nwpw_scf_mixing_reset(rho_in);
      }

      /* Johnson mixing */
      else if (algorithm==2)
      {
         rho_list = new (std::nothrow) double[nsize* (5+(2*max_m))]();
         nwpw_scf_mixing_reset(rho_in);
      }

      /* Anderson mixing */
      else if (algorithm==3)
      {
         rho_list = new (std::nothrow) double[nsize*4]();
         nwpw_scf_mixing_reset(rho_in);
      }

      /* local Thomas-Fermi mixing */
      else if (algorithm==4)
      {
         nwpw_scf_mixing_reset(rho_in);
      }
      std::memcpy(rho_list, rho_in, nsize*sizeof(double));
     

   }



   /* destructor */
   /*******************************************
    *                                         *
    *    nwpw_scf_mixing::~nwpw_scf_mixing    *
    *                                         *
    *******************************************/
   ~nwpw_scf_mixing() 
   {
      m = 0;
      delete[] rho_list;
   }
   /*******************************************
    *                                         *
    *           nwpw_scf_mixing::reset        *
    *                                         *
    *******************************************/
   void nwpw_scf_mixing_reset(double *rho_in)
   {
      m = 1;
      std::memcpy(rho_list,rho_in,nsize*sizeof(double));
   }



   /*******************************************
    *                                         *
    *           nwpw_scf_mixing::mix          *
    *                                         *
    *******************************************/
   void mix(double *vout, double *vnew, const double deltae, double *scf_error0) 
   {
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
         double scf_error = DDOT_PWDFT(nsize,ff,one,ff,one);
         parall->SumAll(1,scf_error);
         scf_error = std::sqrt(scf_error);
         *scf_error0= scf_error;
 
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

            // V1 = V0 + alpha*F0
            std::memcpy(V1,V0,nsize*sizeof(double));
            DAXPY_PWDFT(nsize,alpha,F0,one,V1,one);


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

            // scf_error = sqrt(<F1|F1>)
            double scf_error = DDOT_PWDFT(nsize,F1,one,F1,one);
            parall->SumAll(1,scf_error);
            scf_error = sqrt(scf_error);
            *scf_error0 = scf_error;
           
            // Beta = <F1|F1-F0>/<F1-F0/F1-F0> 
            std::memcpy(Vbar1,F1,nsize*sizeof(double));
            DAXPY_PWDFT(nsize,mrone,F0,one,Vbar1,one);
            double sum0 = DDOT_PWDFT(nsize,F1,one,Vbar1,one);
            double sum1 = DDOT_PWDFT(nsize,Vbar1,one,Vbar1,one);
            parall->SumAll(1,sum0);
            parall->SumAll(1,sum1);
            beta = sum0/sum1;

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
      if (algorithm==2)
      {
         double *V0 = rho_list;
         double *V1 = rho_list + nsize;
         double *F0 = rho_list + 3*nsize;
         double *F1 = rho_list + 4*nsize;
         double *dV = rho_list + 5*nsize;

         if (m==1) 
         {
            // F0 = vout - V0 
            std::memcpy(F0,V0,nsize*sizeof(double));
            DSCAL_PWDFT(nsize,mrone,F0,one);
            DAXPY_PWDFT(nsize,rone,vout,one,F0,one);

            // scf_error = sqrt(<F1|F1>) 
            double scf_error = DDOT_PWDFT(nsize,F0,one,F0,one);
            parall->SumAll(1,scf_error);
            scf_error = sqrt(scf_error);
            *scf_error0 = scf_error;
          
            // V1 = V0 + alpha*F0
            std::memcpy(V1,V0,nsize*sizeof(double));
            DAXPY_PWDFT(nsize,alpha,F0,one,V1,one);  

            // vnew = V1
            std::memcpy(vnew,V1,nsize*sizeof(double));

            std::memset(A,0,max_list*max_list*sizeof(double));

            ++m;
         }
         else //(m.gt.1) 
         {
           // F1 = vout - V1
           std::memcpy(F1,V1,nsize*sizeof(double));
           DSCAL_PWDFT(nsize,mrone,F1,one);
           DAXPY_PWDFT(nsize,rone,vout,one,F1,one);

           // scf_error = sqrt(<F1|F1>) 
           double scf_error = DDOT_PWDFT(nsize,F1,one,F1,one);
           parall->SumAll(1,scf_error);
           scf_error = std::sqrt(scf_error);
           *scf_error0 = scf_error;

           // dF = dF(m-1), U = U(m-1)
           double *dF = rho_list + (5+m)*nsize;
           double *U  = rho_list + (5+max_m+m)*nsize;

           // dF = (F1-F0) 
           std::memcpy(F1,V1,nsize*sizeof(double));
           DAXPY_PWDFT(nsize,mrone,F0,one,dF,one);

           // dV  = (V1-V0) 
           std::memcpy(dV,V1,nsize*sizeof(double));
           DAXPY_PWDFT(nsize,mrone,V0,one,dV,one);

           // U  = alpha*dF + dV
           std::memcpy(U,dV,nsize*sizeof(double));
           DAXPY_PWDFT(nsize,alpha,dF,one,U,one);

           // Define A,c and B
           for (auto i=0; i<(m-1); ++i)
           {
              double *dFi = rho_list + (5+i+1)*nsize;
              double sum0 = DDOT_PWDFT(nsize,dFi,one,dF,one);
              double sum1 = DDOT_PWDFT(nsize,dFi,one,F1,one);
              parall->SumAll(1,sum0);
              parall->SumAll(1,sum1);
              A[i+(m-1)*(m-1)] = sum0;
              A[(m-1)+i*(m-1)] = A[i+(m-1)*(m-1)];
              c[i] = sum1;
           }

           std::memset(B,0,max_list*max_list*sizeof(double));

           double small = 0.0;
           for (auto i=0; i<(m-1); ++i)
           for (auto j=0; j<(m-1); ++j)
           {
              B[i+j*(m-1)] = A[i+j*(m-1)];
              small += std::abs(A[i+j*(m-1)]);
           }
           small /= std::pow(static_cast<double>(m - 1), 2);

           for (auto i=0; i<(m-1); ++i)
           for (auto j=0; j<(m-1); ++j)
           {
              Binv[i+j*(m-1)] = 0.0;
           }
           for (auto i=0; i<(m-1); ++i)
           {
              Binv[i+i*(m-1)] = 1.0*small;
           }

           int mm = m-1;
           DGESV_PWDFT(mm,mm,B,max_list,ipiv,Binv,max_list,ierr);


           // Define d 
           for (auto i=0; i<(m-1); ++i)
           {
              d[i] = 0.0;
              for (auto j=0; j<(m-1); ++j)
                d[i] -= (c[j]/small)*Binv[j+i*(m-1)];
           }


           // F0 = F1,  V0 = V1
           std::memcpy(F0,F1,nsize*sizeof(double));
           std::memcpy(V0,V1,nsize*sizeof(double));

           // V1 = V0 + alpha*F0 - Sum(i=1,m-1) d(i)*U(i)
           DAXPY_PWDFT(nsize,alpha,F0,one,V1,one);

           for (auto i=0; i<m-1; ++i)
           {
              //call nwpw_list_ptr(1,(5+max_m +i),U)
              DAXPY_PWDFT(nsize,(d[i]),U,one,V1,one);
           }
           


           if (m<max_m) 
              ++m;
           else
           {
             // Shift A matrix 
             for (auto j=0; j<m-2; ++j)
             for (auto i=0; i<m-2; ++i)
               A[i+j*(m-1)] = A[(i+1)+(j+1)*(m-1)];

              // Shift dF and U 
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
            parall->SumAll(1,scf_error);
            scf_error = std::sqrt(scf_error);
            *scf_error0 = scf_error;
       
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
            parall->SumAll(1,scf_error);
            scf_error = std::sqrt(scf_error);
            *scf_error0 = scf_error;

            for (auto ms=0; ms<ispin; ++ms)
            {
               int shift = ms*n2ft3d;

               kerker_G(ff + shift);
               kerker_G(tt + shift);

               // generate beta 
               double p00 = DDOT_PWDFT(n2ft3d,ff+shift,one,ff+shift,one); // p00 = <ff_ptr|ff_ptr>
               double p01 = DDOT_PWDFT(n2ft3d,ff+shift,one,tt+shift,one); // p01 = <ff_ptr|tt_ptr>
               double p11 = DDOT_PWDFT(n2ft3d,tt+shift,one,tt+shift,one); // p11 = <tt_ptr|tt_ptr>
               parall->SumAll(1,p00);
               parall->SumAll(1,p01);
               parall->SumAll(1,p11);
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
         double *rr = rho_list;
         double *ss = rho_list + nsize;
         double *tt = rho_list + 2*nsize;
         double *ff = rho_list + 3*nsize;

         std::memcpy(ss,rr,nsize*sizeof(double));
         std::memcpy(rr,vnew,nsize*sizeof(double));
         std::memcpy(tt,vout,nsize*sizeof(double));
      }
   }
};

} // namespace pwdft

#endif

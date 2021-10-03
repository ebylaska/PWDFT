/* nwpw_apc.cpp -
   Author - Eric Bylaska
*/


#include        <iostream>
#include        <cstring>
#include        <cmath>

#include        "nwpw_timing.hpp"
#include        "gdevice.hpp"


#include        "blas.h"

#include        "nwpw_apc.hpp"

namespace pwdft {



/* Constructors */

/*******************************************
 *                                         *
 *            nwpw_apc::nwpw_apc           *
 *                                         *
 *******************************************/
nwpw_apc::nwpw_apc(Ion *myionin, Pneb *mypnebin, Strfac *mystrfacin, Control2& control)
{

   myion    = myionin;
   mypneb   = mypnebin;
   mystrfac = mystrfacin;
   apc_on = control.APC_on();

   if (apc_on)  
   {
      Gc  = control.APC_Gc();
      nga = control.APC_nga();
      if (nga<=0) 
      {
         nga = 3;
         gamma = new double [nga];
         gamma[0] = 0.6;
         gamma[1] = 0.9;
         gamma[2] = 1.35;
      }
      else
      {
         gamma = new double [nga];
         for (auto i=0; i<nga; ++i)
            gamma[i] = control.APC_gamma(i);
      }
      ngs = nga*myion->nion;

      /* allocate APC memory */
      A  = new double [4*ngs*ngs];
      Am = new double [ngs*ngs];
      b  = new double [4*ngs];
      q  = new double [ngs];
      u  = new double [ngs];

      int npack0 = mypneb->npack(0);
      w    = new double[npack0];
      gaus = new double[nga*npack0];



      /* define weight function */
      double gg,xx;
      double fourpi = 16.0*atan(1.0);
      double *Gx = mypneb->Gpackxyz(0,0);
      double *Gy = mypneb->Gpackxyz(0,1);
      double *Gz = mypneb->Gpackxyz(0,2);
      for (auto i=0; i<npack0; ++i)
      {
         gg = Gx[i]*Gx[i] + Gy[i]*Gy[i] + Gz[i]*Gz[i];
         w[i] = 0.0;
         if ((gg>1.0e-6) && (gg<(Gc*Gc))) 
         {
            xx = gg-Gc*Gc;
            w[i] = fourpi*xx*xx/(gg*Gc*Gc);
         }
      }

      /* define Gaussians in G-space */
      double coef = 1.0/mypneb->lattice->omega();
      for (auto n=0; n<nga; ++n) 
      {
         xx = gamma[n]*gamma[n]/4.0;
         for (auto i=0; i<npack0; ++i)
         {
            gg = Gx[i]*Gx[i] + Gy[i]*Gy[i] + Gz[i]*Gz[i];
            gaus[n*npack0 + i] = coef*exp(-xx*gg);

         }
      }

      /* write out APC header */   
      if (mypneb->PGrid::parall->is_master())
      {
         std::cout << std::endl;
         std::cout << " initializing nwpw_APC object" << std::endl;
         std::cout << " ----------------------------" << std::endl;
         std::cout << " nga, ngs:  " << nga << " " << ngs << std::endl;
         std::cout << " Gc      :  " << Gc << std::endl;
         for (auto i=0; i<nga; ++i)
            std::cout << " APC gamma: " << i << " " << gamma[i] << std::endl;
      }
   }
}

/*******************************************
 *                                         *
 *            nwpw_apc::gen_APC            *
 *                                         *
 *******************************************/
void nwpw_apc::gen_APC(double *dng, bool move)
{
   if (apc_on) 
   {
      int npack0 = mypneb->npack(0);
      int ispin  = mypneb->ispin;

      int i,j,indx,indxt;
      double fourpi = 16.0*atan(1.0);
      double *Gx = mypneb->Gpackxyz(0,0);
      double *Gy = mypneb->Gpackxyz(0,1);
      double *Gz = mypneb->Gpackxyz(0,2);
      double omega = mypneb->lattice->omega();

      double exi[2*npack0];
      double exj[2*npack0];
      double gaus_i[2*npack0];
      double gaus_j[2*npack0];
      double xtmp[npack0];

      /* calculate N = dng(G=0)*omega */
      double N = ((double) (mypneb->ne[0] + mypneb->ne[ispin-1]));

      /* calculate c_i = omega*gaus_i(G=0) = 1.0d0 */

      /* calculate b_i = omega*Sum(G) w(G)*Re(dcongj(dng(G))*gaus_i(G)) */
      for (auto ii=0; ii<myion->nion; ++ii)
      {
         mystrfac->strfac_pack(0,ii,exi);

         for (auto iii=0; iii<nga; ++iii)
         {
            i = iii + ii*nga;

            /* gaus_i(G)) */
            /* w(G)*gaus_i(G)) */
            mypneb->tcc_Mul(0,&gaus[npack0*iii],exi,gaus_i);
            mypneb->tc_Mul(0,w,gaus_i);

            /* bi = omega*Sum(G) w(G)*Re(dcongj(dng(G))*gaus_i(G)) */
            b[i] = omega*mypneb->cc_pack_dot(0,dng,gaus_i);

            if (move) 
            {
               mypneb->cct_iconjgMulb(0,dng,gaus_i,xtmp);
               b[i+ngs]   = omega*mypneb->tt_pack_dot(0,Gx,xtmp);
               b[i+2*ngs] = omega*mypneb->tt_pack_dot(0,Gy,xtmp);
               b[i+3*ngs] = omega*mypneb->tt_pack_dot(0,Gz,xtmp);
            }
                   
         }
      }


      /* calculate A_ij = omega*Sum(G) w(G)*dcongj(gaus_i(G))*gaus_j(G)) */
      for (auto ii=0; ii<myion->nion; ++ii)
      {
         mystrfac->strfac_pack(0,ii,exi);

         for (auto jj=ii; jj<myion->nion; ++jj)
         {
            mystrfac->strfac_pack(0,jj,exj);

            for (auto iii=0; iii<nga; ++iii)
            {
               /* gaus_i(G)) */
               /* w(G)*gaus_i(G)) */
               mypneb->tcc_Mul(0,&gaus[npack0*iii],exi,gaus_i);
               mypneb->tc_Mul(0,w,gaus_i);

               for (auto jjj=0; jjj<nga; ++jjj)
               {
                  /* gaus_j(G)) */
                  mypneb->tcc_Mul(0,&gaus[npack0*jjj],exj,gaus_j);

                  /* omega*Sum(G) w(G)*gaus_i(G)*gaus_j(G)) */
                  double e1 = omega*mypneb->cc_pack_dot(0,gaus_i,gaus_j);
                  i = iii + ii*nga;
                  j = jjj + jj*nga;

                  indx  = i + j*ngs;
                  indxt = j + i*ngs;

                  A[indx]  = e1;
                  A[indxt] = e1;

                  if (move)
                  {
                     mypneb->cct_iconjgMulb(0,gaus_i,gaus_j,xtmp);
                     A[indx+ngs*ngs]   = omega*mypneb->tt_pack_dot(0,Gx,xtmp);
                     A[indx+2*ngs*ngs] = omega*mypneb->tt_pack_dot(0,Gy,xtmp);
                     A[indx+3*ngs*ngs] = omega*mypneb->tt_pack_dot(0,Gz,xtmp);
                     if (indx!=indxt)
                     {
                        A[indxt+ngs*ngs]   = -A[indx+ngs*ngs];
                        A[indxt+2*ngs*ngs] = -A[indx+2*ngs*ngs];
                        A[indxt+3*ngs*ngs] = -A[indx+3*ngs*ngs]; 
                     }
                  }
               }
            }
         }
      }

      /* perform matrix operations in serial */
      memset(Am,0,ngs*ngs*sizeof(double));
      memset(q,0,ngs*sizeof(double));

      if (mypneb->PGrid::parall->is_master())
      {
         int ierr,rank;

         int lwork = 5*ngs*ngs;
         double work[lwork];
         double AAA[ngs*ngs];
         double rcond = 1.0e-9;

        for (i=0; i<ngs*ngs; ++i) AAA[i]        = A[i];
        for (i=0; i<ngs;     ++i) Am[i + i*ngs] = 1.0;

        DGELSS_PWDFT(ngs,ngs,ngs,AAA,ngs,Am,ngs,q,rcond,rank,work,lwork,ierr);


        /* calculate q_i */
        double sum  = 0.0;
        double sum1 = 0.0;
        for (i=0; i<ngs; ++i)
        for (j=0; j<ngs; ++j)
        {
           indx = i + j*ngs;
           sum  += Am[indx]*b[j];
           sum1 += Am[indx];
        }
        sum = (sum-N)/sum1;

        memset(q,0,ngs*sizeof(double));
        for (i=0; i<ngs; ++i)
        {
           sum1 = 0.0;
           for (j=0; j<ngs; ++j)
           {
              indx = i + j*ngs;
              sum1 += Am[indx]*(b[j]-sum);
           }
           q[i] = sum1;
        }
      }

      /* synchronization */
      mypneb->PGrid::parall->Vector_SumAll(1,ngs,q);
      mypneb->PGrid::parall->Vector_SumAll(1,ngs*ngs,Am);



   } /*apc_on*/
}

/*******************************************
 *                                         *
 *            nwpw_apc::dngen_APC          *
 *                                         *
 *******************************************/
void nwpw_apc::dngen_APC(double *dn, bool move)
{
}


/*******************************************
 *                                         *
 *            nwpw_apc::Qtot_APC           *
 *                                         *
 *******************************************/
double nwpw_apc::Qtot_APC(const int ii)
{
   double qq = 0.0;
   if (nga>0) for (auto n=0; n<nga; ++n) qq += q[n+ii*nga];
   return qq;
}

/*******************************************
 *                                         *
 *            nwpw_apc::print_APC          *
 *                                         *
 *******************************************/
std::string nwpw_apc::print_APC(const double *zv)
{
   std::stringstream stream;

   stream << std::endl 
          << "*************************************************************" << std::endl 
          << "**                                                         **" << std::endl
          << "**          PSPW Atomic Point Charge (APC) Analysis        **" << std::endl
          << "**                                                         **" << std::endl
          << "**   Point charge analysis based on paper by P.E. Blochl   **" << std::endl
          << "**         (J. Chem. Phys. vol 103, page 7422, 1995)       **" << std::endl
          << "**                                                         **" << std::endl
          << "*************************************************************" << std::endl
          << std::endl
          << " nwpw_APC object" << std::endl
          << " ---------------" << std::endl 
          << " nga = " << std::setw(3) << nga 
          << " ngs = " << std::setw(5) << ngs << std::endl 
          << " Gc  = " << std::fixed   << std::setprecision(5) << std::setw(9) << Gc << std::endl;
   for (auto i=0; i<nga; ++i)
      stream  << " APC gamma: " << i << " " << gamma[i] << std::endl;

   stream << std::endl;
   stream << " charge analysis on each atom" << std::endl 
          << " ----------------------------" << std::endl << std::endl
          << "      no  atom        Qelc        Qion      Qtotal" << std::endl
          << "   -----  ----     -------     -------     -------" << std::endl;
   double sume = 0.0;
   double sumi = 0.0;
   for (auto ii=0; ii<myion->nion; ++ii)
   {
      int    ia  = myion->katm[ii];
      double sum = 0.0;
      for (auto i=0; i<nga; ++i) sum += q[i+ii*nga];
      sume -= sum;
      sumi += zv[ia];
      stream << std::setw(8)  << ii+1
             << std::setw(6)  << myion->symbol(ii)
             << std::setw(12) << std::fixed << std::setprecision(3) << -sum
             << std::setw(12) << std::fixed << std::setprecision(3) << zv[ia]
             << std::setw(12) << std::fixed << std::setprecision(3) << zv[ia]-sum
             << std::endl;
   }
   stream << std::setw(14) << "Total Q"
          << std::setw(12) << std::fixed << std::setprecision(3) << sume
          << std::setw(12) << std::fixed << std::setprecision(3) << sumi
          << std::setw(12) << std::fixed << std::setprecision(3) << sume+sumi
          << std::endl;


   stream << std::endl  << std::endl
          << " gaussian coefficients of model density" << std::endl
          << " --------------------------------------" << std::endl
          << std::endl;
   stream << "      no  atom";
   stream << setw(7) << "g=" << std::fixed << std::setprecision(3) << 0.0;
   for (auto i=0; i<nga; ++i) 
      stream << setw(7) << "g=" << std::fixed << std::setprecision(3) << gamma[i];
   stream << std::endl;
   stream << "   -----  ----";
   for (auto i=0; i<(nga+1); ++i) stream << setw(12) << "-------";
   stream << std::endl;
   for (auto ii=0; ii<myion->nion; ++ii)
   {
      int    ia  = myion->katm[ii];
      stream << std::setw(8)  << ii+1
             << std::setw(6)  << myion->symbol(ii)
             << std::setw(12) << std::fixed << std::setprecision(3) << zv[ia];
      for (auto i=0; i<nga; ++i)
          stream << std::setw(12) << std::fixed << std::setprecision(3) << -q[i+ii*nga];
      stream << std::endl;
   }
   stream << std::endl
          << std::endl;

   return stream.str();
}


}

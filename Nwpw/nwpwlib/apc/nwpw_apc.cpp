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
}

/*******************************************
 *                                         *
 *            nwpw_apc::dngen_APC          *
 *                                         *
 *******************************************/
void nwpw_apc::dngen_APC(double *dn, bool move)
{
}



}

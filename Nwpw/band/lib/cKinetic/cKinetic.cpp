/* cKinetic.C -
   Author - Eric Bylaska
*/

/*


#include        <cmath>
#include        <cstdio>
#include        <cstdlib>
#include        <iostream>
#include        <stdio.h>
#include	<string>
*/

#include "cKinetic.hpp"

namespace pwdft {

/*********************************************
 *                                           *
 *     cKinetic_Operator::cKinetic_Operator  *
 *                                           *
 *********************************************/
cKinetic_Operator::cKinetic_Operator(Cneb *mygrid) 
{
   mycneb = mygrid;

   int nbrillq = mycneb->nbrillq;
   int npack1  = mycneb->npack(1);
   tg = new double[nbrillq*npack1];
   //double *tmp = new double[mygrid->nfft3d];
 
   for (auto nbq=0; nbq<nbrillq; ++nbq)
   {
      double *Gpackx = mycneb->Gpackxyz(nbq,0);
      double *Gpacky = mycneb->Gpackxyz(nbq,1);
      double *Gpackz = mycneb->Gpackxyz(nbq,2);

      double *kv = mycneb->pbrill_kvector(nbq);
      double *tmp_tg = tg + nbq*npack1;
 
      for (auto k=0; k<npack1; ++k) 
      {
         double gx = Gpackx[k] + kv[0];
         double gy = Gpacky[k] + kv[1];
         double gz = Gpackz[k] + kv[2];
         double gg = gx*gx + gy*gy + gz*gz;
         tmp_tg[k] = -0.5 * gg;
      }
      //mycneb->tt_pack_copy(1, tmp, tg + nbq*npack1);
   }
   //delete[] tmp;
}


/*******************************************
 *                                         *
 *         cKinetic_Operator::ke           *
 *                                         *
 *******************************************/
void cKinetic_Operator::ke(const double *psi, double *tpsi) 
{
   int nbqsize = (mycneb->nbrillq);
   int nsize   = (mycneb->neq[0] + mycneb->neq[1]);
   int npack1  = (mycneb->npack(1));

   int k1 = 0;
   int k2 = 1;
   for (auto nbq=0; nbq<nbqsize; ++nbq)
   {
      double *tmp_tg = tg + nbq*npack1;

      for (auto n=0; n<nsize; ++n)
      for (auto k=0; k<npack1; ++k) 
      {
        tpsi[k1] += tmp_tg[k] * psi[k1];
        tpsi[k2] += tmp_tg[k] * psi[k2];
        k1 += 2;
        k2 += 2;
      }
   }
}


/*******************************************
 *                                         *
 *        cKinetic_Operator::ke_ave        *
 *                                         *
 *******************************************/
double cKinetic_Operator::ke_ave(const double *psi) 
{
 
   int nbqsize = (mycneb->nbrillq);
   int nsize   = (mycneb->neq[0] + mycneb->neq[1]);
   int npack1  = (mycneb->npack(1));
 
   double ave = 0.0;
   int k1 = 0;
   int k2 = 1;
   for (auto nbq=0; nbq<nbqsize; ++nbq)
   {
      double weight  = mycneb->pbrill_weight(nbq);
      double *tmp_tg = tg + nbq*npack1;

      for (auto n=0; n<nsize; ++n) 
      for (auto k=0; k<npack1; ++k) 
      {
         ave += weight*tmp_tg[k] * (psi[k1] * psi[k1] + psi[k2] * psi[k2]);
         k1 += 2;
         k2 += 2;
      }
   }

   ave = mycneb->c3db::parall->SumAll(0, ave);

   if (mycneb->ispin == 1)
     ave *= 2.0;
   ave = -ave;
 
   return ave;
}

} // namespace pwdft

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

   int nbrillq    = mycneb->nbrillq;
   int npack1_max = mycneb->npack1_max();
   tg = new double[nbrillq*npack1_max];
   //double *tmp = new double[mygrid->nfft3d];
   std::memset(tg,0.0,nbrillq*npack1_max);
 
   for (auto nbq=0; nbq<nbrillq; ++nbq)
   {
      int npack1     = mycneb->npack(1+nbq);
      double *Gpackx = mycneb->Gpackxyz(1+nbq,0);
      double *Gpacky = mycneb->Gpackxyz(1+nbq,1);
      double *Gpackz = mycneb->Gpackxyz(1+nbq,2);

      double *kv = mycneb->pbrill_kvector(nbq);
      double *tmp_tg = tg + nbq*npack1_max;
 
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
   int nbqsize    = mycneb->nbrillq;
   int nsize      = mycneb->neq[0] + mycneb->neq[1];
   int npack1_max = mycneb->npack1_max();
   int shift1 = 2*npack1_max;
   int indx1n = 0;

   for (auto nbq=0; nbq<nbqsize; ++nbq)
   {
      int npack1  = mycneb->npack(1+nbq);
      double *tmp_tg = tg + nbq*npack1_max;

      for (auto n=0; n<nsize; ++n)
      {
         const double *tmp_psi  = psi  + indx1n;
         double       *tmp_tpsi = tpsi + indx1n;
         int k1 = 0;
         int k2 = 1;
         for (auto k=0; k<npack1; ++k) 
         {
           //tmp_tpsi[k1] += tmp_tg[k]*tmp_psi[k1];
           //tmp_tpsi[k2] += tmp_tg[k]*tmp_psi[k2];
           tmp_tpsi[k1] = tmp_tg[k]*tmp_psi[k1];
           tmp_tpsi[k2] = tmp_tg[k]*tmp_psi[k2];
           k1 += 2;
           k2 += 2;
         }
         indx1n += shift1;
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
 
   int nbqsize    = mycneb->nbrillq;
   int nsize      = mycneb->neq[0] + mycneb->neq[1];
   int npack1_max = mycneb->npack1_max();
   int shift1 = 2*npack1_max;
   int indx1n = 0;
 
   double ave = 0.0;
   for (auto nbq=0; nbq<nbqsize; ++nbq)
   {
      int npack1     = mycneb->npack(1+nbq);
      double weight  = mycneb->pbrill_weight(nbq);
      double *tmp_tg = tg + nbq*npack1_max;

      for (auto n=0; n<nsize; ++n) 
      {
         const double *tmp_psi = psi + indx1n;
         int k1 = 0; 
         int k2 = 1;
         for (auto k=0; k<npack1; ++k) 
         {
            ave += weight*tmp_tg[k]*(tmp_psi[k1]*tmp_psi[k1] + tmp_psi[k2]*tmp_psi[k2]);
            k1 += 2;
            k2 += 2;
         }
         indx1n += shift1;
      }
   }

   ave = mycneb->c3db::parall->SumAll(0, ave);

   if (mycneb->ispin == 1)
      ave *= 2.0;
   ave = -ave;
 
   return ave;
}

/*******************************************
 *                                         *
 *         cKinetic_Operator::ke_orb       *
 *                                         *
 *******************************************/
void cKinetic_Operator::ke_orb(const int nbq1, const double *psi, double *tpsi) 
{
   int npack1_max = mycneb->npack1_max();
   double *tmp_tg = tg + (nbq1-1)*npack1_max;

   int npack1  = mycneb->npack(nbq1);

   int k1 = 0;
   int k2 = 1;
   for (auto k=0; k<npack1; ++k)
   {
     //tmp_tpsi[k1] += tmp_tg[k]*tmp_psi[k1];
     //tmp_tpsi[k2] += tmp_tg[k]*tmp_psi[k2];
     tpsi[k1] = tmp_tg[k]*psi[k1];
     tpsi[k2] = tmp_tg[k]*psi[k2];
     k1 += 2;
     k2 += 2;
   }
}


/***********************************************
 *                                             *
 *       cKinetic_Operator::ke_precondition     *
 *                                             *
 ***********************************************/
// **** My preconditioner ****
void cKinetic_Operator::ke_precondition(const int nbq1, const double Ep, const int neall, double *psi, double *tpsi) 
{
   int npack1  = (mycneb->npack(nbq1));
   int npack2  = 2*npack1;
   int npack1_max = mycneb->npack1_max();
   double *tmp    = new double[npack2];
   double *tg1    = tg + (nbq1-1);
   int k1 = 0;
   int k2 = 1;
   for (auto n=0; n<neall; ++n)
   {
      double *worb = psi + n*npack2;
      mycneb->tcc_pack_Mul(nbq1,tg1,worb,tmp);
      double sum = mycneb->cc_pack_dot(nbq1,tmp,worb);
      for (auto k=0; k<npack1; ++k) 
      {
         double x = tg1[k];
         x =  x*(worb[k1]*worb[k1] + worb[k2]*worb[k2])/sum;
         double cm = 27.00+(18.00+(12.00+8.00*x)*x)*x;
         double dm = (cm + 16.00* x*x*x*x);
         cm = cm/dm;
         tpsi[k1] *= cm;
         tpsi[k2] *= cm;
         k1 += 2;
         k2 += 2;
      }
   }
   delete [] tmp;
}

} // namespace pwdft

/* Kinetic.C -
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

#include "Kinetic.hpp"
#include "PGrid.hpp"

namespace pwdft {

/*******************************************
 *                                         *
 *     Kinetic_Operator::Kinetic_Operator  *
 *                                         *
 *******************************************/
Kinetic_Operator::Kinetic_Operator(Pneb *mygrid) {
  int k;
  double gg;
  double *Gx = mygrid->Gxyz(0);
  double *Gy = mygrid->Gxyz(1);
  double *Gz = mygrid->Gxyz(2);
  tg = new double[mygrid->npack(1)];
  double *tmp = new double[mygrid->nfft3d];

  mypneb = mygrid;

  for (k = 0; k < (mypneb->nfft3d); ++k) {
    gg = Gx[k] * Gx[k] + Gy[k] * Gy[k] + Gz[k] * Gz[k];
    tmp[k] = -0.5 * gg;
  }
  mypneb->t_pack(1, tmp);
  mypneb->tt_pack_copy(1, tmp, tg);

  delete[] tmp;
}

/*******************************************
 *                                         *
 *         Kinetic_Operator::ke            *
 *                                         *
 *******************************************/
void Kinetic_Operator::ke(double *psi, double *tpsi) 
{
   int k, k1, k2, n, nsize, ksize;
 
   nsize = (mypneb->neq[0] + mypneb->neq[1]);
   ksize = (mypneb->npack(1));
   k1 = 0;
   k2 = 1;
   for (n = 0; n < nsize; ++n)
     for (k = 0; k < ksize; ++k) 
     {
        tpsi[k1] += tg[k] * psi[k1];
        tpsi[k2] += tg[k] * psi[k2];
        k1 += 2;
        k2 += 2;
     }
}

/*******************************************
 *                                         *
 *         Kinetic_Operator::ke_orb        *
 *                                         *
 *******************************************/
void Kinetic_Operator::ke_orb(double *orb, double *torb) 
{
   int k1, k2, ksize;
 
   ksize = (mypneb->npack(1));
   k1 = 0;
   k2 = 1;
   for (auto k=0; k<ksize; ++k) 
   {
      torb[k1] += tg[k] * orb[k1];
      torb[k2] += tg[k] * orb[k2];
      k1 += 2;
      k2 += 2;
   }
}

/*******************************************
 *                                         *
 *        Kinetic_Operator::ke_ave         *
 *                                         *
 *******************************************/
double Kinetic_Operator::ke_ave(double *psi) 
{
   int k, k1, k2, n, nsize, ksize1, ksize2;
   double ave;
 
   nsize = (mypneb->neq[0] + mypneb->neq[1]);
   ksize1 = (mypneb->nzero(1));
   ksize2 = (mypneb->npack(1));
 
   ave = 0.0;
   k1 = 0;
   k2 = 1;
   for (n = 0; n < nsize; ++n) {
     for (k = 0; k < ksize1; ++k) {
       ave += tg[k] * (psi[k1] * psi[k1] + psi[k2] * psi[k2]);
       k1 += 2;
       k2 += 2;
     }
     for (k = ksize1; k < ksize2; ++k) {
       ave += 2.0 * tg[k] * (psi[k1] * psi[k1] + psi[k2] * psi[k2]);
       k1 += 2;
       k2 += 2;
     }
   }
   ave = mypneb->d3db::parall->SumAll(0, ave);
   if (mypneb->ispin == 1)
     ave *= 2.0;
   ave = -ave;
 
   return ave;
}

/***********************************************
 *                                             *
 *       Kinetic_Operator::ke_precondition     *
 *                                             *
 ***********************************************/
// **** My preconditioner ****
void Kinetic_Operator::ke_precondition(const double Ep, const int neall, double *psi, double *tpsi) 
{
   int npack1  = (mypneb->npack(1));
   int npack2  = 2*npack1;
   double *tmp = new double[npack2];
   int k1 = 0;
   int k2 = 1;
   for (auto n=0; n<neall; ++n)
   {
      double *worb = psi + n*npack2;
      mypneb->tcc_pack_Mul(1,tg,worb,tmp);
      double sum = mypneb->cc_pack_dot(1,tmp,worb);
      for (auto k=0; k<npack1; ++k) 
      {
         double x = tg[k];
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

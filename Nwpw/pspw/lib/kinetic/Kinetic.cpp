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

} // namespace pwdft

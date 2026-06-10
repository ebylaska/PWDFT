/* Ewald.C -
   Author - Eric Bylaska
*/

#include <cstring>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "units.hpp"
#include "Control2.hpp"
#include "Ewald.hpp"
//#include	"Pseudopotential.hpp"

namespace pwdft {

/*************************************
 *                                   *
 *          mandelung_get            *
 *                                   *
 *************************************/

static double mandelung_get(Lattice *lattice) 
{
   int n1, n2, n3;
   double ax, ay, az, gx, gy, gz, gg;
   double alpha, alpha1, alpha2, sum, ea;
   double rc, rs, epsilon;

   constexpr int N = 40;
 
   //pi = 4.0 * atan(1.0);
   constexpr double pi     = units::PI;

   rs = pow((3.0 * lattice->omega() / (4.0 * pi)), (1.0 / 3.0));
   rc = rs;
   epsilon = 1.0 / rc;
 
   sum = 0.0;
   for (n1 = (-N + 1); n1 <= (N - 1); ++n1)
     for (n2 = (-N + 1); n2 <= (N - 1); ++n2)
       for (n3 = (-N + 1); n3 <= (N - 1); ++n3)
         if ((n1 != 0) || (n2 != 0) || (n3 != 0)) {
           ax = n1 * lattice->unita(0, 0) + n2 * lattice->unita(0, 1) +
                n3 * lattice->unita(0, 2);
           ay = n1 * lattice->unita(1, 0) + n2 * lattice->unita(1, 1) +
                n3 * lattice->unita(1, 2);
           az = n1 * lattice->unita(2, 0) + n2 * lattice->unita(2, 1) +
                n3 * lattice->unita(2, 2);
           ea = sqrt(ax * ax + ay * ay + az * az);
           sum += erfc(epsilon * ea) / ea;
         }
   alpha1 = sum;
 
   /* calculate alpha2 */
   sum = 0.0;
   for (n1 = (-N + 1); n1 <= (N - 1); ++n1)
     for (n2 = (-N + 1); n2 <= (N - 1); ++n2)
       for (n3 = (-N + 1); n3 <= (N - 1); ++n3)
         if ((n1 != 0) || (n2 != 0) || (n3 != 0)) {
           gx = n1 * lattice->unitg(0, 0) + n2 * lattice->unitg(0, 1) +
                n3 * lattice->unitg(0, 2);
           gy = n1 * lattice->unitg(1, 0) + n2 * lattice->unitg(1, 1) +
                n3 * lattice->unitg(1, 2);
           gz = n1 * lattice->unitg(2, 0) + n2 * lattice->unitg(2, 1) +
                n3 * lattice->unitg(2, 2);
           gg = gx * gx + gy * gy + gz * gz;
           sum += (4.0 * pi / gg) * exp(-gg * rc * rc / 4.0);
         }
   alpha2 = sum / (lattice->omega());
 
   sum = alpha1 + alpha2 - pi * rc * rc / (lattice->omega()) -
         2.0 * epsilon / sqrt(pi);
   alpha = -sum * rs;
   return alpha;
}

/* Constructors */

/*********************************
 *                               *
 *          Ewald::Ewald         *
 *                               *
 *********************************/
// Ewald::Ewald(Parallel *inparall, Ion *inion, Pseudopotential *inpsp)
Ewald::Ewald(Parallel *inparall, Ion *inion, Lattice *inlattice, Control2 &control, double *inzv) 
{
   int i, j, k, l, k1, k2, k3;
   int enxh, enyh, enzh, enpack0;
   int tnp, tid, dutask;
   double g1, g2, g3, gg1, gg2, gg3, gg, ggcut;
   double pi4, rs, w, term;
   double q, z, zz;
   double eps = 1.0e-12;
 
   ewaldparall = inparall;
   ewaldion = inion;
   ewaldlattice = inlattice;
   tnp = ewaldparall->np();
   tid = ewaldparall->taskid();
 
   for (j=0; j<3; ++j)
   for (i=0; i<3; ++i) 
   {
      unitg[i + j*3] = ewaldlattice->unitg(i, j);
      unita[i + j*3] = ewaldlattice->unita(i, j);
   }
 
   enx = control.ewald_ngrid(0);
   eny = control.ewald_ngrid(1);
   enz = control.ewald_ngrid(2);
   enxh = enx / 2;
   enyh = eny / 2;
   enzh = enz / 2;
 
   /* determine ggcut */
   g1 = unitg[0] * (enxh);
   g2 = unitg[1] * (enxh);
   g3 = unitg[2] * (enxh);
   gg1 = g1 * g1 + g2 * g2 + g3 * g3;
 
   g1 = unitg[3] * (enyh);
   g2 = unitg[4] * (enyh);
   g3 = unitg[5] * (enyh);
   gg2 = g1 * g1 + g2 * g2 + g3 * g3;
 
   g1 = unitg[6] * (enzh);
   g2 = unitg[7] * (enzh);
   g3 = unitg[8] * (enzh);
   gg3 = g1 * g1 + g2 * g2 + g3 * g3;
 
   ggcut = gg1;
   if (gg2 < ggcut)
     ggcut = gg2;
   if (gg3 < ggcut)
     ggcut = gg3;
   if ((2.0 * control.ecut()) < ggcut)
     ggcut = 2.0 * control.ecut();
   eecut = 0.5 * ggcut;
 
   /* determine enpack */
   dutask = 0;
   enpack = 0;
   enida = 0;
   k1 = 0;
   k2 = 0;
   k3 = 0;
   g1 = k1 * unitg[0] + k2 * unitg[3] + k3 * unitg[6];
   g2 = k1 * unitg[1] + k2 * unitg[4] + k3 * unitg[7];
   g3 = k1 * unitg[2] + k2 * unitg[5] + k3 * unitg[8];
   gg = g1 * g1 + g2 * g2 + g3 * g3;
   if ((gg - ggcut) < (-eps)) {
     if (dutask == tid) {
       enpack++;
       enida++;
     }
     dutask = (dutask + 1) % tnp;
   }
   for (k = 1; k < enzh; ++k) 
   {
      k1 = 0;
      k2 = 0;
      k3 = k;
      g1 = k1 * unitg[0] + k2 * unitg[3] + k3 * unitg[6];
      g2 = k1 * unitg[1] + k2 * unitg[4] + k3 * unitg[7];
      g3 = k1 * unitg[2] + k2 * unitg[5] + k3 * unitg[8];
      gg = g1 * g1 + g2 * g2 + g3 * g3;
      if ((gg - ggcut) < (-eps)) 
      {
         if (dutask == tid)
            enpack++;
         dutask = (dutask + 1) % tnp;
      }
   }
 
   for (k = (-enzh + 1); k < enzh; ++k)
     for (j = 1; j < enyh; ++j) {
       k1 = 0;
       k2 = j;
       k3 = k;
       g1 = k1 * unitg[0] + k2 * unitg[3] + k3 * unitg[6];
       g2 = k1 * unitg[1] + k2 * unitg[4] + k3 * unitg[7];
       g3 = k1 * unitg[2] + k2 * unitg[5] + k3 * unitg[8];
       gg = g1 * g1 + g2 * g2 + g3 * g3;
       if ((gg - ggcut) < (-eps)) {
         if (dutask == tid)
           enpack++;
         dutask = (dutask + 1) % tnp;
       }
     }
 
   for (k = (-enzh + 1); k < enzh; ++k)
     for (j = (-enyh + 1); j < enyh; ++j)
       for (i = 1; i < enxh; ++i) {
         k1 = i;
         k2 = j;
         k3 = k;
         g1 = k1 * unitg[0] + k2 * unitg[3] + k3 * unitg[6];
         g2 = k1 * unitg[1] + k2 * unitg[4] + k3 * unitg[7];
         g3 = k1 * unitg[2] + k2 * unitg[5] + k3 * unitg[8];
         gg = g1 * g1 + g2 * g2 + g3 * g3;
 
         if ((gg - ggcut) < (-eps)) {
           if (dutask == tid)
             enpack++;
           dutask = (dutask + 1) % tnp;
         }
       }
   enpack_all = ewaldparall->ISumAll(0, enpack);
 
   encut = control.ewald_ncut();
   enshl3d = (2 * encut + 1) * (2 * encut + 1) * (2 * encut + 1);
   ercut = control.ewald_rcut();
   
   constexpr double pi     = units::PI;

   pi4 = 4.00 * pi;
   if (encut <= 0)
     encut = 1;
   if (ercut <= 0.00) {
     rs = unita[0] * unita[0] + unita[1] * unita[1] + unita[2] * unita[2];
     rs = sqrt(rs);
     ercut = rs / pi;
 
     rs = unita[3] * unita[3] + unita[4] * unita[4] + unita[5] * unita[5];
     rs = sqrt(rs);
     w = rs / pi;
     if (w < ercut)
       ercut = w;
 
     rs = unita[6] * unita[6] + unita[7] * unita[7] + unita[8] * unita[8];
     rs = sqrt(rs);
     w = rs / pi;
     if (w < ercut)
       ercut = w;
   }
   w = 0.25 * ercut * ercut;
 
   /* allocate memory */
   eG = new double[3 * enpack];
   vg = new double[enpack];
   ss = new double[2*enpack];
   exi = new double[2*enpack];
   tmp3 = new double[enpack];
   ftmp = new double[3 * (ewaldion->nion)];
   vcx = new double[enpack];
   rcell = new double[3 * enshl3d];
   ewx1 = new double[2 * (ewaldion->nion) * enx];
   ewy1 = new double[2 * (ewaldion->nion) * eny];
   ewz1 = new double[2 * (ewaldion->nion) * enz];
   zv = new double[ewaldion->nkatm];
   i_indx = new int[enpack];
   j_indx = new int[enpack];
   k_indx = new int[enpack];
 
   /* determine eG */
   for (i = 0; i < (3 * enpack); ++i)
     eG[i] = 0.0;
   dutask = 0;
   enpack0 = 0;
   k1 = 0;
   k2 = 0;
   k3 = 0;
   g1 = k1 * unitg[0] + k2 * unitg[3] + k3 * unitg[6];
   g2 = k1 * unitg[1] + k2 * unitg[4] + k3 * unitg[7];
   g3 = k1 * unitg[2] + k2 * unitg[5] + k3 * unitg[8];
   gg = g1 * g1 + g2 * g2 + g3 * g3;
   if ((gg - ggcut) < (-eps)) {
     if (dutask == tid) {
       eG[enpack0] = g1;
       eG[enpack0 + enpack] = g2;
       eG[enpack0 + 2 * enpack] = g3;
       i = k1;
       if (i < 0)
         i += enx;
       j = k2;
       if (j < 0)
         j += eny;
       k = k3;
       if (k < 0)
         k += enz;
       i_indx[enpack0] = i;
       j_indx[enpack0] = j;
       k_indx[enpack0] = k;
       enpack0++;
     }
     dutask = (dutask + 1) % tnp;
   }
   k1 = 0;
   k2 = 0;
   for (k3 = 1; k3 < enzh; ++k3) {
     g1 = k1 * unitg[0] + k2 * unitg[3] + k3 * unitg[6];
     g2 = k1 * unitg[1] + k2 * unitg[4] + k3 * unitg[7];
     g3 = k1 * unitg[2] + k2 * unitg[5] + k3 * unitg[8];
     gg = g1 * g1 + g2 * g2 + g3 * g3;
     if ((gg - ggcut) < (-eps)) {
       if (dutask == tid) {
         eG[enpack0] = g1;
         eG[enpack0 + enpack] = g2;
         eG[enpack0 + 2 * enpack] = g3;
         i = k1;
         if (i < 0)
           i += enx;
         j = k2;
         if (j < 0)
           j += eny;
         k = k3;
         if (k < 0)
           k += enz;
         i_indx[enpack0] = i;
         j_indx[enpack0] = j;
         k_indx[enpack0] = k;
         enpack0++;
       }
       dutask = (dutask + 1) % tnp;
     }
   }
 
   k1 = 0;
   for (k3 = (-enzh + 1); k3 < enzh; ++k3)
     for (k2 = 1; k2 < enyh; ++k2) {
       g1 = k1 * unitg[0] + k2 * unitg[3] + k3 * unitg[6];
       g2 = k1 * unitg[1] + k2 * unitg[4] + k3 * unitg[7];
       g3 = k1 * unitg[2] + k2 * unitg[5] + k3 * unitg[8];
       gg = g1 * g1 + g2 * g2 + g3 * g3;
       if ((gg - ggcut) < (-eps)) {
         if (dutask == tid) {
           eG[enpack0] = g1;
           eG[enpack0 + enpack] = g2;
           eG[enpack0 + 2 * enpack] = g3;
           i = k1;
           if (i < 0)
             i += enx;
           j = k2;
           if (j < 0)
             j += eny;
           k = k3;
           if (k < 0)
             k += enz;
           i_indx[enpack0] = i;
           j_indx[enpack0] = j;
           k_indx[enpack0] = k;
           enpack0++;
         }
         dutask = (dutask + 1) % tnp;
       }
     }
 
   for (k3 = (-enzh + 1); k3 < enzh; ++k3)
     for (k2 = (-enyh + 1); k2 < enyh; ++k2)
       for (k1 = 1; k1 < enxh; ++k1) {
         g1 = k1 * unitg[0] + k2 * unitg[3] + k3 * unitg[6];
         g2 = k1 * unitg[1] + k2 * unitg[4] + k3 * unitg[7];
         g3 = k1 * unitg[2] + k2 * unitg[5] + k3 * unitg[8];
         gg = g1 * g1 + g2 * g2 + g3 * g3;
 
         if ((gg - ggcut) < (-eps)) {
           if (dutask == tid) {
             eG[enpack0] = g1;
             eG[enpack0 + enpack] = g2;
             eG[enpack0 + 2 * enpack] = g3;
             i = k1;
             if (i < 0)
               i += enx;
             j = k2;
             if (j < 0)
               j += eny;
             k = k3;
             if (k < 0)
               k += enz;
             i_indx[enpack0] = i;
             j_indx[enpack0] = j;
             k_indx[enpack0] = k;
             enpack0++;
           }
           dutask = (dutask + 1) % tnp;
         }
       }
 
   /* find vg and vcx */
   for (i = 0; i < enpack; ++i)
     vg[i] = 0.0;
   for (i = 0; i < enpack; ++i)
     vcx[i] = 0.0;
 
   for (k = enida; k < enpack; ++k) {
     g1 = eG[k];
     g2 = eG[k + enpack];
     g3 = eG[k + 2 * enpack];
     gg = g1 * g1 + g2 * g2 + g3 * g3;
     term = pi4 / gg;
     vcx[k] = term;
     vg[k] = term * exp(-w * gg);
   }
 
   /* set the mandelung constant */
   alpha = mandelung_get(ewaldlattice);
 
   /* set the ion charges */
   for (i = 0; i < (ewaldion->nkatm); ++i)
     zv[i] = inzv[i];
   // zv[i] = inpsp->zv[i];
 
   /* ewald summation */
   rs = pow(3.0 * ewaldlattice->omega() / pi4, 1.0 / 3.0);
   zz = 0.0;
   z = 0.0;
   for (i = 0; i < (ewaldion->nion); ++i) {
     q = zv[ewaldion->katm[i]];
     zz += q * q;
     z += q;
   }
   cewald = 0.0;
   for (i = 0; i < enpack; ++i)
     cewald += vg[i];
   cewald *= 2.0;
   if (tnp > 1)
     cewald = ewaldparall->SumAll(0, cewald);
 
   cewald = -0.50 * zz * (alpha / rs + cewald / ewaldlattice->omega()) -
            0.50 * (z * z - zz) * ercut * ercut * pi / ewaldlattice->omega();
 
   /* set rcell */
   l = 0;
   rcell[l] = 0.0;
   rcell[l + enshl3d] = 0.0;
   rcell[l + 2 * enshl3d] = 0.0;
   for (k = -encut; k <= encut; ++k)
     for (j = -encut; j <= encut; ++j)
       for (i = -encut; i <= encut; ++i)
         if (!((i == 0) && (j == 0) && (k == 0))) {
           ++l;
           rcell[l] = i * unita[0] + j * unita[3] + k * unita[6];
           rcell[l + enshl3d] = i * unita[1] + j * unita[4] + k * unita[7];
           rcell[l + 2 * enshl3d] = i * unita[2] + j * unita[5] + k * unita[8];
         }
}

/*********************************
 *                               *
 *          Ewald::phafac        *
 *                               *
 *********************************/
void Ewald::phafac() {
  int i, k, enxh, enyh, enzh;
  double a, b, sw1, sw2, sw3;
  double cw1x, cw2x, cw3x;
  double cw1y, cw2y, cw3y;
  
  constexpr double pi     = units::PI;
  
  enxh = enx / 2;
  enyh = eny / 2;
  enzh = enz / 2;

  for (i = 0; i < (ewaldion->nion); ++i) {
    sw1 = unitg[0] * ewaldion->rion1[0 + 3 * i] +
          unitg[1] * ewaldion->rion1[1 + 3 * i] +
          unitg[2] * ewaldion->rion1[2 + 3 * i] + pi;
    sw2 = unitg[3] * ewaldion->rion1[0 + 3 * i] +
          unitg[4] * ewaldion->rion1[1 + 3 * i] +
          unitg[5] * ewaldion->rion1[2 + 3 * i] + pi;
    sw3 = unitg[6] * ewaldion->rion1[0 + 3 * i] +
          unitg[7] * ewaldion->rion1[1 + 3 * i] +
          unitg[8] * ewaldion->rion1[2 + 3 * i] + pi;

    cw1x = cos(sw1);
    cw1y = -sin(sw1);
    cw2x = cos(sw2);
    cw2y = -sin(sw2);
    cw3x = cos(sw3);
    cw3y = -sin(sw3);

    ewx1[2 * i * enx] = 1.0;
    ewx1[2 * i * enx + 1] = 0.0;
    ewy1[2 * i * eny] = 1.0;
    ewy1[2 * i * eny + 1] = 0.0;
    ewz1[2 * i * enz] = 1.0;
    ewz1[2 * i * enz + 1] = 0.0;
    for (k = 1; k <= enxh; ++k) {
      a = ewx1[2 * (k - 1 + i * enx)];
      b = ewx1[2 * (k - 1 + i * enx) + 1];
      ewx1[2 * (k + i * enx)] = a * cw1x - b * cw1y;
      ewx1[2 * (k + i * enx) + 1] = a * cw1y + b * cw1x;
      ewx1[2 * (enx - k + i * enx)] = ewx1[2 * (k + i * enx)];
      ewx1[2 * (enx - k + i * enx) + 1] = -ewx1[2 * (k + i * enx) + 1];
    }
    for (k = 1; k <= enyh; ++k) {
      a = ewy1[2 * (k - 1 + i * eny)];
      b = ewy1[2 * (k - 1 + i * eny) + 1];
      ewy1[2 * (k + i * eny)] = a * cw2x - b * cw2y;
      ewy1[2 * (k + i * eny) + 1] = a * cw2y + b * cw2x;
      ewy1[2 * (eny - k + i * eny)] = ewy1[2 * (k + i * eny)];
      ewy1[2 * (eny - k + i * eny) + 1] = -ewy1[2 * (k + i * eny) + 1];
    }
    for (k = 1; k <= enzh; ++k) {
      a = ewz1[2 * (k - 1 + i * enz)];
      b = ewz1[2 * (k - 1 + i * enz) + 1];
      ewz1[2 * (k + i * enz)] = a * cw3x - b * cw3y;
      ewz1[2 * (k + i * enz) + 1] = a * cw3y + b * cw3x;
      ewz1[2 * (enz - k + i * enz)] = ewz1[2 * (k + i * enz)];
      ewz1[2 * (enz - k + i * enz) + 1] = -ewz1[2 * (k + i * enz) + 1];
    }

    ewx1[2 * (enxh + i * enx)] = 0.0;
    ewx1[2 * (enxh + i * enx) + 1] = 0.0;
    ewy1[2 * (enyh + i * eny)] = 0.0;
    ewy1[2 * (enyh + i * eny) + 1] = 0.0;
    ewz1[2 * (enzh + i * enz)] = 0.0;
    ewz1[2 * (enzh + i * enz) + 1] = 0.0;
  }
}

void ewald_strfac_add_sub(const int npack, const int indxi[], const int indxj[],
                          const int indxk[], const double exi[],
                          const double exj[], const double exk[],
                          const double alpha, double strx[]) 
{
   double ai, aj, ak, c, d;
   double bi, bj, bk;
   for (int i=0; i<npack; ++i) 
   {
      ai = exi[2*indxi[i]];
      bi = exi[2*indxi[i]+1];
      aj = exj[2*indxj[i]];
      bj = exj[2*indxj[i]+1];
      ak = exk[2*indxk[i]];
      bk = exk[2*indxk[i]+1];
      c = aj * ak - bj * bk;
      d = aj * bk + ak * bj;
      strx[2*i]   += alpha * (ai * c - bi * d);
      strx[2*i+1] += alpha * (ai * d + bi * c);
   }
}

/*********************************
 *                               *
 *          Ewald::energy        *
 *                               *
 *********************************/
double Ewald::energy() 
{
   int i, j, k, l, nion, tnp, tid, dutask;
   double x, y, z, dx, dy, dz, zz, r, w;
   double etmp1, etmp2, eall;
 
   tnp = ewaldparall->np();
   tid = ewaldparall->taskid();
   nion = ewaldion->nion;
 
   for (k = 0; k < (2 * enpack); ++k)
     ss[k] = 0.0;
 
   for (i = 0; i < nion; ++i) 
   {
      ewald_strfac_add_sub(enpack, i_indx, j_indx, k_indx, &ewx1[2 * i * enx],
                           &ewy1[2 * i * eny], &ewz1[2 * i * enz],
                           zv[ewaldion->katm[i]], ss);
   }
   etmp1 = 0.0;
   for (k = 0; k < (enpack); ++k) {
     x = ss[2 * k] * ss[2 * k];
     y = ss[2 * k + 1] * ss[2 * k + 1];
     etmp1 += (x + y) * vg[k];
   }
   if (tnp > 1)
     etmp1 = ewaldparall->SumAll(0, etmp1);
   etmp1 = etmp1 / ewaldlattice->omega() + cewald;
 
   dutask = 0;
   etmp2 = 0.0;
   for (i = 0; i < (nion - 1); ++i)
     for (j = i + 1; j < nion; ++j) {
       if (dutask == tid) {
         dx = ewaldion->rion1[3 * i] - ewaldion->rion1[3 * j];
         dy = ewaldion->rion1[3 * i + 1] - ewaldion->rion1[3 * j + 1];
         dz = ewaldion->rion1[3 * i + 2] - ewaldion->rion1[3 * j + 2];
         zz = zv[ewaldion->katm[i]] * zv[ewaldion->katm[j]];
         for (l = 0; l < enshl3d; ++l) 
         {
            x = rcell[l] + dx;
            y = rcell[l + enshl3d] + dy;
            z = rcell[l + 2 * enshl3d] + dz;
            r = sqrt(x * x + y * y + z * z);
            w = r / ercut;
            etmp2 += zz * erfc(w) / r;
         }
       }
       dutask = (dutask + 1) % tnp;
     }
   if (tnp > 1)
     etmp2 = ewaldparall->SumAll(0, etmp2);
 
   eall = etmp1 + etmp2;
   return eall;
}


void ewald_strfac_sub(const int npack, const int indxi[], const int indxj[],
                      const int indxk[], const double exi[], const double exj[],
                      const double exk[], double strx[]) 
{
   double ai, aj, ak, c, d;
   double bi, bj, bk;
   for (int i=0; i<npack; ++i) 
   {
      ai = exi[2*indxi[i]];
      bi = exi[2*indxi[i] + 1];
      aj = exj[2*indxj[i]];
      bj = exj[2*indxj[i] + 1];
      ak = exk[2*indxk[i]];
      bk = exk[2*indxk[i] + 1];
      c = aj * ak - bj * bk;
      d = aj * bk + ak * bj;
      strx[2 * i] = (ai * c - bi * d);
      strx[2 * i + 1] = (ai * d + bi * c);
   }
}

void ewald_f_tmp3_sub(const int n, const double e[], const double s[],
                      const double v[], double t[]) 
{
   for (int i = 0; i < n; ++i) 
   {
      t[i] = v[i]*(e[2*i]*s[2*i+1] - e[2*i+1]*s[2*i]);
   }
}

double ewald_f_ddot_sub(const int n, const double a[], const double b[]) 
{
   double sum = 0.0;
   for (int i = 0; i<n; ++i)
      sum += a[i] * b[i];

   return sum;
}

/*********************************
 *                               *
 *          Ewald::force         *
 *                               *
 *********************************/
void Ewald::force(double *fion) 
{
   int i, j, k, l, nion, tnp, tid, dutask;
   double x, y, z, dx, dy, dz, zz, r, w, zi, f;
   double scal2, sw1, sw2, sw3;
   double cerfc = 1.128379167;
 
   scal2 = 1.0 / ewaldlattice->omega();
   tnp = ewaldparall->np();
   tid = ewaldparall->taskid();
   nion = ewaldion->nion;
 
   for (k=0; k<(2*enpack); ++k)
      ss[k] = 0.0;
   for (i=0; i<nion; ++i) 
   {
      ewald_strfac_add_sub(enpack, i_indx, j_indx, k_indx, &ewx1[2 * i * enx],
                           &ewy1[2 * i * eny], &ewz1[2 * i * enz],
                           zv[ewaldion->katm[i]], ss);
   }
   for (i=0; i<nion; ++i) 
   {
      zi = zv[ewaldion->katm[i]];
      ewald_strfac_sub(enpack, i_indx, j_indx, k_indx, &ewx1[2 * i * enx],
                       &ewy1[2 * i * eny], &ewz1[2 * i * enz], exi);
      ewald_f_tmp3_sub(enpack, exi, ss, vg, tmp3);
  
      ftmp[3*i]   = 2.0 * scal2 * zi * ewald_f_ddot_sub(enpack, eG, tmp3);
      ftmp[3*i+1] = 2.0 * scal2 * zi * ewald_f_ddot_sub(enpack, &eG[enpack], tmp3);
      ftmp[3*i+2] = 2.0 * scal2 * zi * ewald_f_ddot_sub(enpack, &eG[2 * enpack], tmp3);
   }
 
   dutask = 0;
   for (i=0; i<(nion-1); ++i)
      for (j=i+1; j<nion; ++j) 
      {
         if (dutask == tid) 
         {
            dx = ewaldion->rion1[3 * i] - ewaldion->rion1[3 * j];
            dy = ewaldion->rion1[3 * i + 1] - ewaldion->rion1[3 * j + 1];
            dz = ewaldion->rion1[3 * i + 2] - ewaldion->rion1[3 * j + 2];
            zz = zv[ewaldion->katm[i]] * zv[ewaldion->katm[j]];
            sw1 = 0.0;
            sw2 = 0.0;
            sw3 = 0.0;
            for (l = 0; l < enshl3d; ++l) 
            {
               x = rcell[l] + dx;
               y = rcell[l + enshl3d] + dy;
               z = rcell[l + 2 * enshl3d] + dz;
               r = sqrt(x * x + y * y + z * z);
               w = r / ercut;
               f = zz * (erfc(w) + cerfc * w * exp(-w * w)) / (r * r * r);
               sw1 += x * f;
               sw2 += y * f;
               sw3 += z * f;
            }
            ftmp[3 * i] += sw1;
            ftmp[3 * i + 1] += sw2;
            ftmp[3 * i + 2] += sw3;
            ftmp[3 * j] -= sw1;
            ftmp[3 * j + 1] -= sw2;
            ftmp[3 * j + 2] -= sw3;
         }
         dutask = (dutask + 1) % tnp;
      }
   if (tnp > 1) ewaldparall->Vector_SumAll(0, 3 * nion, ftmp);
   for (i=0; i < 3 * nion; ++i)
      fion[i] += ftmp[i];
}


/*********************************
 *                               *
 *          Ewald::stress        *
 *                               *
 *********************************/
void Ewald::stress(double *stress) 
{
   constexpr int N = 40;
   constexpr double third  = 1.0 / 3.0;
   constexpr double pi     = units::PI;
   constexpr double fourpi = 4.0*pi;
   constexpr double scal   = 1.0/(2.0*pi);

   int tnp = ewaldparall->np();
   int tid = ewaldparall->taskid();
   int dutask = 0;

   double omega = ewaldlattice->omega();
   double scal2  = 1.0/omega;
   int nion = ewaldion->nion;

   for (int j=0; j<3; ++j)
   for (int i=0; i<3; ++i) 
   {
      unitg[i + j*3] = ewaldlattice->unitg(i, j);
      unita[i + j*3] = ewaldlattice->unita(i, j);
   }


   double hm[9];
   for (size_t i=0; i<3; ++i)
   for (size_t j=0; j<3; ++j)
      hm[i+3*j] = scal*ewaldlattice->unitg(i,j);


   double zz = 0.0;
   double z  = 0.0;
   for (auto i=0; i<nion; ++i)
   {
      double zi = zv[ewaldion->katm[i]];
      zz += zi*zi;
      z  += zi;
   }

   // Miscellaneous contributions - stress from cewald term
   for (size_t v=0; v<3; ++v)
   for (size_t u=0; u<3; ++u)
   {
      stress[u+3*v] =  0.5*z*z*pi*ercut*ercut * hm[u+3*v]*scal2;
   }


   // G-space contributions 

   // --- 1. Setup and Initialization ---
   double etmp2 = 0.0;
   const double ss_inv_omega = 1.0 / omega; // Pre-calculate for speed

   // Use vectors for automatic memory management
   std::vector<double> H_vec(enpack);
   std::vector<double> ex1_vec(2*enpack, 0.0);
   std::vector<double> strf_vec(2*enpack, 0.0);

   // F is a pointer to an existing buffer (e.g., from the physics engine)
   const double* F = vg; 

   // --- 2. G-space Structure Factor Contribution --- Process ions
   for (int i = 0; i < nion; ++i) 
   {
      double zi = zv[ewaldion->katm[i]];
      
      // Legacy call: requires raw pointers to the start of segments
      ewald_strfac_sub(enpack, i_indx, j_indx, k_indx,
                       &ewx1[2*i*enx], &ewy1[2*i*eny], &ewz1[2*i*enz],
                       ex1_vec.data());
      
      // Vectorized accumulation
      for (size_t k = 0; k < 2 * enpack; ++k) 
      {
          strf_vec[k] += ex1_vec[k] * zi;
      }
   }

   // --- 3. Compute Squared Magnitudes (H-array) ---
   for (size_t k = 0; k < enpack; ++k) 
   {
      double a = strf_vec[2*k];
      double b = strf_vec[2*k + 1];
      H_vec[k] = (a * a) + (b * b);
   }


   // Perform the dot product - Assuming ewald_f_ddot_sub is your legacy function
   etmp2 = ewald_f_ddot_sub(enpack, F, H_vec.data());

   // Parallel Reduction
   if (tnp > 1) etmp2 = ewaldparall->SumAll(0, etmp2);

   for (size_t v=0; v<3; ++v)
   for (size_t u=0; u<3; ++u)
   {
      stress[u+3*v] +=  etmp2 * hm[u+3*v];
   }

   // tmp2(G) = F(G)*H(G)/G**2 + F(G)*H(G)*rcut*rcut/4 

   // --- 5. Calculate Scaled Scaling Factors (tmp1 and tmp2) ---
   std::vector<double> tmp1(enpack);
   std::vector<double> tmp2(enpack, 0.0); 


   if (enpack>0)
   {
      // Pre-calculate constants outside the loop
      const double ss1 = 0.25 * ercut * ercut;
      const double ss2 = 1.0 / fourpi;
     
      // LOOP FUSION: Perform all 5 operations in ONE single pass.
      // This is much faster because we only load F[i], H[i], and vcx[i] into the CPU once.
      for (size_t k=0; k<enpack; ++k)
      {
         // Step A: The base product (F * H)
         double fh_product = F[k] * H_vec[k];
     
         // Step B: Calculate tmp2 component: (F*H) * (rcut^2 / 4)
         tmp2[k] = fh_product * ss1;
     
         // Step C: Calculate tmp1 component: (F*H) * (1/4pi) * vcx
         // We do it all in one line to keep the intermediate value in a CPU register
         tmp1[k] = (fh_product * ss2) * vcx[k];
     
         // Step D: Accumulate tmp1 into tmp2
         // This completes the formula: tmp2 = [F*H*ss1] + [F*H*ss2*vcx]
         tmp2[k] += tmp1[k];
      }
   }

   // --- 6. Calculate Cus (Reciprocal Lattice Stress Contribution) ---
   std::array<double, 9> Cus = {0.0}; // Fixed size symmetric tensor

   // We use an array of pointers to represent the 3 segments of eG.
   // This is compatible with C++11/14/17 and avoids std::span errors.
   const double* g_segments[3]; 

   g_segments[0] = eG;                  // Start of segment 1
   g_segments[1] = eG + enpack;         // Start of segment 2
   g_segments[2] = eG + (2 * enpack);   // Start of segment 3

   if (enpack > 0)
   {
    for (size_t u = 0; u < 3; ++u)
    {
        for (size_t s = u; s < 3; ++s)
        {
            // We need a temporary buffer to hold the product of the two segments.
            // Using a vector here ensures automatic memory cleanup.
            std::vector<double> g_prod(enpack);

            // Step A: Element-wise multiplication (G_u * G_s)
            for (size_t i = 0; i < enpack; ++i)
            {
                g_prod[i] = g_segments[u][i] * g_segments[s][i];
            }

            // Step B: Weighted dot product using the legacy function.
            // We pass the pointer to the start of our vector.
            double sum_val = ewald_f_ddot_sub(enpack, g_prod.data(), tmp2.data());

            // Step C: Assign to the symmetric tensor index
            Cus[u + 3 * s] = ss_inv_omega * sum_val;
        }
    }
   }



   // R-space contributions

   // calculate alpha1 - stress from cewald term
   double rs = std::pow((3.0*omega/fourpi),third);
   double epsilon = 1.0/ercut;
   double zz_part = -0.5*zz;
   dutask = 0;
   for (int n1=-N+1; n1<N; ++n1)
   for (int n2=-N+1; n2<N; ++n2)
   for (int n3=-N+1; n3<N; ++n3)
   {
      if (dutask == tid) 
      {
         if ( n1==0  &&  n2==0  &&  n3==0 ) continue;

         double ax = n1*unita[0] + n2*unita[3] + n3*unita[6];
         double ay = n1*unita[1] + n2*unita[4] + n3*unita[7];
         double az = n1*unita[2] + n2*unita[5] + n3*unita[8];
 
         //double ea = std::sqrt(ax*ax + ay*ay + az*az);
         double ea = std::hypot(ax, ay, az); 
         double w  = ea * epsilon;
         double ss = std::erfc(w)/ea + (2.0 * epsilon / std::sqrt(pi)) * std::exp(-w*w);
         ss = (zz_part*ss) / (ea*ea);

         // Update symmetric tensor components
         Cus[0] += ss * ax*ax; // 1,1
         Cus[3] += ss * ax*ay; // 1,2
         Cus[6] += ss * ax*az; // 1,3
         Cus[4] += ss * ay*ay; // 2,2
         Cus[7] += ss * ay*az; // 2,3
         Cus[8] += ss * az*az; // 3,3
      }
      dutask = (dutask + 1) % tnp;
   }

   // Calculate erfc contribution 
   dutask=0;
   for (auto i=0; i<(nion-1); ++i)
   for (auto j=i+1; j<nion; ++j) 
   {
      if (dutask == tid) 
      {
         double dx = ewaldion->rion1[3*i]   - ewaldion->rion1[3*j];
         double dy = ewaldion->rion1[3*i+1] - ewaldion->rion1[3*j+1];
         double dz = ewaldion->rion1[3*i+2] - ewaldion->rion1[3*j+2];
         double zizj = zv[ewaldion->katm[i]] * zv[ewaldion->katm[j]];

         for (auto l=0; l<enshl3d; ++l) 
         {
            double ax = rcell[l]           + dx;
            double ay = rcell[l+enshl3d]   + dy;
            double az = rcell[l+2*enshl3d] + dz;

            double ea = std::hypot(ax, ay, az); 
            if (ea>1.0e-6)
            {
               double w  = ea * epsilon;
               double ss = -std::erfc(w)/ea - (2.0 * epsilon / std::sqrt(pi)) * std::exp(-w*w);
               ss = ss/(ea*ea);
               Cus[0] +=  ss * ax*ax * zizj;
               Cus[3] +=  ss * ax*ay * zizj;
               Cus[6] +=  ss * ax*az * zizj;
               Cus[4] +=  ss * ay*ay * zizj;
               Cus[7] +=  ss * ay*az * zizj;
               Cus[8] +=  ss * az*az * zizj;
            }
         }
      }
      dutask = (dutask + 1) % tnp;
   }
   if (tnp > 1) ewaldparall->Vector_SumAll(0, 9, Cus.data());

   for (size_t u=0; u<3; ++u)
   for (size_t s=u+1; s<3; ++s)
       Cus[s+3*u] = Cus[u+3*s];

   // Calculate stress = Sum(s) hm(s,v)*Cus(u,s)
   for (size_t v=0; v<3; ++v)
   for (size_t u=0; u<3; ++u)
      for (size_t s=0; s<3; ++s)
         stress[u+3*v] += Cus[u+3*s]*hm[s+3*v];

}

} // namespace pwdft

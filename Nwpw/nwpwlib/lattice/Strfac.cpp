/* Strfac.C -
   Author - Eric Bylaska
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>

#include "Strfac.hpp"

namespace pwdft {

/* Constructors */

/*********************************
 *                               *
 *         Strfac::Strfac        *
 *                               *
 *********************************/

Strfac::Strfac(Ion *inion, PGrid *ingrid) {
  int i, j, k, l, k1, k2, k3;
  int nx, ny, nz, nxh, nb, p, indx;
  int tnp, tid;
  int *ii_indx, *jj_indx, *kk_indx;

  mygrid = ingrid;
  myion = inion;
  tnp = mygrid->parall->np_i();
  tid = mygrid->parall->taskid_i();
  Lattice *lattice = mygrid->lattice;

  nx = mygrid->nx;
  ny = mygrid->ny;
  nz = mygrid->nz;
  nxh = nx / 2;

  for (j = 0; j < 3; ++j)
    for (i = 0; i < 3; ++i) {
      unitg[i + j * 3] = lattice->unitg(i, j);
      unita[i + j * 3] = lattice->unita(i, j);
    }

  /* allocate memory */
  wx1 = new double[2 * (myion->nion) * (mygrid->nx)];
  wy1 = new double[2 * (myion->nion) * (mygrid->ny)];
  wz1 = new double[2 * (myion->nion) * (mygrid->nz)];

  i_indx[0] = new int[mygrid->npack(0)];
  j_indx[0] = new int[mygrid->npack(0)];
  k_indx[0] = new int[mygrid->npack(0)];
  i_indx[1] = new int[mygrid->npack(1)];
  j_indx[1] = new int[mygrid->npack(1)];
  k_indx[1] = new int[mygrid->npack(1)];

  ii_indx = new int[mygrid->nfft3d];
  jj_indx = new int[mygrid->nfft3d];
  kk_indx = new int[mygrid->nfft3d];
  for (nb = 0; nb < 2; ++nb) {
    for (k = 0; k < nz; ++k)
      for (j = 0; j < ny; ++j)
        for (i = 0; i <= nxh; ++i) {
          p = mygrid->ijktop(i, j, k);
          if (p == tid) {
            indx = mygrid->ijktoindex(i, j, k);
            ii_indx[indx] = i;
            jj_indx[indx] = j;
            kk_indx[indx] = k;
          }
        }
    mygrid->i_pack(nb, ii_indx);
    mygrid->i_pack(nb, jj_indx);
    mygrid->i_pack(nb, kk_indx);
    mygrid->ii_pack_copy(nb, ii_indx, i_indx[nb]);
    mygrid->ii_pack_copy(nb, jj_indx, j_indx[nb]);
    mygrid->ii_pack_copy(nb, kk_indx, k_indx[nb]);
  }

  delete[] kk_indx;
  delete[] jj_indx;
  delete[] ii_indx;
}

/*********************************
 *                               *
 *          Strfac::phafac       *
 *                               *
 *********************************/
void Strfac::phafac()
{
   int i, k, nxh, nyh, nzh, nx, ny, nz;
   double a, b, sw1, sw2, sw3, pi;
   double cw1x, cw2x, cw3x;
   double cw1y, cw2y, cw3y;
 
   pi = 4.00 * std::atan(1.0);
 
   nx = (mygrid->nx);
   ny = (mygrid->ny);
   nz = (mygrid->nz);
   nxh = nx / 2;
   nyh = ny / 2;
   nzh = nz / 2;
 
   for (i = 0; i < (myion->nion); ++i) 
   {
      sw1 = unitg[0]*myion->rion1[0+3*i] +
            unitg[1]*myion->rion1[1+3*i] +
            unitg[2]*myion->rion1[2+3*i] + pi;
      sw2 = unitg[3]*myion->rion1[0+3*i] +
            unitg[4]*myion->rion1[1+3*i] +
            unitg[5]*myion->rion1[2+3*i] + pi;
      sw3 = unitg[6]*myion->rion1[0+3*i] +
            unitg[7]*myion->rion1[1+3*i] +
            unitg[8]*myion->rion1[2+3*i] + pi;
     
      cw1x =  std::cos(sw1);
      cw1y = -std::sin(sw1);
      cw2x =  std::cos(sw2);
      cw2y = -std::sin(sw2);
      cw3x =  std::cos(sw3);
      cw3y = -std::sin(sw3);
      wx1[2*i*nx]   = 1.0;
      wx1[2*i*nx+1] = 0.0;
      wy1[2*i*ny]   = 1.0;
      wy1[2*i*ny+1] = 0.0;
      wz1[2*i*nz]   = 1.0;
      wz1[2*i*nz+1] = 0.0;
      for (k=1; k<=nxh; ++k) 
      {
         a = wx1[2 * (k - 1 + i * nx)];
         b = wx1[2 * (k - 1 + i * nx) + 1];
         wx1[2 * (k + i * nx)] = a * cw1x - b * cw1y;
         wx1[2 * (k + i * nx) + 1] = a * cw1y + b * cw1x;
         wx1[2 * (nx - k + i * nx)] = wx1[2 * (k + i * nx)];
         wx1[2 * (nx - k + i * nx) + 1] = -wx1[2 * (k + i * nx) + 1];
      }
     
      for (k=1; k<=nyh; ++k)
      {
         a = wy1[2 * (k - 1 + i * ny)];
         b = wy1[2 * (k - 1 + i * ny) + 1];
         wy1[2 * (k + i * ny)] = a * cw2x - b * cw2y;
         wy1[2 * (k + i * ny) + 1] = a * cw2y + b * cw2x;
         wy1[2 * (ny - k + i * ny)] = wy1[2 * (k + i * ny)];
         wy1[2 * (ny - k + i * ny) + 1] = -wy1[2 * (k + i * ny) + 1];
      }
      for (k=1; k<=nzh; ++k)
      {
         a = wz1[2 * (k - 1 + i * nz)];
         b = wz1[2 * (k - 1 + i * nz) + 1];
         wz1[2 * (k + i * nz)] = a * cw3x - b * cw3y;
         wz1[2 * (k + i * nz) + 1] = a * cw3y + b * cw3x;
         wz1[2 * (nz - k + i * nz)] = wz1[2 * (k + i * nz)];
         wz1[2 * (nz - k + i * nz) + 1] = -wz1[2 * (k + i * nz) + 1];
      }
     
      wx1[2*(nxh+i*nx)]     = 0.0;
      wx1[2*(nxh+i*nx) + 1] = 0.0;
      wy1[2*(nyh+i*ny)]     = 0.0;
      wy1[2*(nyh+i*ny) + 1] = 0.0;
      wz1[2*(nzh+i*nz)]     = 0.0;
      wz1[2*(nzh+i*nz) + 1] = 0.0;
   }
}

/*********************************
 *                               *
 *       Strfac::strfac_pack     *
 *                               *
 *********************************/
void Strfac::strfac_pack(const int nb, const int ii, double *strx) {
   int npack, nx, ny, nz;
   npack = mygrid->npack(nb);
   nx = mygrid->nx;
   ny = mygrid->ny;
   nz = mygrid->nz;
 
   const int *indxi = i_indx[nb];
   const int *indxj = j_indx[nb];
   const int *indxk = k_indx[nb];
 
   const double *exi = &wx1[2*ii*nx];
   const double *exj = &wy1[2*ii*ny];
   const double *exk = &wz1[2*ii*nz];
 
   double ai, aj, ak, bi, bj, bk;
   double c, d;
   for (int i = 0; i < npack; ++i)
   {
      ai = exi[2*indxi[i]];
      bi = exi[2*indxi[i] + 1];
      aj = exj[2*indxj[i]];
      bj = exj[2*indxj[i] + 1];
      ak = exk[2*indxk[i]];
      bk = exk[2*indxk[i] + 1];
      c = aj*ak - bj*bk;
      d = aj*bk + ak*bj;
      strx[2*i]   = (ai*c - bi*d);
      strx[2*i+1] = (ai*d + bi*c);
   }
}

} // namespace pwdft

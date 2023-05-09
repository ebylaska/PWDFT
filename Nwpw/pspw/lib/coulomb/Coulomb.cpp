/* Coulomb.cpp -
   Author - Eric Bylaska
*/

/**
 * class: Coulomb_Operator
 *
 * brief: Class for computing the Coulomb operator in electronic structure calculations.
 *
 * The Coulomb operator is an integral operator that appears in many electronic
 * structure methods. This class provides methods to compute the Coulomb operator
 * and its associated energy, given the density in reciprocal space. The Coulomb
 * operator is defined as:
 *
 *   V(r) = ∫ d³r' ρ(r')/|r - r'|,
 *
 * where ρ(r') is the electron density at position r'.
 *
 * This class requires a Pneb object to define the grid on which the density is
 * represented, and a Control2 object to provide some options for the computation.
 *
 * author: Eric J. Bylaska
 * date: 5/9/2023
 */

/*
#include        <cstdio>
#include        <cstdlib>
#include        <iostream>
#include        <stdio.h>
#include	<string>
*/
#include "Coulomb.hpp"
#include <cmath>

namespace pwdft {

/*******************************************
 *                                         *
 *     Coulomb_Operator::Coulomb_Operator  *
 *                                         *
 *******************************************/
Coulomb_Operator::Coulomb_Operator(Pneb *mygrid, Control2 &control) {
  int k, pzero, zero, taskid;
  double gg;
  double *Gx = mygrid->Gxyz(0);
  double *Gy = mygrid->Gxyz(1);
  double *Gz = mygrid->Gxyz(2);
  vg = new double[mygrid->npack(0)];
  double *tmp = new double[mygrid->nfft3d];
  double fourpi = 16.0 * atan(1.0);
  // double kk0 = 0.0; // where k0 = the uniform background potential

  mypneb = mygrid;

  taskid = mypneb->d3db::parall->taskid_i();
  pzero = mypneb->ijktop(0, 0, 0);
  zero = mypneb->ijktoindex(0, 0, 0);

  for (k = 0; k < (mypneb->nfft3d); ++k) {
    gg = Gx[k] * Gx[k] + Gy[k] * Gy[k] + Gz[k] * Gz[k];
    if ((pzero == taskid) && (k == zero))
      tmp[k] = 0.0;
    else
      tmp[k] = fourpi / gg;
  }
  mypneb->t_pack(0, tmp);
  mypneb->tt_pack_copy(0, tmp, vg);

  delete[] tmp;
}

/*******************************************
 *                                         *
 *        Coulomb_Operator::vcoulomb       *
 *                                         *
 *******************************************/
void Coulomb_Operator::vcoulomb(const double *dng, double *vcout) {
  int k, k1, ksize;

  ksize = (mypneb->npack(0));
  k1 = 0;
  for (k = 0; k < ksize; ++k) {
    vcout[k1] = vg[k] * dng[k1];
    vcout[k1 + 1] = vg[k] * dng[k1 + 1];
    k1 += 2;
  }
}

/*******************************************
 *                                         *
 *        Coulomb_Operator::ecoulomb       *
 *                                         *
 *******************************************/
double Coulomb_Operator::ecoulomb(const double *dng) {
  int k, k1, k2, n, nsize, ksize1, ksize2;
  double ave;

  ksize1 = (mypneb->nzero(0));
  ksize2 = (mypneb->npack(0));
  ave = 0.0;
  k1 = 0;
  for (k = 0; k < ksize1; ++k) {
    ave += vg[k] * (dng[k1] * dng[k1] + dng[k1 + 1] * dng[k1 + 1]);
    k1 += 2;
  }
  for (k = ksize1; k < ksize2; ++k) {
    ave += 2.0 * vg[k] * (dng[k1] * dng[k1] + dng[k1 + 1] * dng[k1 + 1]);
    k1 += 2;
  }
  ave = mypneb->d3db::parall->SumAll(1, ave);
  // ave *= 0.5*lattice_omega();
  ave *= 0.5 * (mypneb->lattice->omega());

  return ave;
}

} // namespace pwdft

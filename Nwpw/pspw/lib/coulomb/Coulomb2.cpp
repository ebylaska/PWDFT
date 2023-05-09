/* Coulomb2.cpp -
   Author - Eric Bylaska
*/
/**
 * The Coulomb2_Operator class computes the Coulomb interaction energy between two charge distributions in three dimensions
 * using the aperiodic convolution algorithm. The class takes in the positions and charges of the two distributions and
 * the box size (i.e. the size of the periodic cell) as inputs. The class provides a public method, vcoulomb(),
 * which returns the Coulomb potential energy of the two distributions.
 *
 * Author: Eric J. Bylaska
 * Date: 5/9/2023
 */

/*
#include        <cstdio>
#include        <cstdlib>
#include        <iostream>
#include        <stdio.h>
#include	<string>
*/

#include "Control2.hpp"
#include "Coulomb2.hpp"
#include "Pneb.hpp"
#include "Strfac.hpp"
#include <cmath>

namespace pwdft {

/*******************************************
 *                                         *
 *   Coulomb2_Operator::Coulomb2_Operator  *
 *                                         *
 *******************************************/
Coulomb2_Operator::Coulomb2_Operator(Pneb *mygrid, Control2 &control) {
  mypneb = mygrid;

  double pi = 4.0 * std::atan(1.0);
  double twopi = 2.0 * pi;
  double fourpi = 4.0 * pi;

  int taskid = mypneb->d3db::parall->taskid_i();

  // allocated double grid d3db object
  int dnx = 2 * mypneb->nx;
  int dny = 2 * mypneb->ny;
  int dnz = 2 * mypneb->nz;
  int dmaptype = -mypneb->maptype;

  myd3db2 =
      new (std::nothrow) d3db(mypneb->d3db::parall, dmaptype, dnx, dny, dnz);

  dnfft3d = myd3db2->nfft3d;
  dn2ft3d = myd3db2->n2ft3d;
  dscale = 1.0 / ((double)(dnx * dny * dnz));

  /* By default, when the new operator is used to attempt to allocate memory and
     the handling function is unable to do so, a bad_alloc exception is thrown.
      But when nothrow is used as argument for new, it returns a null pointer
     instead.  */
  gk = new (std::nothrow) double[dnfft3d]();
  double *glr = new (std::nothrow) double[dn2ft3d]();
  tmpx = glr;

  double *tmp = new (std::nothrow) double[dnfft3d]();

  // define lattice on expanded grid
  for (auto i = 0; i < 9; ++i)
    dunita[i] = 2 * mypneb->PGrid::lattice->unita1d(i);

  // reciprical vectors for expanded grid
  dunitg[0] = dunita[4] * dunita[8] - dunita[5] * dunita[7];
  dunitg[1] = dunita[5] * dunita[6] - dunita[3] * dunita[8];
  dunitg[2] = dunita[3] * dunita[7] - dunita[4] * dunita[6];
  dunitg[3] = dunita[7] * dunita[2] - dunita[8] * dunita[1];
  dunitg[4] = dunita[8] * dunita[0] - dunita[6] * dunita[2];
  dunitg[5] = dunita[6] * dunita[1] - dunita[7] * dunita[0];
  dunitg[6] = dunita[1] * dunita[5] - dunita[2] * dunita[4];
  dunitg[7] = dunita[2] * dunita[3] - dunita[0] * dunita[5];
  dunitg[8] = dunita[0] * dunita[4] - dunita[1] * dunita[3];
  double volume =
      dunita[0] * dunitg[0] + dunita[1] * dunitg[1] + dunita[2] * dunitg[2];
  for (auto i = 0; i < 9; ++i)
    dunitg[i] *= twopi / volume;
  double domega = fabs(volume);

  // *******************************
  // Set up the Green's function ***
  // *******************************
  int dnxh = dnx / 2;
  int dnyh = dny / 2;
  int dnzh = dnz / 2;

  double fdnx = ((double)dnx);
  double fdny = ((double)dny);
  double fdnz = ((double)dnz);
  double dV = domega / (fdnx * fdny * fdnz);
  double temp;

  // short-range part of Greens function
  std::memset(gk, 0, dnfft3d * sizeof(double));
  for (auto k = 0; k < dnz; ++k)
    for (auto j = 0; j < dny; ++j)
      for (auto i = 0; i <= dnxh; ++i) {
        auto i1 = i;
        auto j1 = j;
        auto k1 = k;

        auto indx = myd3db2->ijktoindex(i, j, k);
        auto p = myd3db2->ijktop(i, j, k);
        if (p == taskid) {
          if ((i1 + j1 + k1) > 0) {
            if (j1 > dnyh)
              j1 -= dny;
            if (k1 > dnzh)
              k1 -= dnz;

            auto gx = i1 * dunitg[0] + j1 * dunitg[3] + k1 * dunitg[6];
            auto gy = i1 * dunitg[1] + j1 * dunitg[4] + k1 * dunitg[7];
            auto gz = i1 * dunitg[2] + j1 * dunitg[5] + k1 * dunitg[8];
            auto gg = gx * gx + gy * gy + gz * gz;

            temp = (fourpi / gg) *
                   (1.0 - std::exp(-gg / (4.00 * EPSILON * EPSILON)));
          } else
            temp = pi / (EPSILON * EPSILON);

          gk[indx] = temp;
        }
      }

  // long-range part of Greens function
  std::memset(glr, 0, dn2ft3d * sizeof(double));

  if (mypneb->PGrid::lattice->fast_erf()) {
    // Error function parameters
    double xerf, yerf;
    double c1 = 0.07052307840;
    double c2 = 0.04228201230;
    double c3 = 0.00927052720;
    double c4 = 0.00015201430;
    double c5 = 0.00027656720;
    double c6 = 0.00004306380;

    for (auto k = 0; k < dnz; ++k)
      for (auto j = 0; j < dny; ++j)
        for (auto i = 0; i < dnx; ++i) {
          auto k1 = k;
          auto j1 = j;
          auto i1 = i;
          if (k1 >= dnz / 2)
            k1 -= dnz;
          if (j1 >= dny / 2)
            j1 -= dny;
          if (i1 >= dnx / 2)
            i1 -= dnx;

          // call D3dB_ijktoindex2p(2,i,j,k,index1,p)
          auto indx = myd3db2->ijktoindex2(i, j, k);
          auto p = myd3db2->ijktop2(i, j, k);
          if (p == taskid) {
            auto x = i1 * dunita[0] / fdnx + j1 * dunita[3] / fdny +
                     k1 * dunita[6] / fdnz;
            auto y = i1 * dunita[1] / fdnx + j1 * dunita[4] / fdny +
                     k1 * dunita[7] / fdnz;
            auto z = i1 * dunita[2] / fdnx + j1 * dunita[5] / fdny +
                     k1 * dunita[8] / fdnz;
            auto temp = std::sqrt(x * x + y * y + z * z);
            if (temp > 1.0e-15) {
              auto xerf = EPSILON * temp;
              auto yerf1 =
                  (1.00 +
                   xerf * (c1 +
                           xerf * (c2 +
                                   xerf * (c3 +
                                           xerf * (c4 +
                                                   xerf * (c5 + xerf * c6))))));
              auto yerf4 = yerf1 * yerf1 * yerf1 * yerf1;
              auto yerf16 = yerf4 * yerf4 * yerf4 * yerf4;
              yerf = (1.0 - 1.0 / yerf16);
              temp = yerf / temp;
            } else
              temp = 2.0 * EPSILON / std::sqrt(pi);

            glr[indx] = temp * dV;
          }
        }
  } else {
    for (auto k = 0; k < dnz; ++k)
      for (auto j = 0; j < dny; ++j)
        for (auto i = 0; i < dnx; ++i) {
          auto k1 = k;
          auto j1 = j;
          auto i1 = i;
          if (k1 >= dnz / 2)
            k1 -= dnz;
          if (j1 >= dny / 2)
            j1 -= dny;
          if (i1 >= dnx / 2)
            i1 -= dnx;

          // call D3dB_ijktoindex2p(2,i,j,k,index1,p)
          auto indx = myd3db2->ijktoindex2(i, j, k);
          auto p = myd3db2->ijktop2(i, j, k);
          if (p == taskid) {
            auto x = i1 * dunita[0] / fdnx + j1 * dunita[3] / fdny +
                     k1 * dunita[6] / fdnz;
            auto y = i1 * dunita[1] / fdnx + j1 * dunita[4] / fdny +
                     k1 * dunita[7] / fdnz;
            auto z = i1 * dunita[2] / fdnx + j1 * dunita[5] / fdny +
                     k1 * dunita[8] / fdnz;
            auto temp = std::sqrt(x * x + y * y + z * z);
            if (temp > 1.0e-15)
              temp = std::erf(EPSILON * temp) / temp;
            else
              temp = 2.0 * EPSILON / std::sqrt(pi);

            glr[indx] = temp * dV;
          }
        }
  }

  // Put glr in k-space
  myd3db2->rc_fft3d(glr);

  // add long-range part to short-range part
  // note that only real parts of tranformed grl are used
  for (auto k = 0; k < dnfft3d; ++k) {
    gk[k] += glr[2 * k];
  }

  delete[] tmp;
}

/*******************************************
 *                                         *
 *       Coulomb2_Operator::vcoulomb       *
 *                                         *
 *******************************************/
/*  This routine calculates Poisson's equation for infinite
   space boundry conditions.

        Laplacian(vh) = -4*pi*dn

        vh(r-->infinity) = 0


    Entry:
          dn --- the density of the region in real-space
    Exit:
          vh --- the solution to Poisson's equation in real-space
*/
void Coulomb2_Operator::vcoulomb(const double *dn, double *vcout) {
  std::memset(tmpx, 0, dn2ft3d * sizeof(double));

  // Expand the density
  myd3db2->hr2r_expand(dn, tmpx);

  // Convolution g*dn
  myd3db2->rc_fft3d(tmpx);
  myd3db2->tc_Mul(gk, tmpx);
  myd3db2->cr_fft3d(tmpx);

  // Contract tmpx to extract vcout
  myd3db2->r2hr_contract(tmpx, vcout);
  mypneb->r_SMul(dscale, vcout);
}

/*
Abramowitz and Stegun give several approximations of varying accuracy
(equations 7.1.25–28). This allows one to choose the fastest approximation
suitable for a given application. In order of increasing accuracy, they are:

erf(x) \approx 1-{\frac
{1}{\left(1+a_{1}x+a_{2}x^{2}+a_{3}x^{3}+a_{4}x^{4}\right)^{4}}}, for x>=0
(maximum error: 5×10−4)
where a1 = 0.278393, a2 = 0.230389, a3 = 0.000972, a4 = 0.078108

erf(x) \approx 1-\left(a_{1}t+a_{2}t^{2}+a_{3}t^{3}\right)e^{-x^{2}},\quad
t={\frac {1}{1+px}}, for x>=0. (maximum error: 2.5×10−5) where p = 0.47047, a1 =
0.3480242, a2 = −0.0958798, a3 = 0.7478556

Approximation used by Ryoichi:
erf(x) \approx 1-{\frac {1}{\left(1+a_{1}x+a_{2}x^{2}+\cdots
+a_{6}x^{6}\right)^{16}}}, for x>= 0 (maximum error: 3×10−7) where a1 =
0.0705230784, a2 = 0.0422820123, a3 = 0.0092705272, a4 = 0.0001520143, a5 =
0.0002765672, a6 = 0.0000430638

erf(x) \approx 1-\left(a_{1}t+a_{2}t^{2}+\cdots
+a_{5}t^{5}\right)e^{-x^{2}},\quad t={\frac {1}{1+px}}} (maximum
error: 1.5×10−7) where p = 0.3275911, a1 = 0.254829592, a2 = −0.284496736, a3
= 1.421413741, a4 = −1.453152027, a5 = 1.061405429

All of these approximations are valid for x ≥ 0. To use these approximations for
negative x, use the fact that erf(x) is an odd function, so erf(x) = −erf(−x).

*/

} // namespace pwdft

/* Lattice.cpp
   Author - Eric Bylaska

*/

#include "Lattice.hpp"
#include "Control2.hpp"
#include <cmath>

namespace pwdft {

static void get_cube(double *unita, double *unitg, double *omega) {
  double twopi = 8.0 * atan(1.0);
  unitg[0] = unita[4] * unita[8] - unita[5] * unita[7];
  unitg[1] = unita[5] * unita[6] - unita[3] * unita[8];
  unitg[2] = unita[3] * unita[7] - unita[4] * unita[6];
  unitg[3] = unita[7] * unita[2] - unita[8] * unita[1];
  unitg[4] = unita[8] * unita[0] - unita[6] * unita[2];
  unitg[5] = unita[6] * unita[1] - unita[7] * unita[0];
  unitg[6] = unita[1] * unita[5] - unita[2] * unita[4];
  unitg[7] = unita[2] * unita[3] - unita[0] * unita[5];
  unitg[8] = unita[0] * unita[4] - unita[1] * unita[3];
  double volume =
      unita[0] * unitg[0] + unita[1] * unitg[1] + unita[2] * unitg[2];
  for (int i = 0; i < 9; ++i)
    unitg[i] *= twopi / volume;
  *omega = fabs(volume);
}

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

Lattice::Lattice(Control2 &control) {
  int nx, ny, nz, nxh, nyh, nzh;
  double gx, gy, gz, gg, gg1, gg2, gg3, ecut0, wcut0;

  ecut0 = control.ecut();
  wcut0 = control.wcut();
  punita[0] = control.unita(0, 0);
  punita[1] = control.unita(1, 0);
  punita[2] = control.unita(2, 0);

  punita[3] = control.unita(0, 1);
  punita[4] = control.unita(1, 1);
  punita[5] = control.unita(2, 1);

  punita[6] = control.unita(0, 2);
  punita[7] = control.unita(1, 2);
  punita[8] = control.unita(2, 2);
  get_cube(punita, punitg, &pomega);

  nx = control.ngrid(0);
  ny = control.ngrid(1);
  nz = control.ngrid(2);
  nxh = nx / 2;
  nyh = ny / 2;
  nzh = nz / 2;

  gx = punitg[0] * ((double)nxh);
  gy = punitg[1] * ((double)nxh);
  gz = punitg[2] * ((double)nxh);
  gg1 = gx * gx + gy * gy + gz * gz;

  gx = punitg[3] * ((double)nyh);
  gy = punitg[4] * ((double)nyh);
  gz = punitg[5] * ((double)nyh);
  gg2 = gx * gx + gy * gy + gz * gz;

  gx = punitg[6] * ((double)nzh);
  gy = punitg[7] * ((double)nzh);
  gz = punitg[8] * ((double)nzh);
  gg3 = gx * gx + gy * gy + gz * gz;

  gg = gg1;
  if (gg2 < gg)
    gg = gg2;
  if (gg3 < gg)
    gg = gg3;

  pecut = 0.50 * gg;
  if (ecut0 < pecut)
    pecut = ecut0;
  pwcut = pecut;
  if (wcut0 < pwcut)
    pwcut = wcut0;

  pfast_erf = control.fast_erf();
  paperiodic = (control.version == 4);
}

/********************************
 *                              *
 *           abc_abg            *
 *                              *
 ********************************/
void Lattice::abc_abg(double *a1, double *b1, double *c1, double *alpha1,
                      double *beta1, double *gamma1) {
  double pi = 4.00 * atan(1.00);
  double a = sqrt(pow(unita(0, 0), 2.0) + pow(unita(1, 0), 2.0) +
                  pow(unita(2, 0), 2.0));
  double b = sqrt(pow(unita(0, 1), 2.0) + pow(unita(1, 1), 2.0) +
                  pow(unita(2, 1), 2.0));
  double c = sqrt(pow(unita(0, 2), 2.0) + pow(unita(1, 2), 2.0) +
                  pow(unita(2, 2), 2.0));

  double d2 = pow((unita(0, 1) - unita(0, 2)), 2.0) +
              pow((unita(1, 1) - unita(1, 2)), 2.0) +
              pow((unita(2, 1) - unita(2, 2)), 2.0);
  double alpha = (b * b + c * c - d2) / (2.00 * b * c);
  alpha = acos(alpha) * 180.00 / pi;

  d2 = pow((unita(0, 2) - unita(0, 0)), 2.0) +
       pow((unita(1, 2) - unita(1, 0)), 2.0) +
       pow((unita(2, 2) - unita(2, 0)), 2.0);
  double beta = (c * c + a * a - d2) / (2.00 * c * a);
  beta = acos(beta) * 180.00 / pi;

  d2 = pow((unita(0, 0) - unita(0, 1)), 2.0) +
       pow((unita(1, 0) - unita(1, 1)), 2.0) +
       pow((unita(2, 0) - unita(2, 1)), 2.0);
  double gamma = (a * a + b * b - d2) / (2.00 * a * b);
  gamma = acos(gamma) * 180.00 / pi;

  /* return values */
  *a1 = a;
  *b1 = b;
  *c1 = c;
  *alpha1 = alpha;
  *beta1 = beta;
  *gamma1 = gamma;
}

} // namespace pwdft


#include "util_radial_switching.hpp"
#include <cmath>

namespace pwdft {
using namespace pwdft;

/**************************************
 *                                    *
 *        util_radial_switching       *
 *                                    *
 **************************************/

/* Cubic spline switching function

   entry - Rin,Rout,r: inner,outer,current radiii
   exit - return a value between 0 and 1, r<Rin: 0, r>Rout: 1
*/

double util_radial_switching(const double Rin, const double Rout,
                             const double r) {
  double s = 0.0;
  if (r <= Rin)
    s = 0.0;
  else if (r >= Rout)
    s = 1.0;
  else {
    double dR = Rout - Rin;
    double c1 = dR * dR * dR;
    double c2 = 3.0 * Rout - Rin;
    s = dR * dR * (c2 - 2.00 * r) / c1;
  }

  return s;
}

/**************************************
 *                                    *
 *        util_radial_dswitching      *
 *                                    *
 **************************************/

/* Derivative of hubic spline switching function

   entry - Rin,Rout,r: inner,outer,current radiii
   exit - return a value between 0 and 1, r<Rin: 0, r>Rout: 1 and its derivative
*/
void util_radial_dswitching(const double Rin, const double Rout, const double r,
                            double *S, double *dS) {
  *S = 0.0;
  *dS = 0.0;
  if (r <= Rin) {
    *S = 0.0;
    *dS = 0.0;
  } else if (r >= Rout) {
    *S = 1.0;
    *dS = 0.0;
  } else {
    double dR = Rout - Rin;
    double du = r - Rin;
    double c1 = dR * dR * dR;
    double c2 = 3.0 * Rout - Rin;
    *S = dR * dR * (c2 - 2.00 * r) / c1;
    *dS = 2.0 * du * ((c2 - 2.0 * r) / c1) - 2.0 * (du * du) / c1;
  }
}

} // namespace pwdft

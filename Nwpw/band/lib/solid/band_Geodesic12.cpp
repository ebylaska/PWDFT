/* band_Geodesic12.cpp -
   Author - Eric Bylaska
*/

#include "band_Geodesic12.hpp"

namespace pwdft {

/*******************************************
 *                                         *
 *    band_Geodesic12::band_Geodesic12     *
 *                                         *
 *******************************************/
band_Geodesic12::band_Geodesic12(int minimizer0, Solid *mysolid0, Control2 &control) {
  int minimizer = control.minimizer();
  if ((minimizer == 1) || (minimizer == 2)) {
    has_geodesic1 = true;
    mygeodesic1 = new (std::nothrow) band_Geodesic(minimizer0, mysolid0);
  }

  if ((minimizer == 4) || (minimizer == 7)) {
    has_geodesic2 = true;
    mygeodesic2 = new (std::nothrow) band_Geodesic2(minimizer0, mysolid0);
  }
}

} // namespace pwdft

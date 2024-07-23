#pragma once

#include "Control2.hpp"
#include "band_Geodesic.hpp"
#include "band_Geodesic2.hpp"
#include "Cneb.hpp"

namespace pwdft {

class band_Geodesic12 {

public:
  bool has_geodesic1 = false;
  band_Geodesic *mygeodesic1;

  bool has_geodesic2 = false;
  band_Geodesic2 *mygeodesic2;

  /* Constructors */
  band_Geodesic12(int, Solid *, Control2 &);

  /* destructor */
  ~band_Geodesic12() {
    if (has_geodesic1)
      delete mygeodesic1;
    if (has_geodesic2)
      delete mygeodesic2;
  }
};

} // namespace pwdft

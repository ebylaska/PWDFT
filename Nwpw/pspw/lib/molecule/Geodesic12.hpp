#pragma once

#include "Control2.hpp"
#include "Geodesic.hpp"
#include "Geodesic2.hpp"
#include "Pneb.hpp"

namespace pwdft {

class Geodesic12 {

public:
  bool has_geodesic1 = false;
  Geodesic *mygeodesic1;

  bool has_geodesic2 = false;
  Geodesic2 *mygeodesic2;

  /* Constructors */
  Geodesic12(int, Molecule *, Control2 &);

  /* destructor */
  ~Geodesic12() {
    if (has_geodesic1)
      delete mygeodesic1;
    if (has_geodesic2)
      delete mygeodesic2;
  }
};

} // namespace pwdft

#ifndef _BRILLOUIN_HPP_
#define _BRILLOUIN_HPP_

#include <cmath>
#include <cstdlib>
#include <cstring>

#include "Control2.hpp"
#include "Lattice.hpp"
#include "Brillouin.hpp"

namespace pwdft {

class Brillouin {


public:
   int nbrillouin;
   double *weight;
   double *kvector, *ksvector;

   /* Constructors */
   // Ion(RTDB&, Control2&);
   Brillouin(std::string, Lattice *, Control2 &);
 
   /* destructor */
   ~Brillouin() {
      delete [] kvector;
      delete [] ksvector;
      delete [] weight;
   }
 
   /* functions */
 
   std::string print_zone();
   std::string print_zone_point(const int);


};
} // namespace pwdft

#endif

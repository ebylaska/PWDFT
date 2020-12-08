#ifndef _LATTICE_HPP_
#define _LATTICE_HPP_
/* Lattice.hpp
   Author - Eric Bylaska

*/

#include	"Control2.hpp"

class Lattice {

   float punita[9],punitg[9],pecut,pwcut,pomega;

public:

   /* constructor */
   Lattice(Control2&);

   float unita1d(const int i) { return punita[i]; }
   float unita(const int i, const int j) { return punita[i+j*3]; }
   float unitg(const int i, const int j) { return punitg[i+j*3]; }
   float ecut() { return pecut; }
   float wcut() { return pwcut; }
   float omega() { return pomega; }
   float eggcut() { return 2*pecut; }
   float wggcut() { return 2*pwcut; }
};

#endif

#ifndef _LATTICE_HPP_
#define _LATTICE_HPP_

#pragma once

/* Lattice.hpp
   Author - Eric Bylaska

*/

#include	"Control2.hpp"

namespace pwdft {


class Lattice {

   bool   pfast_erf,paperiodic;
   double punita[9],punitg[9],pecut,pwcut,pomega;

public:

   /* constructor */
   Lattice(Control2&);

   double unita1d(const int i) { return punita[i]; }
   double unita(const int i, const int j) { return punita[i+j*3]; }
   double unitg(const int i, const int j) { return punitg[i+j*3]; }
   double ecut() { return pecut; }
   double wcut() { return pwcut; }
   double omega() { return pomega; }
   double eggcut() { return 2*pecut; }
   double wggcut() { return 2*pwcut; }

   double *unita_ptr() { return punita; } 
   void abc_abg(double *, double *, double *, double *, double *, double *);

   bool fast_erf() { return pfast_erf;}
   bool aperiodic() { return paperiodic; }


};
}

#endif

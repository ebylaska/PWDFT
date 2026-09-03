#ifndef _LATTICE_HPP_
#define _LATTICE_HPP_

#pragma once

/* Lattice.hpp
   Author - Eric Bylaska

*/

#include "Control2.hpp"

namespace pwdft {

/**
 * @class Lattice
 * @brief Represents a lattice in a DFT simulation.
 */
class Lattice {

   bool pfast_erf, paperiodic;
   double punita[9], punitg[9], pub[9], pecut, pwcut, pomega;
   double punita_frozen[9], punitg_frozen[9], pub_frozen[9],  pecut_frozen, pwcut_frozen, pomega_frozen;

public:
   /* constructor */
   Lattice(Control2 &);
 

   //lattice operations
   double unita1d(const int i) { return punita[i]; }
   double unitg1d(const int i) { return punitg[i]; }
   double unita(const int i, const int j) { return punita[i+j*3]; }
   double unitg(const int i, const int j) { return punitg[i+j*3]; }
   double ub(const int i, const int j)    { return pub[i+j*3]; }
   double ecut() { return pecut; }
   double wcut() { return pwcut; }
   double omega() { return pomega; }
   double eggcut() { return 2 * pecut; }
   double wggcut() { return 2 * pwcut; }
 
   double *unita_ptr() { return punita; }
   double *unitg_ptr() { return punitg; }
   double *ub_ptr() { return pub; }

   void abc_abg(double *, double *, double *, double *, double *, double *);
   void min_diff_xyz(double *, double *, double *);
   void min_diff(double *);

   double unita_frozen1d(const int i) { return punita_frozen[i]; }
   double unitg_frozen1d(const int i) { return punitg_frozen[i]; }
   double unita_frozen(const int i, const int j) { return punita_frozen[i+j*3]; }
   double unitg_frozen(const int i, const int j) { return punitg_frozen[i+j*3]; }
   double ub_frozen(const int i, const int j)    { return pub_frozen[i+j*3]; }
   double ecut_frozen() { return pecut_frozen; }
   double wcut_frozen() { return pwcut_frozen; }
   double omega_frozen() { return pomega_frozen; }
   double eggcut_frozen() { return 2 * pecut_frozen; }
   double wggcut_frozen() { return 2 * pwcut_frozen; }

   double *unita_frozen_ptr() { return punita_frozen; }
   double *unitg_frozen_ptr() { return punitg_frozen; }
   double *ub_frozen_ptr() { return pub_frozen; }
 
   bool fast_erf() { return pfast_erf; }
   bool aperiodic() { return paperiodic; }

   // update lattice
   void update_unita_keep_basis(const double unita_new[9]);
   

};
} // namespace pwdft

#endif

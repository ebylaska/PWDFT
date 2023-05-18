#ifndef _ION_BOND_HPP_
#define _ION_BOND_HPP_

#pragma once

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "iofmt.hpp"
#include "Control2.hpp"

namespace pwdft {

class ion_bond
{

   bool bond_exists=false;;
   bool periodic;
   int nhb;
   int    *i0,*j0;
   double *K0;
   double *R0;
  
   double ua[9],ub[9];
   double *rion;


public:

  /* Constructors */
  ion_bond(double *, Control2 &);


  /* destructor */
  ~ion_bond() {
      if (bond_exists) {
         delete [] i0;
         delete [] j0;
         delete [] K0;
         delete [] R0;
      }
  }

  /* functions */
  bool has_bond() {return bond_exists;}

  void min_diff_xyz(double *, double *, double *);
  double spring_energy(const int) ;
  double energy();
  double energyfion(double *);
  std::string print_all(const int);
};
} // namespace pwdft

#endif

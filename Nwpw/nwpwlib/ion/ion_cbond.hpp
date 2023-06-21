#ifndef _ION_CBOND_HPP_
#define _ION_CBOND_HPP_

#pragma once

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "iofmt.hpp"
#include "Control2.hpp"

namespace pwdft {

class ion_cbond
{

   bool cbond_exists=false;;
   bool periodic;
   int nhcb;
   int    *i0,*j0,*k0,*l0;
   double *Kspring0;
   double *Rij0, *Rkl0;
  
   double ua[9],ub[9];
   double *rion;


public:

  /* Constructors */
  ion_cbond(double *, Control2 &);

  /* destructor */
  ~ion_cbond() {
      if (cbond_exists) {
         delete [] i0;
         delete [] j0;
         delete [] k0;
         delete [] l0;
         delete [] Kspring0;
         delete [] Rij0;
         delete [] Rkl0;
      }
  }

  /* functions */
  bool has_cbond() {return cbond_exists;}

  void min_diff_xyz(double *, double *, double *);
  double spring_energy(const int) ;
  double energy();
  double energyfion(double *);
  std::string print_all(const int);
};
} // namespace pwdft

#endif

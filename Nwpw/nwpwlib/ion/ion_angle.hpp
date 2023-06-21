#ifndef _ION_ANGLE_HPP_
#define _ION_ANGLE_HPP_

#pragma once

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "iofmt.hpp"
#include "Control2.hpp"

namespace pwdft {

class ion_angle
{

   bool angle_exists=false;;
   bool periodic;
   int nha;
   int    *i0,*j0,*k0;
   double *Kspring0;
   double *Theta0;
  
   double ua[9],ub[9];
   double *rion;


public:

  /* Constructors */
  ion_angle(double *, Control2 &);


  /* destructor */
  ~ion_angle() {
      if (angle_exists) {
         delete [] i0;
         delete [] j0;
         delete [] k0;
         delete [] Kspring0;
         delete [] Theta0;
      }
  }

  /* functions */
  bool has_angle() {return angle_exists;}

  void min_diff_xyz(double *, double *, double *);
  double spring_energy(const int) ;
  double energy();
  double energyfion(double *);
  std::string print_all(const int);
};
} // namespace pwdft

#endif

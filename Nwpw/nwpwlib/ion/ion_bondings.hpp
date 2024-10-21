#ifndef _ION_BONDINGS_HPP_
#define _ION_BONDINGS_HPP_

#pragma once

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>

#include "iofmt.hpp"
#include "Control2.hpp"

namespace pwdft {

class ion_bondings 
{

   bool bondings_exists=false;
   bool periodic;
   int nhc;
   int    *n0;
   double *Kspring0;
   double *gamma0;

   int    **indx;
   double **coef;
   double ua[9],ub[9];
   double *rion;


public:

  /* Constructors */
  ion_bondings(double *, Control2 &);


  /* destructor */
  ~ion_bondings() {
      if (bondings_exists) {
         for (auto i=0; i<nhc; ++i) {
            delete [] coef[i];
            delete [] indx[i];
         }
         delete [] indx;
         delete [] coef;
         delete [] Kspring0;
         delete [] gamma0;
         delete [] n0;
      }
  }

  /* functions */
  bool has_bondings() {return bondings_exists;}

  void min_diff_xyz(double *, double *, double *);
  double ij_spring_distance(const int, const int);
  double spring_gamma(const int);
  double spring_energy(const int) ;
  double energy();
  double energyfion(double *);
  std::string print_all(const int);
};
} // namespace pwdft

#endif

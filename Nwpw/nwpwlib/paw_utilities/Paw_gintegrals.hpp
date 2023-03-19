#ifndef _PAW_GINTEGRALS_HPP_
#define _PAW_GINTEGRALS_HPP_

#pragma once

#include        "Parallel.hpp"
#include        "Ion.hpp"
#include        "Ewald.hpp"
#include        "Control2.hpp"

namespace pwdft {



class Paw_gintegrals {
   int nthr=1;
   int ngauss_max;

   Parallel *myparall;
   Ion      *myion;
   Ewald    *myewald;

   int lm_size_max;
   int nion_paw,*katm_paw,*mult_l, *ion_pawtoion;
   double *sigma_paw, sigma_smooth;

public:

   bool periodic;
   int ngauss;

   int *tgauss, *tgauss_shift;                          /* used for threading */
   int *lm1_gauss, *lm2_gauss, *iii1_gauss, *iii2_gauss;   /* indexing used for Gaussian integrals */
   double *e_gauss, *f_gauss;                           /* Gaussian integrals */

   /* Constructors */
   Paw_gintegrals(Parallel *, Ion *, Ewald *, const bool,
                  const int, const int, int *, int *, int *, double *, const double);

   /* destructor */
   ~Paw_gintegrals() {
      delete [] tgauss;
      delete [] tgauss_shift;
      delete [] lm1_gauss;
      delete [] lm2_gauss;
      delete [] iii1_gauss;
      delete [] iii2_gauss;
      delete [] e_gauss;
      delete [] f_gauss;
   }

   /* sets the gaussian integrals */
   void set(const bool);

};

}

#endif


#ifndef _COULOMB2_HPP_
#define _COULOMB2_HPP_

#pragma once

#include "Control2.hpp"
#include "Dielectric.hpp"
#include	"Pneb.hpp"

namespace pwdft {

class	Coulomb2_Operator {

   double *gk;
   double *tmpx;
   double dunita[9],dunitg[9],dscale;
   int    dnfft3d,dn2ft3d;

   Pneb   *mypneb;
   d3db   *myd3db2;

   /* cutoff constant */
   double EPSILON=1.0;

   bool has_dielec = false;
   Dielectric_Operator *mydielec;


public:

   /* Constructors */
   Coulomb2_Operator(Pneb *, Control2&);

   /* destructor */
   ~Coulomb2_Operator() {
         delete [] gk;
         delete [] tmpx;
         delete myd3db2;

         if (has_dielec) delete mydielec;
    }

    void   vcoulomb(const double *, double *);
};

}

#endif

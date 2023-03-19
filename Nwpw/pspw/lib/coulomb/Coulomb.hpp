#ifndef _COULOMB_HPP_
#define _COULOMB_HPP_

#pragma once

#include "Control2.hpp"
#include	"Pneb.hpp"
#include "Dielectric.hpp"

namespace pwdft {

class	Coulomb_Operator {

   double *vg;
   Pneb   *mypneb;

   bool has_dielec = false;
   Dielectric_Operator *mydielec;

public:

   /* Constructors */
   Coulomb_Operator(Pneb *, Control2&);

   /* destructor */
   ~Coulomb_Operator() {
         delete [] vg;
         if (has_dielec) delete mydielec;
    }

    void   vcoulomb(double *, double *);
    double ecoulomb(double *);
};

}

#endif

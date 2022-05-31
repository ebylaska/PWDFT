#pragma once

#include	"Pneb.hpp"

namespace pwdft {

class	Coulomb_Operator {

   double *vg;
   Pneb   *mypneb;

public:

   /* Constructors */
   Coulomb_Operator(Pneb *);

   /* destructor */
   ~Coulomb_Operator() {
         delete [] vg;
    }

    void   vcoulomb(double *, double *);
    double ecoulomb(double *);
};

}

#ifndef _DIELECTRIC_HPP_
#define _DIELECTRIC_HPP_

#pragma once

#include "Control2.hpp"
#include	"Pneb.hpp"

namespace pwdft {

class	Dielectric_Operator {

   Pneb   *mypneb;
   int    n2ft3d;
   double dielec, rho0, beta;

public:
   double *epsilon, *sw, *depsilondx, *depsilondy, *depsilondz;
   double *r_grid;

   /* Constructors */
   Dielectric_Operator(Pneb *, Control2&);

   /* destructor */
   ~Dielectric_Operator() {
         mypneb->r_dealloc(epsilon);
         mypneb->r_dealloc(sw);
         mypneb->r_dealloc(depsilondx);
         mypneb->r_dealloc(depsilondy);
         mypneb->r_dealloc(depsilondz);

         mypneb->r_dealloc(r_grid);
    }

    void   rho_generate(double *);
};

}

#endif

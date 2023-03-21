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
   double *r_grid;

public:
   double *epsilon, *sw, *p;

   /* Constructors */
   Dielectric_Operator(Pneb *, Control2&);

   /* destructor */
   ~Dielectric_Operator() {
       mypneb->r_dealloc(epsilon);
       mypneb->r_dealloc(sw);
       mypneb->r_dealloc(p);
    }

    void generate_dielec(const double *);
    void generate_scaled(double *);
};

}

#endif

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
   double *epsilon_x, *epsilon_y, *epsilon_z;
   double *w_x, *w_y, *w_z;

   /* Constructors */
   Dielectric_Operator(Pneb *, Control2&);

   /* destructor */
   ~Dielectric_Operator() {
       mypneb->r_dealloc(epsilon);
       mypneb->r_dealloc(sw);
       mypneb->r_dealloc(p);
       mypneb->r_dealloc(epsilon_x);
       mypneb->r_dealloc(epsilon_y);
       mypneb->r_dealloc(epsilon_z);
       mypneb->r_dealloc(w_x);
       mypneb->r_dealloc(w_y);
       mypneb->r_dealloc(w_z);
    }

    void generate_dielec(const double *);
    void generate_scaled(double *);
    void generate_over_epsilon(double *);
    void generate_dpotential(const double *, double *);
};

}

#endif

#ifndef _NWPW_DIPOLE_HPP_
#define _NWPW_DIPOLE_HPP_

// ********************************************************************
// *                                                                  *
// *       nwpw_dipole: used to generate aperiodic and periodic       *
// *                  dipoles from a plane-wave density.              *
// *                                                                  *
// *   The algorithms used in this module are based on the work of    *
// *                                                                  *
// *                                                                  *
// ********************************************************************
#include        <iostream>
#include        "blas.h"

#include        "gdevice.hpp"
#include        "Control2.hpp"
#include        "Ion.hpp"
#include        "Pneb.hpp"
#include        "Strfac.hpp"


namespace pwdft {


class nwpw_dipole {

   Pneb   *mypneb;
   Ion    *myion;
   Strfac *mystrfac;

public:

   bool dipole_on;

   double *gamma,*A,*Am,*b,*q,*u,*utmp,*w,*gaus,*vtmp;
   double *qion, *uion;
   double Eapc,Papc;

   /* constructor */
   nwpw_dipole(Ion *, Pneb *, Strfac *, Control2&, std::ostream&);

   /* destructor */
   ~nwpw_dipole() {
      if (dipole_on) {
         delete [] gamma;
         delete [] A;
         delete [] Am;
         delete [] b;
         delete [] q;
         delete [] u;
         delete [] w;
         delete [] gaus;
         delete [] vtmp;
         delete [] qion;
         delete [] uion;
      }
   }

   void gen_molecular_dipole(double *, double *);
   void gen_aperiodic_dipole(double *, double *);
   void gen_periodic_dipole(double *, double *);

   double Qtot_APC(const int);
   std::string shortprint_APC();
   std::string print_APC(const double *);


};

}

#endif



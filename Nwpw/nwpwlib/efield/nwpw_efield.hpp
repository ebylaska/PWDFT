#ifndef _NWPW_EFIELD_HPP_
#define _NWPW_EFIELD_HPP_

// ********************************************************************
// *                                                                  *
// *       nwpw_efield: used to generate aperiodic and periodic       *
// *                  efield  potentials                              *
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


class nwpw_efield {

   Pneb   *mypneb;
   Ion    *myion;
   Strfac *mystrfac;
   double dv;
   int n2ft3d,ispin,*ne;

public:

   double autoDebye = 2.5416;

   double mdipole[3],mcdv1[3],mcdv2[3],mcdv3[3],mqv1[3];
   bool dipole_on;


   /* constructor */
   nwpw_efield(Ion *, Pneb *, Strfac *, Control2&);

   /* destructor */
   ~nwpw_efield() {
   }

   //void gen_dipole(const double *);
   //void gen_Resta_dipole(const double *, double *);

   //double Qtot_APC(const int);
   //std::string print_APC(const double *);

   std::string shortprint_efield();

};

}

#endif



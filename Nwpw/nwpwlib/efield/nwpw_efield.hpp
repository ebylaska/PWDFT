#ifndef _NWPW_EFIELD_HPP_
#define _NWPW_EFIELD_HPP_

#pragma once

// ********************************************************************
// *                                                                  *
// *       nwpw_efield: used to generate aperiodic and periodic       *
// *                  efield  potentials                              *
// *                                                                  *
// *   The algorithms used in this module are based on the work of    *
// *                                                                  *
// *                                                                  *
// ********************************************************************
#include "blas.h"
#include <iostream>

#include "Control2.hpp"
#include "Ion.hpp"
#include "Pneb.hpp"
#include "Strfac.hpp"
#include "gdevice.hpp"

namespace pwdft {

class nwpw_efield {

  Pneb *mypneb;
  Ion *myion;
  Strfac *mystrfac;
  double dv;
  int n2ft3d, ispin, *ne;

public:
  double autoDebye = 2.5416;

  bool efield_on;
  int efield_type = 2; // default efield_type = rgrid
  double efield_vector[3] = {0.0, 0.0, 0.0};
  double efield_center[3] = {0.0, 0.0, 0.0};

  double *v_field;

  /* constructor */
  nwpw_efield(Ion *, Pneb *, Strfac *, Control2 &, std::ostream &);

  /* destructor */
  ~nwpw_efield() {
    if (efield_on)
      mypneb->r_dealloc(v_field);
  }

  // void gen_dipole(const double *);
  // void gen_Resta_dipole(const double *, double *);

  // double Qtot_APC(const int);
  // std::string print_APC(const double *);

  /**********************************
   *                                *
   *        efield_ion_energy       *
   *                                *
   **********************************/
  // calculates the energy between the QM ions and efield. Note the ions are
  // positive.
  double efield_ion_energy() {
    double eetmp = 0.0;
    if (efield_on) {
       for (auto ii=0; ii<myion->nion; ++ii) {
          double qii = myion->zv_psp[myion->katm[ii]];
          eetmp += qii*efield_vector[0]*(myion->rion1[3*ii]   - efield_center[0]) 
                 + qii*efield_vector[1]*(myion->rion1[3*ii+1] - efield_center[1])
                 + qii*efield_vector[2]*(myion->rion1[3*ii+2] - efield_center[2]);
       }
    }
    return eetmp;
  }

  /**********************************
   *                                *
   *         efield_ion_fion        *
   *                                *
   **********************************/
  // calculates the forces between the QM ions and efield. Note the ions are
  // positive.
  void efield_ion_fion(double *fion) {
     if (efield_on && ((efield_type==0)||(efield_type==2)) ) {
        for (auto ii=0; ii<myion->nion; ++ii) {
           double qii = myion->zv_psp[myion->katm[ii]];
           fion[3*ii]   -= qii*efield_vector[0];
           fion[3*ii+1] -= qii*efield_vector[1];
           fion[3*ii+2] -= qii*efield_vector[2];
        }
     }
  }

  /**********************************
   *                                *
   *         shortprint_efield      *
   *                                *
   **********************************/
  std::string shortprint_efield();
};

} // namespace pwdft

#endif

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
#include "blas.h"
#include <iostream>

#include "Control2.hpp"
#include "Ion.hpp"
#include "Pneb.hpp"
#include "Strfac.hpp"
//#include "gdevice.hpp"

namespace pwdft {

class nwpw_dipole {

  Pneb *mypneb;
  Ion *myion;
  Strfac *mystrfac;
  double dv;
  int n2ft3d, ispin, *ne;

public:
  double autoDebye = 2.5416;

  double mdipole[3], mcdv1[3], mcdv2[3], mcdv3[3], mqv1[3];
  bool dipole_on;

  /* constructor */
  nwpw_dipole(Ion *, Pneb *, Strfac *, Control2 &);

  /* destructor */
  ~nwpw_dipole() {}

  void gen_dipole(const double *);
  void gen_Resta_dipole(const double *, double *);

  void gen_molecular_dipole(const double *, double *);

  // double Qtot_APC(const int);
  // std::string print_APC(const double *);

  std::string shortprint_dipole();
};

} // namespace pwdft

#endif

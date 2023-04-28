#ifndef _COULOMB12_HPP_
#define _COULOMB12_HPP_

#pragma once

#include "Control2.hpp"
#include "Coulomb.hpp"
#include "Coulomb2.hpp"
#include "Ion.hpp"
#include "Pneb.hpp"
#include "Strfac.hpp"

#include "nwpw_dplot.hpp"

namespace pwdft {

class Coulomb12_Operator {
  Pneb *mypneb;


  // Dielectric variables
  bool has_dielec = false;
  double dielec, rhomin, rhomax, beta, rho0;
  double *epsilon, *depsilon, *ddepsilon, *sw, *p;
  double *epsilon_x, *epsilon_y, *epsilon_z, *epsilon_lap;
  double *w_x, *w_y, *w_z;
  double *rho_ind0, *rho_ind1;
  double tole_pol = 1.0e-7;
  double alpha_pol;
  int    model_pol, maxit_pol;

  double *rho_ion, *dng_ion, *v_ion, rcut_ion;

  Ion    *myion;
  Strfac *mystrfac;

public:
  bool has_coulomb1 = false;
  Coulomb_Operator *mycoulomb1;

  bool has_coulomb2 = false;
  Coulomb2_Operator *mycoulomb2;

  /* Constructors */
  Coulomb12_Operator(Pneb *, Control2 &);

  /* destructor */
  ~Coulomb12_Operator() {
    if (has_coulomb1)
      delete mycoulomb1;
    if (has_coulomb2)
      delete mycoulomb2;
    if (has_dielec) {
      mypneb->r_dealloc(epsilon);
      mypneb->r_dealloc(depsilon);
      mypneb->r_dealloc(ddepsilon);
      mypneb->r_dealloc(epsilon_x);
      mypneb->r_dealloc(epsilon_y);
      mypneb->r_dealloc(epsilon_z);
      mypneb->r_dealloc(epsilon_lap);
      mypneb->r_dealloc(sw);
      mypneb->r_dealloc(p);
      mypneb->r_dealloc(w_x);
      mypneb->r_dealloc(w_y);
      mypneb->r_dealloc(w_z);
      mypneb->r_dealloc(rho_ind0);
      mypneb->r_dealloc(rho_ind1);
      mypneb->r_dealloc(rho_ion);
    }
  }

  void initialize_dielectric(Ion *, Strfac *);

  double v_dielectric(const double *, const double *, const double *,
                      const double *, double *);

  void v_dielectric_aperiodic(const double *, const double *, const double *,
                              double *, double *, double *, bool, double *);


  double v_dielectric2_aperiodic(const double *, const double *, const double *,
                                 const double *, const double *, const double *,
                                 const bool, double *, double *, nwpw_dplot *);

  void generate_dng_ion(Pneb *, Ion *, Strfac *, double, double *);
};

} // namespace pwdft

#endif

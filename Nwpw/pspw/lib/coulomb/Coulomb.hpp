#ifndef _COULOMB_HPP_
#define _COULOMB_HPP_

#pragma once

#include "Control2.hpp"
#include "Pneb.hpp"

namespace pwdft {

class Coulomb_Operator {

  double *vg;
  Pneb *mypneb;

public:
  /* Constructors */
  Coulomb_Operator(Pneb *, Control2 &);

  /* destructor */
  ~Coulomb_Operator() { delete[] vg; }

  void vcoulomb(const double *, double *);
  double ecoulomb(const double *);
  // void   vcoulomb_dielec(const double *, double *);
  // void   vcoulomb_dielec2(const double *, double *, double *);
};

} // namespace pwdft

#endif

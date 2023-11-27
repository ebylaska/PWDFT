#ifndef _cCOULOMB_HPP_
#define _cCOULOMB_HPP_

#pragma once

#include "Control2.hpp"
#include "Cneb.hpp"

namespace pwdft {

class cCoulomb_Operator {

  double *vg;
  Cneb *mycneb;

public:
  /* Constructors */
  cCoulomb_Operator(Cneb *, Control2 &);

  /* destructor */
  ~cCoulomb_Operator() { delete[] vg; }

  void vcoulomb(const double *, double *);
  double ecoulomb(const double *);
  // void   vcoulomb_dielec(const double *, double *);
  // void   vcoulomb_dielec2(const double *, double *, double *);
};

} // namespace pwdft

#endif

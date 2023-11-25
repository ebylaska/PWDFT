#ifndef _cKINETIC_HPP_
#define _cKINETIC_HPP_

#include "Cneb.hpp"

namespace pwdft {

class cKinetic_Operator {

  double *tg;
  Cneb *mycneb;

public:
  /* Constructors */
  cKinetic_Operator(Cneb *);

  /* destructor */
  ~cKinetic_Operator() { delete[] tg; }

  void ke(const double *, double *);
  double ke_ave(const double *);
};

} // namespace pwdft

#endif

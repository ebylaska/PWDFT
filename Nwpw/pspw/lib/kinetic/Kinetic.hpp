#ifndef _KINETIC_HPP_
#define _KINETIC_HPP_

#include "Pneb.hpp"

namespace pwdft {

class Kinetic_Operator {

  double *tg;
  Pneb *mypneb;

public:
  /* Constructors */
  Kinetic_Operator(Pneb *);

  /* destructor */
  ~Kinetic_Operator() { delete[] tg; }

  void ke(double *, double *);
  double ke_ave(double *);
};

} // namespace pwdft

#endif

#ifndef _KINETIC_HPP_
#define _KINETIC_HPP_

#include "Pneb.hpp"

namespace pwdft {

class Kinetic_Operator {

  double *tg;
  Pneb *mypneb;

  void rebuild_kinetic_coefficients();

public:
  /* Constructors */
  explicit Kinetic_Operator(Pneb *);

  /* destructor */
  ~Kinetic_Operator() { delete[] tg; }

  void ke(double *, double *);
  void ke_orb(double *, double *);
  //double ke_ave(double *);
  //double ke_ave(double *, double *);
  double ke_ave(double *psi, double *occ = nullptr);

  void ke_euv(double *psi, double *stress, double *occ = nullptr);
 
  void ke_precondition(const double, const int, double *, double*);

  void update_lattice_keep_basis();

};

} // namespace pwdft

#endif

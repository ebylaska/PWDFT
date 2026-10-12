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
   double ke_ave(const double *, const double *);
   void ke_orb(const int, const double *, double *);

   void ke_precondition(const int, const double, const int, double *, double*);

   // Kinetic contribution to the lattice stress tensor.
   // euv(u,v) = -sum_s Aus(u,s) * hm(s,v)
   // where hm = (1/2π) * unitg and Aus is the k- and G-weighted
   // second moment of the density.
   void euv(const double *psi, double euv[3][3]);
   void euv(const double *psi, const double *occ, double euv[3][3]);

};

} // namespace pwdft

#endif

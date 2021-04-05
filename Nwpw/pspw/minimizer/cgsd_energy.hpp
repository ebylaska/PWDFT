#ifndef _CGSD_ENERGY_HPP_
#define _CGSD_ENERGY_HPP_

#include        "Molecule.hpp"

extern double cgsd_noit_energy(Molecule&);
extern double cgsd_energy(Control2&, Molecule&);
extern void   cgsd_energy_gradient(Molecule&, double *);


#endif

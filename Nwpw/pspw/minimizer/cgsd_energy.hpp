#ifndef _CGSD_ENERGY_HPP_
#define _CGSD_ENERGY_HPP_
namespace pwdft {

#include        "Molecule.hpp"

extern double cgsd_noit_energy(Molecule&, bool, std::ostream&);
extern double cgsd_energy(Control2&, Molecule&, bool, std::ostream&);
extern void   cgsd_energy_gradient(Molecule&, double *);


}
#endif

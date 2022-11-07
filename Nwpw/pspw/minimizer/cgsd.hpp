#ifndef _CGSD_HPP_
#define _CGSD_HPP_
namespace pwdft {

#include        "Molecule.hpp"
#include        "Geodesic.hpp"
#include        "pspw_lmbfgs.hpp"

extern double cgsd_cgminimize(Molecule&, Geodesic&, double *, double *, double *, int, int, double, double);
extern double cgsd_bfgsminimize(Molecule&, Geodesic&, pspw_lmbfgs&, double *, double *, double *, int, int, double, double);

}
#endif

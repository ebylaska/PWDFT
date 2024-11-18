#ifndef _CGSD_HPP_
#define _CGSD_HPP_

#pragma once

namespace pwdft {

#include "Geodesic.hpp"
#include "Geodesic2.hpp"
#include "Molecule.hpp"
#include "pspw_lmbfgs.hpp"
#include "nwpw_scf_mixing.hpp"

extern double cgsd_cgminimize(Molecule &, Geodesic *, double *, double *,
                              double *, int, int, double, double);
extern double cgsd_bfgsminimize(Molecule &, Geodesic *, pspw_lmbfgs &, double *,
                                double *, double *, int, int, double, double);

extern double cgsd_cgminimize2(Molecule &, Geodesic2 *, double *, double *,
                               double *, int, int, double, double);
extern double cgsd_bfgsminimize2(Molecule &, Geodesic2 *, pspw_lmbfgs2 &,
                                 double *, double *, double *, int, int, double,
                                 double);

extern double cgsd_bybminimize(Molecule &, Geodesic *, double *, double *,
                              double *, int, int,
                              int, double, double, int,  
                              double, double);

extern double cgsd_bybminimize2(Molecule &, Geodesic *, double *, double *,
                              double *, int, int,int,
                              double, double);

extern double cgsd_bybminimize0(Molecule &, Geodesic *, double *, double *,
                              double *, int, int,int,
                              double, double);

} // namespace pwdft
#endif

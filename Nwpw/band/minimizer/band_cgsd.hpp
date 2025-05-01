#ifndef _BAND_CGSD_HPP_
#define _BAND_CGSD_HPP_

#pragma once

namespace pwdft {

#include "band_Geodesic.hpp"
#include "band_Geodesic2.hpp"
#include "Solid.hpp"
#include "band_lmbfgs.hpp"

extern double band_cgsd_cgminimize(Solid &, band_Geodesic *, double *, double *,
                                   double *, int, int, double, double);

extern double band_cgsd_cgksminimize(Solid &, band_Geodesic *, double *, double *,
                                     double *, int, int, double, double);

extern double band_cgsd_bfgsminimize(Solid &, band_Geodesic *, band_lmbfgs &, double *,
                                     double *, double *, int, int, double, double);

extern double band_cgsd_cgminimize2(Solid &, band_Geodesic2 *, double *, double *,
                                    double *, int, int, double, double);
extern double band_cgsd_bfgsminimize2(Solid &, band_Geodesic2 *, band_lmbfgs2 &,
                                      double *, double *, double *, int, int, double, double);

extern double band_cgsd_bybminimize2(Solid &, band_Geodesic *, double *, double *,
                                     double *, int, int, int, double, double);

} // namespace pwdft
#endif

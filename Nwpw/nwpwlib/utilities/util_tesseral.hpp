#ifndef _UTIL_TESSERAL_HPP_
#define _UTIL_TESSERAL_HPP_

#pragma once

#include <complex>

namespace pwdft {

extern std::complex<double> util_rSphericalHarmonic3(const int, const int,
                                                     const double, const double,
                                                     const double);

extern void util_Tesseral3_vector_lm(const int, const int, const int,
                                     const double *, const double *,
                                     const double *, double *);

extern void util_Tesseral3_rgrid_lm(const int, const int, const int,
                                    const double *, double *);

extern void util_Tesseral3_rgrid_lm_rl(const int, const int, const int,
                                       const double *, double *);

extern void util_Tesseral3_rgrid_lm_over_rl(const int, const int, const int,
                                            const double *, double *);

extern double util_Tesseral3_lm(const int, const int, const double,
                                const double, const double);

extern double util_Tesseral3_lm_over_rl(const int, const int, const double,
                                        const double, const double);

extern void util_dTesseral3_lm(const int, const int, const double, const double,
                               const double, double &, double &, double &);

extern void util_dTesseral3_lm_over_rl(const int, const int, const double,
                                       const double, const double, double &,
                                       double &, double &);

extern void dTesseral3_rgrid_lm_rl(const int, const int, const int,
                                   const double *, double *);

extern double util_Tesseral_lm(const int, const int, const double,
                               const double);

extern std::complex<double> util_YSpherical_lm(const int, const int,
                                               const double, const double);

extern void util_dTesseral_lm(const int, const int, const double, const double,
                              double &, double &, double &);

extern void util_dYspherical_lm(const int, const int, const double,
                                const double, std::complex<double> &,
                                std::complex<double> &, std::complex<double> &);

} // namespace pwdft
#endif

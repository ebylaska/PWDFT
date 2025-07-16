#ifndef _UTIL_HPP_
#define _UTIL_HPP_

#pragma once

#include "Parallel.hpp"
#include "util_gamma.hpp"
#include "util_legendre.hpp"
#include "util_log_integrate.hpp"
#include "util_tesseral.hpp"
//#include	"util_wgaussian.hpp"
#include "util_linesearch.hpp"
#include <cstdlib>
#include <string>

inline std::string get_initial_wavefunction_guess() {
    const char* env = std::getenv("PWDFT_INITIAL_WAVEFUNCTION_GUESS");
    if (env) return std::string(env);
    return "random";
}

namespace pwdft {

extern void transpose2DArray(const int, const int, double *, double *);

extern void c_aindexcopy(const int, const int *, double *, double *);
extern void c_aindexcopy_stride(const int, const int, const int *, double *, double *);

extern void c_bindexcopy(const int, const int *, double *, double *);
extern void c_bindexcopy_stride(const int, const int, const int *, double *, double *);

extern void c_bindexcopy_conjg(const int, const int *, double *, double *);
extern void c_bindexcopy_conjg_stride(const int, const int, const int *, double *, double *);

extern void c_bindexzero(const int, const int *, double *);

extern void t_aindexcopy(const int, const int *, double *, double *);
extern void t_bindexcopy(const int, const int *, double *, double *);

extern void i_aindexcopy(const int, const int *, int *, int *);

extern void eigsrt(double *, double *, int);
extern double util_random(const int);

extern void util_getfilling(int, int *, int *, double *);

extern bool util_filefind(Parallel *, char *);

extern void util_spline(const double *, const double *, const int, const double,
                        const double, double *, double *);
extern double util_splint(const double *, const double *, const double *,
                          const int, const int, const double);

extern void util_filter(int, double *, double, double *);

extern void util_fattebert_dielec(const int, const double, const double,
                                  const double, const double *, double *);
extern void util_andreussi_dielec(const int, const double, const double,
                                  const double, const double *, double *);
extern void util_andreussi2_dielec(const int, const double, const double,
                                   const double, const double *, double *);
extern void util_dfattebert_dielec(const int, const double, const double,
                                   const double, const double *, double *);
extern void util_dandreussi_dielec(const int, const double, const double,
                                   const double, const double *, double *);
extern void util_dandreussi2_dielec(const int, const double, const double,
                                    const double, const double *, double *);
extern void util_ddandreussi_dielec(const int, const double, const double,
                                    const double, const double *, double *);
extern void util_ddandreussi2_dielec(const int, const double, const double,
                                     const double, const double *, double *);
extern void util_sphere_dielec(const int, const double *, const double *,
                               const double, const double, const double,
                               double *);
extern void util_sphere_gradient_dielec(const int, const double *, const double *,
                                        const double, const double, const double,
                                        double *, double *, double *);

extern void util_weighted_fattebert_dielec(const int, const double,
                                           const double, const double,
                                           const double *, const double *,
                                           double *);
extern double util_switching_function(const double, const double, const double);
extern double util_dswitching_function(const double, const double,
                                       const double);

extern double util_kiril_coulomb_transform(const int, const double, const double, const double);
extern double util_kiril_coulomb_transform0(const int, const double, const double);

extern double util_occupation_distribution(const int, const double);
extern double util_smearcorrection(const int, const double, const double, const double, const double);

extern void util_print_elapsed_time(const double);

//#define _MATPRINT_
#ifdef _MATPRINT_
#include <iostream>
extern void util_matprint(std::string, int, double *);
#endif


extern std::complex<double> util_zdotc(int n,
                                       const double* x, int incx,
                                       const double* y, int incy);


} // namespace pwdft

#endif

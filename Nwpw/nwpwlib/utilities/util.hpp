#ifndef _UTIL_H_
#define _UTIL_H_

#include        "Parallel.hpp"

namespace pwdft {
using namespace pwdft;

extern void c_aindexcopy(const int, const int *, double *, double *);
extern void c_bindexcopy(const int, const int *, double *, double *);
extern void c_bindexcopy_conjg(const int, const int *, double *, double *);

extern void t_aindexcopy(const int, const int *, double *, double *);
extern void t_bindexcopy(const int, const int *, double *, double *);

extern void i_aindexcopy(const int, const int *, int *, int *);

extern void eigsrt(double *, double *, int);
extern double util_random(const int);

extern double util_double_factorial(const int);

extern void util_compcharge_gen_rgaussian(const int, const double, const int, const double *, double *);

extern double util_log_integrate_def(const int, const double *,
                                     const int, const double *,
                                     const double, const int);
extern double util_log_integrate_def0(const int,
                                      const double *,
                                      const double *,
                                      const double,
                                      const int);
extern void util_log_integrate_indef(const int, const double *,
                                     const int, const double *,
                                     const double,
                                     const int, double *);
extern double util_log_multipole_energy(const int, const int, const double *,
                                        const int, const double *,
                                        const int, const double *,
                                        const double);
 
extern double util_log_r2integrate_eric(const int, const double, const double *, const double *);
extern double util_log_corrector_iF(const int, const double *);
extern double util_log_coulomb0_energy(const double *, const double,
                                       const double *, const int,
                                       const double, const double);
extern double util_log_coulomb_energy(const double *, const double,
                                      const double *, const int,
                                      const double);

extern void util_getfilling(int, int *, int *, double *);

extern bool util_filefind(Parallel *, char *);

extern void   util_spline(double *, double *, int, double, double, double *, double *);
extern double util_splint(double *, double *, double *, int, int, double);

extern void util_filter(int, double *, double, double *);

extern double util_ln_gamma(const double);
extern double util_gamma(const double);
extern double util_gser(const double, const double);
extern double util_gcf(const double, const double);
extern double util_gammp(const double, const double);

extern void util_gauss_weights(const double, const double,
                               double *, double *, const int);

extern double util_legendre_lm(const int, const int, const double);
extern double util_rlegendre_lm(const int, const int, const double);
extern double util_theta_lm(const int, const int, const double);


extern double util_linesearch(double, double (*)(double), double (*)(double), double, double, double, int);

//#define _MATPRINT_
#ifdef _MATPRINT_
#include <iostream>
extern void util_matprint(std::string, int, double *);
#endif

}

#endif

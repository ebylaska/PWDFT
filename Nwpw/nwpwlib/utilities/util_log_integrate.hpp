#ifndef _UTIL_LOG_GRID_HPP_
#define _UTIL_LOG_GRID_HPP_


namespace pwdft {
using namespace pwdft;

extern double util_double_factorial(const int);
extern void util_gauss_weights(const double, const double, double *, double *, const int);

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

extern double util_SpecialKummer(const int, const int, const double);
extern double util_GaussBessel(const int, const int, const double, const double);
extern double util_dGaussBessel(const int, const int, const double, const double);

}

#endif


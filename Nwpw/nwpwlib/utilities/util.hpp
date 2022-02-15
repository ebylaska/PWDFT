#ifndef _UTIL_H_
#define _UTIL_H_

#include        "Parallel.hpp"
#include	"util_legendre.hpp"
#include	"util_gamma.hpp"
#include	"util_log_integrate.hpp"
#include	"util_tesseral.hpp"
//#include	"util_wgaussian.hpp"
#include	"util_linesearch.hpp"

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

extern void util_getfilling(int, int *, int *, double *);

extern bool util_filefind(Parallel *, char *);

extern void   util_spline(const double *, const double *, const int, const double, const double, double *, double *);
extern double util_splint(const double *, const double *, const double *, const int, const int, const double);

extern void util_filter(int, double *, double, double *);


//#define _MATPRINT_
#ifdef _MATPRINT_
#include <iostream>
extern void util_matprint(std::string, int, double *);
#endif

}

#endif

#ifndef _UTIL_H_
#define _UTIL_H_

#include        "Parallel.hpp"

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

extern void   util_spline(double *, double *, int, double, double, double *, double *);
extern double util_splint(double *, double *, double *, int, int, double);

extern void util_filter(int, double *, double, double *);

#endif

#ifndef _UTIL_H_
#define _UTIL_H_

#include        "Parallel.hpp"

extern void c_aindexcopy(const int, const int *, float *, float *);
extern void c_bindexcopy(const int, const int *, float *, float *);
extern void c_bindexcopy_conjg(const int, const int *, float *, float *);

extern void t_aindexcopy(const int, const int *, float *, float *);
extern void t_bindexcopy(const int, const int *, float *, float *);

extern void i_aindexcopy(const int, const int *, int *, int *);

extern void eigsrt(float *, float *, int);
extern float util_random(const int);

extern void util_getfilling(int, int *, int *, float *);

extern bool util_filefind(Parallel *, char *);

extern void   util_spline(float *, float *, int, float, float, float *, float *);
extern float util_splint(float *, float *, float *, int, int, float);

extern void util_filter(int, float *, float, float *);

#endif

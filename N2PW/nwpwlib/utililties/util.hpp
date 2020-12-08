#ifndef _UTIL_H_
#define _UTIL_H_

extern void c_aindexcopy(const int, const int *, float *, float *);
extern void c_bindexcopy(const int, const int *, float *, float *);
extern void c_bindexcopy_conjg(const int, const int *, float *, float *);

extern void t_aindexcopy(const int, const int *, float *, float *);
extern void t_bindexcopy(const int, const int *, float *, float *);

extern void i_aindexcopy(const int, const int *, int *, int *);

extern void eigsrt(float *, float *, int);

#endif

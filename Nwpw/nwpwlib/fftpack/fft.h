#ifndef _FFT_H_
#define _FFT_H_

extern "C" void drffti_(const int *, double *);
extern "C" void dcffti_(const int *, double *);

extern "C" void drfftb_(const int *, double *, double*);
extern "C" void dcfftb_(const int *, double *, double*);

extern "C" void drfftf_(const int *, double *, double*);
extern "C" void dcfftf_(const int *, double *, double*);

#endif


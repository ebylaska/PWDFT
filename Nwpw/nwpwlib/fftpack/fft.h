#ifndef _FFT_H_
#define _FFT_H_

extern "C" void drffti_(const int *, float *);
extern "C" void dcffti_(const int *, float *);

extern "C" void drfftb_(const int *, float *, float*);
extern "C" void dcfftb_(const int *, float *, float*);

extern "C" void drfftf_(const int *, float *, float*);
extern "C" void dcfftf_(const int *, float *, float*);

#endif


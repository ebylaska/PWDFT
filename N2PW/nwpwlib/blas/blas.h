#ifndef _BLAS_H_
#define _BLAS_H_

extern "C" void dcopy_(int *, float *, int *, float *, int *);
extern "C" float ddot_(int *, float *, int *, float *, int *);
extern "C" void daxpy_(int *, float *, float *, int *, float *, int *);
extern "C" void dscal_(int *, float *, float *, int *);
extern "C" void dgemm_(char *, char *, int *, int *, int *,
                       float *, 
                       float *, int *,
                       float *, int *,
                       float *,
                       float *, int *);

//extern "C" void eigen_(int *, int *, float *, float *, float *, int *);

extern "C" void dsyev_(char *, char *, int *,
                       float *, int *,
                       float *,
                       float *, int *, int*);
extern "C" int  idamax_(int *, float *, int *);


#endif


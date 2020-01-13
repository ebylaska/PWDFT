#ifndef _BLAS_H_
#define _BLAS_H_

extern "C" void dcopy_(int *, double *, int *, double *, int *);
extern "C" double ddot_(int *, double *, int *, double *, int *);
extern "C" void daxpy_(int *, double *, double *, int *, double *, int *);
extern "C" void dscal_(int *, double *, double *, int *);
extern "C" void dgemm_(char *, char *, int *, int *, int *,
                       double *, 
                       double *, int *,
                       double *, int *,
                       double *,
                       double *, int *);

//extern "C" void eigen_(int *, int *, double *, double *, double *, int *);

extern "C" void dsyev_(char *, char *, int *,
                       double *, int *,
                       double *,
                       double *, int *, int*);
extern "C" int  idamax_(int *, double *, int *);


#endif


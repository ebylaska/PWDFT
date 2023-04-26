#ifndef _BLAS_H_
#define _BLAS_H_

#if defined(NWPW_INTEL_MKL)
#include "mkl.h"

#define DSCAL_PWDFT(n, alpha, a, ida) cblas_dscal(n, alpha, a, ida);
#define DCOPY_PWDFT(n, a, ida, b, idb) cblas_dcopy(n, a, ida, b, idb)
#define DAXPY_PWDFT(n, alpha, a, ida, b, idb)                                  \
  cblas_daxpy(n, alpha, a, ida, b, idb)

//#define TRANSCONV(a)    ( (std::strcmp(a, "N")) ?  CblasNoTrans : CblasTrans )
#define TRANSCONV(a) ((a == "N") ? CblasNoTrans : CblasTrans)
#define DGEMM_PWDFT(s1, s2, n, m, k, alpha, a, ida, b, idb, beta, c, idc)      \
  cblas_dgemm(CblasColMajor, TRANSCONV(s1), TRANSCONV(s2), n, m, k, alpha, a,  \
              ida, b, idb, beta, c, idc)

#define IDAMAX_PWDFT(nn, hml, one) cblas_idamax(nn, hml, one)

#define EIGEN_PWDFT(n, hml, eig, xtmp, nn, ierr)                               \
  ierr = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, hml, n, eig)

#define DDOT_PWDFT(n, a, ida, b, idb) cblas_ddot(n, a, ida, b, idb);

#define DGELSS_PWDFT(m, n, nrhs, a, ida, b, idb, s1, rcond, rank, work, iwork, \
                     ierr)                                                     \
  ierr = LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', m, n, nrhs, a, ida, b, idb);

#define ZGEMM_PWDFT(s1, s2, n, m, k, alpha, a, ida, b, idb, beta, c, idc)      \
  cblas_zgemm(CblasColMajor, TRANSCONV(s1), TRANSCONV(s2), n, m, k, alpha, a,  \
              ida, b, idb, beta, c, idc)

#else

extern "C" void dcopy_(int *, double *, int *, double *, int *);
extern "C" double ddot_(int *, double *, int *, double *, int *);
extern "C" void daxpy_(int *, double *, double *, int *, double *, int *);
extern "C" void dscal_(int *, double *, double *, int *);
extern "C" void dgemm_(char *, char *, int *, int *, int *, double *, double *,
                       int *, double *, int *, double *, double *, int *);

// extern "C" void eigen_(int *, int *, double *, double *, double *, int *);

extern "C" void dsyev_(char *, char *, int *, double *, int *, double *,
                       double *, int *, int *);
extern "C" int idamax_(int *, double *, int *);

extern "C" void dlacpy_(char *, int *, int *, double *, int *, double *, int *);

extern "C" void dgelss_(int *, int *, int *, double *, int *, double *, int *,
                        double *, double *, int *, double *, int *, int *);

extern "C" void zgemm_(char *, char *, int *, int *, int *, double *, double *,
                       int *, double *, int *, double *, double *, int *);

#define DSCAL_PWDFT(n, alpha, a, ida) dscal_(&(n), &(alpha), a, &(ida))
#define DCOPY_PWDFT(n, a, ida, b, idb) dcopy_(&(n), a, &(ida), b, &(idb))
#define DAXPY_PWDFT(n, alpha, a, ida, b, idb)                                  \
  daxpy_(&(n), &(alpha), a, &(ida), b, &(idb))
#define DGEMM_PWDFT(s1, s2, n, m, k, alpha, a, ida, b, idb, beta, c, idc)      \
  dgemm_(s1, s2, &(n), &(m), &(k), &(alpha), a, &(ida), b, &(idb), &(beta), c, \
         &(idc))

#define IDAMAX_PWDFT(nn, hml, one) idamax_(&(nn), hml, &(one))

#define EIGEN_PWDFT(n, hml, eig, xtmp, nn, ierr)                               \
  dsyev_((char *)"V", (char *)"U", &(n), hml, &(n), eig, xtmp, &(nn), &ierr)

#define DDOT_PWDFT(n, a, ida, b, idb) ddot_(&(n), a, &ida, b, &(idb))
#define DLACPY_PWDFT(s1, m, n, a, ida, b, idb)                                 \
  dlacpy_(s1, &(m), &(n), a, &(ida), b, &(idb))

#define DGELSS_PWDFT(m, n, nrhs, a, ida, b, idb, s1, rcond, rank, work, iwork, \
                     ierr)                                                     \
  dgelss_(&(m), &(n), &(nrhs), a, &(ida), b, &(idb), s1, &(rcond), &(rank),    \
          work, &(iwork), &(ierr))

#define ZGEMM_PWDFT(s1, s2, n, m, k, alpha, a, ida, b, idb, beta, c, idc)      \
  zgemm_(s1, s2, &(n), &(m), &(k), &(alpha), a, &(ida), b, &(idb), &(beta), c, \
         &(idc))

#endif

extern "C" void factor_skew_(int *, double *, double *, double *, double *);

#define FACTOR_SKEW_PWDFT(n, k, v, w, s) factor_skew_(&(n), k, v, w, s)

#endif

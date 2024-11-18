#ifndef _BLAS_H_
#define _BLAS_H_

#if defined(NWPW_INTEL_MKL)
#include "mkl.h"
#include "mkl_lapacke.h"
#include <cstring> // Required for strcmp

#define DSCAL_PWDFT(n, alpha, a, ida) cblas_dscal(n, alpha, a, ida);
#define DCOPY_PWDFT(n, a, ida, b, idb) cblas_dcopy(n, a, ida, b, idb)
#define DAXPY_PWDFT(n, alpha, a, ida, b, idb)                                  \
  cblas_daxpy(n, alpha, a, ida, b, idb)

//#define TRANSCONV(a)    ( (std::strcmp(a, "N")) ?  CblasNoTrans : CblasTrans )
//#define TRANSCONV(a) ((a == "N") ? CblasNoTrans : CblasTrans)
//#define CTRANSCONV(a) ((a == "N") ? CblasNoTrans : CblasConjTrans)

#define TRANSCONV(a) ((strcmp(a, "N") == 0) ? CblasNoTrans : CblasTrans)
#define CTRANSCONV(a) ((strcmp(a, "N") == 0) ? CblasNoTrans : CblasConjTrans)


#define DGEMM_PWDFT(s1, s2, n, m, k, alpha, a, ida, b, idb, beta, c, idc)      \
  cblas_dgemm(CblasColMajor, TRANSCONV(s1), TRANSCONV(s2), n, m, k, alpha, a,  \
              ida, b, idb, beta, c, idc)

#define IDAMAX_PWDFT(nn, hml, one) cblas_idamax(nn, hml, one)

#define EIGEN_PWDFT(n, hml, eig, xtmp, nn, ierr)                               \
  ierr = LAPACKE_dsyev(LAPACK_COL_MAJOR, 'V', 'U', n, hml, n, eig)


#define DDOT_PWDFT(n, a, ida, b, idb) cblas_ddot(n, a, ida, b, idb);

#define DLACPY_PWDFT(s1, m, n, a, ida, b, idb)                                 \
  auto ierr0 = LAPACKE_dlacpy(LAPACK_COL_MAJOR, (s1)[0], m, n, a, ida, b, idb)

#define DGELSS_PWDFT(m, n, nrhs, a, ida, b, idb, s1, rcond, rank, work, iwork, \
                     ierr)                                                     \
  ierr = LAPACKE_dgels(LAPACK_COL_MAJOR, 'N', m, n, nrhs, a, ida, b, idb);


#define ZSCAL_PWDFT(n, alpha, a, ida) cblas_zscal(n, alpha, a, ida);

#define ZDOTC_PWDFT(n, a, ida, b, idb) ([](int _n, const double* _a, int _ida, const double* _b, int _idb) -> std::complex<double> { \
    std::vector<std::complex<double>> temp_a(_n), temp_b(_n); \
    for (int i = 0; i < _n; ++i) { \
        temp_a[i] = std::complex<double>(_a[2*i], _a[2*i + 1]); \
        temp_b[i] = std::complex<double>(_b[2*i], _b[2*i + 1]); \
    } \
    MKL_Complex16 result; \
    cblas_zdotc_sub(_n, reinterpret_cast<const MKL_Complex16*>(temp_a.data()), _ida, reinterpret_cast<const MKL_Complex16*>(temp_b.data()), _idb, &result); \
    return std::complex<double>(result.real, result.imag); \
}((n), (a), (ida), (b), (idb)))



#define ZAXPY_PWDFT(n, alpha, a, ida, b, idb)                                  \
  cblas_zaxpy(n, alpha, a, ida, b, idb)

#define ZGEMM_PWDFT(s1, s2, n, m, k, alpha, a, ida, b, idb, beta, c, idc)      \
  cblas_zgemm(CblasColMajor, CTRANSCONV(s1), CTRANSCONV(s2), n, m, k, alpha, a,  \
              ida, b, idb, beta, c, idc)

#define IZAMAX_PWDFT(nn, hml, one) cblas_izamax(nn, hml, one)

#define ZEIGEN_PWDFT(n, hml, eig, xtmp, nn, rtmp,ierr)                         \
  ierr = LAPACKE_zheev(LAPACK_COL_MAJOR, 'V', 'L', n, reinterpret_cast<MKL_Complex16*> (hml), n, eig)

#define ZLACPY_PWDFT(s1, m, n, a, ida, b, idb)                                 \
  auto ierr0 = LAPACKE_dlacpy(LAPACK_COL_MAJOR, (s1)[0], m, n, a, ida, b, idb)

#define DGESV_PWDFT(n, nrhs, a, lda, ipiv, b, ldb, ierr)                       \
  ierr = LAPACKE_dgesv(LAPACK_COL_MAJOR, n, nrhs, a, lda, ipiv, b, ldb)

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
extern "C" void dgesv_(int *, int *,  double *, int *, int *, double *, int *, int *);


extern "C" void zscal_(int *, double *, double *, int *);
extern "C" double zdotc_(int *, double *, int *, double *, int *);
extern "C" void zaxpy_(int *, double *, double *, int *, double *, int *);
extern "C" void zgemm_(char *, char *, int *, int *, int *, double *, double *,
                       int *, double *, int *, double *, double *, int *);

extern "C" int izamax_(int *, double *, int *);

extern "C" void zheev_(char *, char *, int *, double *, int *, double *,
                       double *, int *, double *, int *);

extern "C" void zlacpy_(char *, int *, int *, double *, int *, double *, int *);

#define DSCAL_PWDFT(n, alpha, a, ida) dscal_(&(n), &(alpha), a, &(ida))
#define DCOPY_PWDFT(n, a, ida, b, idb) dcopy_(&(n), a, &(ida), b, &(idb))
#define DAXPY_PWDFT(n, alpha, a, ida, b, idb)                                  \
  daxpy_(&(n), &(alpha), a, &(ida), b, &(idb))
#define DGEMM_PWDFT(s1, s2, n, m, k, alpha, a, ida, b, idb, beta, c, idc)      \
  dgemm_(s1, s2, &(n), &(m), &(k), &(alpha), (a), &(ida), (b), &(idb), &(beta), (c), \
         &(idc))

#define IDAMAX_PWDFT(nn, hml, one) idamax_(&(nn), hml, &(one))


#define EIGEN_PWDFT(n, hml, eig, xtmp, nn, ierr)                               \
  dsyev_((char *)"V", (char *)"U", &(n), hml, &(n), eig, xtmp, &(nn), &ierr)

#define DDOT_PWDFT(n, a, ida, b, idb) ddot_(&(n), (a), &ida, (b), &(idb))

#define DLACPY_PWDFT(s1, m, n, a, ida, b, idb)                                 \
  dlacpy_(s1, &(m), &(n), a, &(ida), b, &(idb))

#define DGELSS_PWDFT(m, n, nrhs, a, ida, b, idb, s1, rcond, rank, work, iwork, \
                     ierr)                                                     \
  dgelss_(&(m), &(n), &(nrhs), a, &(ida), b, &(idb), s1, &(rcond), &(rank),    \
          work, &(iwork), &(ierr))



#define ZSCAL_PWDFT(n, alpha, a, ida) zscal_(&(n), alpha, a, &(ida))
#define ZDOTC_PWDFT(n, a, ida, b, idb) zdotc_(&(n), a, &ida, b, &(idb))

#define ZAXPY_PWDFT(n, alpha, a, ida, b, idb)                                  \
  zaxpy_(&(n), alpha, a, &(ida), b, &(idb))

#define ZGEMM_PWDFT(s1, s2, n, m, k, alpha, a, ida, b, idb, beta, c, idc)      \
  zgemm_(s1, s2, &(n), &(m), &(k), alpha, a, &(ida), b, &(idb), beta, c, \
         &(idc))

#define IZAMAX_PWDFT(nn, hml, one) izamax_(&(nn), hml, &(one))

#define ZEIGEN_PWDFT(n, hml, eig, xtmp, nn, rtmp, ierr)                               \
  zheev_((char *)"V", (char *)"L", &(n), hml, &(n), eig, xtmp, &(nn), rtmp, &ierr)

#define ZLACPY_PWDFT(s1, m, n, a, ida, b, idb)                                 \
  zlacpy_(s1, &(m), &(n), a, &(ida), b, &(idb))

#define DGESV_PWDFT(n, nrhs, a, lda, ipiv, b, ldb, ierr) \
  dgesv_(&(n), &(nrhs), a, &(lda), ipiv, b, &(ldb), &(ierr))

#endif


extern "C" void factor_skew_(int *, double *, double *, double *, double *);


#define FACTOR_SKEW_PWDFT(n, k, v, w, s) factor_skew_(&(n), k, v, w, s)

#endif

#include "gdevice2.hpp"


// using namespace pwdft;
namespace pwdft {

static Gdevices mygdevice0;


/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
gdevice2::gdevice2() {mygdevice2 = &mygdevice0;}


/********************************
 *                              *
 *         Access functions     *
 *                              *
 ********************************/

void gdevice2::TN4_dgemm(int npack, int ne, double alpha, double *a, double *b,
                       double beta, double *caa, double *cab, double *cba,
                       double *cbb) {
   mygdevice2->TN4_dgemm(npack, ne, alpha, a, b, beta, caa, cab, cba, cbb);
}

void gdevice2::TN3_dgemm(int npack, int ne, double alpha, double *a, double *b,
                       double beta, double *caa, double *cab, double *cbb) {
   mygdevice2->TN3_dgemm(npack, ne, alpha, a, b, beta, caa, cab, cbb);
}

void gdevice2::TN1_dgemm(int npack, int ne, double alpha, double *a, double *b,
                       double beta, double *c) {
   mygdevice2->TN1_dgemm(npack, ne, alpha, a, b, beta, c);
}

void gdevice2::TN_dgemm(int npack, int ne, int nprj, double alpha, double *a,
                      double *b, double beta, double *c) {
   mygdevice2->TN_dgemm(npack, ne, nprj, alpha, a, b, beta, c);
}

void gdevice2::T_free() { mygdevice2->T_free(); }

void gdevice2::NN_dgemm(int npack, int ne, double alpha, double *a, double *b,
                      double beta, double *c) {
   mygdevice2->NN_dgemm(npack, ne, alpha, a, b, beta, c);
}

void gdevice2::NT_dgemm(int npack, int ne, int nprj, double alpha, double *a,
                      double *b, double beta, double *c) {
   mygdevice2->NT_dgemm(npack, ne, nprj, alpha, a, b, beta, c);
}

void gdevice2::MM6_dgemm(int ne, double *s21, double *s12, double *s11,
                       double *sa0, double *sa1, double *st1) {
   mygdevice2->MM6_dgemm(ne, s21, s12, s11, sa0, sa1, st1);
}

void gdevice2::NN_eigensolver(int ispin, int ne[], double *a, double *w) {
   mygdevice2->NN_eigensolver(ispin, ne, a, w);
}

void gdevice2::psi_alloc(int npack, int ne, int tfac0 = 1) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->psi_alloc(npack, ne, tfac0);
#endif
}

void gdevice2::psi_dealloc() {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->psi_dealloc();
#endif
}

void gdevice2::psi_copy_host2gpu(int npack, int ne, double *psi) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->psi_copy_host2gpu(npack, ne, psi);
#endif
}

void gdevice2::psi_copy_gpu2host(int npack, int ne, double *psi) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->psi_copy_gpu2host(npack, ne, psi);
#endif
}

void gdevice2::hpsi_copy_host2gpu(int npack, int ne, double *psi) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->hpsi_copy_host2gpu(npack, ne, psi);
#endif
}

void gdevice2::hpsi_copy_gpu2host(int npack, int ne, double *psi) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->hpsi_copy_gpu2host(npack, ne, psi);
#endif
}

/* fft functions*/
void gdevice2::batch_fft_init(int nx, int ny, int nz, int nq1, int nq2, int nq3) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_fft_init(nx, ny, nz, nq1, nq2, nq3);
#endif
}
void gdevice2::batch_fft_end() {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_fft_end();
#endif
}

void gdevice2::batch_cfftx(bool forward, int nx, int nq, int n2ft3d, double *a) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cfftx(forward, nx, nq, n2ft3d, a);
#endif
}
void gdevice2::batch_cfftx_tmpx(bool forward, int nx, int nq, int n2ft3d,
                              double *a, double *tmpx) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cfftx(forward, nx, nq, n2ft3d, a);
#else
   mygdevice2->batch_cfftx_tmpx(forward, nx, nq, n2ft3d, a, tmpx);
#endif
}

void gdevice2::batch_cffty(bool forward, int ny, int nq, int n2ft3d, double *a) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cffty(forward, ny, nq, n2ft3d, a);
#endif
}

void gdevice2::batch_cffty_tmpy(bool forward, int ny, int nq, int n2ft3d,
                              double *a, double *tmpy) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cffty(forward, ny, nq, n2ft3d, a);
#else
   mygdevice2->batch_cffty_tmpy(forward, ny, nq, n2ft3d, a, tmpy);
#endif
}


void gdevice2::batch_cffty_tmpy_zero(bool forward, int ny, int nq, int n2ft3d,
                                   double *a, double *tmpy, bool *zero) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cffty(forward, ny, nq, n2ft3d, a);
#else
   mygdevice2->batch_cffty_tmpy_zero(forward, ny, nq, n2ft3d, a, tmpy, zero);
#endif
}

void gdevice2::batch_cfftz(bool forward, int nz, int nq, int n2ft3d, double *a) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cfftz(forward, nz, nq, n2ft3d, a);
#endif
}

void gdevice2::batch_cfftz_tmpz(bool forward, int nz, int nq, int n2ft3d,
                              double *a, double *tmpz) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cfftz(forward, nz, nq, n2ft3d, a);
#else
   mygdevice2->batch_cfftz_tmpz(forward, nz, nq, n2ft3d, a, tmpz);
#endif
}
void gdevice2::batch_cfftz_tmpz_zero(bool forward, int nz, int nq, int n2ft3d,
                                   double *a, double *tmpz, bool *zero) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cfftz(forward, nz, nq, n2ft3d, a);
#else
   mygdevice2->batch_cfftz_tmpz_zero(forward, nz, nq, n2ft3d, a, tmpz, zero);
#endif
}

} // namespace pwdft

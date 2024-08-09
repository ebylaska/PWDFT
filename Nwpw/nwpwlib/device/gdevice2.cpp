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
void gdevice2::TN_dgemm2c(int n, int m, int npack2, int nida2,
                          double *a, double *b, double *c) {
   mygdevice2->TN_dgemm2c(n, m, npack2, nida2, a, b, c);
}

void gdevice2::T_free() { mygdevice2->T_free(); }

void gdevice2::NN_dgemm(int npack, int ne, double alpha, double *a, double *b,
                      double beta, double *c) {
   mygdevice2->NN_dgemm(npack, ne, alpha, a, b, beta, c);
}

void gdevice2::NN_dgemm1(int n, int m, int k, double alpha, double *a, int lda, double *b, int ldb, double beta, double *c, int ldc) {
   mygdevice2->NN_dgemm1(n,m,k,alpha,a,lda,b,ldb,beta,c,ldc);
}
void gdevice2::TN_dgemm2(int n, int m, int k, double alpha, double *a, int lda, double *b, int ldb, double beta, double *c, int ldc) {
   mygdevice2->TN_dgemm2(n,m,k,alpha,a,lda,b,ldb,beta,c,ldc);
}
void gdevice2::NT_dgemm3(int n, int m, int k, double alpha, double *a, int lda, double *b, int ldb, double beta, double *c, int ldc) {
   mygdevice2->NT_dgemm3(n,m,k,alpha,a,lda,b,ldb,beta,c,ldc);
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

void gdevice2::NN1_zgemm(int npack1_max, int npack, int ne, double *alpha, double *a, double *b,
                         double *beta, double *c) {
   mygdevice2->NN1_zgemm(npack1_max, npack, ne, alpha, a, b, beta, c);
}

void gdevice2::CN1_zgemm(int npack1_max, int npack, int ne, double *alpha, double *a, double *b,
                         double *beta, double *c) {
   mygdevice2->CN1_zgemm(npack1_max, npack, ne, alpha, a, b, beta, c);
}

void gdevice2::CN2_zgemm(int ne, int nprj, int npack, int npack1_max, double *alpha, double *a, double *b,
                         double *beta, double *c) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->CN2_stride_zgemm(ne,nprj,npack,npack1_max, alpha, a, b, beta, c);

   else
#endif
      mygdevice2->CN2_zgemm(ne, nprj, npack, npack1_max, alpha, a, b, beta, c);
}

void gdevice2::NC2_zgemm(int npack1_max, int npack, int ne, int nprj, double *alpha, double *a, double *b,
                         double *beta, double *c) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->NC2_stride_zgemm(npack1_max, npack, ne, nprj, alpha, a, b, beta, c);
   else
#endif
      mygdevice2->NC2_zgemm(npack1_max, npack, ne, nprj, alpha, a, b, beta, c);
}


void gdevice2::NN_zgemm(int n, int m, int k, double *alpha, double *a, int lda, double *b, int ldb, double *beta, double *c, int ldc) {
   mygdevice2->NN_zgemm(n,m,k,alpha,a,lda,b,ldb,beta,c,ldc);
}



void gdevice2::NC_zgemm(int n, int m, int k, double *alpha, double *a, int lda, double *b, int ldb, double *beta, double *c, int ldc) {
   mygdevice2->NC_zgemm(n,m,k,alpha,a,lda,b,ldb,beta,c,ldc);
}

void gdevice2::CN_zgemm(int n, int m, int k, double *alpha, double *a, int lda, double *b, int ldb, double *beta, double *c, int ldc) {
   mygdevice2->CN_zgemm(n,m,k,alpha,a,lda,b,ldb,beta,c,ldc);
}

void gdevice2::CN4_zgemm(int npack1_max, int npack, int ne, double *alpha, double *a, double *b,
                         double *beta, double *caa, double *cab, double *cba,
                         double *cbb) {
   mygdevice2->CN4_zgemm(npack1_max, npack, ne, alpha, a, b, beta, caa, cab, cba, cbb);
}

void gdevice2::CN3_zgemm(int npack1_max, int npack, int ne, double *alpha, double *a, double *b,
                         double *beta, double *caa, double *cab, double *cbb) {
   mygdevice2->CN3_zgemm(npack1_max, npack, ne, alpha, a, b, beta, caa, cab, cbb);
}

/*
void gdevice2::computeTrans3_Mult(const int ne, const int nprj, 
                                  const double *psi, const double *prj, int ng, int ng0,
                                  double *Gx, double *Gy, double *Gz, double *xtmp1,
                                  double *sum3) {
   if (!mygdevice2->hasgpu)
      mygdevice2->computeTrans3_Mult(ne,nprj,psi,prj,ng,ng0,Gx,Gy,Gz,xtmp1,sum3);
}
*/

/*
void gdevice2::computeTrans_Mult(int ne, int nprj, double alpha, double alpha1, int ng, int ng0, 
                                 double *psi, double *prj, double beta, double beta1, double *sum1) 
{
   if (!mygdevice2->hasgpu)
      mygdevice2->computeTrans_Mult(ne,nprj,alpha,alpha1,ng,ng0,psi,prj,beta,beta1,sum1);
}
*/



void gdevice2::WW6_zgemm(int ne, double *s21, double *s12, double *s11,
                       double *sa0, double *sa1, double *st1) {
   mygdevice2->WW6_zgemm(ne, s21, s12, s11, sa0, sa1, st1);
}



void gdevice2::WW_eigensolver(int ispin, int ne[], double *a, double *w) {
   mygdevice2->WW_eigensolver(ispin, ne, a, w);
}





void gdevice2::psi_alloc(int npack, int ne, int tfac0 = 1) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->psi_alloc(npack, ne, tfac0);
#endif
}

void gdevice2::psi_dealloc() {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->psi_dealloc();
#endif
}

void gdevice2::psi_copy_host2gpu(int npack, int ne, double *psi) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->psi_copy_host2gpu(npack, ne, psi);
#endif
}

void gdevice2::psi_copy_gpu2host(int npack, int ne, double *psi) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->psi_copy_gpu2host(npack, ne, psi);
#endif
}

void gdevice2::hpsi_copy_host2gpu(int npack, int ne, double *psi) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->hpsi_copy_host2gpu(npack, ne, psi);
#endif
}

void gdevice2::hpsi_copy_gpu2host(int npack, int ne, double *psi) {
#if defined(NWPW_SYCL) || defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->hpsi_copy_gpu2host(npack, ne, psi);
#endif
}

/* fft functions - SYCL ffts removed */
int  gdevice2::batch_fft_init(int nx, int ny, int nz, int nq1, int nq2, int nq3) {
   int tag = -1;
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      tag = mygdevice2->batch_fft_init(nx, ny, nz, nq1, nq2, nq3);
#endif
   return tag;
}
void gdevice2::batch_fft_end(const int tag) {
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_fft_end(tag);
#endif
}



void gdevice2::batch_fft_pipeline_mem_init(const int nstages, const int n2ft3d) {
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_fft_pipeline_mem_init(nstages,n2ft3d);
#endif
}


                                   
void gdevice2::batch_rfftx_tmpx(const int tag,bool forward, int nx, int nq, int n2ft3d,
                                double *a, double *tmpx) {
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_rfftx(tag,forward, nx, nq, n2ft3d, a);
#else
   mygdevice2->batch_rfftx_tmpx(forward, nx, nq, n2ft3d, a, tmpx);
#endif
}
 
void gdevice2::batch_rfftx_stages_tmpx(const int stage, const int tag,bool forward, int nx, int nq, int n2ft3d,
                                       double *a, double *tmpx, int da) { 
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_rfftx_stages(stage,tag,forward, nx, nq, n2ft3d, a,da);
#endif
}



void gdevice2::batch_cfftx_tmpx(const int tag,bool forward, int nx, int nq, int n2ft3d,
                                double *a, double *tmpx) {
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cfftx(tag,forward, nx, nq, n2ft3d, a);
#else
   mygdevice2->batch_cfftx_tmpx(forward, nx, nq, n2ft3d, a, tmpx);
#endif
}

void gdevice2::batch_cfftx_stages_tmpx(const int stage, const int tag,bool forward, int nx, int nq, int n2ft3d,
                                       double *a, double *tmpx, int da) {
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cfftx_stages(stage,tag,forward, nx, nq, n2ft3d, a,da);
#endif
}



void gdevice2::batch_cffty_tmpy(const int tag, bool forward, int ny, int nq, int n2ft3d, double *a, double *tmpy) 
{
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cffty(tag,forward, ny, nq, n2ft3d, a);
#else
   mygdevice2->batch_cffty_tmpy(forward, ny, nq, n2ft3d, a, tmpy);
#endif
}

void gdevice2::batch_cffty_tmpy_zero(const int tag, bool forward, int ny, int nq, int nffts, int n2ft3d,
                                     double *a, double *tmpy, bool *zero) {
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cffty(tag,forward,ny,nq,n2ft3d,a);
#else
   mygdevice2->batch_cffty_tmpy_zero(forward, ny, nq, nffts, n2ft3d, a, tmpy, zero);
   //mygdevice2->batch_cffty_tmpy_zero(forward, ny, nq, n2ft3d, a, tmpy, zero);
#endif
}

void gdevice2::batch_cffty_stages_tmpy(const int stage,const int tag,bool forward, int ny, int nq, int n2ft3d,
                                       double *a, double *tmpy, int da) {
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cffty_stages(stage,tag,forward,ny,nq,1,n2ft3d,a,da);
#endif
}

void gdevice2::batch_cffty_stages_tmpy_zero(const int stage, const int tag, bool forward, int ny, int nq, int nffts, int n2ft3d,
                                            double *a, double *tmpy, bool *zero, int da) {
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cffty_stages(stage,tag,forward,ny,nq,nffts,n2ft3d,a,da);
#endif
}





void gdevice2::batch_cfftz_tmpz(const int tag, bool forward, int nz, int nq, int n2ft3d,
                              double *a, double *tmpz) {
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cfftz(tag,forward,nz,nq,n2ft3d,a);
#else
   mygdevice2->batch_cfftz_tmpz(forward,nz,nq,n2ft3d,a,tmpz);
#endif
}

void gdevice2::batch_cfftz_tmpz_zero(const int tag, bool forward, int nz, int nq, int nffts, int n2ft3d,
                                     double *a, double *tmpz, bool *zero) {
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cfftz(tag,forward,nz,nq,n2ft3d,a);
#else
   mygdevice2->batch_cfftz_tmpz_zero(forward,nz,nq,nffts,n2ft3d,a,tmpz,zero);
   //mygdevice2->batch_cfftz_tmpz_zero(forward,nz,nq,n2ft3d,a,tmpz,zero);
#endif
}

void gdevice2::batch_cfftz_stages_tmpz(const int stage, const int tag, bool forward, int nz, int nq, int n2ft3d,
                                       double *a, double *tmpz, const int da) {
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cfftz_stages(stage, tag, forward, nz, nq, 1, n2ft3d, a,da);
#endif
}

void gdevice2::batch_cfftz_stages_tmpz_zero(const int stage, const int tag, bool forward, int nz, int nq, int nffts, int n2ft3d,
                                            double *a, double *tmpz, bool *zero, int da) {
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
      mygdevice2->batch_cfftz_stages(stage, tag, forward, nz, nq, nffts, n2ft3d, a,da);
#endif
}





int gdevice2::size_fft_twiddle(const int n) {
   return mygdevice2->size_fft_twiddle(n);
}

void gdevice2::set_fft_twiddle(const int isgn, const int n, double *twiddle) {
   mygdevice2->set_fft_twiddle(isgn,n,twiddle);
}





void gdevice2::batch_cfft(const int tag, const bool forward, const int nz, const int nq, const int nffts, const int n2ft3d, double *a, 
                          const double *twiddle, const double *tmpz, const int xyz_gpu) {
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
   {
      if (xyz_gpu==0) mygdevice2->batch_cfftx(tag, forward, nz, nq, n2ft3d, a);
      if (xyz_gpu==1) mygdevice2->batch_cffty(tag, forward, nz, nq, n2ft3d, a);
      if (xyz_gpu==2) mygdevice2->batch_cfftz(tag, forward, nz, nq, n2ft3d, a);
   }
#else
   mygdevice2->batch_cfft(forward, nz, nq, nffts, n2ft3d, a, twiddle, tmpz);
#endif
}

void gdevice2::batch_cfft_zero(const int tag, const bool forward, const int nz, const int nq, const int nffts, const int n2ft3d, double *a, 
                               const double *twiddle, const double *tmpz, const bool *zero, const int xyz_gpu) {
#if defined(NWPW_CUDA) || defined(NWPW_HIP)
   if (mygdevice2->hasgpu)
   {
      if (xyz_gpu==0) mygdevice2->batch_cfftx(tag,forward,nz,nq,n2ft3d,a);
      if (xyz_gpu==1) mygdevice2->batch_cffty(tag,forward,nz,nq,n2ft3d,a);
      if (xyz_gpu==2) mygdevice2->batch_cfftz(tag,forward,nz,nq,n2ft3d,a);
   }
#else
   mygdevice2->batch_cfft_zero(forward,nz,nq,nffts,n2ft3d,a,twiddle,tmpz,zero);
   //mygdevice2->batch_cffty_tmpy_zero(forward, ny, nq, n2ft3d, a, tmpy, zero);

#endif
}








} // namespace pwdft

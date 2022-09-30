#include "gdevices.hpp"

//using namespace pwdft;
namespace pwdft {

static Gdevices mygdevice;

#ifdef NWPW_SYCL
sycl::queue* get_syclQue() {
  return mygdevice.device_queue;
}
#endif // NWPW_SYCL

void gdevice_TN3_dgemm(int npack, int ne, double alpha, double *a, double *b, double beta, double *caa, double *cab, double *cbb)
{
  mygdevice.TN3_dgemm(npack,ne,alpha,a,b,beta,caa,cab,cbb);
}

void gdevice_TN_dgemm(int npack, int ne, int nprj, double alpha, double *a, double *b, double beta, double *c)
{
  mygdevice.TN_dgemm(npack,ne,nprj,alpha,a,b,beta,c);
}

void gdevice_T_free() { 
   mygdevice.T_free(); 
}

void gdevice_NN_dgemm(int npack, int ne, double alpha, double *a, double *b, double beta, double *c)
{
  mygdevice.NN_dgemm(npack,ne,alpha,a,b,beta,c);
}

void gdevice_NT_dgemm(int npack, int ne, int nprj, double alpha, double *a, double *b, double beta, double *c)
{
  mygdevice.NT_dgemm(npack,ne,nprj,alpha,a,b,beta,c);
}


void gdevice_psi_alloc(int npack, int ne)
{ 
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
  if (mygdevice.hasgpu) mygdevice.psi_alloc(npack,ne);
#endif
}

void gdevice_psi_dealloc()
{
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
  if (mygdevice.hasgpu) mygdevice.psi_dealloc();
#endif
}

void gdevice_psi_copy_host2gpu(int npack , int ne, double *psi)
{ 
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
  if (mygdevice.hasgpu) mygdevice.psi_copy_host2gpu(npack, ne, psi);
#endif
}

void gdevice_psi_copy_gpu2host(int npack, int ne, double *psi)
{ 
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
  if (mygdevice.hasgpu) mygdevice.psi_copy_gpu2host(npack,ne,psi);
#endif
}

void gdevice_hpsi_copy_host2gpu(int npack , int ne, double *psi)
{
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
  if (mygdevice.hasgpu) mygdevice.hpsi_copy_host2gpu(npack, ne, psi);
#endif
}

void gdevice_hpsi_copy_gpu2host(int npack, int ne, double *psi)
{
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
  if (mygdevice.hasgpu) mygdevice.hpsi_copy_gpu2host(npack,ne,psi);
#endif
}

/* fft functions*/
void gdevice_batch_fft_init(int nx,int ny,int nz, int nq1, int nq2, int nq3) 
{
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
  if (mygdevice.hasgpu) mygdevice.batch_fft_init(nx,ny,nz,nq1,nq2,nq3);
#endif
}

void gdevice_batch_cfftx(bool forward,int nx,int nq,int n2ft3d,double *a)
{
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
  if (mygdevice.hasgpu) mygdevice.batch_cfftx(forward,nx,nq,n2ft3d,a);
#endif
}

void gdevice_batch_cffty(bool forward,int ny,int nq,int n2ft3d,double *a)
{
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
  if (mygdevice.hasgpu) mygdevice.batch_cffty(forward,ny,nq,n2ft3d,a);
#endif
}

void gdevice_batch_cfftz(bool forward,int nz,int nq,int n2ft3d,double *a)
{
#if defined(NWPW_SYCL) || defined(NWPW_CUDA)
  if (mygdevice.hasgpu) mygdevice.batch_cfftz(forward,nz,nq,n2ft3d,a);
#endif
}

}

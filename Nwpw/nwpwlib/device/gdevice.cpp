#include "gdevices.hpp"

static Gdevices mygdevice;

#ifdef NWPW_SYCL
#include        <cstdio>
#include        <iostream>
#include        <limits>
#include        <CL/sycl.hpp>
#include        <oneapi/mkl.hpp>


cl::sycl::queue* get_syclQue() {
  return mygdevice.device_queue;
}
double* get_sycl_mem(const size_t mem_in_bytes) {
  return mygdevice.getGpuMem(mem_in_bytes);
}
void free_sycl_mem(double* ptr) {
  mygdevice.freeGpuMem(ptr);
}
double* get_host_mem(const size_t mem_in_bytes) {
  return mygdevice.getHostMem(mem_in_bytes);
}
void free_host_mem(double* ptr) {
  mygdevice.freeHostMem(ptr);
}
#endif

void gdevice_TN3_dgemm(int npack, int ne, double alpha, double *a, double *b, double beta, double *caa, double *cab, double *cbb)
{
  mygdevice.TN3_dgemm(npack,ne,alpha,a,b,beta,caa,cab,cbb);
}

void gdevice_TN_dgemm(int npack, int ne, int nprj, double alpha, double *a, double *b, double beta, double *c)
{
  mygdevice.TN_dgemm(npack,ne,nprj,alpha,a,b,beta,c);
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
#ifdef NWPW_SYCL
  if (mygdevice.hasgpu) mygdevice.psi_alloc(npack,ne);
#endif
}

void gdevice_psi_dealloc()
{
#ifdef NWPW_SYCL
  if (mygdevice.hasgpu) mygdevice.psi_dealloc();
#endif
}

void gdevice_psi_copy_host2gpu(int npack , int ne, double *psi)
{ 
#ifdef NWPW_SYCL
  if (mygdevice.hasgpu) mygdevice.psi_copy_host2gpu(npack, ne, psi);
#endif
}

void gdevice_psi_copy_gpu2host(int npack, int ne, double *psi)
{ 
#ifdef NWPW_SYCL
  if (mygdevice.hasgpu) mygdevice.psi_copy_gpu2host(npack,ne,psi);
#endif
}

void gdevice_hpsi_copy_host2gpu(int npack , int ne, double *psi)
{
#ifdef NWPW_SYCL
  if (mygdevice.hasgpu) mygdevice.hpsi_copy_host2gpu(npack, ne, psi);
#endif
}

void gdevice_hpsi_copy_gpu2host(int npack, int ne, double *psi)
{
#ifdef NWPW_SYCL
  if (mygdevice.hasgpu) mygdevice.hpsi_copy_gpu2host(npack,ne,psi);
#endif
}



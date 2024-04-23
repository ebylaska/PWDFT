// NWPW_CUDA Routines

#pragma once

#include "blas.h"

#define NDEV_MAX  39
#define DEBUG_IO  false

//#include        "gdevice.hpp"

#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cusolverDn.h>
//#include "cusolver_utils.h"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

#include <cuComplex.h>
#include <cublas_api.h>
#include <cuda_runtime_api.h>
#include <cusolverDn.h>
#include <library_types.h>

//
#ifndef cusolver_int_t
#define cusolver_int_t int
#endif

// cusolver API error checking
#define CUSOLVER_CHECK(err)                                                    \
  do {                                                                         \
    cusolverStatus_t err_ = (err);                                             \
    if (err_ != CUSOLVER_STATUS_SUCCESS) {                                     \
      printf("cusolver error %d at %s:%d\n", err_, __FILE__, __LINE__);        \
      throw std::runtime_error("cusolver error");                              \
    }                                                                          \
  } while (0)

class cuda_exception : public std::exception {

  std::string file_;
  int line_;
  cudaError_t err_code_;

  const char *what() const noexcept override {
    std::stringstream ss;
    ss << "CUDA Exception, " << cudaGetErrorString(err_code_) << " at "
       << std::endl
       << file_ << " : " << line_ << std::endl;
    auto msg = ss.str();
    return strdup(msg.c_str());
  }

public:
  cuda_exception(std::string file, int line, cudaError_t err)
      : file_(file), line_(line), err_code_(err) {}
};

class cufft_exception : public std::exception {

  std::string file_;
  int line_;
  cufftResult err_code_;

  static const char *_cudaGetErrorEnum(cufftResult error) {
    switch (error) {
    case CUFFT_SUCCESS:
      return "CUFFT_SUCCESS";

    case CUFFT_INVALID_PLAN:
      return "CUFFT_INVALID_PLAN";

    case CUFFT_ALLOC_FAILED:
      return "CUFFT_ALLOC_FAILED";

    case CUFFT_INVALID_TYPE:
      return "CUFFT_INVALID_TYPE";

    case CUFFT_INVALID_VALUE:
      return "CUFFT_INVALID_VALUE";

    case CUFFT_INTERNAL_ERROR:
      return "CUFFT_INTERNAL_ERROR";

    case CUFFT_EXEC_FAILED:
      return "CUFFT_EXEC_FAILED";

    case CUFFT_SETUP_FAILED:
      return "CUFFT_SETUP_FAILED";

    case CUFFT_INVALID_SIZE:
      return "CUFFT_INVALID_SIZE";

    case CUFFT_UNALIGNED_DATA:
      return "CUFFT_UNALIGNED_DATA";
    }
    return "<unknown>";
  }

  const char *what() const noexcept override {
    std::stringstream ss;
    ss << "CUFFT Exception, "
       << " Error Code: " << _cudaGetErrorEnum(err_code_) << std::endl
       << " at " << file_ << " : " << line_ << std::endl;
    auto msg = ss.str();
    return strdup(msg.c_str());
  }

public:
  cufft_exception(std::string file, int line, cufftResult err)
      : file_(file), line_(line), err_code_(err) {}
};

class cublas_exception : public std::exception {

  std::string file_;
  int line_;
  cublasStatus_t err_code_;

  const char *what() const noexcept override {
    std::stringstream ss;
    ss << "CUBLAS Exception, "
       << " Error Code: " << cublasGetStatusString(err_code_) << std::endl
       << " at " << file_ << " : " << line_ << std::endl;
    auto msg = ss.str();
    return strdup(msg.c_str());
  }

public:
  cublas_exception(std::string file, int line, cublasStatus_t err)
      : file_(file), line_(line), err_code_(err) {}
};

#define NWPW_CUDA_ERROR(ERR)                                                   \
  if (ERR != cudaSuccess)                                                      \
    throw cuda_exception(__FILE__, __LINE__, ERR);

#define NWPW_CUFFT_ERROR(ERR)                                                  \
  if (ERR != CUFFT_SUCCESS)                                                    \
    throw cufft_exception(__FILE__, __LINE__, ERR);

#define NWPW_CUBLAS_ERROR(ERR)                                                 \
  if (ERR != CUBLAS_STATUS_SUCCESS)                                            \
    throw cublas_exception(__FILE__, __LINE__, ERR);

/* Gdevices (CUDA) object -

   Operations:
      blas functions (using cublas):
      TN3_dgemm - computes <A|A>, <B|B>, and <A|B> overlaps
      TN1_dgemm - computes  <B|B> overlap

      TN_dgemm
      NT_dgemm
      psi_alloc
      psi_dealloc

      fft functions (uses cuFFT)


*/

class Gdevices {

   int fftcount = 0;
   int nxfft[2],nyfft[2], nzfft[2];
   cufftHandle forward_plan_x[2]  = {0,0}, plan_x[2]={0,0}, plan_y[2] = {0,0}, plan_z[2] = {0,0};
   cufftHandle backward_plan_x[2] = {0,0};
   int ifft_dev[15];
   int ifft_n;
 
   cublasHandle_t master_handle = 0;
   cublasOperation_t matC = CUBLAS_OP_C;
   cublasOperation_t matT = CUBLAS_OP_T;
   cublasOperation_t matN = CUBLAS_OP_N;
 
   cusolverEigMode_t jobz =
       CUSOLVER_EIG_MODE_VECTOR; // compute eigenvalues and eigenvectors.
   cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
   cusolverDnHandle_t cusolverH = NULL;
   cudaStream_t solverstream = NULL;
 
   cudaStream_t stream[12];

public:
   int  typegpu= 1;
   bool hasgpu = true;
 
   /* device memory */
   int ndev_mem = 0;
   bool inuse[NDEV_MAX] = {false};
   size_t ndsize_mem[NDEV_MAX];
   double *dev_mem[NDEV_MAX];
   int tile_fac = 1;
   int tile_npack2_max; // ,tile_npack1_max;
   int tile_npack2[19], tile_start2[19];
   int tile_npack1[19], tile_start1[19];

   double *a_psi, *a_hpsi, *b_prj;
   int ia_psi[2], ia_hpsi[2], ib_prj[2];
 
   int *d_info[2];
   int lwork = 0;
   double *d_work;
 
   /* constructor */
   /**************************************
    *                                    *
    *              Gdevices              *
    *                                    *
    **************************************/
   Gdevices() {
     ndev_mem = 0;
 
     if (DEBUG_IO) std::cout << "Into cublasCreate" << std::endl;
 
     NWPW_CUBLAS_ERROR(cublasCreate(&master_handle));
 
     // allocate cuda streams
     for (auto i=0; i<12; ++i)
       NWPW_CUDA_ERROR(cudaStreamCreate(&stream[i]));
 
     // create cusolver handle, bind a stream
     CUSOLVER_CHECK(cusolverDnCreate(&cusolverH));
     NWPW_CUDA_ERROR( cudaStreamCreateWithFlags(&solverstream, cudaStreamNonBlocking));
     CUSOLVER_CHECK(cusolverDnSetStream(cusolverH, solverstream));
 
     // query working space of syevd
     for (int i=0; i<2; ++i)
       NWPW_CUDA_ERROR(cudaMalloc(reinterpret_cast<void **>(&d_info[i]), sizeof(int)));
   }
 
   /* deconstructor */
   /**************************************
    *                                    *
    *              ~Gdevices             *
    *                                    *
    **************************************/
   ~Gdevices() noexcept(false) {
      // free dev_mem
      for (auto i=0; i<ndev_mem; ++i)
      {
         inuse[i] = false;
         NWPW_CUDA_ERROR(cudaFree(dev_mem[i]));
      }
      ndev_mem = 0;
     
      // free cuda streams
      for (auto i=0; i<12; ++i)
         NWPW_CUDA_ERROR(cudaStreamDestroy(stream[i]));
     
      NWPW_CUBLAS_ERROR(cublasDestroy(master_handle));
     
      /* free cusolver resources */
      for (auto i=0; i<2; ++i)
         NWPW_CUDA_ERROR(cudaFree(d_info[i]));
      NWPW_CUDA_ERROR(cudaFree(d_work));
      CUSOLVER_CHECK(cusolverDnDestroy(cusolverH));
     
      NWPW_CUDA_ERROR(cudaStreamDestroy(solverstream));
      cudaDeviceReset();
     
      // free fft descriptors
      // cufftDestroy(forward_plan_x);
      // cufftDestroy(plan_y);
      // cufftDestroy(plan_z);
      // cufftDestroy(backward_plan_x);
   }
 

   /**************************************
    *                                    *
    *           fetch_dev_mem_indx       *
    *                                    *
    **************************************/
   int fetch_dev_mem_indx(const size_t ndsize) 
   {
      int ii = 0;
      while ((((ndsize != ndsize_mem[ii]) || inuse[ii])) && (ii < ndev_mem))
         ++ii;
     
      if (ii < ndev_mem) 
      {
         inuse[ii] = true;
      } 
      else 
      {
         ii = ndev_mem;
         inuse[ii] = true;
         ndsize_mem[ii] = ndsize;
         NWPW_CUDA_ERROR(cudaMalloc((void **)&(dev_mem[ii]), ndsize * sizeof(double)));
         ndev_mem += 1;
         if (ndev_mem>NDEV_MAX) std::cout << "ERROR: ndev_mem > NDEV_MAX" << std::endl;
      }
     
      NWPW_CUDA_ERROR(cudaMemset(dev_mem[ii], 0, ndsize * sizeof(double)));
      return ii;
   }
 
   /**************************************
    *                                    *
    *              TN4_dgemm             *
    *                                    *
    **************************************/
   /* This function computes <host_a|host_a>, <host_a|host_b>, and
      <host_b|host_b> overlap matrices.
 
       host_caa = beta*host_caa + alpha*host_a'*host_a
       host_cab = beta*host_cab + alpha*host_a'*host_b
       host_cba = beta*host_cba + alpha*host_b'*host_a
       host_cbb = beta*host_cbb + alpha*host_b'*host_b
 
      Entry - npack2,ne: matrix size
              alpha, beta: standard dgemm parameters
              host_a: (npack2xne) matrix
              host_b: (npack2xne) matrix
      Exit - host_caa,host_cab,host_cba,host_cbb: (nexne) matrices
      Uses - device memory for (npack2xne) matrices ia_psi, and ia_hpsi allocated
      previously with psi_alloc
           - temporary device memory for (nexne) matrices ic11, ic12, ic21, and
      ic22.
   */
   void TN4_dgemm(int npack2, int ne, double alpha, double *host_a,
                  double *host_b, double beta, double *host_caa,
                  double *host_cab, double *host_cba, double *host_cbb) {
     int ic11 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
     int ic12 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
     int ic21 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
     int ic22 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
 
     if (std::fabs(beta) > 0.0) {
       NWPW_CUDA_ERROR(cudaMemcpy(dev_mem[ic11], host_caa, ne * ne * sizeof(double), cudaMemcpyHostToDevice));
       NWPW_CUDA_ERROR(cudaMemcpy(dev_mem[ic12], host_cab, ne * ne * sizeof(double), cudaMemcpyHostToDevice));
       NWPW_CUDA_ERROR(cudaMemcpy(dev_mem[ic21], host_cba, ne * ne * sizeof(double), cudaMemcpyHostToDevice));
       NWPW_CUDA_ERROR(cudaMemcpy(dev_mem[ic22], host_cbb, ne * ne * sizeof(double), cudaMemcpyHostToDevice));
     }
 
     // copy host_a,host_b --> dev_mem
     NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[0], ne, sizeof(double), &host_a[tile_start2[0]], npack2,
                                            dev_mem[ia_psi[0]], tile_npack2[0], stream[0]));
     NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[0], ne, sizeof(double), &host_b[tile_start2[0]], npack2,
                                            dev_mem[ia_hpsi[0]], tile_npack2[0], stream[0]));
 
     double beta0 = beta;
     for (auto tt = 0; tt < tile_fac; ++tt) {
       int ttp1 = tt + 1;
       if (ttp1 < tile_fac) {
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double), &host_a[tile_start2[ttp1]], npack2, 
                                                dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1], stream[ttp1 % 2]));
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double), &host_b[tile_start2[ttp1]], npack2, 
                                                dev_mem[ia_hpsi[ttp1 % 2]], tile_npack2[ttp1], stream[ttp1 % 2]));
       }
       NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[tt % 2]));
       NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
                                     dev_mem[ia_psi[tt % 2]], tile_npack2[tt], 
                                     dev_mem[ia_psi[tt % 2]], tile_npack2[tt], 
                                     &beta0, 
                                     dev_mem[ic11], ne));
       NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
                                     dev_mem[ia_psi[tt % 2]], tile_npack2[tt], 
                                     dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], 
                                     &beta0, 
                                     dev_mem[ic12], ne));
       NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
                                     dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], 
                                     dev_mem[ia_psi[tt % 2]], tile_npack2[tt], 
                                     &beta0, 
                                     dev_mem[ic21], ne));
       NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
                                     dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], 
                                     dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], 
                                     &beta0, 
                                     dev_mem[ic22], ne));
       beta0 = 1.0;
     }
 
     NWPW_CUDA_ERROR(cudaMemcpy(host_caa, dev_mem[ic11], ne * ne * sizeof(double), cudaMemcpyDeviceToHost));
     NWPW_CUDA_ERROR(cudaMemcpy(host_cab, dev_mem[ic12], ne * ne * sizeof(double), cudaMemcpyDeviceToHost));
     NWPW_CUDA_ERROR(cudaMemcpy(host_cba, dev_mem[ic21], ne * ne * sizeof(double), cudaMemcpyDeviceToHost));
     NWPW_CUDA_ERROR(cudaMemcpy(host_cbb, dev_mem[ic22], ne * ne * sizeof(double), cudaMemcpyDeviceToHost));
 
     inuse[ic11] = false;
     inuse[ic12] = false;
     inuse[ic21] = false;
     inuse[ic22] = false;
   }
 
   /**************************************
    *                                    *
    *              TN3_dgemm             *
    *                                    *
    **************************************/
   /* This function computes <host_a|host_a>, <host_a|host_b>, and
      <host_b|host_b> overlap matrices.
 
       host_caa = beta*host_caa + alpha*host_a'*host_a
       host_cab = beta*host_cab + alpha*host_a'*host_b
       host_cbb = beta*host_cbb + alpha*host_b'*host_b
 
      Entry - npack2,ne: matrix size
              alpha, beta: standard dgemm parameters
              host_a: (npack2xne) matrix
              host_b: (npack2xne) matrix
      Exit - host_caa,host_cab,host_cbb: (nexne) matrices
      Uses - device memory for (npack2xne) matrices ia_psi, and ia_hpsi allocated
      previously with psi_alloc
           - temporary device memory for (nexne) matrices ic11, ic12, and ic22.
   */
   void TN3_dgemm(int npack2, int ne, double alpha, double *host_a,
                  double *host_b, double beta, double *host_caa,
                  double *host_cab, double *host_cbb) 
   {
      int ic11 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
      int ic12 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
      int ic22 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
     
      if (std::fabs(beta) > 0.0) 
      {
         NWPW_CUDA_ERROR(cudaMemcpy(dev_mem[ic11], host_caa, ne * ne * sizeof(double), cudaMemcpyHostToDevice));
         NWPW_CUDA_ERROR(cudaMemcpy(dev_mem[ic12], host_cab, ne * ne * sizeof(double), cudaMemcpyHostToDevice));
         NWPW_CUDA_ERROR(cudaMemcpy(dev_mem[ic22], host_cbb, ne * ne * sizeof(double), cudaMemcpyHostToDevice));
      }
     
      // copy host_a,host_b --> dev_mem
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[0], ne, sizeof(double), &host_a[tile_start2[0]], npack2,
                                             dev_mem[ia_psi[0]], tile_npack2[0], stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[0], ne, sizeof(double), &host_b[tile_start2[0]], npack2,
                                             dev_mem[ia_hpsi[0]], tile_npack2[0], stream[0]));
     
      double beta0 = beta;
      for (auto tt = 0; tt < tile_fac; ++tt) 
      {
         int ttp1 = tt + 1;
         if (ttp1 < tile_fac) 
         {
            NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double), &host_a[tile_start2[ttp1]], npack2, 
                                                   dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1], stream[ttp1 % 2]));
            NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double), &host_b[tile_start2[ttp1]], npack2, 
                                                   dev_mem[ia_hpsi[ttp1 % 2]], tile_npack2[ttp1], stream[ttp1 % 2]));
         }
         NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[tt % 2]));
         NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
                                       dev_mem[ia_psi[tt % 2]], tile_npack2[tt], 
                                       dev_mem[ia_psi[tt % 2]], tile_npack2[tt], 
                                       &beta0, 
                                       dev_mem[ic11], ne));
         NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
                                       dev_mem[ia_psi[tt % 2]], tile_npack2[tt], 
                                       dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], 
                                       &beta0, 
                                       dev_mem[ic12], ne));
         NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
                                       dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], 
                                       dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], 
                                       &beta0, 
                                       dev_mem[ic22], ne));
         beta0 = 1.0;
      }
     
      NWPW_CUDA_ERROR(cudaMemcpy(host_caa, dev_mem[ic11], ne * ne * sizeof(double), cudaMemcpyDeviceToHost));
      NWPW_CUDA_ERROR(cudaMemcpy(host_cab, dev_mem[ic12], ne * ne * sizeof(double), cudaMemcpyDeviceToHost));
      NWPW_CUDA_ERROR(cudaMemcpy(host_cbb, dev_mem[ic22], ne * ne * sizeof(double), cudaMemcpyDeviceToHost));
     
      inuse[ic11] = false;
      inuse[ic12] = false;
      inuse[ic22] = false;
   }
 
   /**************************************
    *                                    *
    *              TN1_dgemm             *
    *                                    *
    **************************************/
   /* This function computes  <host_a|host_b> overlap matrix.
 
      host_cab = beta*host_cab + alpha*host_a'*host_b
 
      Entry - npack2,ne: matrix size
              alpha, beta: standard dgemm parameters
              host_a:  (npack2xne) matrix
              host_b:  (npack2xne) matrix
      Exit - host_cab: (nexne) matrices
      Uses - device memory for (npack2xne) matrices ia_psi, and ia_hpsi allocated
      previously with psi_alloc
           - temporary device memory for (nexne) matrix, ic12.
   */
   void TN1_dgemm(int npack2, int ne, double alpha, double *host_a,
                  double *host_b, double beta, double *host_cab) 
   {
      int ic12 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
     
      if (std::fabs(beta) > 0.0) 
      {
         NWPW_CUDA_ERROR(cudaMemcpy(dev_mem[ic12], host_cab, ne * ne * sizeof(double), cudaMemcpyHostToDevice));
      }
     
      // copy host_a,host_b --> dev_mem
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[0], ne, sizeof(double), &host_a[tile_start2[0]], npack2,
                                             dev_mem[ia_psi[0]], tile_npack2[0], stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[0], ne, sizeof(double), &host_b[tile_start2[0]], npack2,
                                             dev_mem[ia_hpsi[0]], tile_npack2[0], stream[0]));
     
      double beta0 = beta;
      for (auto tt = 0; tt < tile_fac; ++tt) 
      {
         int ttp1 = tt + 1;
         if (ttp1 < tile_fac) 
         {
            NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double), &host_a[tile_start2[ttp1]], npack2, 
                                                   dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1], stream[ttp1 % 2]));
            NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double), &host_b[tile_start2[ttp1]], npack2, 
                                                   dev_mem[ia_hpsi[ttp1 % 2]], tile_npack2[ttp1], stream[ttp1 % 2]));
         }
         NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[tt % 2]));
         NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
                                       dev_mem[ia_psi[tt % 2]], tile_npack2[tt], 
                                       dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], 
                                       &beta0, 
                                       dev_mem[ic12], ne));
         beta0 = 1.0;
      }
     
      NWPW_CUDA_ERROR(cudaMemcpy(host_cab, dev_mem[ic12], ne * ne * sizeof(double), cudaMemcpyDeviceToHost));
     
      inuse[ic12] = false;
   }
 
   /**************************************
    *                                    *
    *              TN_dgemm2c            *
    *                                    *
    **************************************/
    void TN_dgemm2c(int n, int m, int npack2, int nida2, double *host_a, double *host_b, double *host_c) 
    {
       double rtwo  = 2.0;
       double rone  = 1.0;
       double rmone = -1.0;
       double rzero = 0.0;
 
       int ia = fetch_dev_mem_indx(static_cast<size_t>(npack2) * n);
       int ib = fetch_dev_mem_indx(static_cast<size_t>(npack2) * m);
       int ic = fetch_dev_mem_indx(static_cast<size_t>(n) * m);
 
       //syclSetMatrixAsync(npack2,n,sizeof(double),host_a,npack2,dev_mem[ia],npack2,stream[0]);
       //syclSetMatrixAsync(npack2,m,sizeof(double),host_b,npack2,dev_mem[ib],npack2,stream[0]);
       NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(npack2,n,sizeof(double),host_a,npack2,dev_mem[ia],npack2,stream[0]));
       NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(npack2,m,sizeof(double),host_b,npack2,dev_mem[ib],npack2,stream[0]));
 
       NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
       NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matT, matN, n, m, npack2, &rtwo,
                                     dev_mem[ia], npack2, 
                                     dev_mem[ib], npack2, 
                                     &rzero, 
                                     dev_mem[ic],n));
 
       NWPW_CUDA_ERROR(cudaMemcpy(host_c,dev_mem[ic],n*m*sizeof(double),cudaMemcpyDeviceToHost));
 
       if (nida2 > 0) {
          DGEMM_PWDFT((char *) "T", (char *) "N", n, m, nida2, rmone, host_a, npack2, host_b, npack2, rzero, host_c, n);
       }
 
       inuse[ia] = false;
       inuse[ib] = false;
       inuse[ic] = false;
    }
 
 
 
 
   /**************************************
    *                                    *
    *              NN_dgemm              *
    *                                    *
    **************************************/
   void NN_dgemm(int npack2, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_c) 
   {
      // DGEMM_PWDFT((char *) "N",(char *)
      // "N",npack2,ne,ne,alpha,host_a,npack2,host_b,ne,beta,host_c,npack2);
     
      int ib = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
     
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(ne, ne, sizeof(double), host_b, ne, dev_mem[ib], ne, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[0], ne, sizeof(double), host_a+tile_start2[0], npack2,
                                             dev_mem[ia_psi[0]], tile_npack2[0], stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[0], ne, sizeof(double), host_c+tile_start2[0], npack2,
                                             dev_mem[ia_hpsi[0]], tile_npack2[0], stream[0]));
     
      // double beta0 = beta;
      for (auto tt = 0; tt < tile_fac; ++tt) 
      {
         int ttp1 = tt + 1;
         if (ttp1 < tile_fac) 
         {
            NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double), host_a+tile_start2[ttp1], npack2, 
                                                   dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1], stream[ttp1 % 2]));
            NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double), host_c+tile_start2[ttp1], npack2, 
                                                   dev_mem[ia_hpsi[ttp1 % 2]], tile_npack2[ttp1], stream[ttp1 % 2]));
         }
         NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[tt % 2]));
         NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matN, matN, tile_npack2[tt], ne, ne, &alpha, 
                                       dev_mem[ia_psi[tt % 2]],   tile_npack2[tt], 
                                       dev_mem[ib], ne, 
                                       &beta,
                                       dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt]));
         NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(tile_npack2[tt], ne, sizeof(double), dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], 
                                                host_c+tile_start2[tt], npack2, stream[tt % 2]));
      }
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[(tile_fac - 1) % 2]));
      // NWPW_CUDA_ERROR(
      // cudaMemcpy(host_c,dev_mem[ia_hpsi[0]],npack2*ne*sizeof(double),cudaMemcpyDeviceToHost)
      // );
     
      inuse[ib] = false;
   }
 
   /**************************************
    *                                    *
    *              NN_dgemm1             *
    *                                    *
    **************************************/
   void NN_dgemm1(int m, int n, int k,
                  double alpha,
                  double *host_a, int lda,
                  double *host_b, int ldb,
                  double beta,
                  double *host_c,int ldc) 
   {
      int ia = fetch_dev_mem_indx(((size_t)lda) * ((size_t)k));
      int ib = fetch_dev_mem_indx(((size_t)ldb) * ((size_t)n));
      int ic = fetch_dev_mem_indx(((size_t)ldc) * ((size_t)n));
     
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(lda, k, sizeof(double), host_a, lda, dev_mem[ia], lda, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(ldb, n, sizeof(double), host_b, ldb, dev_mem[ib], ldb, stream[0]));
     
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      NWPW_CUBLAS_ERROR(cublasDgemm(master_handle,matN,matN,m,n,k, 
                                    &alpha, 
                                    dev_mem[ia],lda,
                                    dev_mem[ib],ldb,
                                    &beta,
                                    dev_mem[ic],ldc));
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(ldc,n,sizeof(double),dev_mem[ic],ldc,host_c,ldc,stream[0]));
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
     
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }

 
   /**************************************
    *                                    *
    *              TN_dgemm2             *
    *                                    *
    **************************************/
   void TN_dgemm2(int m, int n, int k,
                  double alpha,
                  double *host_a, int lda,
                  double *host_b, int ldb,
                  double beta,
                  double *host_c,int ldc) 
   {
      int ia = fetch_dev_mem_indx(((size_t)lda) * ((size_t)m));
      int ib = fetch_dev_mem_indx(((size_t)ldb) * ((size_t)n));
      int ic = fetch_dev_mem_indx(((size_t)ldc) * ((size_t)n));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(lda, m, sizeof(double), host_a, lda, dev_mem[ia], lda, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(ldb, n, sizeof(double), host_b, ldb, dev_mem[ib], ldb, stream[0]));
     
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matT, matN, m, n, k, 
                                    &alpha, 
                                    dev_mem[ia], lda,
                                    dev_mem[ib], ldb,
                                    &beta,
                                    dev_mem[ic],ldc));
     
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(ldc,n,sizeof(double),dev_mem[ic],ldc,host_c,ldc,stream[0]));
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }
 

   /**************************************
    *                                    *
    *              NT_dgemm3             *
    *                                    *
    **************************************/
   void NT_dgemm3(int m, int n, int k,
                  double alpha,
                  double *host_a, int lda,
                  double *host_b, int ldb,
                  double beta,
                  double *host_c,int ldc) 
   {
      int ia = fetch_dev_mem_indx(((size_t)lda) * ((size_t)k));
      int ib = fetch_dev_mem_indx(((size_t)ldb) * ((size_t)k));
      int ic = fetch_dev_mem_indx(((size_t)ldc) * ((size_t)n));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(lda,k,sizeof(double),host_a,lda,dev_mem[ia],lda,stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(ldb,k,sizeof(double),host_b,ldb,dev_mem[ib],ldb,stream[0]));
     
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matN, matT, m, n, k,
                                    &alpha, 
                                    dev_mem[ia], lda,
                                    dev_mem[ib], ldb,
                                    &beta,
                                    dev_mem[ic],ldc));
     
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(ldc,n,sizeof(double),dev_mem[ic],ldc,host_c,ldc,stream[0]));
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }
 
 
   /**************************************
    *                                    *
    *              TN_dgemm              *
    *                                    *
    **************************************/
   void TN_dgemm(int ne, int nprj, int npack2, double alpha, double *host_a,
                 double *host_b, double beta, double *host_c) 
   {
      // DGEMM_PWDFT((char *) "T",(char *)
      // "N",ne,nprj,npack2,alpha,host_a,npack2,host_b,npack2,beta,host_c,ne);
     
      // gdevice_TN_dgemm(nn,nprj,ng,rtwo,a,b,rzero,sum);
     
      // int ia = fetch_dev_mem_indx(((size_t) npack2) * ((size_t) ne));
      // int ib = fetch_dev_mem_indx(((size_t) npack2) * ((size_t) nprj));
      b_prj = host_b;
      ib_prj[0] = fetch_dev_mem_indx(((size_t)tile_npack2_max) * ((size_t)nprj));
      if (tile_fac > 1)
         ib_prj[1] = fetch_dev_mem_indx(((size_t)tile_npack2_max) * ((size_t)nprj));
      int ic = fetch_dev_mem_indx(((size_t)ne) * ((size_t)nprj));
     
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(ne, nprj, sizeof(double), host_c, ne,
                                             dev_mem[ic], ne, stream[0]));
     
      if (tile_fac > 1)
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[0], ne, sizeof(double), a_psi+tile_start2[0], npack2,
                                                dev_mem[ia_psi[0]], tile_npack2[0], stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[0], nprj, sizeof(double), b_prj+tile_start2[0], npack2,
                                             dev_mem[ib_prj[0]], tile_npack2[0], stream[0]));
     
      double beta0 = beta;
      for (auto tt = 0; tt < tile_fac; ++tt) {
         int ttp1 = tt + 1;
         if (ttp1 < tile_fac) {
            NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double), a_psi+tile_start2[ttp1], npack2, 
                                                   dev_mem[ia_psi[ttp1%2]], tile_npack2[ttp1], stream[ttp1 % 2]));
            NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttp1], nprj, sizeof(double), b_prj+tile_start2[ttp1], npack2, 
                                                   dev_mem[ib_prj[ttp1%2]], tile_npack2[ttp1], stream[ttp1 % 2]));
         }
         NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[tt % 2]));
         NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matT, matN, ne, nprj, tile_npack2[tt], &alpha,
                                       dev_mem[ia_psi[tt % 2]], tile_npack2[tt], 
                                       dev_mem[ib_prj[tt % 2]], tile_npack2[tt], 
                                       &beta0, 
                                       dev_mem[ic], ne));
         beta0 = 1.0;
      }
      cudaMemcpy(host_c, dev_mem[ic], ne * nprj * sizeof(double), cudaMemcpyDeviceToHost);
     
      // inuse[ia] = false;
      // inuse[ib_prj[0]] = false;
      // if (tile_fac>1) inuse[ib_prj[1]] = false;
      inuse[ic] = false;
   }
 
 
   
   /**************************************
    *                                    *
    *              T_free                *
    *                                    *
    **************************************/
   void T_free() 
   {
      inuse[ib_prj[0]] = false;
      if (tile_fac > 1)
         inuse[ib_prj[1]] = false;
   }
 
   /**************************************
    *                                    *
    *             NT_dgemm               *
    *                                    *
    **************************************/
   void NT_dgemm(int npack2, int ne, int nprj, double alpha, double *host_a,
                 double *host_b, double beta, double *host_c) 
   {
      // DGEMM_PWDFT((char *) "N",(char *)
      // "T",npack2,ne,nprj,alpha,host_a,npack2,host_b,ne,beta,host_c,npack2);
     
      int ib = fetch_dev_mem_indx(((size_t)ne) * ((size_t)nprj));
     
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(ne, nprj, sizeof(double), host_b, ne,
                                             dev_mem[ib], ne, stream[(tile_fac - 1) % 2]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[tile_fac - 1], ne, sizeof(double), &host_c[tile_start2[tile_fac - 1]], npack2,
                                             dev_mem[ia_hpsi[(tile_fac - 1) % 2]], tile_npack2[tile_fac - 1], stream[(tile_fac - 1) % 2]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[tile_fac - 1], nprj, sizeof(double), &host_a[tile_start2[tile_fac - 1]], npack2,
                                             dev_mem[ib_prj[(tile_fac - 1) % 2]], tile_npack2[tile_fac - 1], stream[(tile_fac - 1) % 2]));
      for (auto tt = tile_fac - 1; tt >= 0; --tt) {
         int ttm1 = tt - 1;
         if (ttm1 >= 0) {
            NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttm1], ne, sizeof(double), &host_c[tile_start2[ttm1]], npack2, 
                                                   dev_mem[ia_hpsi[ttm1 % 2]], tile_npack2[ttm1], stream[ttm1 % 2]));
            NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttm1], nprj, sizeof(double), &host_a[tile_start2[ttm1]], npack2, 
                                                   dev_mem[ib_prj[ttm1 % 2]], tile_npack2[ttm1], stream[ttm1 % 2]));
         }
         NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[tt % 2]));
         NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matN, matT, tile_npack2[tt],
                                       ne, nprj, &alpha, 
                                       dev_mem[ib_prj[tt % 2]], tile_npack2[tt], 
                                       dev_mem[ib], ne, 
                                       &beta,
                                       dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt]));
         NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(tile_npack2[tt], ne, sizeof(double), dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], 
                                                &host_c[tile_start2[tt]], npack2, stream[tt % 2]));
      }
     
      inuse[ib] = false;
      inuse[ib_prj[0]] = false;
      if (tile_fac > 1)
         inuse[ib_prj[1]] = false;
   }
 
   /**************************************
    *                                    *
    *              MM6_dgemm             *
    *                                    *
    **************************************/
   void MM6_dgemm(int ne, double *host_s21, double *host_s12, double *host_s11,
                  double *host_sa0, double *host_sa1, double *host_st1) 
   {
      double rzero = 0.0;
      double rone = 1.0;
      int i_s21 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input
      int i_s12 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input
      int i_s11 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input
      int i_sa0 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input
      int i_st1 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // tmp
      int i_sa1 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input-output
     
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(ne, ne, sizeof(double), host_s12, ne, dev_mem[i_s12], ne, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(ne, ne, sizeof(double), host_s11, ne, dev_mem[i_s11], ne, stream[1]));
     
      NWPW_CUBLAS_ERROR(cublasSetMatrix(ne, ne, sizeof(double), host_s21, ne, dev_mem[i_s21], ne));
      NWPW_CUBLAS_ERROR(cublasSetMatrix(ne, ne, sizeof(double), host_sa0, ne, dev_mem[i_sa0], ne));
      NWPW_CUBLAS_ERROR(cublasSetMatrix(ne, ne, sizeof(double), host_s11, ne, dev_mem[i_s11], ne));
      NWPW_CUBLAS_ERROR(cublasSetMatrix(ne, ne, sizeof(double), host_sa1, ne, dev_mem[i_sa1], ne));
     
      // mmm_Multiply(ms, s21, sa0, 1.0, sa1, 1.0);
      NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matN, matN, ne, ne, ne, &rone,
                                    dev_mem[i_s21], ne, 
                                    dev_mem[i_sa0], ne, 
                                    &rone,
                                    dev_mem[i_sa1], ne));
     
      // mmm_Multiply(ms, sa0, s12, 1.0, sa1, 1.0);
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matN, matN, ne, ne, ne, &rone,
                                    dev_mem[i_sa0], ne, 
                                    dev_mem[i_s12], ne, 
                                    &rone,
                                    dev_mem[i_sa1], ne));
     
      // mmm_Multiply(ms, s11, sa0, 1.0, st1, 0.0);
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[1]));
      NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matN, matN, ne, ne, ne, &rone,
                                    dev_mem[i_s11], ne, 
                                    dev_mem[i_sa0], ne,
                                    &rzero, 
                                    dev_mem[i_st1], ne));
     
      // mmm_Multiply(ms, sa0, st1, 1.0, sa1, 1.0);
      NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matN, matN, ne, ne, ne, &rone,
                                    dev_mem[i_sa0], ne, 
                                    dev_mem[i_st1], ne, 
                                    &rone,
                                    dev_mem[i_sa1], ne));
     
      NWPW_CUBLAS_ERROR(cublasGetMatrix(ne, ne, sizeof(double), dev_mem[i_sa1], ne, host_sa1, ne));
     
      inuse[i_s21] = false;
      inuse[i_s12] = false;
      inuse[i_s11] = false;
      inuse[i_sa0] = false;
      inuse[i_st1] = false;
      inuse[i_sa1] = false;
   }


   /**************************************
    *                                    *
    *             isCmplxZero            *
    *                                    *
    **************************************/
   bool isCmplxZero(const double *beta) 
   {
      const double EPSILON = 1e-9; // Threshold for comparison
      return std::abs(beta[0]) < EPSILON && std::abs(beta[1]) < EPSILON;
   }


   /**************************************
    *                                    *
    *             NN1_zgemm              *
    *                                    *
    **************************************/
   void NN1_zgemm(int npack1_max, int npack, int ne, double *alpha, double *host_a, double *host_b,
                 double *beta, double *host_c) { 
      // Assuming fetch_dev_mem_indx, NWPW_CUBLAS_ERROR, NWPW_CUDA_ERROR, master_handle, matN,
      // dev_mem, stream, and inuse are properly defined and initialized elsewhere.
      // ZGEMM_PWDFT((char *)"N", (char *)"N", npack, ne, ne, alpha, host_a, npack1,
      // host_b, ne, beta, host_c, npack1);

      int ia = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int ib = fetch_dev_mem_indx(((size_t) 2*ne)     * ((size_t) ne));
      int ic = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
     
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_a, 2*npack1_max, dev_mem[ia], 2*npack1_max, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne,         ne, sizeof(double), host_b, 2*ne,         dev_mem[ib], 2*ne,         stream[0]));
      if (!isCmplxZero(beta))
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_c, 2*npack1_max, dev_mem[ic], 2*npack1_max, stream[0]));
      
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matN,matN,npack,ne,ne,
                                    reinterpret_cast<const cuDoubleComplex*>(alpha),
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia]),npack1_max,
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]),ne,
                                    reinterpret_cast<const cuDoubleComplex*>(beta),
                                    reinterpret_cast<cuDoubleComplex*>(dev_mem[ic]),npack1_max));
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*npack1_max,ne,sizeof(double),dev_mem[ic],2*npack1_max,host_c,2*npack1_max,stream[0]));
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
     
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }



   /**************************************
    *                                    *
    *              CN1_zgemm             *
    *                                    *
    **************************************/
   void CN1_zgemm(int npack1_max, int npack, int ne, double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) {
       // Assuming fetch_dev_mem_indx, NWPW_CUBLAS_ERROR, NWPW_CUDA_ERROR, master_handle, matN,
       // dev_mem, stream, and inuse are properly defined and initialized elsewhere.
       // ZGEMM_PWDFT((char *)"C", (char *)"N", ne, ne, npack, alpha, host_a, npack1,
       //             host_b, npack1, beta, host_c, ne);
     
      int ia = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int ib = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int ic = fetch_dev_mem_indx(((size_t) 2*ne)     * ((size_t) ne));
     
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_a, 2*npack1_max, dev_mem[ia], 2*npack1_max, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_b, 2*npack1_max, dev_mem[ib], 2*npack1_max, stream[0]));
      if (!isCmplxZero(beta))
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, ne, sizeof(double), host_c, 2*ne, dev_mem[ic], 2*ne, stream[0]));
     
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matC,matN,ne,ne,npack,
                                    reinterpret_cast<const cuDoubleComplex*>(alpha),
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia]),npack1_max,
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]),npack1_max,
                                    reinterpret_cast<const cuDoubleComplex*>(beta),
                                    reinterpret_cast<cuDoubleComplex*>(dev_mem[ic]),ne));
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*ne,ne,sizeof(double),dev_mem[ic],2*ne,host_c,2*ne,stream[0]));
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
     
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }


   /**************************************
    *                                    *
    *              CN2_zgemm             *
    *                                    *
    **************************************/
   void CN2_zgemm(int ne, int nprj, int npack, int npack1_max, double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) {
      // Assuming fetch_dev_mem_indx, NWPW_CUBLAS_ERROR, NWPW_CUDA_ERROR, master_handle, matN,
      // dev_mem, stream, and inuse are properly defined and initialized elsewhere.
      //  ZGEMM_PWDFT((char *)"C", (char *)"N", ne, nprj, npack, alpha, host_a, npack1,
      //              host_b, npack1, beta, host_c, ne);

      int ia = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int ib = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) nprj));
      int ic = fetch_dev_mem_indx(((size_t) 2*ne)     * ((size_t) nprj));
     
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*npack1_max, ne,   sizeof(double), host_a, 2*npack1_max, dev_mem[ia], 2*npack1_max, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*npack1_max, nprj, sizeof(double), host_b, 2*npack1_max, dev_mem[ib], 2*npack1_max, stream[0]));
      if (!isCmplxZero(beta))
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, nprj, sizeof(double), host_c, 2*ne, dev_mem[ic], 2*ne, stream[0]));
     
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matC,matN,ne,nprj,npack,
                                    reinterpret_cast<const cuDoubleComplex*>(alpha),
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia]),npack1_max,
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]),npack1_max,
                                    reinterpret_cast<const cuDoubleComplex*>(beta),
                                    reinterpret_cast<cuDoubleComplex*>(dev_mem[ic]),ne));
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*ne,nprj,sizeof(double),dev_mem[ic],2*ne,host_c,2*ne,stream[0]));
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
     
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }

   /**************************************
    *                                    *
    *            CN2_stride_zgemm        *
    *                                    *
    **************************************/
   void CN2_stride_zgemm(int ne, int nprj, int npack, int npack1_max, double *alpha, double *host_a,
                        double *host_b, double *beta, double *host_c) {

      b_prj = host_b;
      ib_prj[0] = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t)nprj));
      if (tile_fac > 1)
         ib_prj[1] = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t)nprj));
      int ic = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t)nprj));

      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, nprj, sizeof(double), host_c, 2*ne, dev_mem[ic], 2*ne, stream[0]));

      if (tile_fac > 1)
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[0], ne, sizeof(double), a_psi+tile_start2[0], 2*npack1_max, dev_mem[ia_psi[0]], tile_npack2[0], stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[0],  nprj, sizeof(double), b_prj+tile_start2[0], 2*npack1_max, dev_mem[ib_prj[0]], tile_npack2[0], stream[0]));

      double beta0[2] = {beta[0],beta[1]};
      for (auto tt = 0; tt < tile_fac; ++tt) 
       {
         int ttp1 = tt + 1;
         if (ttp1 < tile_fac) 
         {
            NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttp1], ne, sizeof(double), a_psi+tile_start2[ttp1], 2*npack1_max,
                                                   dev_mem[ia_psi[ttp1%2]], tile_npack2[ttp1], stream[ttp1 % 2]));
            NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttp1], nprj, sizeof(double), b_prj+tile_start2[ttp1], 2*npack1_max,
                                                   dev_mem[ib_prj[ttp1%2]], tile_npack2[ttp1], stream[ttp1 % 2]));
         }
         NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[tt % 2]));
         NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matC,matN,ne,nprj,tile_npack1[tt],
                                       reinterpret_cast<const cuDoubleComplex*>(alpha),
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia_psi[tt%2]]),tile_npack1[tt],
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib_prj[tt%2]]),tile_npack1[tt],
                                       reinterpret_cast<const cuDoubleComplex*>(beta0),
                                       reinterpret_cast<cuDoubleComplex*>(dev_mem[ic]),ne));

         //NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matT, matN, ne, nprj, tile_npack2[tt], &alpha,
         //                              dev_mem[ia_psi[tt % 2]], tile_npack2[tt],
         //                              dev_mem[ib_prj[tt % 2]], tile_npack2[tt],
         //                              &beta0,
         //                              dev_mem[ic], ne));
         beta0[0] = 1.0; 
         beta0[1] = 0.0;
      }
      cudaMemcpy(host_c, dev_mem[ic], 2*ne*nprj*sizeof(double), cudaMemcpyDeviceToHost);

      inuse[ic] = false;
   }



   /**************************************
    *                                    *
    *             NC2_zgemm              *
    *                                    *
    **************************************/
   void NC2_zgemm(int npack1_max, int npack, int ne, int nprj,  double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) {
      // Assuming fetch_dev_mem_indx, NWPW_CUBLAS_ERROR, NWPW_CUDA_ERROR, master_handle, matN,
      // dev_mem, stream, and inuse are properly defined and initialized elsewhere.
      // ZGEMM_PWDFT((char *)"N", (char *)"C", npack,ne,nprj, alpha, host_a, npack1,
      //             host_b, ne, beta, host_c, npack1);

      int ia = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) nprj));
      int ib = fetch_dev_mem_indx(((size_t) 2*ne)         * ((size_t) nprj));
      int ic = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));

      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*npack1_max,  nprj, sizeof(double), host_a, 2*npack1_max, dev_mem[ia], 2*npack1_max, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne,          nprj, sizeof(double), host_b, 2*ne,         dev_mem[ib], 2*ne,         stream[0]));
      if (!isCmplxZero(beta))
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_c, 2*npack1_max, dev_mem[ic], 2*npack1_max, stream[0]));
                
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matN,matC,npack,ne,nprj,
                                    reinterpret_cast<const cuDoubleComplex*>(alpha),
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia]),npack1_max,
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]),ne,
                                    reinterpret_cast<const cuDoubleComplex*>(beta),
                                    reinterpret_cast<cuDoubleComplex*>(dev_mem[ic]),npack1_max));
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*npack1_max,ne,sizeof(double),dev_mem[ic],2*npack1_max,host_c,2*npack1_max,stream[0]));
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));

      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }   



   /**************************************
    *                                    *
    *          NC2_stride_zgemm          *
    *                                    *
    **************************************/
   void NC2_stride_zgemm(int npack1_max, int npack, int ne, int nprj,  double *alpha, double *host_a,
                         double *host_b, double *beta, double *host_c) {
      // Assuming fetch_dev_mem_indx, NWPW_CUBLAS_ERROR, NWPW_CUDA_ERROR, master_handle, matN,
      // dev_mem, stream, and inuse are properly defined and initialized elsewhere.
      // ZGEMM_PWDFT((char *)"N", (char *)"C", npack,ne,nprj, alpha, host_a, npack1,
      //             host_b, ne, beta, host_c, npack1);

      int ib = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t) nprj));

      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, nprj, sizeof(double), host_b, 2*ne, dev_mem[ib], 2*ne, stream[(tile_fac-1) % 2]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[tile_fac-1], ne, sizeof(double), host_c + tile_start2[tile_fac-1], 2*npack1_max,
                                             dev_mem[ia_hpsi[(tile_fac-1) % 2]], tile_npack2[tile_fac-1], stream[(tile_fac - 1) % 2]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[tile_fac-1], nprj, sizeof(double), host_a+tile_start2[tile_fac-1], 2*npack1_max,
                                             dev_mem[ib_prj[(tile_fac-1) % 2]], tile_npack2[tile_fac-1], stream[(tile_fac-1) % 2]));
      for (auto tt = tile_fac - 1; tt >= 0; --tt) 
      {
         int ttm1 = tt - 1;
         if (ttm1 >= 0) {
            NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttm1], ne, sizeof(double), host_c+tile_start2[ttm1], 2*npack1_max,
                                                   dev_mem[ia_hpsi[ttm1%2]], tile_npack2[ttm1], stream[ttm1 % 2]));
            NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[ttm1], nprj, sizeof(double), host_a+tile_start2[ttm1], 2*npack1_max,
                                                   dev_mem[ib_prj[ttm1%2]], tile_npack2[ttm1], stream[ttm1 % 2]));
         }
         NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[tt%2]));

         NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matN,matC,tile_npack1[tt],ne,nprj,
                                    reinterpret_cast<const cuDoubleComplex*>(alpha),
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib_prj[tt%2]]),tile_npack1[tt],
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]),ne,
                                    reinterpret_cast<const cuDoubleComplex*>(beta),
                                    reinterpret_cast<cuDoubleComplex*>(dev_mem[ia_hpsi[tt]]),tile_npack1[tt]));
         //NWPW_CUBLAS_ERROR(cublasDgemm(master_handle, matN, matT, tile_npack2[tt],
         //                              ne, nprj, &alpha,
         //                              dev_mem[ib_prj[tt % 2]], tile_npack2[tt],
         //                              dev_mem[ib], ne,
         //                              &beta,
         //                              dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt]));
         NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(tile_npack2[tt], ne, sizeof(double), dev_mem[ia_hpsi[tt%2]], tile_npack2[tt],
                                                host_c+tile_start2[tt], 2*npack1_max, stream[tt%2]));
      }

      inuse[ib] = false;
      inuse[ib_prj[0]] = false;
      if (tile_fac > 1)
         inuse[ib_prj[1]] = false;
   }

      
                 
   /**************************************
    *                                    *
    *              NN_zgemm              *
    *                                    *
    **************************************/
   void NN_zgemm(int m, int n, int k,
                 double *alpha,
                 double *host_a, int lda,
                 double *host_b, int ldb,
                 double *beta,
                 double *host_c,int ldc) {
      // ZGEMM_PWDFT((char *)"N", (char *)"N", m, n, k, alpha, host_a, lda, host_b, ldb, beta, host_c, ldc);

      int ia = fetch_dev_mem_indx(((size_t) 2*lda) * ((size_t) k));
      int ib = fetch_dev_mem_indx(((size_t) 2*ldb) * ((size_t) n));
      int ic = fetch_dev_mem_indx(((size_t) 2*ldc) * ((size_t) n));

      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*lda, k, sizeof(double), host_a, 2*lda, dev_mem[ia], 2*lda, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ldb, n, sizeof(double), host_b, 2*ldb, dev_mem[ib], 2*ldb, stream[0]));
      if (!isCmplxZero(beta))
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ldc, n, sizeof(double), host_c, 2*ldc, dev_mem[ic], 2*ldc, stream[0]));

      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matN,matN,m,n,k,
                                    reinterpret_cast<const cuDoubleComplex*>(alpha),
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia]),lda,
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]),ldb,
                                    reinterpret_cast<const cuDoubleComplex*>(beta),
                                    reinterpret_cast<cuDoubleComplex*>(dev_mem[ic]),ldc));
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*ldc,n,sizeof(double),dev_mem[ic],2*ldc,host_c,2*ldc,stream[0]));
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));

      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }    


   /**************************************
    *                                    *
    *              CN_zgemm              *
    *                                    *
    **************************************/
   void CN_zgemm(int m, int n, int k,
                 double *alpha,
                 double *host_a, int lda,
                 double *host_b, int ldb,
                 double *beta,
                 double *host_c,int ldc) {
      // ZGEMM_PWDFT((char *)"C", (char *)"N", m, n, k, alpha, host_a, lda, host_b, ldb, beta, host_c, ldc);
     
      int ia = fetch_dev_mem_indx(((size_t) 2*lda) * ((size_t) m));
      int ib = fetch_dev_mem_indx(((size_t) 2*ldb) * ((size_t) n));
      int ic = fetch_dev_mem_indx(((size_t) 2*ldc) * ((size_t) n));
      
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*lda, m, sizeof(double), host_a, 2*lda, dev_mem[ia], 2*lda, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ldb, n, sizeof(double), host_b, 2*ldb, dev_mem[ib], 2*ldb, stream[0]));
      if (!isCmplxZero(beta))
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ldc, n, sizeof(double), host_c, 2*ldc, dev_mem[ic], 2*ldc, stream[0]));
      
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matC,matN,m,n,k,
                                    reinterpret_cast<const cuDoubleComplex*>(alpha),
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia]),lda,
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]),ldb,
                                    reinterpret_cast<const cuDoubleComplex*>(beta),
                                    reinterpret_cast<cuDoubleComplex*>(dev_mem[ic]),ldc));
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*ldc,n,sizeof(double),dev_mem[ic],2*ldc,host_c,2*ldc,stream[0]));
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }


   /**************************************
    *                                    *
    *              NC_zgemm              *
    *                                    *
    **************************************/
   void NC_zgemm(int m, int n, int k,
                 double *alpha,
                 double *host_a, int lda,
                 double *host_b, int ldb,
                 double *beta,
                 double *host_c,int ldc) {
      // ZGEMM_PWDFT((char *)"N", (char *)"C", m, n, k, alpha, host_a, lda, host_b, ldb, beta, host_c, ldc);
     
      int ia = fetch_dev_mem_indx(((size_t) 2*lda) * ((size_t) k));
      int ib = fetch_dev_mem_indx(((size_t) 2*ldb) * ((size_t) n));
      int ic = fetch_dev_mem_indx(((size_t) 2*ldc) * ((size_t) n));
      
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*lda, k, sizeof(double), host_a, 2*lda, dev_mem[ia], 2*lda, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ldb, n, sizeof(double), host_b, 2*ldb, dev_mem[ib], 2*ldb, stream[0]));
      if (!isCmplxZero(beta))
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ldc, n, sizeof(double), host_c, 2*ldc, dev_mem[ic], 2*ldc, stream[0]));
      
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matN,matC,m,n,k,
                                    reinterpret_cast<const cuDoubleComplex*>(alpha),
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia]),lda,
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]),ldb,
                                    reinterpret_cast<const cuDoubleComplex*>(beta),
                                    reinterpret_cast<cuDoubleComplex*>(dev_mem[ic]),ldc));
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*ldc,n,sizeof(double),dev_mem[ic],2*ldc,host_c,2*ldc,stream[0]));
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[ic] = false;
   }



   /**************************************
    *                                    *
    *              CN3_zgemm             *
    *                                    *
    **************************************/
   void CN3_zgemm(int npack1_max, int npack, int ne, double *alpha, double *host_a,
                 double *host_b, double *beta, double *host_caa,
                 double *host_cab, double *host_cbb) {
      int one = 1;
      int shift1 = 0;
      int mshift1 = 0;

      int ia = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int ib = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int icaa = fetch_dev_mem_indx(((size_t) 2*ne)   * ((size_t) ne));
      int icab = fetch_dev_mem_indx(((size_t) 2*ne)   * ((size_t) ne));
      int icbb = fetch_dev_mem_indx(((size_t) 2*ne)   * ((size_t) ne));

      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_a, 2*npack1_max, dev_mem[ia], 2*npack1_max, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_b, 2*npack1_max, dev_mem[ib], 2*npack1_max, stream[1]));
      if (!isCmplxZero(beta))
      {
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, ne, sizeof(double), host_caa, 2*ne, dev_mem[icaa], 2*ne, stream[0]));
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, ne, sizeof(double), host_cab, 2*ne, dev_mem[icab], 2*ne, stream[1]));
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, ne, sizeof(double), host_cbb, 2*ne, dev_mem[icbb], 2*ne, stream[1]));
      }
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));

      for (auto k=1; k<=ne; ++k)
      {
         /*
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack,
                     alpha,
                     host_a, npack1,
                     host_a + shift1, npack1,
                     beta,
                     host_caa + mshift1, k);
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack,
                     alpha,
                     host_a, npack1,
                     host_b + shift1, npack1,
                     beta,
                     host_cab + mshift1, k);
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack,
                     alpha,
                     host_b, npack1,
                     host_b + shift1, npack1,
                     beta,
                     host_cbb + mshift1, k);
         */
         NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matC,matN,k,one,npack,
                                       reinterpret_cast<const cuDoubleComplex*>(alpha),
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia]),npack1_max,
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia]+shift1),npack1_max,
                                       reinterpret_cast<const cuDoubleComplex*>(beta),
                                       reinterpret_cast<cuDoubleComplex*>(dev_mem[icaa]+mshift1),k));
         shift1 += 2*npack1_max;
         mshift1 += 2*ne;
      }
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*ne,ne,sizeof(double),dev_mem[icaa],2*ne,host_caa,2*ne,stream[0]));


      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[1]));
      for (auto k=1; k<=ne; ++k)
      {
         NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matC,matN,k,one,npack,
                                       reinterpret_cast<const cuDoubleComplex*>(alpha),
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia]),npack1_max,
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]+shift1),npack1_max,
                                       reinterpret_cast<const cuDoubleComplex*>(beta),
                                       reinterpret_cast<cuDoubleComplex*>(dev_mem[icab]+mshift1),k));
         NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matC,matN,k,one,npack,
                                       reinterpret_cast<const cuDoubleComplex*>(alpha),
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]),npack1_max,
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]+shift1),npack1_max,
                                       reinterpret_cast<const cuDoubleComplex*>(beta),
                                       reinterpret_cast<cuDoubleComplex*>(dev_mem[icbb]+mshift1),k));

         shift1 += 2*npack1_max;
         mshift1 += 2*ne;
      }
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));

      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*ne,ne,sizeof(double),dev_mem[icab],2*ne,host_cab,2*ne,stream[1]));
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*ne,ne,sizeof(double),dev_mem[icbb],2*ne,host_cbb,2*ne,stream[1]));
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[1]));
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[icaa] = false;
      inuse[icab] = false;
      inuse[icbb] = false;
   }

   /**************************************
    *                                    *
    *              CN4_zgemm             *
    *                                    *
    **************************************/
   void CN4_zgemm(int npack1_max, int npack, int ne, double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_caa,
                  double *host_cab, double *host_cba, double *host_cbb) {
      int one = 1;
      int shift1 = 0;
      int mshift1 = 0;

      int ia = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int ib = fetch_dev_mem_indx(((size_t) 2*npack1_max) * ((size_t) ne));
      int icaa = fetch_dev_mem_indx(((size_t) 2*ne)   * ((size_t) ne));
      int icab = fetch_dev_mem_indx(((size_t) 2*ne)   * ((size_t) ne));
      int icba = fetch_dev_mem_indx(((size_t) 2*ne)   * ((size_t) ne));
      int icbb = fetch_dev_mem_indx(((size_t) 2*ne)   * ((size_t) ne));

      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_a, 2*npack1_max, dev_mem[ia], 2*npack1_max, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*npack1_max, ne, sizeof(double), host_b, 2*npack1_max, dev_mem[ib], 2*npack1_max, stream[1]));
      if (!isCmplxZero(beta))
      {
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, ne, sizeof(double), host_caa, 2*ne, dev_mem[icaa], 2*ne, stream[0]));
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, ne, sizeof(double), host_cab, 2*ne, dev_mem[icab], 2*ne, stream[1]));
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, ne, sizeof(double), host_cba, 2*ne, dev_mem[icba], 2*ne, stream[1]));
         NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, ne, sizeof(double), host_cbb, 2*ne, dev_mem[icbb], 2*ne, stream[1]));
      }
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));

      for (auto k=1; k<=ne; ++k)
      {
         /*
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack,
                     alpha,
                     host_a, npack1,
                     host_a + shift1, npack1,
                     beta,
                     host_caa + mshift1, k);
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack,
                     alpha,
                     host_a, npack1,
                     host_b + shift1, npack1,
                     beta,
                     host_cab + mshift1, k);
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack,
                     alpha,
                     host_b, npack1,
                     host_a + shift1, npack1,
                     beta,
                     host_cba + mshift1, k);
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack,
                     alpha,
                     host_b, npack1,
                     host_b + shift1, npack1,
                     beta,
                     host_cbb + mshift1, k);
         */
         NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matC,matN,k,one,npack,
                                       reinterpret_cast<const cuDoubleComplex*>(alpha),
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia]),npack1_max,
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia]+shift1),npack1_max,
                                       reinterpret_cast<const cuDoubleComplex*>(beta),
                                       reinterpret_cast<cuDoubleComplex*>(dev_mem[icaa]+mshift1),k));
         shift1 += 2*npack1_max;
         mshift1 += 2*ne;
      }
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*ne,ne,sizeof(double),dev_mem[icaa],2*ne,host_caa,2*ne,stream[0]));

      shift1 = 0;
      mshift1 = 0;
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[1]));
      for (auto k=1; k<=ne; ++k)
      {
         NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matC,matN,k,one,npack,
                                       reinterpret_cast<const cuDoubleComplex*>(alpha),
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia]),npack1_max,
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]+shift1),npack1_max,
                                       reinterpret_cast<const cuDoubleComplex*>(beta),
                                       reinterpret_cast<cuDoubleComplex*>(dev_mem[icab]+mshift1),k));
         NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matC,matN,k,one,npack,
                                       reinterpret_cast<const cuDoubleComplex*>(alpha),
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]),npack1_max,
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ia]+shift1),npack1_max,
                                       reinterpret_cast<const cuDoubleComplex*>(beta),
                                       reinterpret_cast<cuDoubleComplex*>(dev_mem[icba]+mshift1),k));
         NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matC,matN,k,one,npack,
                                       reinterpret_cast<const cuDoubleComplex*>(alpha),
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]),npack1_max,
                                       reinterpret_cast<const cuDoubleComplex*>(dev_mem[ib]+shift1),npack1_max,
                                       reinterpret_cast<const cuDoubleComplex*>(beta),
                                       reinterpret_cast<cuDoubleComplex*>(dev_mem[icbb]+mshift1),k));
         shift1 += 2*npack1_max;
         mshift1 += 2*ne;
      }
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*ne,ne,sizeof(double),dev_mem[icab],2*ne,host_cab,2*ne,stream[1]));
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*ne,ne,sizeof(double),dev_mem[icba],2*ne,host_cba,2*ne,stream[1]));
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*ne,ne,sizeof(double),dev_mem[icbb],2*ne,host_cbb,2*ne,stream[1]));
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[1]));
      
      inuse[ia] = false;
      inuse[ib] = false;
      inuse[icaa] = false;
      inuse[icab] = false;
      inuse[icba] = false;
      inuse[icbb] = false;
   }


   /**************************************
    *                                    *
    *              WW6_zgemm             *
    *                                    *
    **************************************/
   void WW6_zgemm(int ne, double *host_s21, double *host_s12, double *host_s11,
                 double *host_sa0, double *host_sa1, double *host_st1) {
      double rone[2]  = {1.0,0.0};
      double rzero[2] = {0.0,0.0};

      int is11 = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t) ne));
      int is12 = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t) ne));
      int is21 = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t) ne));
      int isa0 = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t) ne));
      int isa1 = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t) ne));
      int ist1 = fetch_dev_mem_indx(((size_t) 2*ne) * ((size_t) ne));

      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, ne, sizeof(double), host_s21, 2*ne, dev_mem[is21], 2*ne, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, ne, sizeof(double), host_sa0, 2*ne, dev_mem[isa0], 2*ne, stream[0]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, ne, sizeof(double), host_sa1, 2*ne, dev_mem[isa1], 2*ne, stream[0]));

      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, ne, sizeof(double), host_s12, 2*ne, dev_mem[is12], 2*ne, stream[1]));
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(2*ne, ne, sizeof(double), host_s11, 2*ne, dev_mem[is11], 2*ne, stream[2]));

      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));

      // www_Multiply1(ms, s21, sa0, 1.0, sa1, 1.0);
      //ZGEMM_PWDFT((char *)"N", (char *)"N", ne, ne, ne, rone, host_s21, ne,
      //            host_sa0, ne, rone, host_sa1, ne);
      NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matN,matN,ne,ne,ne,
                                    reinterpret_cast<const cuDoubleComplex*>(rone),
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[is21]),ne,
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[isa0]),ne,
                                    reinterpret_cast<const cuDoubleComplex*>(rone),
                                    reinterpret_cast<cuDoubleComplex*>(dev_mem[isa1]),ne));

      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[1]));
      // www_Multiply2(ms, sa0, s12, 1.0, sa1, 1.0);
      //ZGEMM_PWDFT((char *)"C", (char *)"N", ne, ne, ne, rone, host_sa0, ne,
      //            host_s12, ne, rone, host_sa1, ne);
      
      NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matN,matN,ne,ne,ne,
                                    reinterpret_cast<const cuDoubleComplex*>(rone),
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[isa0]),ne,
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[is12]),ne,
                                    reinterpret_cast<const cuDoubleComplex*>(rone),
                                    reinterpret_cast<cuDoubleComplex*>(dev_mem[isa1]),ne));
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[2]));

      // www_Multiply3(ms, s11, sa0, 1.0, st1, 0.0);
      //ZGEMM_PWDFT((char *)"N", (char *)"C", ne, ne, ne, rone, host_s11, ne,
      //            host_sa0, ne, rzero, host_st1, ne);
      
      NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matN,matN,ne,ne,ne,
                                    reinterpret_cast<const cuDoubleComplex*>(rone),
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[is11]),ne,
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[isa0]),ne,
                                    reinterpret_cast<const cuDoubleComplex*>(rzero),
                                    reinterpret_cast<cuDoubleComplex*>(dev_mem[ist1]),ne));

      // www_Multiply1(ms, sa0, st1, 1.0, sa1, 1.0);
      //ZGEMM_PWDFT((char *)"N", (char *)"N", ne, ne, ne, rone, host_sa0, ne,
      //            host_st1, ne, rone, host_sa1, ne);
      NWPW_CUBLAS_ERROR(cublasZgemm(master_handle,matN,matN,ne,ne,ne,
                                    reinterpret_cast<const cuDoubleComplex*>(rone),
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[isa0]),ne,
                                    reinterpret_cast<const cuDoubleComplex*>(dev_mem[ist1]),ne,
                                    reinterpret_cast<const cuDoubleComplex*>(rone),
                                    reinterpret_cast<cuDoubleComplex*>(dev_mem[isa1]),ne));

      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*ne,ne,sizeof(double),dev_mem[ist1],2*ne,host_st1,2*ne,stream[0]));
      NWPW_CUBLAS_ERROR(cublasGetMatrixAsync(2*ne,ne,sizeof(double),dev_mem[isa1],2*ne,host_sa1,2*ne,stream[0]));
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));

      inuse[is11] = false;
      inuse[is12] = false;
      inuse[is21] = false;
      inuse[isa0] = false;
      inuse[isa1] = false;
      inuse[ist1] = false;
  }     

   /**************************************
    *                                    *
    *          WW_eigensolver            *
    *                                    *
    **************************************/
   void WW_eigensolver(int ispin, int ne[], double *host_hml, double *host_eig)
   {     
      int n, ierr;
      int nn = ne[0] * ne[0] + 14; 
      double xmp1[nn];
      double rmp1[nn];
      // double *xmp1 = new (std::nothrow) double[nn]();
            
      int shift1 = 0;
      int shift2 = 0;
      for (int ms=0; ms<ispin; ++ms)
      {
         n = ne[ms];
       
         // eigen_(&n,&n,&hml[shift2],&eig[shift1],xmp1,&ierr);
         //  d3db::parall->Barrier();
         ZEIGEN_PWDFT(n, host_hml + shift2, host_eig + shift1, xmp1, nn, rmp1, ierr);
         // if (ierr != 0) throw std::runtime_error(std::string("NWPW Error:
         // EIGEN_PWDFT failed!"));
            
         //eigsrt_device(host_eig + shift1, host_hml + shift2, n);
         shift1 += ne[0];
         shift2 += 2*ne[0]*ne[0]; 
      }
   }




 
   /********************/
   /* psi_dev functions*/
   /********************/
 
   /**************************************
    *                                    *
    *              psi_alloc             *
    *                                    *
    **************************************/
   void psi_alloc(int npack1, int ne, int tile_fac0 = 1) 
   {
      tile_fac = tile_fac0;
     
      tile_npack2_max = ( (  (2*npack1) % tile_fac) == 0)
                           ? (2*npack1) / tile_fac
                           : (2*npack1) / tile_fac + 1;
      //tile_npack1_max = tile_npack2_max / 2;
      // for (auto i=0; i<tile_fac; ++i) tile_npack2[i] =
      // (i<((2*npack1)%tile_fac)) ? (2*npack1)/tile_fac+1 : (2*npack1)/tile_fac;
      for (auto i=0; i<tile_fac; ++i)
      {
         tile_npack1[i] =   (npack1) / tile_fac;
         tile_npack2[i] = (2*npack1) / tile_fac;
      }
      for (auto i=0; i<((npack1) % tile_fac); ++i)
         tile_npack1[i] += 1;
      for (auto i=0; i<((2*npack1) % tile_fac); ++i)
         tile_npack2[i] += 1;
     
      tile_start1[0] = 0;
      tile_start2[0] = 0;
      for (auto i=1; i<tile_fac; ++i)
      {
         tile_start1[i] = tile_start1[i-1] + tile_npack1[i-1];
         tile_start2[i] = tile_start2[i-1] + tile_npack2[i-1];
      }
     
      ia_psi[0]  = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) ne));
      ia_hpsi[0] = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) ne));
     
      if (tile_fac > 1) {
         ia_psi[1]  = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) ne));
         ia_hpsi[1] = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) ne));
      }
      if (DEBUG_IO) std::cout << "Into psi_alloc, tile_factor = " << tile_fac << " ndev_mem=" << ndev_mem << std::endl;
   }
 
   /**************************************
    *                                    *
    *              psi_dealloc           *
    *                                    *
    **************************************/
   void psi_dealloc() 
   {
      inuse[ia_psi[0]] = false;
      inuse[ia_hpsi[0]] = false;
      if (tile_fac > 1) 
      {
         inuse[ia_psi[1]] = false;
         inuse[ia_hpsi[1]] = false;
      }
   }
 
   /**************************************
    *                                    *
    *           psi_copy_host2gpu        *
    *                                    *
    **************************************/
   void psi_copy_host2gpu(int npack1, int ne, double *psi) 
   {
      // cudaMemcpy(dev_mem[ia_psi[0]],psi,tile_npack2_max*ne*sizeof(double),cudaMemcpyHostToDevice);
      a_psi = psi;
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[0], ne, sizeof(double), psi, 2*npack1, 
                                             dev_mem[ia_psi[0]], tile_npack2[0], stream[0]));
   }
 
   /**************************************
    *                                    *
    *          hpsi_copy_host2gpu        *
    *                                    *
    **************************************/
   void hpsi_copy_host2gpu(int npack1, int ne, double *hpsi) 
   {
      // cudaMemcpy(dev_mem[ia_hpsi[0]],hpsi,2*npack1*ne*sizeof(double),cudaMemcpyHostToDevice);
      int tt = tile_fac - 1;
      a_hpsi = hpsi;
      NWPW_CUBLAS_ERROR(cublasSetMatrixAsync(tile_npack2[tt], ne, sizeof(double), &hpsi[tile_start2[tt]], 2 * npack1,
                                             dev_mem[ia_hpsi[tt%2]], tile_npack2[tt], stream[tt%2]));
   }
 
   /**************************************
    *                                    *
    *           psi_copy_gpu2host        *
    *                                    *
    **************************************/
   void psi_copy_gpu2host(int npack1, int ne, double *psi) 
   {
      // cudaMemcpy(psi, dev_mem[ia_psi[0]],
      // 2*ne*npack1*sizeof(double),cudaMemcpyDeviceToHost);
      if (tile_fac == 1) 
      {
         NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
         NWPW_CUBLAS_ERROR(cublasGetMatrix(tile_npack2[0], ne, sizeof(double), dev_mem[ia_psi[0]], tile_npack2[0], 
                                           psi, 2*npack1));
      }
   }
 
 
   /**************************************
    *                                    *
    *          hpsi_copy_gpu2host        *
    *                                    *
    **************************************/
   void hpsi_copy_gpu2host(int npack1, int ne, double *hpsi) 
   {
      NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[0]));
      // cudaMemcpy(hpsi, dev_mem[ia_hpsi[0]],
      // 2*ne*npack1*sizeof(double),cudaMemcpyDeviceToHost); if (tile_fac==1) {
      //     NWPW_CUBLAS_ERROR( cublasGetMatrix(tile_npack2[0],ne,sizeof(double),
      //                                        dev_mem[ia_hpsi[0]],tile_npack2[0],
      //                                        hpsi,2*npack1));
      // }
   }
 
 
   /******************************/
   /* fft functions (uses cuFFT) */
   /******************************/
 
   /**************************************
    *                                    *
    *          batch_fft_init            *
    *                                    *
    **************************************/
   int batch_fft_init(int nx, int ny, int nz, int nq1, int nq2, int nq3) 
   {
      if (DEBUG_IO) std::cout << "Into batch_fft_init" << std::endl;
      NWPW_CUFFT_ERROR(cufftPlan1d(&forward_plan_x[fftcount], nx, CUFFT_D2Z, nq1));
      NWPW_CUFFT_ERROR(cufftPlan1d(&backward_plan_x[fftcount], nx, CUFFT_Z2D, nq1));

      int x_inembed[] = {nx};
      int x_onembed[] = {nx};
      NWPW_CUFFT_ERROR(cufftPlanMany(&plan_x[fftcount], 1, &nx, x_inembed, 1, nx, x_onembed, 1, nx, CUFFT_Z2Z, nq1));
      int y_inembed[] = {ny};
      int y_onembed[] = {ny};
      NWPW_CUFFT_ERROR(cufftPlanMany(&plan_y[fftcount], 1, &ny, y_inembed, 1, ny, y_onembed, 1, ny, CUFFT_Z2Z, nq2));
      int z_inembed[] = {nz};
      int z_onembed[] = {nz};
      NWPW_CUFFT_ERROR(cufftPlanMany(&plan_z[fftcount], 1, &nz, z_inembed, 1, nz, z_onembed, 1, nz, CUFFT_Z2Z, nq3));
      nxfft[fftcount] = nx;
      nyfft[fftcount] = ny;
      nzfft[fftcount] = nz;
 
      int tag = fftcount;
      ++fftcount;
 
      return tag;
   }



 
   /**************************************
    *                                    *
    *    batch_fft_pipeline_mem_init     *
    *                                    *
    **************************************/
   void batch_fft_pipeline_mem_init(const int nstages, const int n2ft3d) 
   {
      ifft_n = nstages;
 
       // allocate memory and cuda streams
      for (auto i=0; i<ifft_n; ++i)
         ifft_dev[i] = fetch_dev_mem_indx(((size_t)n2ft3d));
   }
 
 
   /**************************************
    *                                    *
    *          batch_fft_end             *
    *                                    *
    **************************************/
   void batch_fft_end(const int tag) 
   {
      // free fft descriptors
      NWPW_CUFFT_ERROR(cufftDestroy(forward_plan_x[tag]));
      NWPW_CUFFT_ERROR(cufftDestroy(plan_x[tag]));
      NWPW_CUFFT_ERROR(cufftDestroy(plan_y[tag]));
      NWPW_CUFFT_ERROR(cufftDestroy(plan_z[tag]));
      NWPW_CUFFT_ERROR(cufftDestroy(backward_plan_x[tag]));
 
      --fftcount;
 
      // free dev_mem
      for (auto i=0; i<ndev_mem; ++i)
      {
         inuse[i] = false;
         NWPW_CUDA_ERROR(cudaFree(dev_mem[i]));
      }
      ndev_mem = 0;
   }

 
   /**************************************
    *                                    *
    *          batch_rfftx               *
    *                                    *
    **************************************/
   void batch_rfftx(const int fft_indx, bool forward, int nx, int nq, int n2ft3d, double *a) 
   {
      int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
      NWPW_CUDA_ERROR(cudaMemcpy(dev_mem[ia_dev], a, n2ft3d * sizeof(double), cudaMemcpyHostToDevice));
 
      if (forward) {
        NWPW_CUFFT_ERROR(cufftExecD2Z(
            forward_plan_x[fft_indx], reinterpret_cast<cufftDoubleReal *>(dev_mem[ia_dev]),
            reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev])));
      } else {
        NWPW_CUFFT_ERROR(
            cufftExecZ2D(backward_plan_x[fft_indx],
                         reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
                         reinterpret_cast<cufftDoubleReal *>(dev_mem[ia_dev])));
      }
 
      NWPW_CUDA_ERROR(cudaMemcpy(a, dev_mem[ia_dev], n2ft3d * sizeof(double), cudaMemcpyDeviceToHost));
 
      inuse[ia_dev] = false;
   }
 
   /**************************************
    *                                    *
    *          batch_rfftx_stages        *
    *                                    *
    **************************************/
   void batch_rfftx_stages(const int stage, const int fft_indx, bool forward, int nx, int nq, int n2ft3d, double *a, int da)
   {
      //int ia_dev = fetch_dev_mem_indx(((size_t) n2ft3d));
      int ia_dev = ifft_dev[da];
 
      if (stage==0)
      {
         inuse[ia_dev] = true;
         NWPW_CUDA_ERROR(cudaMemcpyAsync(dev_mem[ia_dev],a,n2ft3d*sizeof(double),cudaMemcpyHostToDevice,stream[da]));
      }
      else if (stage==1)
      {
         //NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[da]));
         if (forward)
         {
           NWPW_CUFFT_ERROR(cufftExecD2Z(forward_plan_x[fft_indx],
                            reinterpret_cast<cufftDoubleReal *> (dev_mem[ia_dev]),
                            reinterpret_cast<cufftDoubleComplex *> (dev_mem[ia_dev])));
         }
         else
         {
            NWPW_CUFFT_ERROR(cufftExecZ2D(backward_plan_x[fft_indx],
                             reinterpret_cast<cufftDoubleComplex *> (dev_mem[ia_dev]),
                             reinterpret_cast<cufftDoubleReal *> (dev_mem[ia_dev])));
         }
         NWPW_CUDA_ERROR(cudaMemcpyAsync(a,dev_mem[ia_dev],n2ft3d*sizeof(double),cudaMemcpyDeviceToHost,stream[da]));
      }
      else if (stage==2)
      {
         NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[da]));
         inuse[ia_dev] = false;
      }
   }



   /**************************************
    *                                    *
    *          batch_cfftx               *
    *                                    *
    **************************************/
   void batch_cfftx(const int fft_indx, bool forward, int nx, int nq, int n2ft3d, double *a)
   {
      int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
      NWPW_CUDA_ERROR(cudaMemcpy(dev_mem[ia_dev], a, n2ft3d * sizeof(double), cudaMemcpyHostToDevice));

      if (forward) {
        NWPW_CUFFT_ERROR(cufftExecZ2Z(
            plan_x[fft_indx], reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
            reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
            CUFFT_FORWARD));
      } else {
        NWPW_CUFFT_ERROR(cufftExecZ2Z(
            plan_x[fft_indx], reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
            reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
            CUFFT_INVERSE));
      }

      NWPW_CUDA_ERROR(cudaMemcpy(a, dev_mem[ia_dev], n2ft3d * sizeof(double), cudaMemcpyDeviceToHost));

      inuse[ia_dev] = false;
   }



   /**************************************
    *                                    *
    *          batch_cfftx_stages        *
    *                                    *
    **************************************/
   void batch_cfftx_stages(const int stage, const int fft_indx, bool forward, int nx, int nq, int n2ft3d, double *a, int da)
   {
      //int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
      int ia_dev = ifft_dev[da];
      if (stage==0)
      {
         inuse[ia_dev] = true;
         NWPW_CUDA_ERROR(cudaMemcpyAsync(dev_mem[ia_dev],a,n2ft3d*sizeof(double),cudaMemcpyHostToDevice,stream[da]));
      }
      else if (stage==1)
      {
         //NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[da]));
         if (forward) {
           NWPW_CUFFT_ERROR(cufftExecZ2Z(plan_x[fft_indx],
               reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
               reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
               CUFFT_FORWARD));
         } else {
           NWPW_CUFFT_ERROR(cufftExecZ2Z(plan_x[fft_indx],
               reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
               reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
               CUFFT_INVERSE));
         }
         NWPW_CUDA_ERROR(cudaMemcpyAsync(a,dev_mem[ia_dev],n2ft3d*sizeof(double),cudaMemcpyDeviceToHost,stream[da]));
      }
      else if (stage==2)
      {
         NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[da]));
         inuse[ia_dev] = false;
      }
   }





 
 
   /**************************************
    *                                    *
    *          batch_cffty               *
    *                                    *
    **************************************/
   void batch_cffty(const int fft_indx, bool forward, int ny, int nq, int n2ft3d, double *a) 
   {
      int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
      NWPW_CUDA_ERROR(cudaMemcpy(dev_mem[ia_dev], a, n2ft3d * sizeof(double), cudaMemcpyHostToDevice));
 
      if (forward) {
        NWPW_CUFFT_ERROR(cufftExecZ2Z(
            plan_y[fft_indx], reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
            reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
            CUFFT_FORWARD));
      } else {
        NWPW_CUFFT_ERROR(cufftExecZ2Z(
            plan_y[fft_indx], reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
            reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
            CUFFT_INVERSE));
      }
 
      NWPW_CUDA_ERROR(cudaMemcpy(a, dev_mem[ia_dev], n2ft3d * sizeof(double), cudaMemcpyDeviceToHost));
 
      inuse[ia_dev] = false;
   }
 
 
   /**************************************
    *                                    *
    *          batch_cffty_stages        *
    *                                    *
    **************************************/
   void batch_cffty_stages(const int stage, const int fft_indx, bool forward, int ny, int nq, int n2ft3d, double *a, int da)
   {
      //int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
      int ia_dev = ifft_dev[da];
      if (stage==0)
      {
         inuse[ia_dev] = true;
         NWPW_CUDA_ERROR(cudaMemcpyAsync(dev_mem[ia_dev],a,n2ft3d*sizeof(double),cudaMemcpyHostToDevice,stream[da]));
      }
      else if (stage==1)
      {
         //NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[da]));
         if (forward) {
           NWPW_CUFFT_ERROR(cufftExecZ2Z(plan_y[fft_indx],
               reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
               reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
               CUFFT_FORWARD));
         } else {
           NWPW_CUFFT_ERROR(cufftExecZ2Z(plan_y[fft_indx],
               reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
               reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
               CUFFT_INVERSE));
         }
         NWPW_CUDA_ERROR(cudaMemcpyAsync(a,dev_mem[ia_dev],n2ft3d*sizeof(double),cudaMemcpyDeviceToHost,stream[da]));
      }
      else if (stage==2)
      {
         NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[da]));
         inuse[ia_dev] = false;
      }
   }
 
 
   /**************************************
    *                                    *
    *          batch_cfftz               *
    *                                    *
    **************************************/
   void batch_cfftz(const int fft_indx, bool forward, int nz, int nq, int n2ft3d, double *a) 
   {
      int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
      cudaMemcpy(dev_mem[ia_dev], a, n2ft3d * sizeof(double), cudaMemcpyHostToDevice);
 
      if (forward) {
        NWPW_CUFFT_ERROR(cufftExecZ2Z(
            plan_z[fft_indx], reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
            reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
            CUFFT_FORWARD));
      } else {
        NWPW_CUFFT_ERROR(cufftExecZ2Z(
            plan_z[fft_indx], reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
            reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
            CUFFT_INVERSE));
      }
 
      cudaMemcpy(a, dev_mem[ia_dev], n2ft3d * sizeof(double), cudaMemcpyDeviceToHost);
 
      inuse[ia_dev] = false;
   }
 
 
   /**************************************
    *                                    *
    *          batch_cfftz_stages        *
    *                                    *
    **************************************/
   void batch_cfftz_stages(const int stage,const int fft_indx,bool forward,int nz,int nq,int n2ft3d,double *a,int da)
   {
      //int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
      int ia_dev = ifft_dev[da];
 
      if (stage==0)
      {
         inuse[ia_dev] = true;
         NWPW_CUDA_ERROR(cudaMemcpyAsync(dev_mem[ia_dev],a,n2ft3d*sizeof(double),cudaMemcpyHostToDevice,stream[da]));
      }
      else if (stage==1)
      {
         //NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[da]));
         if (forward) {
           NWPW_CUFFT_ERROR(cufftExecZ2Z(
               plan_z[fft_indx], reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
               reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
               CUFFT_FORWARD));
         } else {
           NWPW_CUFFT_ERROR(cufftExecZ2Z(
               plan_z[fft_indx], reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
               reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
               CUFFT_INVERSE));
         }
         NWPW_CUDA_ERROR(cudaMemcpyAsync(a,dev_mem[ia_dev],n2ft3d*sizeof(double),cudaMemcpyDeviceToHost,stream[da]));
      }
      else if (stage==2)
      {
         NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[da]));
         inuse[ia_dev] = false;
      }
   }
 
 
 
   // routines below need to be made into sycl or removed
 
   /**************************************
    *                                    *
    *          eigsrt_device             *
    *                                    *
    **************************************/
   static void eigsrt_device(double *D, double *V, int n) 
   {
      int i, j, k;
      double p;
     
      for (i=0; i<(n-1); ++i) 
      {
         k = i;
         p = D[i];
         for (j = i + 1; j < n; ++j)
           if (D[j] >= p) 
           {
              k = j;
              p = D[j];
           }
        
         if (k != i) 
         {
            D[k] = D[i];
            D[i] = p;
            for (j = 0; j < n; ++j) 
            {
               p = V[j + i * n];
               V[j + i * n] = V[j + k * n];
               V[j + k * n] = p;
            }
         }
      }
   }
 
   /*
      void NN_eigensolver(int ispin, int ne[], double *host_hml, double
      *host_eig) { int i_a1[ispin],i_w1[ispin],info[ispin]; int shift1 = 0; int
      shift2 = 0; for (auto ms=0; ms<ispin; ++ms)
         {
            int nn = ne[ms]*ne[ms];
            i_a1[ms] = fetch_dev_mem_indx(((size_t) ne[ms]) * ((size_t) ne[ms]));
      //input-output i_w1[ms] = fetch_dev_mem_indx(((size_t) ne[ms]) ); //output
            NWPW_CUDA_ERROR(cudaMemcpyAsync(dev_mem[i_a1[ms]],host_hml+shift2,nn*sizeof(double),cudaMemcpyHostToDevice,stream[ms]));
            //NWPW_CUDA_ERROR(cudaMemcpy(dev_mem[i_a1[ms]],host_hml+shift2,nn*sizeof(double),cudaMemcpyHostToDevice));
            //NWPW_CUDA_ERROR(cudaMemcpyAsync(dev_mem[i_w1[ms]],host_eig+shift1,ne[ms]*sizeof(double),cudaMemcpyHostToDevice));
 
            shift1 += ne[0];
            shift2 += ne[0]*ne[0];
         }
 
         // allocate work space for syevd
         if (lwork==0)
         {
            // query working space of syevd
            CUSOLVER_CHECK(cusolverDnDsyevd_bufferSize(cusolverH,jobz,uplo,ne[0],dev_mem[i_a1[0]],ne[0],dev_mem[i_w1[0]],&lwork));
         }
         NWPW_CUDA_ERROR(cudaMalloc(reinterpret_cast<void
      **>(&d_work),sizeof(double) * lwork));
 
         shift1 = 0;
         shift2 = 0;
         for (auto ms=0; ms<ispin; ++ms)
         {
            NWPW_CUDA_ERROR( cudaStreamSynchronize(stream[ms]) );
 
            // compute spectrum
            int nn = ne[ms]*ne[ms];
            CUSOLVER_CHECK(cusolverDnDsyevd(cusolverH,jobz,uplo,ne[ms],dev_mem[i_a1[ms]],ne[ms],dev_mem[i_w1[ms]],d_work,lwork,d_info[ms]));
 
            NWPW_CUDA_ERROR( cudaStreamSynchronize(solverstream) );
 
            NWPW_CUDA_ERROR(cudaMemcpyAsync(host_hml+shift2,dev_mem[i_a1[ms]],nn*sizeof(double),cudaMemcpyDeviceToHost,stream[ms]));
            NWPW_CUDA_ERROR(cudaMemcpyAsync(host_eig+shift1,dev_mem[i_w1[ms]],ne[ms]*sizeof(double),cudaMemcpyDeviceToHost,stream[ms]));
            NWPW_CUDA_ERROR(cudaMemcpyAsync(info+ms,d_info[ms],sizeof(int),cudaMemcpyDeviceToHost,stream[ms]));
 
            shift1 += ne[0];
            shift2 += ne[0]*ne[0];
         }
 
         shift1 = 0;
         shift2 = 0;
         for (auto ms=0; ms<ispin; ++ms)
         {
            NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[ms]));
            NWPW_CUDA_ERROR( cudaFree(d_work) );
            eigsrt_device(host_eig+shift1,host_hml+shift2,ne[ms]);
            shift1 += ne[0];
            shift2 += ne[0]*ne[0];
         }
         for (auto ms=0; ms<ispin; ++ms)
         {
            inuse[i_a1[ms]] = false;
            inuse[i_w1[ms]] = false;
         }
 
         for (auto ms=0; ms<ispin; ++ms)
         {
            inuse[i_a1[ms]];
            inuse[i_w1[ms]];
         }
      }
   */
 
   /**************************************
    *                                    *
    *          NN_eigensolver            *
    *                                    *
    **************************************/
   void NN_eigensolver(int ispin, int ne[], double *host_hml, double *host_eig) 
   {
      int n, ierr;
      int nn = ne[0] * ne[0] + 14;
      double xmp1[nn];
      // double *xmp1 = new (std::nothrow) double[nn]();
     
      int shift1 = 0;
      int shift2 = 0;
      for (int ms = 0; ms < ispin; ++ms) 
      {
         n = ne[ms];
        
         // eigen_(&n,&n,&hml[shift2],&eig[shift1],xmp1,&ierr);
         //  d3db::parall->Barrier();
         EIGEN_PWDFT(n, host_hml + shift2, host_eig + shift1, xmp1, nn, ierr);
         // if (ierr != 0) throw std::runtime_error(std::string("NWPW Error:
         // EIGEN_PWDFT failed!"));
        
         eigsrt_device(host_eig + shift1, host_hml + shift2, n);
         shift1 += ne[0];
         shift2 += ne[0] * ne[0];
      }
   }

   ////////////////////////// special complex-complex fft ////////////////////////////

typedef std::complex<double> complex_t;
         
   // Define a constants for the radix values
static constexpr int radix_values[] = {17, 16, 11, 9, 8, 7, 6, 5, 4, 3, 2};
//static constexpr int radix_values[] = {17, 11, 9, 8, 7, 6, 5, 4, 3, 2};
//static constexpr int radix_values[] = {17, 11, 9, 8, 7, 6, 5, 4, 3, 2};
            
   /**************************************
    *                                    *
    *             fft_radix              *
    *                                    *
    **************************************/
   inline static void fft_radix(const int n, const int s, bool eo, const int radix, const complex_t* twiddle, complex_t* x, complex_t *y)
   {
      if (n == 1) 
      { 
         // Use std::copy to efficiently copy elements from x to y
         if (eo) std::copy(x, x + s, y); 
         //if (eo) std::memcpy(y, x, s * sizeof(double));
         return;
      }

      const int m = n/radix;
      complex_t Atwiddle[radix*radix];
 
      // Precompute twiddle factors for matrix multiplication
      for (int r2=0; r2<radix; ++r2)
      for (int r1=0; r1<radix; ++r1)
         Atwiddle[r1+r2*radix] = twiddle[(r1*r2)%radix];
 
      for (int p=0; p<m; ++p) 
      {
         for (int q=0; q<s; ++q) 
         {
            complex_t* x1 = x+q+s*p;
            complex_t* y1 = y+q+s*radix*p;
 
            // Initialize the output y1 vector
            for (int r1=0; r1<radix; ++r1) 
               y1[s * r1] = 0.0;
 
            // Matrix-vector multiplication
            for (int r2 = 0; r2 < radix; ++r2)
            for (int r1 = 0; r1 < radix; ++r1)
               y1[s*r1] += x1[s*m*r2] * Atwiddle[r1 + r2*radix];
               //y1[s*r1] += x1[s*m*r2] * twiddle[(r1*r2)%radix];
 
            // Apply phase factor to each result
            for (int r1 = 0; r1 < radix; ++r1)
            {
               //complex_t fac = std::pow(twiddle[radix + r1], p);
               complex_t fac = complex_t(1.0, 0.0);
               for (int ps = 0; ps < p; ++ps)
                  fac *= twiddle[radix + r1];
 
               y1[s * r1] *= fac;
            }
         }
      }
   }

   /**************************************
    *                                    *
    *            fft_twiddle             *
    *                                    *
    **************************************/
   inline static void fft_twiddle(const int n, const complex_t* twiddle, complex_t* x) // Fourier transform
   {
      //complex_t* y = new complex_t[n];
      complex_t y[n];
      int eo = 0;
      int s  = 1;
      int nn = n;
      int nsize = 0;
      while (s<=n) 
      {
         // Identify the largest radix applicable for current nn
         int radix = 2;  // Default to radix-2
         for (int r : radix_values) {
            if (nn % r == 0) {
               radix = r;
               break;
            }
         }
 
         // Perform FFT with the determined radix
         if (eo)
            fft_radix(nn, s, eo, radix, twiddle + nsize, y, x);
         else
            fft_radix(nn, s, eo, radix, twiddle + nsize, x, y);

         nsize += 2*radix;
         nn /= radix;
         s *= radix;
         eo = !eo;  // Toggle the 'even-odd' flag
      }
   }


   /**************************************
    *                                    *
    *         set_sub_fft_twiddle        *
    *                                    *
    **************************************/
   static void set_sub_fft_twiddle(const int isgn, const int n, const int radix, complex_t *twiddle)
   {
      const double theta0 = 2*M_PI/((double) n);
      const double theta_radix = 2*M_PI/((double) radix);
 
      // Calculate radix-specific twiddle factors
      for (int r=0; r<radix; ++r)
         twiddle[r] = complex_t(cos(r*theta_radix), isgn*sin(r*theta_radix));
 
      // Calculate the main twiddle factors for the FFT
      for (int r=0; r<radix; ++r)
         twiddle[radix+r] = complex_t(cos(r*theta0), isgn*sin(r*theta0));
   }



   /**************************************
    *                                    *
    *            set_fft_twiddle         *
    *                                    *
    **************************************/
   void set_fft_twiddle(const int isgn, const int n, double *twiddle) 
   {
      complex_t* complex_twiddle = reinterpret_cast<complex_t*>(twiddle);
      int nsize = 0;
      int s = 1;
      int nn = n;
 
      while (s <= n) 
      {
         bool found = false;
         // Loop through possible radix values to find the largest factor of nn
         for (int radix : radix_values) 
         {
            if (nn % radix == 0) 
            {
               set_sub_fft_twiddle(isgn, nn, radix, complex_twiddle + nsize);
               nsize += 2*radix;
               nn /= radix;
               s *= radix;
               found = true;
               break;
            }
         }
         if (!found) break;
      }
   }

   /**************************************
    *                                    *
    *            size_fft_twiddle        *
    *                                    *
    **************************************/
   int size_fft_twiddle(const int n) 
   {
      int nsize = 0;
      int s = 1;
      int nn = n;
 
      while (s <= n) 
      {
         bool found = false;
         // Loop through possible radix values to find the largest factor of nn
         for (int radix : radix_values) 
         {
            if (nn % radix == 0) 
            {
               nsize += 2 * radix;
               nn /= radix;
               s *= radix;
               found = true;
               break;
            }
         }
         if (!found) break;
      }
 
      return nsize;
   }

   /**************************************
    *                                    *
    *              batch_cfft            *
    *                                    *
    **************************************/
   void batch_cfft(const bool forward, const int nz, const int nq, const int nfft3d, double *a, const double *twiddle, const double *tmpz)
   {
      // Ensure the function processes the right type of data
      // If twiddle is indeed of type complex_t, this cast is necessary
      //complex_t*       complex_a       = reinterpret_cast<complex_t*>(a);
      //const complex_t* complex_twiddle = reinterpret_cast<const complex_t*>(twiddle);

      //int shift = nfft3d;
      int shift = 2*nfft3d;
      int indx = 0;

      // Process each FFT batch
      for (auto q=0; q<nq; ++q)
      {
         (forward ? dcfftf_(&nz, a + indx, tmpz) : dcfftb_(&nz, a + indx, tmpz));
         //fft_twiddle(nz, complex_twiddle, complex_a+indx);

         indx += (shift);
      }
   }

   /**************************************
    *                                    *
    *          batch_cfft_zero           *
    *                                    *
    **************************************/
   void batch_cfft_zero(const bool forward, int nz, int nq, int nfft3d, double *a, const double *twiddle, const double *tmpz, const bool *zero) 
   {
      // Ensure the function processes the right type of data
      // If twiddle is indeed of type complex_t, this cast is necessary
      //complex_t*       complex_a       = reinterpret_cast<complex_t*>(a);
      //const complex_t* complex_twiddle = reinterpret_cast<const complex_t*>(twiddle);

      //int shift = nfft3d;
      int shift = 2*nfft3d;
      int indx = 0;

      // Process each FFT batch
      for (auto q=0; q<nq; ++q)
      {
         if (!zero[q])
            (forward ? dcfftf_(&nz, a + indx, tmpz) : dcfftb_(&nz, a + indx, tmpz));
            //fft_twiddle(nz, complex_twiddle, complex_a+indx);

         indx += (shift);
      }
   }


   ////////////////////////// special complex-complex fft ////////////////////////////



};

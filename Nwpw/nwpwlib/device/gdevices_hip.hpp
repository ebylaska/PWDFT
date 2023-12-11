// NWPW_HIP Routines

#pragma once

#include "blas.h"

// #include        "gdevice.hpp"

#include <hip/hip_runtime.h>
#include <rocblas.h>
#include <rocfft.h>

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

class hip_exception : public std::exception {

  std::string file_;
  int line_;
  hipError_t err_code_;

  const char *what() const noexcept override {
    std::stringstream ss;
    ss << "HIP Exception, " << hipGetErrorString(err_code_) << " at "
       << std::endl
       << file_ << " : " << line_ << std::endl;
    auto msg = ss.str();
    return strdup(msg.c_str());
  }

public:
  hip_exception(std::string file, int line, hipError_t err)
      : file_(file), line_(line), err_code_(err) {}
};

class rocfft_exception : public std::exception {

  std::string file_;
  int line_;
  rocfft_status err_code_;

  static const char *_rocfftGetErrorEnum(rocfft_status error) {
    switch (error) {
    case rocfft_status_success:
      return "ROCFFT_STATUS_SUCCESS";

    case rocfft_status_failure:
      return "ROCFFT_STATUS_FAILURE";

    case rocfft_status_invalid_arg_value:
      return "ROCFFT_STATUS_INVALID_ARG_VALUE";

    case rocfft_status_invalid_dimensions:
      return "ROCFFT_STATUS_INVALID_DIMENSIONS";

    case rocfft_status_invalid_array_type:
      return "ROCFFT_STATUS_INVALID_ARRAY_TYPE";

    case rocfft_status_invalid_strides:
      return "ROCFFT_STATUS_INVALID_STRIDES";

    case rocfft_status_invalid_distance:
      return "ROCFFT_STATUS_INVALID_DISTANCE";

    case rocfft_status_invalid_offset:
      return "ROCFFT_STATUS_INVALID_OFFSET";

    case rocfft_status_invalid_work_buffer:
      return "ROCFFT_STATUS_INVALID_WORK_BUFFER";
    }
    return "<unknown>";
  }

  const char *what() const noexcept override {
    std::stringstream ss;
    ss << "ROCFFT Exception, "
       << " Error Code: " << _rocfftGetErrorEnum(err_code_) << std::endl
       << " at " << file_ << " : " << line_ << std::endl;
    auto msg = ss.str();
    return strdup(msg.c_str());
  }

public:
  rocfft_exception(std::string file, int line, rocfft_status err)
      : file_(file), line_(line), err_code_(err) {}
};

class rocblas_exception : public std::exception {

  std::string file_;
  int line_;
  rocblas_status err_code_;

  const char *what() const noexcept override {
    std::stringstream ss;
    ss << "rocBLAS Exception, "
       << " Error Code: " << rocblas_status_to_string(err_code_) << std::endl
       << " at " << file_ << " : " << line_ << std::endl;
    auto msg = ss.str();
    return strdup(msg.c_str());
  }

public:
  rocblas_exception(std::string file, int line, rocblas_status err)
      : file_(file), line_(line), err_code_(err) {}
};

#define NWPW_HIP_ERROR(CALL)                                                   \
  do {                                                                         \
    hipError_t err = CALL;                                                     \
    if (err != hipSuccess) {                                                   \
      throw hip_exception(__FILE__, __LINE__, err);                            \
    }                                                                          \
  } while (0)

#define NWPW_ROCFFT_ERROR(CALL)                                                \
  do {                                                                         \
    rocfft_status err = CALL;                                                  \
    if (err != rocfft_status_success) {                                        \
      throw rocfft_exception(__FILE__, __LINE__, err);                         \
    }                                                                          \
  } while (0)

#define NWPW_ROCBLAS_ERROR(CALL)                                               \
  do {                                                                         \
    rocblas_status err = CALL;                                                 \
    if (err != rocblas_status_success) {                                       \
      throw rocblas_exception(__FILE__, __LINE__, err);                        \
    }                                                                          \
  } while (0)

/* Gdevices (HIP) object -

   Operations:
      blas functions (using rocblas):
      TN3_dgemm - computes <A|A>, <B|B>, and <A|B> overlaps
      TN1_dgemm - computes  <B|B> overlap

      TN_dgemm
      NT_dgemm
      psi_alloc
      psi_dealloc

      fft functions (uses rocFFT)


*/

class Gdevices {

  int fftcount = 0;
  int nxfft[2], nyfft[2], nzfft[2];
  rocfft_plan forward_plan_x[2],  forward_plan_y[2],  forward_plan_z[2], forward_plan_rx[2];
  rocfft_plan backward_plan_x[2], backward_plan_y[2], backward_plan_z[2], backward_plan_rx[2];;
  int ifft_dev[15];
  int ifft_n;

  rocblas_handle master_handle = 0;
  rocblas_operation matT = rocblas_operation_transpose;
  rocblas_operation matN = rocblas_operation_none;

  hipStream_t stream[12];

public:
  bool hasgpu = true;

  /* device memory */
  int ndev_mem = 0;
  bool inuse[29] = {false};
  size_t ndsize_mem[29];
  double *dev_mem[29];
  int tile_fac = 1;
  int tile_npack2_max;
  int tile_npack2[19], tile_start2[19];
  double *a_psi, *a_hpsi, *b_prj;
  int ia_psi[2], ia_hpsi[2], ib_prj[2];

  /* constructor */
  Gdevices() {
    ndev_mem = 0;

    NWPW_ROCBLAS_ERROR(rocblas_create_handle(&master_handle));

    // allocate hip streams
    for (auto i = 0; i < 12; ++i)
      NWPW_HIP_ERROR(hipStreamCreate(&stream[i]));

    NWPW_ROCFFT_ERROR(rocfft_setup());
  }

  /* deconstructor */
  ~Gdevices() {
    // free dev_mem
    for (auto i = 0; i < ndev_mem; ++i) {
      inuse[i] = false;
      hipFree(dev_mem[i]);
    }
    ndev_mem = 0;

    // free hip streams
    for (auto i = 0; i < 12; ++i)
      hipStreamDestroy(stream[i]);

    rocblas_destroy_handle(master_handle);

    // free fft descriptors
    rocfft_cleanup();
  }

  int fetch_dev_mem_indx(const size_t ndsize) {
    int ii = 0;
    while ((((ndsize != ndsize_mem[ii]) || inuse[ii])) && (ii < ndev_mem))
      ++ii;

    if (ii < ndev_mem) {
      inuse[ii] = true;
    } else {
      ii = ndev_mem;
      inuse[ii] = true;
      ndsize_mem[ii] = ndsize;
      NWPW_HIP_ERROR(
          hipMalloc((void **)&(dev_mem[ii]), ndsize * sizeof(double)));
      ndev_mem += 1;
    }

    NWPW_HIP_ERROR(hipMemset(dev_mem[ii], 0, ndsize * sizeof(double)));
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
      NWPW_HIP_ERROR(hipMemcpy(dev_mem[ic11], host_caa,
                               ne * ne * sizeof(double),
                               hipMemcpyHostToDevice));
      NWPW_HIP_ERROR(hipMemcpy(dev_mem[ic12], host_cab,
                               ne * ne * sizeof(double),
                               hipMemcpyHostToDevice));
      NWPW_HIP_ERROR(hipMemcpy(dev_mem[ic21], host_cba,
                               ne * ne * sizeof(double),
                               hipMemcpyHostToDevice));
      NWPW_HIP_ERROR(hipMemcpy(dev_mem[ic22], host_cbb,
                               ne * ne * sizeof(double),
                               hipMemcpyHostToDevice));
    }

    // copy host_a,host_b --> dev_mem
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        tile_npack2[0], ne, sizeof(double), &host_a[tile_start2[0]], npack2,
        dev_mem[ia_psi[0]], tile_npack2[0], stream[0]));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        tile_npack2[0], ne, sizeof(double), &host_b[tile_start2[0]], npack2,
        dev_mem[ia_hpsi[0]], tile_npack2[0], stream[0]));

    double beta0 = beta;
    for (auto tt = 0; tt < tile_fac; ++tt) {
      int ttp1 = tt + 1;
      if (ttp1 < tile_fac) {
        NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
            tile_npack2[ttp1], ne, sizeof(double), &host_a[tile_start2[ttp1]],
            npack2, dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1],
            stream[ttp1 % 2]));
        NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
            tile_npack2[ttp1], ne, sizeof(double), &host_b[tile_start2[ttp1]],
            npack2, dev_mem[ia_hpsi[ttp1 % 2]], tile_npack2[ttp1],
            stream[ttp1 % 2]));
      }
      NWPW_HIP_ERROR(hipStreamSynchronize(stream[tt % 2]));
      NWPW_ROCBLAS_ERROR(rocblas_dgemm(
          master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
          dev_mem[ia_psi[tt % 2]], tile_npack2[tt], dev_mem[ia_psi[tt % 2]],
          tile_npack2[tt], &beta0, dev_mem[ic11], ne));
      NWPW_ROCBLAS_ERROR(rocblas_dgemm(
          master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
          dev_mem[ia_psi[tt % 2]], tile_npack2[tt], dev_mem[ia_hpsi[tt % 2]],
          tile_npack2[tt], &beta0, dev_mem[ic12], ne));
      NWPW_ROCBLAS_ERROR(rocblas_dgemm(
          master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
          dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], dev_mem[ia_psi[tt % 2]],
          tile_npack2[tt], &beta0, dev_mem[ic21], ne));
      NWPW_ROCBLAS_ERROR(rocblas_dgemm(
          master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
          dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], dev_mem[ia_hpsi[tt % 2]],
          tile_npack2[tt], &beta0, dev_mem[ic22], ne));
      beta0 = 1.0;
    }

    NWPW_HIP_ERROR(hipMemcpy(host_caa, dev_mem[ic11], ne * ne * sizeof(double),
                             hipMemcpyDeviceToHost));
    NWPW_HIP_ERROR(hipMemcpy(host_cab, dev_mem[ic12], ne * ne * sizeof(double),
                             hipMemcpyDeviceToHost));
    NWPW_HIP_ERROR(hipMemcpy(host_cba, dev_mem[ic21], ne * ne * sizeof(double),
                             hipMemcpyDeviceToHost));
    NWPW_HIP_ERROR(hipMemcpy(host_cbb, dev_mem[ic22], ne * ne * sizeof(double),
                             hipMemcpyDeviceToHost));

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
                 double *host_cab, double *host_cbb) {
    int ic11 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
    int ic12 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));
    int ic22 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));

    if (std::fabs(beta) > 0.0) {
      NWPW_HIP_ERROR(hipMemcpy(dev_mem[ic11], host_caa,
                               ne * ne * sizeof(double),
                               hipMemcpyHostToDevice));
      NWPW_HIP_ERROR(hipMemcpy(dev_mem[ic12], host_cab,
                               ne * ne * sizeof(double),
                               hipMemcpyHostToDevice));
      NWPW_HIP_ERROR(hipMemcpy(dev_mem[ic22], host_cbb,
                               ne * ne * sizeof(double),
                               hipMemcpyHostToDevice));
    }

    // copy host_a,host_b --> dev_mem
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        tile_npack2[0], ne, sizeof(double), &host_a[tile_start2[0]], npack2,
        dev_mem[ia_psi[0]], tile_npack2[0], stream[0]));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        tile_npack2[0], ne, sizeof(double), &host_b[tile_start2[0]], npack2,
        dev_mem[ia_hpsi[0]], tile_npack2[0], stream[0]));

    double beta0 = beta;
    for (auto tt = 0; tt < tile_fac; ++tt) {
      int ttp1 = tt + 1;
      if (ttp1 < tile_fac) {
        NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
            tile_npack2[ttp1], ne, sizeof(double), &host_a[tile_start2[ttp1]],
            npack2, dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1],
            stream[ttp1 % 2]));
        NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
            tile_npack2[ttp1], ne, sizeof(double), &host_b[tile_start2[ttp1]],
            npack2, dev_mem[ia_hpsi[ttp1 % 2]], tile_npack2[ttp1],
            stream[ttp1 % 2]));
      }
      NWPW_HIP_ERROR(hipStreamSynchronize(stream[tt % 2]));
      NWPW_ROCBLAS_ERROR(rocblas_dgemm(
          master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
          dev_mem[ia_psi[tt % 2]], tile_npack2[tt], dev_mem[ia_psi[tt % 2]],
          tile_npack2[tt], &beta0, dev_mem[ic11], ne));
      NWPW_ROCBLAS_ERROR(rocblas_dgemm(
          master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
          dev_mem[ia_psi[tt % 2]], tile_npack2[tt], dev_mem[ia_hpsi[tt % 2]],
          tile_npack2[tt], &beta0, dev_mem[ic12], ne));
      NWPW_ROCBLAS_ERROR(rocblas_dgemm(
          master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
          dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], dev_mem[ia_hpsi[tt % 2]],
          tile_npack2[tt], &beta0, dev_mem[ic22], ne));
      beta0 = 1.0;
    }

    NWPW_HIP_ERROR(hipMemcpy(host_caa, dev_mem[ic11], ne * ne * sizeof(double),
                             hipMemcpyDeviceToHost));
    NWPW_HIP_ERROR(hipMemcpy(host_cab, dev_mem[ic12], ne * ne * sizeof(double),
                             hipMemcpyDeviceToHost));
    NWPW_HIP_ERROR(hipMemcpy(host_cbb, dev_mem[ic22], ne * ne * sizeof(double),
                             hipMemcpyDeviceToHost));

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
                 double *host_b, double beta, double *host_cab) {
    int ic12 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));

    if (std::fabs(beta) > 0.0) {
      NWPW_HIP_ERROR(hipMemcpy(dev_mem[ic12], host_cab,
                               ne * ne * sizeof(double),
                               hipMemcpyHostToDevice));
    }

    // copy host_a,host_b --> dev_mem
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        tile_npack2[0], ne, sizeof(double), &host_a[tile_start2[0]], npack2,
        dev_mem[ia_psi[0]], tile_npack2[0], stream[0]));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        tile_npack2[0], ne, sizeof(double), &host_b[tile_start2[0]], npack2,
        dev_mem[ia_hpsi[0]], tile_npack2[0], stream[0]));

    double beta0 = beta;
    for (auto tt = 0; tt < tile_fac; ++tt) {
      int ttp1 = tt + 1;
      if (ttp1 < tile_fac) {
        NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
            tile_npack2[ttp1], ne, sizeof(double), &host_a[tile_start2[ttp1]],
            npack2, dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1],
            stream[ttp1 % 2]));
        NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
            tile_npack2[ttp1], ne, sizeof(double), &host_b[tile_start2[ttp1]],
            npack2, dev_mem[ia_hpsi[ttp1 % 2]], tile_npack2[ttp1],
            stream[ttp1 % 2]));
      }
      NWPW_HIP_ERROR(hipStreamSynchronize(stream[tt % 2]));
      NWPW_ROCBLAS_ERROR(rocblas_dgemm(
          master_handle, matT, matN, ne, ne, tile_npack2[tt], &alpha,
          dev_mem[ia_psi[tt % 2]], tile_npack2[tt], dev_mem[ia_hpsi[tt % 2]],
          tile_npack2[tt], &beta0, dev_mem[ic12], ne));
      beta0 = 1.0;
    }

    NWPW_HIP_ERROR(hipMemcpy(host_cab, dev_mem[ic12], ne * ne * sizeof(double),
                             hipMemcpyDeviceToHost));

    inuse[ic12] = false;
  }

  /**************************************
   *                                    *
   *              NN_dgemm              *
   *                                    *
   **************************************/
  void NN_dgemm(int npack2, int ne, double alpha, double *host_a,
                double *host_b, double beta, double *host_c) {
    // DGEMM_PWDFT((char *) "N",(char *)
    // "N",npack2,ne,ne,alpha,host_a,npack2,host_b,ne,beta,host_c,npack2);

    int ib = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne));

    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        ne, ne, sizeof(double), host_b, ne, dev_mem[ib], ne, stream[0]));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        tile_npack2[0], ne, sizeof(double), &host_a[tile_start2[0]], npack2,
        dev_mem[ia_psi[0]], tile_npack2[0], stream[0]));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        tile_npack2[0], ne, sizeof(double), &host_c[tile_start2[0]], npack2,
        dev_mem[ia_hpsi[0]], tile_npack2[0], stream[0]));

    // double beta0 = beta;
    for (auto tt = 0; tt < tile_fac; ++tt) {
      int ttp1 = tt + 1;
      if (ttp1 < tile_fac) {
        NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
            tile_npack2[ttp1], ne, sizeof(double), &host_a[tile_start2[ttp1]],
            npack2, dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1],
            stream[ttp1 % 2]));
        NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
            tile_npack2[ttp1], ne, sizeof(double), &host_c[tile_start2[ttp1]],
            npack2, dev_mem[ia_hpsi[ttp1 % 2]], tile_npack2[ttp1],
            stream[ttp1 % 2]));
      }
      NWPW_HIP_ERROR(hipStreamSynchronize(stream[tt % 2]));
      NWPW_ROCBLAS_ERROR(rocblas_dgemm(
          master_handle, matN, matN, tile_npack2[tt], ne, ne, &alpha,
          dev_mem[ia_psi[tt % 2]], tile_npack2[tt], dev_mem[ib], ne, &beta,
          dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt]));
      NWPW_ROCBLAS_ERROR(rocblas_get_matrix_async(
          tile_npack2[tt], ne, sizeof(double), dev_mem[ia_hpsi[tt % 2]],
          tile_npack2[tt], &host_c[tile_start2[tt]], npack2, stream[tt % 2]));
    }
    NWPW_HIP_ERROR(hipStreamSynchronize(stream[(tile_fac - 1) % 2]));
    // NWPW_HIP_ERROR(
    // hipMemcpy(host_c,dev_mem[ia_hpsi[0]],npack2*ne*sizeof(double),hipMemcpyDeviceToHost)
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
                 double *host_c,int ldc) {
    int ia = fetch_dev_mem_indx(((size_t)lda) * ((size_t)k));
    int ib = fetch_dev_mem_indx(((size_t)ldb) * ((size_t)n));
    int ic = fetch_dev_mem_indx(((size_t)ldc) * ((size_t)n));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(lda,k,sizeof(double),host_a,lda,dev_mem[ia],lda,stream[0]));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(ldb,n,sizeof(double),host_b,ldb,dev_mem[ib],ldb,stream[0]));

    NWPW_HIP_ERROR(hipStreamSynchronize(stream[0]));
    NWPW_ROCBLAS_ERROR(roclas_dgemm(master_handle,matN,matN,m,n,k,
                                    &alpha,
                                    dev_mem[ia],lda,
                                    dev_mem[ib],ldb,
                                    &beta,
                                    dev_mem[ic],ldc);
    NWPW_ROCBLAS_ERROR(rocblas_get_matrix_async(ldc,n,sizeof(double),dev_mem[ic],ldc,host_c,ldc,stream[0]);
    NWPW_HIP_ERROR(hipStreamSynchronize(stream[0]));

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
                 double *host_c,int ldc) {
    int ia = fetch_dev_mem_indx(((size_t)lda) * ((size_t)m));
    int ib = fetch_dev_mem_indx(((size_t)ldb) * ((size_t)n));
    int ic = fetch_dev_mem_indx(((size_t)ldc) * ((size_t)n));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(lda,m,sizeof(double),host_a,lda,dev_mem[ia],lda,stream[0]));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(ldb,n,sizeof(double),host_b,ldb,dev_mem[ib],ldb,stream[0]));
   
    NWPW_HIP_ERROR(hipStreamSynchronize(stream[0]));
    NWPW_ROCBLAS_ERROR(roclas_dgemm(master_handle,matT,matN,m,n,k,
                                    &alpha,
                                    dev_mem[ia],lda,
                                    dev_mem[ib],ldb,
                                    &beta,
                                    dev_mem[ic],ldc);
    NWPW_ROCBLAS_ERROR(rocblas_get_matrix_async(ldc,n,sizeof(double),dev_mem[ic],ldc,host_c,ldc,stream[0]);
    NWPW_HIP_ERROR(hipStreamSynchronize(stream[0]));

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
                 double *host_c,int ldc) {
    int ia = fetch_dev_mem_indx(((size_t)lda) * ((size_t)k));
    int ib = fetch_dev_mem_indx(((size_t)ldb) * ((size_t)k));
    int ic = fetch_dev_mem_indx(((size_t)ldc) * ((size_t)n));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(lda,k,sizeof(double),host_a,lda,dev_mem[ia],lda,stream[0]));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(ldb,k,sizeof(double),host_b,ldb,dev_mem[ib],ldb,stream[0]));

    NWPW_HIP_ERROR(hipStreamSynchronize(stream[0]));
    NWPW_ROCBLAS_ERROR(roclas_dgemm(master_handle,matN,matT,m,n,k,
                                    &alpha,
                                    dev_mem[ia],lda,
                                    dev_mem[ib],ldb,
                                    &beta,
                                    dev_mem[ic],ldc);
    NWPW_ROCBLAS_ERROR(rocblas_get_matrix_async(ldc,n,sizeof(double),dev_mem[ic],ldc,host_c,ldc,stream[0]);
    NWPW_HIP_ERROR(hipStreamSynchronize(stream[0]));

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
                double *host_b, double beta, double *host_c) {
    // DGEMM_PWDFT((char *) "T",(char *)
    // "N",ne,nprj,npack2,alpha,host_a,npack2,host_b,npack2,beta,host_c,ne);

    // gdevice_TN_dgemm(nn,nprj,ng,rtwo,a,b,rzero,sum);

    // int ia = fetch_dev_mem_indx(((size_t) npack2) * ((size_t) ne));
    // int ib = fetch_dev_mem_indx(((size_t) npack2) * ((size_t) nprj));
    b_prj = host_b;
    ib_prj[0] = fetch_dev_mem_indx(((size_t)tile_npack2_max) * ((size_t)nprj));
    if (tile_fac > 1)
      ib_prj[1] =
          fetch_dev_mem_indx(((size_t)tile_npack2_max) * ((size_t)nprj));
    int ic = fetch_dev_mem_indx(((size_t)ne) * ((size_t)nprj));

    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        ne, nprj, sizeof(double), host_c, ne, dev_mem[ic], ne, stream[0]));

    if (tile_fac > 1)
      NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
          tile_npack2[0], ne, sizeof(double), &a_psi[tile_start2[0]], npack2,
          dev_mem[ia_psi[0]], tile_npack2[0], stream[0]));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        tile_npack2[0], nprj, sizeof(double), &b_prj[tile_start2[0]], npack2,
        dev_mem[ib_prj[0]], tile_npack2[0], stream[0]));

    double beta0 = beta;
    for (auto tt = 0; tt < tile_fac; ++tt) {
      int ttp1 = tt + 1;
      if (ttp1 < tile_fac) {
        NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
            tile_npack2[ttp1], ne, sizeof(double), &a_psi[tile_start2[ttp1]],
            npack2, dev_mem[ia_psi[ttp1 % 2]], tile_npack2[ttp1],
            stream[ttp1 % 2]));
        NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
            tile_npack2[ttp1], nprj, sizeof(double), &b_prj[tile_start2[ttp1]],
            npack2, dev_mem[ib_prj[ttp1 % 2]], tile_npack2[ttp1],
            stream[ttp1 % 2]));
      }
      NWPW_HIP_ERROR(hipStreamSynchronize(stream[tt % 2]));
      NWPW_ROCBLAS_ERROR(rocblas_dgemm(
          master_handle, matT, matN, ne, nprj, tile_npack2[tt], &alpha,
          dev_mem[ia_psi[tt % 2]], tile_npack2[tt], dev_mem[ib_prj[tt % 2]],
          tile_npack2[tt], &beta0, dev_mem[ic], ne));
      beta0 = 1.0;
    }
    NWPW_HIP_ERROR(hipMemcpy(host_c, dev_mem[ic], ne * nprj * sizeof(double),
                             hipMemcpyDeviceToHost));

    // inuse[ia] = false;
    // inuse[ib_prj[0]] = false;
    // if (tile_fac>1) inuse[ib_prj[1]] = false;
    inuse[ic] = false;
  }

  void T_free() {
    inuse[ib_prj[0]] = false;
    if (tile_fac > 1)
      inuse[ib_prj[1]] = false;
  }

  void NT_dgemm(int npack2, int ne, int nprj, double alpha, double *host_a,
                double *host_b, double beta, double *host_c) {
    // DGEMM_PWDFT((char *) "N",(char *)
    // "T",npack2,ne,nprj,alpha,host_a,npack2,host_b,ne,beta,host_c,npack2);

    int ib = fetch_dev_mem_indx(((size_t)ne) * ((size_t)nprj));

    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(ne, nprj, sizeof(double),
                                                host_b, ne, dev_mem[ib], ne,
                                                stream[(tile_fac - 1) % 2]));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        tile_npack2[tile_fac - 1], ne, sizeof(double),
        &host_c[tile_start2[tile_fac - 1]], npack2,
        dev_mem[ia_hpsi[(tile_fac - 1) % 2]], tile_npack2[tile_fac - 1],
        stream[(tile_fac - 1) % 2]));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        tile_npack2[tile_fac - 1], nprj, sizeof(double),
        &host_a[tile_start2[tile_fac - 1]], npack2,
        dev_mem[ib_prj[(tile_fac - 1) % 2]], tile_npack2[tile_fac - 1],
        stream[(tile_fac - 1) % 2]));
    for (auto tt = tile_fac - 1; tt >= 0; --tt) {
      int ttm1 = tt - 1;
      if (ttm1 >= 0) {
        NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
            tile_npack2[ttm1], ne, sizeof(double), &host_c[tile_start2[ttm1]],
            npack2, dev_mem[ia_hpsi[ttm1 % 2]], tile_npack2[ttm1],
            stream[ttm1 % 2]));
        NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
            tile_npack2[ttm1], nprj, sizeof(double), &host_a[tile_start2[ttm1]],
            npack2, dev_mem[ib_prj[ttm1 % 2]], tile_npack2[ttm1],
            stream[ttm1 % 2]));
      }
      NWPW_HIP_ERROR(hipStreamSynchronize(stream[tt % 2]));
      NWPW_ROCBLAS_ERROR(rocblas_dgemm(
          master_handle, matN, matT, tile_npack2[tt], ne, nprj, &alpha,
          dev_mem[ib_prj[tt % 2]], tile_npack2[tt], dev_mem[ib], ne, &beta,
          dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt]));
      NWPW_ROCBLAS_ERROR(rocblas_get_matrix_async(
          tile_npack2[tt], ne, sizeof(double), dev_mem[ia_hpsi[tt % 2]],
          tile_npack2[tt], &host_c[tile_start2[tt]], npack2, stream[tt % 2]));
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
                 double *host_sa0, double *host_sa1, double *host_st1) {
    double rzero = 0.0;
    double rone = 1.0;
    int i_s21 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input
    int i_s12 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input
    int i_s11 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input
    int i_sa0 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input
    int i_st1 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // tmp
    int i_sa1 = fetch_dev_mem_indx(((size_t)ne) * ((size_t)ne)); // input-output

    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        ne, ne, sizeof(double), host_s12, ne, dev_mem[i_s12], ne, stream[0]));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        ne, ne, sizeof(double), host_s11, ne, dev_mem[i_s11], ne, stream[1]));

    NWPW_ROCBLAS_ERROR(rocblas_set_matrix(ne, ne, sizeof(double), host_s21, ne,
                                          dev_mem[i_s21], ne));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix(ne, ne, sizeof(double), host_sa0, ne,
                                          dev_mem[i_sa0], ne));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix(ne, ne, sizeof(double), host_s11, ne,
                                          dev_mem[i_s11], ne));
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix(ne, ne, sizeof(double), host_sa1, ne,
                                          dev_mem[i_sa1], ne));

    // mmm_Multiply(ms, s21, sa0, 1.0, sa1, 1.0);
    NWPW_ROCBLAS_ERROR(rocblas_dgemm(master_handle, matN, matN, ne, ne, ne,
                                     &rone, dev_mem[i_s21], ne, dev_mem[i_sa0],
                                     ne, &rone, dev_mem[i_sa1], ne));

    // mmm_Multiply(ms, sa0, s12, 1.0, sa1, 1.0);
    NWPW_HIP_ERROR(hipStreamSynchronize(stream[0]));
    NWPW_ROCBLAS_ERROR(rocblas_dgemm(master_handle, matN, matN, ne, ne, ne,
                                     &rone, dev_mem[i_sa0], ne, dev_mem[i_s12],
                                     ne, &rone, dev_mem[i_sa1], ne));

    // mmm_Multiply(ms, s11, sa0, 1.0, st1, 0.0);
    NWPW_HIP_ERROR(hipStreamSynchronize(stream[1]));
    NWPW_ROCBLAS_ERROR(rocblas_dgemm(master_handle, matN, matN, ne, ne, ne,
                                     &rone, dev_mem[i_s11], ne, dev_mem[i_sa0],
                                     ne, &rzero, dev_mem[i_st1], ne));

    // mmm_Multiply(ms, sa0, st1, 1.0, sa1, 1.0);
    NWPW_ROCBLAS_ERROR(rocblas_dgemm(master_handle, matN, matN, ne, ne, ne,
                                     &rone, dev_mem[i_sa0], ne, dev_mem[i_st1],
                                     ne, &rone, dev_mem[i_sa1], ne));

    NWPW_ROCBLAS_ERROR(rocblas_get_matrix(ne, ne, sizeof(double),
                                          dev_mem[i_sa1], ne, host_sa1, ne));

    inuse[i_s21] = false;
    inuse[i_s12] = false;
    inuse[i_s11] = false;
    inuse[i_sa0] = false;
    inuse[i_st1] = false;
    inuse[i_sa1] = false;
  }

  /********************/
  /* psi_dev functions*/
  /********************/
  void psi_alloc(int npack1, int ne, int tile_fac0 = 1) {
    tile_fac = tile_fac0;

    tile_npack2_max = (((2 * npack1) % tile_fac) == 0)
                          ? (2 * npack1) / tile_fac
                          : (2 * npack1) / tile_fac + 1;
    // for (auto i=0; i<tile_fac; ++i) tile_npack2[i] =
    // (i<((2*npack1)%tile_fac)) ? (2*npack1)/tile_fac+1 : (2*npack1)/tile_fac;
    for (auto i = 0; i < tile_fac; ++i)
      tile_npack2[i] = (2 * npack1) / tile_fac;
    for (auto i = 0; i < ((2 * npack1) % tile_fac); ++i)
      tile_npack2[i] += 1;

    tile_start2[0] = 0;
    for (auto i = 1; i < tile_fac; ++i)
      tile_start2[i] = tile_start2[i - 1] + tile_npack2[i - 1];

    ia_psi[0] = fetch_dev_mem_indx(((size_t)tile_npack2_max) * ((size_t)ne));
    ia_hpsi[0] = fetch_dev_mem_indx(((size_t)tile_npack2_max) * ((size_t)ne));

    if (tile_fac > 1) {
      ia_psi[1] = fetch_dev_mem_indx(((size_t)tile_npack2_max) * ((size_t)ne));
      ia_hpsi[1] = fetch_dev_mem_indx(((size_t)tile_npack2_max) * ((size_t)ne));
    }
    std::cout << "Into psi_alloc, tile_factor = " << tile_fac
              << " ndev_mem=" << ndev_mem << std::endl;
  }
  void psi_dealloc() {
    inuse[ia_psi[0]] = false;
    inuse[ia_hpsi[0]] = false;
    if (tile_fac > 1) {
      inuse[ia_psi[1]] = false;
      inuse[ia_hpsi[1]] = false;
    }
  }
  void psi_copy_host2gpu(int npack1, int ne, double *psi) {
    // hipMemcpy(dev_mem[ia_psi[0]],psi,tile_npack2_max*ne*sizeof(double),hipMemcpyHostToDevice);
    a_psi = psi;
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        tile_npack2[0], ne, sizeof(double), psi, 2 * npack1, dev_mem[ia_psi[0]],
        tile_npack2[0], stream[0]));
  }
  void hpsi_copy_host2gpu(int npack1, int ne, double *hpsi) {
    // hipMemcpy(dev_mem[ia_hpsi[0]],hpsi,2*npack1*ne*sizeof(double),hipMemcpyHostToDevice);
    int tt = tile_fac - 1;
    a_hpsi = hpsi;
    NWPW_ROCBLAS_ERROR(rocblas_set_matrix_async(
        tile_npack2[tt], ne, sizeof(double), &hpsi[tile_start2[tt]], 2 * npack1,
        dev_mem[ia_hpsi[tt % 2]], tile_npack2[tt], stream[tt % 2]));
  }
  void psi_copy_gpu2host(int npack1, int ne, double *psi) {
    // hipMemcpy(psi, dev_mem[ia_psi[0]],
    // 2*ne*npack1*sizeof(double),hipMemcpyDeviceToHost);
    if (tile_fac == 1) {
      NWPW_HIP_ERROR(hipStreamSynchronize(stream[0]));
      NWPW_ROCBLAS_ERROR(rocblas_get_matrix(tile_npack2[0], ne, sizeof(double),
                                            dev_mem[ia_psi[0]], tile_npack2[0],
                                            psi, 2 * npack1));
    }
  }
  void hpsi_copy_gpu2host(int npack1, int ne, double *hpsi) {
    NWPW_HIP_ERROR(hipStreamSynchronize(stream[0]));
    // hipMemcpy(hpsi, dev_mem[ia_hpsi[0]],
    // 2*ne*npack1*sizeof(double),hipMemcpyDeviceToHost); if (tile_fac==1) {
    //     NWPW_ROCBLAS_ERROR(
    //     rocblas_get_matrix(tile_npack2[0],ne,sizeof(double),
    //                                        dev_mem[ia_hpsi[0]],tile_npack2[0],
    //                                        hpsi,2*npack1));
    // }
  }

  /*******************************/
  /* fft functions (uses rocFFT) */
  /*******************************/
  int batch_fft_init(int nx, int ny, int nz, int nq1, int nq2, int nq3) {
    size_t length_nx = (size_t)nx;
    size_t length_ny = (size_t)ny;
    size_t length_nz = (size_t)nz;

    NWPW_ROCFFT_ERROR(rocfft_plan_create(
        &forward_plan_rx[fftcount], rocfft_placement_inplace,
        rocfft_transform_type_real_forward, rocfft_precision_double, (size_t)1,
        &length_nx, (size_t)nq1, nullptr));
    NWPW_ROCFFT_ERROR(rocfft_plan_create(
        &backward_plan_rx[fftcount], rocfft_placement_inplace,
        rocfft_transform_type_real_inverse, rocfft_precision_double, (size_t)1,
        &length_nx, (size_t)nq1, nullptr));

    NWPW_ROCFFT_ERROR(rocfft_plan_create(
        &forward_plan_x[fftcount], rocfft_placement_inplace,
        rocfft_transform_type_complex_forward, rocfft_precision_double,
        (size_t)1, &length_nx, (size_t)nq1, nullptr));
    NWPW_ROCFFT_ERROR(rocfft_plan_create(
        &backward_plan_x[fftcount], rocfft_placement_inplace,
        rocfft_transform_type_complex_inverse, rocfft_precision_double,
        (size_t)1, &length_nx, (size_t)nq1, nullptr));

    NWPW_ROCFFT_ERROR(rocfft_plan_create(
        &forward_plan_y[fftcount], rocfft_placement_inplace,
        rocfft_transform_type_complex_forward, rocfft_precision_double,
        (size_t)1, &length_ny, (size_t)nq2, nullptr));
    NWPW_ROCFFT_ERROR(rocfft_plan_create(
        &backward_plan_y[fftcount], rocfft_placement_inplace,
        rocfft_transform_type_complex_inverse, rocfft_precision_double,
        (size_t)1, &length_ny, (size_t)nq2, nullptr));

    NWPW_ROCFFT_ERROR(rocfft_plan_create(
        &forward_plan_z[fftcount], rocfft_placement_inplace,
        rocfft_transform_type_complex_forward, rocfft_precision_double,
        (size_t)1, &length_nz, (size_t)nq3, nullptr));
    NWPW_ROCFFT_ERROR(rocfft_plan_create(
        &backward_plan_z[fftcount], rocfft_placement_inplace,
        rocfft_transform_type_complex_inverse, rocfft_precision_double,
        (size_t)1, &length_nz, (size_t)nq3, nullptr));

    nxfft[fftcount] = nx;
    nyfft[fftcount] = ny;
    nzfft[fftcount] = nz;

    int tag = fftcount;
    ++fftcount;

    return tag;
  }

  void batch_fft_pipeline_mem_init(const int nstages, const int n2ft3d) {
    ifft_n = nstages;

    // allocate memory and cuda streams
    for (auto i = 0; i < ifft_n; ++i)
      ifft_dev[i] = fetch_dev_mem_indx(((size_t)n2ft3d));
  }

  void batch_fft_end(const int tag) {
    NWPW_ROCFFT_ERROR(rocfft_plan_destroy(forward_plan_x[tag]));
    NWPW_ROCFFT_ERROR(rocfft_plan_destroy(forward_plan_y[tag]));
    NWPW_ROCFFT_ERROR(rocfft_plan_destroy(forward_plan_z[tag]));

    NWPW_ROCFFT_ERROR(rocfft_plan_destroy(backward_plan_x[tag]));
    NWPW_ROCFFT_ERROR(rocfft_plan_destroy(backward_plan_y[tag]));
    NWPW_ROCFFT_ERROR(rocfft_plan_destroy(backward_plan_z[tag]));

    --fftcount;

    // free dev_mem
    for (auto i = 0; i < ndev_mem; ++i) {
      inuse[i] = false;
      NWPW_HIP_ERROR(hipFree(dev_mem[i]));
    }
    ndev_mem = 0;
  }

  void batch_rfftx(const int fft_indx, bool forward, int nx, int nq, int n2ft3d, double *a) 
  {
     int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
     NWPW_HIP_ERROR(hipMemcpy(dev_mem[ia_dev], a, n2ft3d * sizeof(double), hipMemcpyHostToDevice));
    
     if (forward) {
       NWPW_ROCFFT_ERROR(rocfft_execute(
           forward_plan_rx[fft_indx],
           reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
     } else {
       NWPW_ROCFFT_ERROR(rocfft_execute(
           backward_plan_rx[fft_indx],
           reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
     }
    
     NWPW_HIP_ERROR(hipMemcpy(a, dev_mem[ia_dev], n2ft3d * sizeof(double), hipMemcpyDeviceToHost));
    
     inuse[ia_dev] = false;
  }

  void batch_rfftx_stages(const int stage, const int fft_indx, bool forward,
                          int nx, int nq, int n2ft3d, double *a, int da) 
  {
     // int ia_dev = fetch_dev_mem_indx(((size_t) n2ft3d));
     int ia_dev = ifft_dev[da];
    
     if (stage == 0) {
       inuse[ia_dev] = true;
       NWPW_HIP_ERROR(hipMemcpyAsync(dev_mem[ia_dev], a, n2ft3d * sizeof(double),
                                     hipMemcpyHostToDevice, stream[da]));
     } else if (stage == 1) {
       // NWPW_HIP_ERROR(hipStreamSynchronize(stream[da]));
       if (forward) {
         NWPW_ROCFFT_ERROR(rocfft_execute(
             forward_plan_rx[fft_indx],
             reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
       } else {
         NWPW_ROCFFT_ERROR(rocfft_execute(
             backward_plan_rx[fft_indx],
             reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
       }
       NWPW_HIP_ERROR(hipMemcpyAsync(a, dev_mem[ia_dev], n2ft3d * sizeof(double),
                                     hipMemcpyDeviceToHost, stream[da]));
     } else if (stage == 2) {
       NWPW_HIP_ERROR(hipStreamSynchronize(stream[da]));
       inuse[ia_dev] = false;
     }
  }


  void batch_cfftx(const int fft_indx, bool forward, int nx, int nq, int n2ft3d, double *a) 
  {
     int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
     NWPW_HIP_ERROR(hipMemcpy(dev_mem[ia_dev], a, n2ft3d * sizeof(double), hipMemcpyHostToDevice));
   
     if (forward) {
       NWPW_ROCFFT_ERROR(rocfft_execute(
           forward_plan_x[fft_indx],
           reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
     } else {
       NWPW_ROCFFT_ERROR(rocfft_execute(
           backward_plan_x[fft_indx],
           reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
     }
   
     NWPW_HIP_ERROR(hipMemcpy(a, dev_mem[ia_dev], n2ft3d * sizeof(double), hipMemcpyDeviceToHost));
   
     inuse[ia_dev] = false;
  } 
      
  void batch_cfftx_stages(const int stage, const int fft_indx, bool forward,
                          int nx, int nq, int n2ft3d, double *a, int da) 
  {
     // int ia_dev = fetch_dev_mem_indx(((size_t) n2ft3d));
     int ia_dev = ifft_dev[da];
     
     if (stage == 0) {
       inuse[ia_dev] = true;
       NWPW_HIP_ERROR(hipMemcpyAsync(dev_mem[ia_dev], a, n2ft3d * sizeof(double), hipMemcpyHostToDevice, stream[da]));
     } else if (stage == 1) {
       // NWPW_HIP_ERROR(hipStreamSynchronize(stream[da]));
       if (forward) {
         NWPW_ROCFFT_ERROR(rocfft_execute(
             forward_plan_x[fft_indx],
             reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
       } else {
         NWPW_ROCFFT_ERROR(rocfft_execute(
             backward_plan_x[fft_indx],
             reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
       }
       NWPW_HIP_ERROR(hipMemcpyAsync(a, dev_mem[ia_dev], n2ft3d * sizeof(double),
                                     hipMemcpyDeviceToHost, stream[da]));
     } else if (stage == 2) {
       NWPW_HIP_ERROR(hipStreamSynchronize(stream[da]));
       inuse[ia_dev] = false;
     }
  }




  void batch_cffty(const int fft_indx, bool forward, int ny, int nq, int n2ft3d,
                   double *a) {
    int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
    NWPW_HIP_ERROR(hipMemcpy(dev_mem[ia_dev], a, n2ft3d * sizeof(double),
                             hipMemcpyHostToDevice));

    if (forward) {
      NWPW_ROCFFT_ERROR(rocfft_execute(
          forward_plan_y[fft_indx],
          reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
    } else {
      NWPW_ROCFFT_ERROR(rocfft_execute(
          backward_plan_y[fft_indx],
          reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
    }

    NWPW_HIP_ERROR(hipMemcpy(a, dev_mem[ia_dev], n2ft3d * sizeof(double),
                             hipMemcpyDeviceToHost));

    inuse[ia_dev] = false;
  }

  void batch_cffty_stages(const int stage, const int fft_indx, bool forward,
                          int ny, int nq, int n2ft3d, double *a, int da) {
    // int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
    int ia_dev = ifft_dev[da];
    if (stage == 0) {
      inuse[ia_dev] = true;
      NWPW_HIP_ERROR(hipMemcpyAsync(dev_mem[ia_dev], a, n2ft3d * sizeof(double),
                                    hipMemcpyHostToDevice, stream[da]));
    } else if (stage == 1) {
      // NWPW_HIP_ERROR(hipStreamSynchronize(stream[da]));
      if (forward) {
        NWPW_ROCFFT_ERROR(rocfft_execute(
            forward_plan_y[fft_indx],
            reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
      } else {
        NWPW_ROCFFT_ERROR(rocfft_execute(
            backward_plan_y[fft_indx],
            reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
      }
      NWPW_HIP_ERROR(hipMemcpyAsync(a, dev_mem[ia_dev], n2ft3d * sizeof(double),
                                    hipMemcpyDeviceToHost, stream[da]));
    } else if (stage == 2) {
      NWPW_HIP_ERROR(hipStreamSynchronize(stream[da]));
      inuse[ia_dev] = false;
    }
  }

  void batch_cfftz(const int fft_indx, bool forward, int nz, int nq, int n2ft3d,
                   double *a) {
    int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
    NWPW_HIP_ERROR(hipMemcpy(dev_mem[ia_dev], a, n2ft3d * sizeof(double),
                             hipMemcpyHostToDevice));

    if (forward) {
      NWPW_ROCFFT_ERROR(rocfft_execute(
          forward_plan_z[fft_indx],
          reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
    } else {
      NWPW_ROCFFT_ERROR(rocfft_execute(
          backward_plan_z[fft_indx],
          reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
    }

    NWPW_HIP_ERROR(hipMemcpy(a, dev_mem[ia_dev], n2ft3d * sizeof(double),
                             hipMemcpyDeviceToHost));

    inuse[ia_dev] = false;
  }

  void batch_cfftz_stages(const int stage, const int fft_indx, bool forward,
                          int nz, int nq, int n2ft3d, double *a, int da) {
    // int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
    int ia_dev = ifft_dev[da];

    if (stage == 0) {
      inuse[ia_dev] = true;
      NWPW_HIP_ERROR(hipMemcpyAsync(dev_mem[ia_dev], a, n2ft3d * sizeof(double),
                                    hipMemcpyHostToDevice, stream[da]));
    } else if (stage == 1) {
      // NWPW_HIP_ERROR(hipStreamSynchronize(stream[da]));
      if (forward) {
        NWPW_ROCFFT_ERROR(rocfft_execute(
            forward_plan_z[fft_indx],
            reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
      } else {
        NWPW_ROCFFT_ERROR(rocfft_execute(
            backward_plan_z[fft_indx],
            reinterpret_cast<void **>(&(dev_mem[ia_dev])), nullptr, nullptr));
      }
      NWPW_HIP_ERROR(hipMemcpyAsync(a, dev_mem[ia_dev], n2ft3d * sizeof(double),
                                    hipMemcpyDeviceToHost, stream[da]));
    } else if (stage == 2) {
      NWPW_HIP_ERROR(hipStreamSynchronize(stream[da]));
      inuse[ia_dev] = false;
    }
  }

  static void eigsrt_device(double *D, double *V, int n) {
    int i, j, k;
    double p;

    for (i = 0; i < (n - 1); ++i) {
      k = i;
      p = D[i];
      for (j = i + 1; j < n; ++j)
        if (D[j] >= p) {
          k = j;
          p = D[j];
        }

      if (k != i) {
        D[k] = D[i];
        D[i] = p;
        for (j = 0; j < n; ++j) {
          p = V[j + i * n];
          V[j + i * n] = V[j + k * n];
          V[j + k * n] = p;
        }
      }
    }
  }

  void NN_eigensolver(int ispin, int ne[], double *host_hml, double *host_eig) {
    int n, ierr;
    int nn = ne[0] * ne[0] + 14;
    double xmp1[nn];
    // double *xmp1 = new (std::nothrow) double[nn]();

    int shift1 = 0;
    int shift2 = 0;
    for (int ms = 0; ms < ispin; ++ms) {
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
        shift1 += 2*ne[0];
        shift2 += 4*ne[0]*ne[0];
     } 
  }

};

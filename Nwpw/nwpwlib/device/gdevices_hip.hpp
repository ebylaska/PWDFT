// NWPW_HIP Routines

#pragma once

#include "blas.h"

//#include        "gdevice.hpp"

#include <hip_runtime.h>
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

  static const char *_rocfftGetErrorEnum(rocfftResult error) {
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
  rocfft_exception(std::string file, int line, rocfftResult err)
      : file_(file), line_(line), err_code_(err) {}
};

class rocblas_exception : public std::exception {

  std::string file_;
  int line_;
  rocblas_status err_code_;

  static const char *_rocblasGetErrorEnum(rocblas_status error) {
    switch (error) {
    case rocblas_status_success:
      return "ROCBLAS_STATUS_SUCCESS";

    case rocblas_status_invalid_handle:
      return "ROCBLAS_STATUS_INVALID_HANDLE";

    case rocblas_status_not_implemented:
      return "ROCBLAS_STATUS_NOT_IMPLEMENTED";

    case rocblas_status_invalid_pointer:
      return "ROCBLAS_STATUS_INVALID_POINTER";

    case rocblas_status_invalid_size:
      return "ROCBLAS_STATUS_INVALID_SIZE";

    case rocblas_status_memory_error:
      return "ROCBLAS_STATUS_MEMORY_ERROR";

    case rocblas_status_internal_error:
      return "ROCBLAS_STATUS_INTERNAL_ERROR";

    case rocblas_status_perf_degraded:
      return "ROCBLAS_STATUS_PERF_DEGRADED";

    case rocblas_status_size_query_mismatch:
      return "ROCBLAS_STATUS_SIZE_QUERY_MISMATCH";

    case rocblas_status_size_increased:
      return "ROCBLAS_STATUS_SIZE_INCREASED";

    case rocblas_status_size_unchanged:
      return "ROCBLAS_STATUS_SIZE_UNCHANGED";

    case rocblas_status_invalid_value:
      return "ROCBLAS_STATUS_INVALID_VALUE";

    case rocblas_status_continue:
      return "ROBLAS_STATUS_CONTINUE";

    case rocblas_status_check_numerics_fail:
      return "ROCBLAS_STATUS_CHECK_NUMERICS_FAIL";
    }

    return "<unknown>";
  }

  const char *what() const noexcept override {
    std::stringstream ss;
    ss << "rocBLAS Exception, "
       << " Error Code: " << _rocblasGetErrorEnum(err_code_) << std::endl
       << " at " << file_ << " : " << line_ << std::endl;
    auto msg = ss.str();
    return strdup(msg.c_str());
  }

public:
  rocblas_exception(std::string file, int line, rocblas_status err)
      : file_(file), line_(line), err_code_(err) {}
};

#define NWPW_HIP_ERROR(ERR)                                                    \
  if (ERR != hipSuccess)                                                       \
    throw hip_exception(__FILE__, __LINE__, ERR);

#define NWPW_ROCFFT_ERROR(ERR)                                                 \
  if (ERR != ROCFFT_SUCCESS)                                                   \
    throw rocfft_exception(__FILE__, __LINE__, ERR);

#define NWPW_ROCBLAS_ERROR(ERR)                                                \
  if (ERR != rocblas_status_success)                                           \
    throw rocblas_exception(__FILE__, __LINE__, ERR);

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

  rocfft_execution_info pExecInfo = nullptr;
  rocfft_plan forward_plan_x = nullptr, plan_y = nullptr, plan_z = nullptr;
  rocfft_plan backward_plan_x = nullptr;

  rocblas_handle master_handle = 0;
  rocblas_operation matT = rocblas_operation_transpose;
  rocblas_operation matN = rocblas_operation_none;

  hipStream_t stream[2];

public:
  bool hasgpu = true;

  /* device memory */
  int ndev_mem = 0;
  bool inuse[19] = {false};
  size_t ndsize_mem[19];
  double *dev_mem[19];
  int tile_fac = 1;
  int tile_npack2_max;
  int tile_npack2[19], tile_start2[19];
  double *a_psi, *a_hpsi, *b_prj;
  int ia_psi[2], ia_hpsi[2], ib_prj[2];

  /* constructor */
  Gdevices() {
    ndev_mem = 0;

    std::cout << "Into rocblasCreate" << std::endl;

    NWPW_ROCBLAS_ERROR(rocblas_create_handle(&master_handle));

    // allocate hip streams
    for (auto i = 0; i < 2; ++i)
      NWPW_HIP_ERROR(hipStreamCreate(&stream[i]));

    NWPW_ROCFFT_ERROR(rocfft_execution_info_create(&pExecInfo));
    NWPW_ROCFFT_ERROR(rocfft_execution_info_set_stream(pExecInfo, stream[0]));
  }

  /* deconstructor */
  ~Gdevices() {

    // free dev_mem
    for (auto i = 0; i < ndev_mem; ++i)
      NWPW_HIP_ERROR(hipFree(dev_mem[i]));
    ndev_mem = 0;

    // free hip streams
    for (auto i = 0; i < 2; ++i)
      NWPW_HIP_ERROR(hipStreamDestroy(stream[i]));

    rocblas_destroy_handle(master_handle);

    // free fft descriptors
    NWPW_ROCFFT_ERROR(rocfft_execution_info_destroy(pExecInfo));
    // cufftDestroy(forward_plan_x);
    // cufftDestroy(plan_y);
    // cufftDestroy(plan_z);
    // cufftDestroy(backward_plan_x);
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
    hipMemcpy(host_c, dev_mem[ic], ne * nprj * sizeof(double),
              hipMemcpyDeviceToHost);

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
  void batch_fft_init(int nx, int ny, int nz, int nq1, int nq2, int nq3) {
    std::cout << "Into batch_fft_init" << std::endl;
    NWPW_ROCFFT_ERROR(cufftPlan1d(&forward_plan_x, nx, CUFFT_D2Z, nq1));
    NWPW_ROCFFT_ERROR(cufftPlan1d(&backward_plan_x, nx, CUFFT_Z2D, nq1));

    int y_inembed[] = {ny};
    int y_onembed[] = {ny};
    NWPW_ROCFFT_ERROR(cufftPlanMany(&plan_y, 1, &ny, y_inembed, 1, ny,
                                    y_onembed, 1, ny, CUFFT_Z2Z, nq2));

    int z_inembed[] = {nz};
    int z_onembed[] = {nz};
    NWPW_ROCFFT_ERROR(cufftPlanMany(&plan_z, 1, &nz, z_inembed, 1, nz,
                                    z_onembed, 1, nz, CUFFT_Z2Z, nq3));
  }

  void batch_fft_end() {
    // free fft descriptors
    cufftDestroy(forward_plan_x);
    cufftDestroy(plan_y);
    cufftDestroy(plan_z);
    cufftDestroy(backward_plan_x);

    // free dev_mem
    for (auto i = 0; i < ndev_mem; ++i)
      NWPW_HIP_ERROR(hipFree(dev_mem[i]));
    ndev_mem = 0;
  }

  void batch_cfftx(bool forward, int nx, int nq, int n2ft3d, double *a) {
    int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
    NWPW_HIP_ERROR(hipMemcpy(dev_mem[ia_dev], a, n2ft3d * sizeof(double),
                             hipMemcpyHostToDevice));

    if (forward) {
      NWPW_ROCFFT_ERROR(rocfft_execute(
          forward_plan_x, reinterpret_cast<cufftDoubleReal *>(dev_mem[ia_dev]),
          reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev])));
    } else {
      NWPW_ROCFFT_ERROR(rocfft_execute(
          backward_plan_x,
          reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
          reinterpret_cast<cufftDoubleReal *>(dev_mem[ia_dev])));
    }

    NWPW_HIP_ERROR(hipMemcpy(a, dev_mem[ia_dev], n2ft3d * sizeof(double),
                             hipMemcpyDeviceToHost));

    inuse[ia_dev] = false;
  }

  void batch_cffty(bool forward, int ny, int nq, int n2ft3d, double *a) {
    int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
    NWPW_HIP_ERROR(hipMemcpy(dev_mem[ia_dev], a, n2ft3d * sizeof(double),
                             hipMemcpyHostToDevice));

    if (forward) {
      NWPW_ROCFFT_ERROR(cufftExecZ2Z(
          plan_y, reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
          reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
          CUFFT_FORWARD));
    } else {
      NWPW_ROCFFT_ERROR(cufftExecZ2Z(
          plan_y, reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
          reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
          CUFFT_INVERSE));
    }

    NWPW_HIP_ERROR(hipMemcpy(a, dev_mem[ia_dev], n2ft3d * sizeof(double),
                             hipMemcpyDeviceToHost));

    inuse[ia_dev] = false;
  }

  void batch_cfftz(bool forward, int nz, int nq, int n2ft3d, double *a) {
    int ia_dev = fetch_dev_mem_indx(((size_t)n2ft3d));
    NWPW_HIP_ERROR(hipMemcpy(dev_mem[ia_dev], a, n2ft3d * sizeof(double),
                             hipMemcpyHostToDevice));

    if (forward) {
      NWPW_ROCFFT_ERROR(cufftExecZ2Z(
          plan_z, reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
          reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
          CUFFT_FORWARD));
    } else {
      NWPW_ROCFFT_ERROR(cufftExecZ2Z(
          plan_z, reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
          reinterpret_cast<cufftDoubleComplex *>(dev_mem[ia_dev]),
          CUFFT_INVERSE));
    }

    NWPW_HIP_ERROR(hipMemcpy(a, dev_mem[ia_dev], n2ft3d * sizeof(double),
                             hipMemcpyDeviceToHost));

    inuse[ia_dev] = false;
  }
};

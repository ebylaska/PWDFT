// NWPW_CUDA Routines

#pragma once

#include        "blas.h"

//#include        "gdevice.hpp"

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cufft.h>
#include <cusolverDn.h>

#include <stdexcept>
#include <string>
#include <sstream>
#include <iostream>

class cuda_exception : public std::exception {

    std::string file_;
    int         line_;
    cudaError_t err_code_;

    const char* what() const noexcept override {
        std::stringstream ss;
        ss << "CUDA Exception, " << cudaGetErrorString( err_code_ ) << " at " << std::endl
           << file_ << " : " << line_ << std::endl;
        auto msg = ss.str();
        return strdup( msg.c_str() );
    }

public:

    cuda_exception( std::string file, int line, cudaError_t err ) :
        file_(file), line_(line), err_code_(err) { }
};

class cufft_exception : public std::exception {

    std::string file_;
    int         line_;
    cufftResult err_code_;

    static const char* _cudaGetErrorEnum(cufftResult error) {
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

    const char* what() const noexcept override {
        std::stringstream ss;
        ss << "CUFFT Exception, " << " Error Code: " << _cudaGetErrorEnum(err_code_) << std::endl
           << " at " << file_ << " : " << line_ << std::endl;
        auto msg = ss.str();
        return strdup( msg.c_str() );
    }

public:

    cufft_exception( std::string file, int line, cufftResult err ) :
        file_(file), line_(line), err_code_(err) { }
};

class cublas_exception : public std::exception {

    std::string file_;
    int         line_;
    cublasStatus_t err_code_;

    static const char* _cudaGetErrorEnum(cublasStatus_t error) {
        switch (error)
        {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";

        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";

        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";

        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";

        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";

        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";

        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";

        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
        }

        return "<unknown>";
    }

    const char* what() const noexcept override {
        std::stringstream ss;
        ss << "CUBLAS Exception, " << " Error Code: " << _cudaGetErrorEnum(err_code_) << std::endl
           << " at " << file_ << " : " << line_ << std::endl;
        auto msg = ss.str();
        return strdup( msg.c_str() );
    }

public:

    cublas_exception( std::string file, int line, cublasStatus_t err ) :
        file_(file), line_(line), err_code_(err) { }
};

#define NWPW_CUDA_ERROR( ERR )                                  \
    if( ERR != cudaSuccess )                                    \
        throw cuda_exception( __FILE__, __LINE__, ERR );

#define NWPW_CUFFT_ERROR( ERR )                                 \
    if( ERR != CUFFT_SUCCESS )                                  \
        throw cufft_exception( __FILE__, __LINE__, ERR );

#define NWPW_CUBLAS_ERROR( ERR )                                \
    if( ERR != CUBLAS_STATUS_SUCCESS )                          \
        throw cublas_exception( __FILE__, __LINE__, ERR );



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

    cufftHandle forward_plan_x = 0 , plan_y = 0, plan_z = 0;
    cufftHandle backward_plan_x = 0;

    cublasHandle_t master_handle = 0;
    cublasOperation_t matT = CUBLAS_OP_T;
    cublasOperation_t matN = CUBLAS_OP_N;

    cusolverEigMode_t jobz = CUSOLVER_EIG_MODE_VECTOR; // compute eigenvalues and eigenvectors.
    cublasFillMode_t uplo = CUBLAS_FILL_MODE_LOWER;
    cusolverDnHandle_t cusolverH = NULL;
    cudaStream_t solverstream = NULL;

    cudaStream_t stream[2];

public:
    bool hasgpu = true;

    /* device memory */
    int    ndev_mem = 0;
    bool   inuse[19] = {false};
    size_t ndsize_mem[19];
    double *dev_mem[19];
    int    tile_fac=1;
    int    tile_npack2_max;
    int    tile_npack2[19],tile_start2[19];
    double *a_psi,*a_hpsi,*b_prj;
    int    ia_psi[2],ia_hpsi[2],ib_prj[2];

    int *d_info[2];
    int lwork=0;
    double *d_work;

    /* constructor */
    Gdevices() {
        ndev_mem = 0;

        std::cout << "Into cublasCreate" << std::endl;

        NWPW_CUBLAS_ERROR( cublasCreate(&master_handle) );

        // allocate cuda streams
        for (auto i=0; i<2; ++i) NWPW_CUDA_ERROR( cudaStreamCreate(&stream[i]) );


        // create cusolver handle, bind a stream 
        CUSOLVER_CHECK(cusolverDnCreate(&cusolverH));
        CUDA_CHECK(cudaStreamCreateWithFlags(&solverstream, cudaStreamNonBlocking));
        CUSOLVER_CHECK(cusolverDnSetStream(cusolverH, solverstream));

        // query working space of syevd
        CUSOLVER_CHECK(cusolverDnDsyevd_bufferSize(cusolverH, jobz, uplo, m, d_A, lda, d_W, &lwork));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&d_work), sizeof(double) * lwork));
        for (int i=0; i<2; ++i) NWPW_CUDA_ERROR(cudaMalloc(reinterpret_cast<void **>(&d_info[i]),sizeof(int)));

    }

    /* deconstructor */
    ~Gdevices() {

        // free dev_mem
        for (auto i=0; i<ndev_mem; ++i)  cudaFree(dev_mem[i]);
        ndev_mem = 0;

        // free cuda streams
        for (auto i=0; i<2; ++i) cudaStreamDestroy(stream[i]);

        cublasDestroy(master_handle);


       /* free cusolver resources */
       for (auto i=0; i<2; ++i) NWPW_CUDA_ERROR(cudaFree(d_info[i]));
       NWPW_CUDA_ERROR(cudaFree(d_work));
       NWPW_CUSOLVER_ERROR(cusolverDnDestroy(cusolverH));

       NWPW_CUDA_ERROR(cudaStreamDestroy(solverstream));
       NWPW_CUDA_ERROR(cudaDeviceReset());

        // free fft descriptors
        //cufftDestroy(forward_plan_x);
        //cufftDestroy(plan_y);
        //cufftDestroy(plan_z);
        //cufftDestroy(backward_plan_x);
    }

    int fetch_dev_mem_indx(const size_t ndsize) {
        int ii = 0;
        while ((((ndsize!=ndsize_mem[ii]) || inuse[ii])) && (ii<ndev_mem))
            ++ii;

        if (ii<ndev_mem) {
            inuse[ii] = true;
        } else {
            ii            = ndev_mem;
            inuse[ii]     = true;
            ndsize_mem[ii] = ndsize;
            NWPW_CUDA_ERROR( cudaMalloc((void**)&(dev_mem[ii]), ndsize*sizeof(double)) );
            ndev_mem += 1;
        }

	NWPW_CUDA_ERROR( cudaMemset(dev_mem[ii], 0, ndsize*sizeof(double)) );
        return ii;
    }

   /**************************************
    *                                    *
    *              TN4_dgemm             *
    *                                    *
    **************************************/
    /* This function computes <host_a|host_a>, <host_a|host_b>, and <host_b|host_b> overlap matrices.

        host_caa = beta*host_caa + alpha*host_a'*host_a
        host_cab = beta*host_cab + alpha*host_a'*host_b
        host_cba = beta*host_cba + alpha*host_b'*host_a
        host_cbb = beta*host_cbb + alpha*host_b'*host_b

       Entry - npack2,ne: matrix size
               alpha, beta: standard dgemm parameters
               host_a: (npack2xne) matrix
               host_b: (npack2xne) matrix
       Exit - host_caa,host_cab,host_cba,host_cbb: (nexne) matrices
       Uses - device memory for (npack2xne) matrices ia_psi, and ia_hpsi allocated previously with psi_alloc
            - temporary device memory for (nexne) matrices ic11, ic12, ic21, and ic22.
    */
    void TN4_dgemm(int npack2, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_caa, double *host_cab, double *host_cba, double *host_cbb) {
       int ic11 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));
       int ic12 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));
       int ic21 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));
       int ic22 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));

       if (std::fabs(beta)>0.0) {
          NWPW_CUDA_ERROR( cudaMemcpy(dev_mem[ic11], host_caa,ne*ne*sizeof(double),cudaMemcpyHostToDevice) );
          NWPW_CUDA_ERROR( cudaMemcpy(dev_mem[ic12], host_cab,ne*ne*sizeof(double),cudaMemcpyHostToDevice) );
          NWPW_CUDA_ERROR( cudaMemcpy(dev_mem[ic21], host_cba,ne*ne*sizeof(double),cudaMemcpyHostToDevice) );
          NWPW_CUDA_ERROR( cudaMemcpy(dev_mem[ic22], host_cbb,ne*ne*sizeof(double),cudaMemcpyHostToDevice) );
       }

      // copy host_a,host_b --> dev_mem
      NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                                              &host_a[tile_start2[0]],npack2,
                                              dev_mem[ia_psi[0]],tile_npack2[0],stream[0]) );
      NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                                              &host_b[tile_start2[0]],npack2,
                                              dev_mem[ia_hpsi[0]],tile_npack2[0],stream[0]) );

      double beta0 = beta;
      for (auto tt=0; tt<tile_fac; ++tt) {
         int ttp1 = tt+1;
         if (ttp1<tile_fac) {
            NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                                                    &host_a[tile_start2[ttp1]],npack2,
                                                    dev_mem[ia_psi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]) );
            NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                                                    &host_b[tile_start2[ttp1]],npack2,
                                                    dev_mem[ia_hpsi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]) );
         }
         NWPW_CUDA_ERROR( cudaStreamSynchronize(stream[tt%2]) );
         NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
                                        matT, matN,
                                        ne,ne,tile_npack2[tt],&alpha,
                                        dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                                        dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                                        &beta0,dev_mem[ic11],ne) );
         NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
                                        matT, matN,
                                        ne,ne,tile_npack2[tt],&alpha,
                                        dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                                        dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                                        &beta0,dev_mem[ic12],ne) );
         NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
                                        matT, matN,
                                        ne,ne,tile_npack2[tt],&alpha,
                                        dev_mem[ia_hpsi[tt%2]], tile_npack2[tt],
                                        dev_mem[ia_psi[tt%2]],  tile_npack2[tt],
                                        &beta0,dev_mem[ic21],ne) );
         NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
                                        matT, matN,
                                        ne,ne,tile_npack2[tt],&alpha,
                                        dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                                        dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                                        &beta0,dev_mem[ic22],ne) );
         beta0 = 1.0;
      }

      NWPW_CUDA_ERROR( cudaMemcpy(host_caa,dev_mem[ic11],ne*ne*sizeof(double),cudaMemcpyDeviceToHost) );
      NWPW_CUDA_ERROR( cudaMemcpy(host_cab,dev_mem[ic12],ne*ne*sizeof(double),cudaMemcpyDeviceToHost) );
      NWPW_CUDA_ERROR( cudaMemcpy(host_cba,dev_mem[ic21],ne*ne*sizeof(double),cudaMemcpyDeviceToHost) );
      NWPW_CUDA_ERROR( cudaMemcpy(host_cbb,dev_mem[ic22],ne*ne*sizeof(double),cudaMemcpyDeviceToHost) );

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
    /* This function computes <host_a|host_a>, <host_a|host_b>, and <host_b|host_b> overlap matrices.

        host_caa = beta*host_caa + alpha*host_a'*host_a
        host_cab = beta*host_cab + alpha*host_a'*host_b
        host_cbb = beta*host_cbb + alpha*host_b'*host_b

       Entry - npack2,ne: matrix size
               alpha, beta: standard dgemm parameters
               host_a: (npack2xne) matrix
               host_b: (npack2xne) matrix
       Exit - host_caa,host_cab,host_cbb: (nexne) matrices
       Uses - device memory for (npack2xne) matrices ia_psi, and ia_hpsi allocated previously with psi_alloc
            - temporary device memory for (nexne) matrices ic11, ic12, and ic22.
    */
    void TN3_dgemm(int npack2, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_caa, double *host_cab, double *host_cbb) {
       int ic11 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));
       int ic12 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));
       int ic22 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));

       if (std::fabs(beta)>0.0) {
          NWPW_CUDA_ERROR( cudaMemcpy(dev_mem[ic11], host_caa,ne*ne*sizeof(double),cudaMemcpyHostToDevice) );
          NWPW_CUDA_ERROR( cudaMemcpy(dev_mem[ic12], host_cab,ne*ne*sizeof(double),cudaMemcpyHostToDevice) );
          NWPW_CUDA_ERROR( cudaMemcpy(dev_mem[ic22], host_cbb,ne*ne*sizeof(double),cudaMemcpyHostToDevice) );
       }

      // copy host_a,host_b --> dev_mem
      NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                                              &host_a[tile_start2[0]],npack2,
                                              dev_mem[ia_psi[0]],tile_npack2[0],stream[0]) );
      NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                                              &host_b[tile_start2[0]],npack2,
                                              dev_mem[ia_hpsi[0]],tile_npack2[0],stream[0]) );

      double beta0 = beta;
      for (auto tt=0; tt<tile_fac; ++tt) {
         int ttp1 = tt+1;
         if (ttp1<tile_fac) {
            NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                                                    &host_a[tile_start2[ttp1]],npack2,
                                                    dev_mem[ia_psi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]) );
            NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                                                    &host_b[tile_start2[ttp1]],npack2,
                                                    dev_mem[ia_hpsi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]) );
         }
         NWPW_CUDA_ERROR( cudaStreamSynchronize(stream[tt%2]) );
         NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
                                        matT, matN,
                                        ne,ne,tile_npack2[tt],&alpha,
                                        dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                                        dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                                        &beta0,dev_mem[ic11],ne) );
         NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
                                        matT, matN,
                                        ne,ne,tile_npack2[tt],&alpha,
                                        dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                                        dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                                        &beta0,dev_mem[ic12],ne) );
         NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
                                        matT, matN,
                                        ne,ne,tile_npack2[tt],&alpha,
                                        dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                                        dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                                        &beta0,dev_mem[ic22],ne) );
         beta0 = 1.0;
      }

      NWPW_CUDA_ERROR( cudaMemcpy(host_caa,dev_mem[ic11],ne*ne*sizeof(double),cudaMemcpyDeviceToHost) );
      NWPW_CUDA_ERROR( cudaMemcpy(host_cab,dev_mem[ic12],ne*ne*sizeof(double),cudaMemcpyDeviceToHost) );
      NWPW_CUDA_ERROR( cudaMemcpy(host_cbb,dev_mem[ic22],ne*ne*sizeof(double),cudaMemcpyDeviceToHost) );

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
       Uses - device memory for (npack2xne) matrices ia_psi, and ia_hpsi allocated previously with psi_alloc
            - temporary device memory for (nexne) matrix, ic12.
    */
    void TN1_dgemm(int npack2, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_cab) {
       int ic12 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));

       if (std::fabs(beta)>0.0) {
          NWPW_CUDA_ERROR( cudaMemcpy(dev_mem[ic12], host_cab,ne*ne*sizeof(double),cudaMemcpyHostToDevice) );
       }

       // copy host_a,host_b --> dev_mem
       NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                                               &host_a[tile_start2[0]],npack2,
                                               dev_mem[ia_psi[0]],tile_npack2[0],stream[0]) );
       NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                                               &host_b[tile_start2[0]],npack2,
                                               dev_mem[ia_hpsi[0]],tile_npack2[0],stream[0]) );

       double beta0 = beta;
       for (auto tt=0; tt<tile_fac; ++tt) {
          int ttp1 = tt+1;
          if (ttp1<tile_fac) {
             NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                                                     &host_a[tile_start2[ttp1]],npack2,
                                                     dev_mem[ia_psi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]) );
             NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                                                     &host_b[tile_start2[ttp1]],npack2,
                                                     dev_mem[ia_hpsi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]) );
          }
          NWPW_CUDA_ERROR( cudaStreamSynchronize(stream[tt%2]) );
          NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
                                         matT, matN,
                                         ne,ne,tile_npack2[tt],&alpha,
                                         dev_mem[ia_psi[tt%2]], tile_npack2[tt],
                                         dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                                         &beta0,dev_mem[ic12],ne) );
          beta0 = 1.0;
       }

       NWPW_CUDA_ERROR( cudaMemcpy(host_cab,dev_mem[ic12],ne*ne*sizeof(double),cudaMemcpyDeviceToHost) );

       inuse[ic12] = false;
    }

   /************************************** 
    *                                    *
    *              NN_dgemm              *
    *                                    *
    **************************************/
    void NN_dgemm(int npack2, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
       //DGEMM_PWDFT((char *) "N",(char *) "N",npack2,ne,ne,alpha,host_a,npack2,host_b,ne,beta,host_c,npack2);
       
       int ib = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne));

       NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(ne,ne,sizeof(double),host_b,ne,dev_mem[ib],ne,stream[0]) );
       NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                                               &host_a[tile_start2[0]],npack2,
                                               dev_mem[ia_psi[0]],tile_npack2[0],stream[0]) );
       NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                                               &host_c[tile_start2[0]],npack2,
                                               dev_mem[ia_hpsi[0]],tile_npack2[0],stream[0]) );

       //double beta0 = beta;
       for (auto tt=0; tt<tile_fac; ++tt) {
          int ttp1 = tt+1;
          if (ttp1<tile_fac) {
             NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                                                     &host_a[tile_start2[ttp1]],npack2,
                                                     dev_mem[ia_psi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]) );
             NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                                                     &host_c[tile_start2[ttp1]],npack2,
                                                     dev_mem[ia_hpsi[ttp1%2]], tile_npack2[ttp1],stream[ttp1%2]) );
          }
          NWPW_CUDA_ERROR( cudaStreamSynchronize(stream[tt%2]) );
          NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
				       matN,matN,
				       tile_npack2[tt],ne,ne,&alpha,
				       dev_mem[ia_psi[tt%2]],tile_npack2[tt],
				       dev_mem[ib],ne,
				       &beta,dev_mem[ia_hpsi[tt%2]],tile_npack2[tt]) );
          NWPW_CUBLAS_ERROR( cublasGetMatrixAsync(tile_npack2[tt],ne,sizeof(double),
                                                   dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                                                   &host_c[tile_start2[tt]],npack2,stream[tt%2]) );

       }
       NWPW_CUDA_ERROR( cudaStreamSynchronize(stream[(tile_fac-1)%2]) );
       //NWPW_CUDA_ERROR( cudaMemcpy(host_c,dev_mem[ia_hpsi[0]],npack2*ne*sizeof(double),cudaMemcpyDeviceToHost) );

       inuse[ib] = false;
    }


   /************************************** 
    *                                    *
    *              TN_dgemm              *
    *                                    *
    **************************************/
    void TN_dgemm(int ne, int nprj, int npack2, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
       //DGEMM_PWDFT((char *) "T",(char *) "N",ne,nprj,npack2,alpha,host_a,npack2,host_b,npack2,beta,host_c,ne);

     
        //gdevice_TN_dgemm(nn,nprj,ng,rtwo,a,b,rzero,sum);

        //int ia = fetch_dev_mem_indx(((size_t) npack2) * ((size_t) ne));
        //int ib = fetch_dev_mem_indx(((size_t) npack2) * ((size_t) nprj));
        b_prj  = host_b;
        ib_prj[0] = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) nprj));
        if (tile_fac>1) ib_prj[1] = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) nprj));
        int ic = fetch_dev_mem_indx(((size_t) ne)         * ((size_t) nprj));

        NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(ne,nprj,sizeof(double),host_c,ne,dev_mem[ic],ne,stream[0]) );

        if (tile_fac>1)
           NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                                                   &a_psi[tile_start2[0]],npack2,
                                                   dev_mem[ia_psi[0]],tile_npack2[0],stream[0]) );
        NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[0],nprj,sizeof(double),
                                                &b_prj[tile_start2[0]],npack2,
                                                dev_mem[ib_prj[0]],tile_npack2[0],stream[0]) );

        double beta0 = beta;
        for (auto tt=0; tt<tile_fac; ++tt) {
           int ttp1 = tt+1;
           if (ttp1<tile_fac) {
              NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[ttp1],ne,sizeof(double),
                                                      &a_psi[tile_start2[ttp1]],npack2,
                                                      dev_mem[ia_psi[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]) );
              NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[ttp1],nprj,sizeof(double),
                                                      &b_prj[tile_start2[ttp1]],npack2,
                                                      dev_mem[ib_prj[ttp1%2]],tile_npack2[ttp1],stream[ttp1%2]) );
           }
           NWPW_CUDA_ERROR( cudaStreamSynchronize(stream[tt%2]) );
           NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
				       matT,matN,
				       ne,nprj,tile_npack2[tt],&alpha,
				       dev_mem[ia_psi[tt%2]],tile_npack2[tt],
				       dev_mem[ib_prj[tt%2]],tile_npack2[tt],
				       &beta0,dev_mem[ic],ne) );
           beta0 = 1.0;
        }
        cudaMemcpy(host_c,dev_mem[ic],ne*nprj*sizeof(double),cudaMemcpyDeviceToHost);

        //inuse[ia] = false;
        //inuse[ib_prj[0]] = false;
        //if (tile_fac>1) inuse[ib_prj[1]] = false;
        inuse[ic] = false;
    
    }

    void T_free() { inuse[ib_prj[0]] = false; if (tile_fac>1) inuse[ib_prj[1]] = false; }

    void NT_dgemm(int npack2, int ne, int nprj, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
       //DGEMM_PWDFT((char *) "N",(char *) "T",npack2,ne,nprj,alpha,host_a,npack2,host_b,ne,beta,host_c,npack2);

       int ib = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) nprj));

       NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(ne,nprj,sizeof(double),
                                               host_b,ne,
                                               dev_mem[ib],ne,stream[(tile_fac-1)%2]) );
       NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[tile_fac-1],ne,sizeof(double),
                                               &host_c[tile_start2[tile_fac-1]],npack2,
                                               dev_mem[ia_hpsi[(tile_fac-1)%2]],tile_npack2[tile_fac-1],stream[(tile_fac-1)%2]) );
       NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[tile_fac-1],nprj,sizeof(double),
                                               &host_a[tile_start2[tile_fac-1]],npack2,
                                               dev_mem[ib_prj[(tile_fac-1)%2]],tile_npack2[tile_fac-1],stream[(tile_fac-1)%2]) );
       for (auto tt=tile_fac-1; tt>=0; --tt) {
          int ttm1 = tt-1;
          if (ttm1>=0) {
             NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[ttm1],ne,sizeof(double),
                                                     &host_c[tile_start2[ttm1]],npack2,
                                                     dev_mem[ia_hpsi[ttm1%2]],tile_npack2[ttm1],stream[ttm1%2]) );
             NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[ttm1],nprj,sizeof(double),
                                                     &host_a[tile_start2[ttm1]],npack2,
                                                     dev_mem[ib_prj[ttm1%2]], tile_npack2[ttm1],stream[ttm1%2]) );
          }
          NWPW_CUDA_ERROR( cudaStreamSynchronize(stream[tt%2]) );
          NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
 		                         matN,matT,
	                                 tile_npack2[tt],ne,nprj, &alpha,
                                         dev_mem[ib_prj[tt%2]],tile_npack2[tt],
                                         dev_mem[ib],ne,
                                         &beta,dev_mem[ia_hpsi[tt%2]],tile_npack2[tt]) );
          NWPW_CUBLAS_ERROR( cublasGetMatrixAsync(tile_npack2[tt],ne,sizeof(double),
                                                  dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],
                                                  &host_c[tile_start2[tt]],npack2,stream[tt%2]) );
       }

       inuse[ib] = false;
       inuse[ib_prj[0]] = false;
       if (tile_fac>1) inuse[ib_prj[1]] = false;
    }

       
   /************************************** 
    *                                    *
    *              MM6_dgemm             *
    *                                    *
    **************************************/
    void MM6_dgemm(int ne, 
                   double *host_s21, double *host_s12, double *host_s11, 
                   double *host_sa0, double *host_sa1, double *host_st1) {
       double rzero=0.0;
       double rone =1.0;
       int i_s21 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne)); //input
       int i_s12 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne)); //input
       int i_s11 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne)); //input
       int i_sa0 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne)); //input
       int i_st1 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne)); //tmp
       int i_sa1 = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne)); //input-output

       NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(ne,ne,sizeof(double),host_s12,ne,dev_mem[i_s12],ne,stream[0]) );
       NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(ne,ne,sizeof(double),host_s11,ne,dev_mem[i_s11],ne,stream[1]) );
     
       NWPW_CUBLAS_ERROR( cublasSetMatrix(ne,ne,sizeof(double),host_s21,ne,dev_mem[i_s21],ne) );
       NWPW_CUBLAS_ERROR( cublasSetMatrix(ne,ne,sizeof(double),host_sa0,ne,dev_mem[i_sa0],ne) );
       NWPW_CUBLAS_ERROR( cublasSetMatrix(ne,ne,sizeof(double),host_s11,ne,dev_mem[i_s11],ne) );
       NWPW_CUBLAS_ERROR( cublasSetMatrix(ne,ne,sizeof(double),host_sa1,ne,dev_mem[i_sa1],ne) );

       //mmm_Multiply(ms, s21, sa0, 1.0, sa1, 1.0);
       NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
 		                      matN,matN,
	                              ne,ne,ne, &rone,
                                      dev_mem[i_s21],ne,
                                      dev_mem[i_sa0],ne,
                                      &rone,dev_mem[i_sa1],ne) );

       //mmm_Multiply(ms, sa0, s12, 1.0, sa1, 1.0);
       NWPW_CUDA_ERROR( cudaStreamSynchronize(stream[0]) );
       NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
 		                      matN,matN,
	                              ne,ne,ne, &rone,
                                      dev_mem[i_sa0],ne,
                                      dev_mem[i_s12],ne,
                                      &rone,dev_mem[i_sa1],ne) );

       //mmm_Multiply(ms, s11, sa0, 1.0, st1, 0.0);
       NWPW_CUDA_ERROR( cudaStreamSynchronize(stream[1]) );
       NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
 		                      matN,matN,
	                              ne,ne,ne, &rone,
                                      dev_mem[i_s11],ne,
                                      dev_mem[i_sa0],ne,
                                      &rzero,dev_mem[i_st1],ne) );

       //mmm_Multiply(ms, sa0, st1, 1.0, sa1, 1.0);
       NWPW_CUBLAS_ERROR( cublasDgemm(master_handle,
 		                      matN,matN,
	                              ne,ne,ne, &rone,
                                      dev_mem[i_sa0],ne,
                                      dev_mem[i_st1],ne,
                                      &rone,dev_mem[i_sa1],ne) );

       NWPW_CUBLAS_ERROR( cublasGetMatrix(ne,ne,sizeof(double),dev_mem[i_sa1],ne,host_sa1,ne) );

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
    void psi_alloc(int npack1, int ne, int tile_fac0=1) {
        tile_fac        = tile_fac0;

        tile_npack2_max = (((2*npack1)%tile_fac)==0) ? (2*npack1)/tile_fac : (2*npack1)/tile_fac + 1;
        //for (auto i=0; i<tile_fac; ++i) tile_npack2[i] = (i<((2*npack1)%tile_fac)) ? (2*npack1)/tile_fac+1 : (2*npack1)/tile_fac;
        for (auto i=0; i<tile_fac; ++i) tile_npack2[i] = (2*npack1)/tile_fac;
        for (auto i=0; i<((2*npack1)%tile_fac); ++i) tile_npack2[i] += 1;

        tile_start2[0] = 0;
        for (auto i=1; i<tile_fac; ++i) tile_start2[i] = tile_start2[i-1] + tile_npack2[i-1];

        ia_psi[0]  = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) ne));
        ia_hpsi[0] = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) ne));

        if (tile_fac>1) {
           ia_psi[1]  = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) ne));
           ia_hpsi[1] = fetch_dev_mem_indx(((size_t) tile_npack2_max) * ((size_t) ne));
        }
        std::cout << "Into psi_alloc, tile_factor = " << tile_fac << " ndev_mem=" << ndev_mem << std::endl;
    }
    void psi_dealloc() {
        inuse[ia_psi[0]]  = false;
        inuse[ia_hpsi[0]] = false;
        if (tile_fac>1) {
           inuse[ia_psi[1]]  = false;
           inuse[ia_hpsi[1]] = false;
        }
    }
    void psi_copy_host2gpu(int npack1, int ne, double *psi) {
        //cudaMemcpy(dev_mem[ia_psi[0]],psi,tile_npack2_max*ne*sizeof(double),cudaMemcpyHostToDevice);
        a_psi = psi;
        NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[0],ne,sizeof(double),
                                                psi,2*npack1,
                                                dev_mem[ia_psi[0]],tile_npack2[0],stream[0]));
    }
    void hpsi_copy_host2gpu(int npack1, int ne, double *hpsi) {
        //cudaMemcpy(dev_mem[ia_hpsi[0]],hpsi,2*npack1*ne*sizeof(double),cudaMemcpyHostToDevice);
        int tt = tile_fac-1;
        a_hpsi = hpsi;
        NWPW_CUBLAS_ERROR( cublasSetMatrixAsync(tile_npack2[tt],ne,sizeof(double),
                                                &hpsi[tile_start2[tt]],2*npack1,
                                                dev_mem[ia_hpsi[tt%2]],tile_npack2[tt],stream[tt%2]) );
    }
    void psi_copy_gpu2host(int npack1, int ne, double *psi) {
       //cudaMemcpy(psi, dev_mem[ia_psi[0]], 2*ne*npack1*sizeof(double),cudaMemcpyDeviceToHost);
       if (tile_fac==1) {
           NWPW_CUDA_ERROR( cudaStreamSynchronize(stream[0]) );
           NWPW_CUBLAS_ERROR( cublasGetMatrix(tile_npack2[0],ne,sizeof(double),
                                              dev_mem[ia_psi[0]],tile_npack2[0],
                                              psi,2*npack1) );
       }
    }
    void hpsi_copy_gpu2host(int npack1, int ne, double *hpsi) {
       NWPW_CUDA_ERROR( cudaStreamSynchronize(stream[0]) );
       //cudaMemcpy(hpsi, dev_mem[ia_hpsi[0]], 2*ne*npack1*sizeof(double),cudaMemcpyDeviceToHost);
       //if (tile_fac==1) {
       //    NWPW_CUBLAS_ERROR( cublasGetMatrix(tile_npack2[0],ne,sizeof(double),
       //                                       dev_mem[ia_hpsi[0]],tile_npack2[0],
       //                                       hpsi,2*npack1));
       //}
    }


    /******************************/
    /* fft functions (uses cuFFT) */
    /******************************/
    void batch_fft_init(int nx, int ny, int nz, int nq1, int nq2, int nq3) {
    std::cout << "Into batch_fft_init" << std::endl;
        NWPW_CUFFT_ERROR( cufftPlan1d(&forward_plan_x, nx, CUFFT_D2Z, nq1) );
        NWPW_CUFFT_ERROR( cufftPlan1d(&backward_plan_x, nx, CUFFT_Z2D, nq1) );

        int y_inembed[] = {ny};
        int y_onembed[] = {ny};
        NWPW_CUFFT_ERROR( cufftPlanMany(&plan_y, 1, &ny, y_inembed, 1, ny, y_onembed, 1, ny, CUFFT_Z2Z, nq2) );

        int z_inembed[] = {nz};
        int z_onembed[] = {nz};
        NWPW_CUFFT_ERROR( cufftPlanMany(&plan_z, 1, &nz, z_inembed, 1, nz, z_onembed, 1, nz, CUFFT_Z2Z, nq3) );
    }

    void batch_fft_end() {
        // free fft descriptors
        cufftDestroy(forward_plan_x);
        cufftDestroy(plan_y);
        cufftDestroy(plan_z);
        cufftDestroy(backward_plan_x);

        // free dev_mem
        for (auto i=0; i<ndev_mem; ++i) NWPW_CUDA_ERROR( cudaFree(dev_mem[i]) );
        ndev_mem = 0;
    }

    void batch_cfftx(bool forward, int nx, int nq, int n2ft3d, double *a) {
        int ia_dev = fetch_dev_mem_indx(((size_t) n2ft3d));
        cudaMemcpy(dev_mem[ia_dev],a,n2ft3d*sizeof(double),cudaMemcpyHostToDevice);

        if (forward) {
            NWPW_CUFFT_ERROR( cufftExecD2Z(forward_plan_x,
                                           reinterpret_cast<cufftDoubleReal*>(dev_mem[ia_dev]),
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev])) );
        }
        else {
            NWPW_CUFFT_ERROR( cufftExecZ2D(backward_plan_x,
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           reinterpret_cast<cufftDoubleReal*>(dev_mem[ia_dev])) );
        }

        cudaMemcpy(a,dev_mem[ia_dev],n2ft3d*sizeof(double),cudaMemcpyDeviceToHost);

        inuse[ia_dev] = false;
    }

    void batch_cffty(bool forward, int ny,int nq,int n2ft3d, double *a) {
        int ia_dev = fetch_dev_mem_indx(((size_t) n2ft3d));
        cudaMemcpy(dev_mem[ia_dev],a,n2ft3d*sizeof(double),cudaMemcpyHostToDevice);

        if (forward) {
            NWPW_CUFFT_ERROR( cufftExecZ2Z(plan_y,
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           CUFFT_FORWARD) );
        }
        else {
            NWPW_CUFFT_ERROR( cufftExecZ2Z(plan_y,
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           CUFFT_INVERSE) );
        }

        cudaMemcpy(a,dev_mem[ia_dev],n2ft3d*sizeof(double),cudaMemcpyDeviceToHost);

        inuse[ia_dev] = false;
    }

    void batch_cfftz(bool forward, int nz,int nq,int n2ft3d, double *a) {
        int ia_dev = fetch_dev_mem_indx(((size_t) n2ft3d));
        cudaMemcpy(dev_mem[ia_dev],a,n2ft3d*sizeof(double),cudaMemcpyHostToDevice);

        if (forward) {
            NWPW_CUFFT_ERROR( cufftExecZ2Z(plan_z,
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           CUFFT_FORWARD) );
        }
        else {
            NWPW_CUFFT_ERROR( cufftExecZ2Z(plan_z,
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           reinterpret_cast<cufftDoubleComplex*>(dev_mem[ia_dev]),
                                           CUFFT_INVERSE) );
        }

        cudaMemcpy(a,dev_mem[ia_dev],n2ft3d*sizeof(double),cudaMemcpyDeviceToHost);

        inuse[ia_dev] = false;
    }


    // routines below need to be made into sycl or removed

static void eigsrt_device(double *D, double *V, int n) {
   int i,j,k;
   double p;

   for (i=0; i<(n-1); ++i)
   {
      k = i;
      p = D[i];
      for(j=i+1; j<n; ++j)
         if (D[j]>=p)
         {
            k = j;
            p = D[j];
         }

      if (k!=i)
      {
         D[k] = D[i];
         D[i] = p;
         for (j=0; j<n; ++j)
         {
            p = V[j+i*n];
            V[j+i*n] = V[j+k*n];
            V[j+k*n] = p;
         }
      }
   }
}


   void NN_eigensolver(int ispin, int ne[], double *host_hml, double *host_eig) {
      int i_a1[ispin],i_w1[ispin],info[ispin];
      int shift1 = 0;
      int shift2 = 0;
      for (auto ms=0; ms<ispin; ++ms)
      {
         int nn = ne[ms]*ne[ms];
         i_a1[ms] = fetch_dev_mem_indx(((size_t) ne[ms]) * ((size_t) ne[ms])); //input-output
         i_w1[ms] = fetch_dev_mem_indx(((size_t) ne[ms]) );                    //output
         NWPW_CUDA_ERROR(cudaMemcpyAsync(dev_mem[i_a1[ms]],host_hml+shift2, nn*sizeof(double),cudaMemcpyHostToDevice,stream[ms]));
         NWPW_CUDA_ERROR(cudaMemcpyAsync(dev_mem[i_w1[ms]],host_eig+shift1, nn*sizeof(double),cudaMemcpyHostToDevice,stream[ms]));
         shift1 += ne[0];
         shift2 += ne[0]*ne[0];
      }

      // allocate work space for syevd
      if (lwork==0)
      {
         // query working space of syevd
         CUSOLVER_CHECK(cusolverDnDsyevd_bufferSize(cusolverH, jobz, uplo, ne[0], dev_mem[i_a1[0]],n,dev_mem[i_w1[0]],&lwork));
         NWPW_CUDA_ERROR(cudaMalloc(reinterpret_cast<void **>(&d_work),sizeof(double) * lwork));
      }

      shift1 = 0;
      shift2 = 0;
      for (auto ms=0; ms<ispin; ++ms)
      {
         NWPW_CUDA_ERROR( cudaStreamSynchronize(stream[ms]) );

         // compute spectrum
         n = ne[ms];
         CUSOLVER_CHECK(cusolverDnDsyevd(cusolverH,jobz,uplo,n,dev_mem[i_a1[ms]],n,dev_mem[i_w1[ms]],d_work,lwork,d_info[ms]));

        NWPW_CUDA_ERROR(cudaMemcpyAsync(host_hml+shift2,dev_mem[i_a1[ms]],nn*sizeof(double),cudaMemcpyDeviceToHost,stream[ms]));
        NWPW_CUDA_ERROR(cudaMemcpyAsync(host_eig+shift2,dev_mem[i_w1[ms]],nn*sizeof(double),cudaMemcpyDeviceToHost,stream[ms]));
        NWPW_CUDA_ERROR(cudaMemcpyAsync(info+ms,d_info[ms],sizeof(int),cudaMemcpyDeviceToHost,stream[ms]));

        shift1 += ne[0];
        shift2 += ne[0]*ne[0];
      }
      for (auto ms=0; ms<ispin; ++ms)
         NWPW_CUDA_ERROR(cudaStreamSynchronize(stream[ms]));
   }


};

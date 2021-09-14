#ifndef _GDEVICE_HPP_
#define _GDEVICE_HPP_

namespace pwdft {
using namespace pwdft;

extern void gdevice_TN3_dgemm(int, int, double, double *, double *, double, double *, double *, double *);
extern void gdevice_TN_dgemm(int, int, int, double, double *, double *, double, double *);
extern void gdevice_NN_dgemm(int, int, double, double *, double *, double, double *);
extern void gdevice_NT_dgemm(int, int, int, double, double *, double *, double, double *);

extern void gdevice_psi_alloc(int, int);
extern void gdevice_psi_dealloc();

extern void gdevice_psi_copy_host2gpu(int, int, double *);
extern void gdevice_hpsi_copy_host2gpu(int, int, double *);

extern void gdevice_psi_copy_gpu2host(int, int, double *);
extern void gdevice_hpsi_copy_gpu2host(int, int, double *);

extern void gdevice_batch_fft_init(int, int, int, int, int, int);
extern void gdevice_batch_cfftx(bool,int,int,int,double *);
extern void gdevice_batch_cffty(bool,int,int,int,double *);
extern void gdevice_batch_cfftz(bool,int,int,int,double *);


#ifdef NWPW_SYCL

#include <CL/sycl.hpp>

#define NWPW_SYCL_ERROR( EXPR )                                         \
    try {                                                               \
        EXPR;                                                           \
    }                                                                   \
    catch (sycl::exception const &ex) {                                 \
        std::stringstream msg;                                          \
        msg << "SYCL Exception at " << __FILE__ << " : " << __LINE__    \
            << std::endl;                                               \
        throw(std::runtime_error( ex.what() ));                         \
    }                                                                   \
    catch(std::runtime_error const& ex) {                               \
        std::stringstream msg;                                          \
        msg << "runtime Exception at " << __FILE__ << " : " << __LINE__ \
            << std::endl;                                               \
        throw(std::runtime_error( ex.what() ));                         \
    }

extern sycl::queue* get_syclQue();

#endif // NWPW_SYCL

#ifdef NWPW_CUDA

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cufft.h>

#include <stdexcept>
#include <string>
#include <sstream>

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

#endif // NWPW_CUDA

}

#endif

#ifndef _GDEVICES_HPP_
#define _GDEVICES_HPP_

#ifdef _NWPW_SYCL_

/* can place sycl mkl code here */
#include        <cstdio>
#include        <iostream>
#include        <limits>
#include        <CL/sycl.hpp>
#include        "mkl_blas_sycl.hpp"


class Gdevices {

   oneapi::mkl::transpose matT = oneapi::mkl::transpose::trans;
   oneapi::mkl::transpose matN = oneapi::mkl::transpose::nontrans;

    auto asyncHandler = [&](cl::sycl::exception_list eL) {
       for (auto& e : eL) {
         try {
           std::rethrow_exception(e);
         } catch (cl::sycl::exception& e) {
           std::cout << e.what() << std::endl;
           std::cout << "fail" << std::endl;
           std::terminate();
         }
       }
    };
    cl::sycl::gpu_selector device_selector;
    cl::sycl::queue device_queue(device_selector,
                                 asyncHandler,
                                 cl::sycl::property_list{cl::sycl::property::queue::in_order{}});
    /* device memory */
    int    ndev_mem;
    bool   inuse[25];
    size_t nsize_mem[25];
    double *dev_mem[25];

public:

    Gdevice() { ndev_mem = 0; }

    ~Gdevice() {
       for (auto i=0; i<ndev_mem; ++i)
          cl::sycl::free(dev_mem[i], device_queue);
     }

     int fetch_dev_mem_indx(const size_t ndsize) {
        ii = 0;
        while ((((ndsize!=ndsize_mem[ii]) || inuse[ii])) && (ii<ndev_mem)) 
          ++ii;

        if (ii<ndev_mem) {
           inuse[ii] = true;
        } else {
           ii            = ndev_mem;
           inuse[ii]     = true;
           nsize_mem[ii] = ndsize;
           dev_mem[ii]   = cl::sycl::malloc_device<double>(ndsize,device_queue);
           ndev_mem += 1
        }

        return ii;
     }

     void TN3_dgemm(const int npack, const int ne, double alpha, const double *host_a, const double *host_b, const double beta,
                    double *host_caa, double *host_cab, double *host_bb) {
        int ia = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne));
        int ib = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne));
        int icaa = fetch_dev_mem_indx(((size_t) ne) * ((size_t) ne));
        int icab = fetch_dev_mem_indx(((size_t) ne) * ((size_t) ne));
        int icbb = fetch_dev_mem_indx(((size_t) ne) * ((size_t) ne));
        try {
           device_queue.submit([&](cl::sycl::handler& cgh) { cgh.memcpy(dev_mem[ia], host_a, npack*ne*sizeof(double)); });
           device_queue.submit([&](cl::sycl::handler& cgh) { cgh.memcpy(dev_mem[ib], host_b, npack*ne*sizeof(double)); });

           oneapi::mkl::blas::gemm(device_queue,matT,matN,ne,ne,npack,alpha,dev_mem[ia],npack,dev_mem[ia],npack,beta,dev_mem[icaa],ne);
           oneapi::mkl::blas::gemm(device_queue,matT,matN,ne,ne,npack,alpha,dev_mem[ia],npack,dev_mem[ib],npack,beta,dev_mem[icab],ne);
           oneapi::mkl::blas::gemm(device_queue,matT,matN,ne,ne,npack,alpha,dev_mem[ib],npack,dev_mem[ib],npack,beta,dev_mem[icbb],ne);
           device_queue.memcpy(host_caa, dev_mem[icaa], ne*ne*sizeof(double));
           device_queue.memcpy(host_cab, dev_mem[icab], ne*ne*sizeof(double));
           device_queue.memcpy(host_cbb, dev_mem[icbb], ne*ne*sizeof(double));
           device_queue.wait();
        }
        catch(cl::sycl::exception const& e) {
           std::cout << "\t\tSYCL exception during GEMM\n" << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
        }

        inuse[ia] = false;
        inuse[ib] = false;
        inuse[icaa] = false;
        inuse[icab] = false;
        inuse[icbb] = false;
     }

     void TN_dgemm(const int npack, const int ne, double alpha, const double *host_a, const double *host_b, const double beta, double *host_c) {
        int ia = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne));
        int ib = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne));
        int ic = fetch_dev_mem_indx(((size_t) ne) * ((size_t) ne));
        try {
           device_queue.submit([&](cl::sycl::handler& cgh) { cgh.memcpy(dev_mem[ia], host_a, npack*ne*sizeof(double)); });
           device_queue.submit([&](cl::sycl::handler& cgh) { cgh.memcpy(dev_mem[ib], host_b, npack*ne*sizeof(double)); });

           oneapi::mkl::blas::gemm(device_queue,matT,matN,ne,ne,npack,alpha,dev_mem[ia],npack,dev_mem[ib],npack,beta,dev_mem[ic],ne);
           device_queue.memcpy(host_c, dev_mem[ic], ne*ne*sizeof(double));
           device_queue.wait();
        }
        catch(cl::sycl::exception const& e) {
           std::cout << "\t\tSYCL exception during GEMM\n" << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
        }

        inuse[ia] = false;
        inuse[ib] = false;
        inuse[ic] = false;
     }

     void NN_dgemm(const int npack, const int ne, double alpha, const double *host_a, const double *host_b, const double beta, double *host_c) {
        int ia = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne));
        int ib = fetch_dev_mem_indx(((size_t) ne) * ((size_t) ne));
        int ic = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne));
        try {
           device_queue.submit([&](cl::sycl::handler& cgh) { cgh.memcpy(dev_mem[ia], host_a, npack*ne*sizeof(double)); });
           device_queue.submit([&](cl::sycl::handler& cgh) { cgh.memcpy(dev_mem[ib], host_b, ne*ne*sizeof(double)); });

           oneapi::mkl::blas::gemm(device_queue, matN, matN, npack,ne,ne, alpha, dev_mem[ia], npack, dev_mem[ib], ne, beta, dev_mem[ic], npack);
           device_queue.memcpy(host_c, dev_mem[ic], npack*ne*sizeof(double));
           device_queue.wait();
        }
        catch(cl::sycl::exception const& e) {
           std::cout << "\t\tSYCL exception during GEMM\n" << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
        }

        inuse[ia] = false;
        inuse[ib] = false;
        inuse[ic] = false;
     }
};


#elif defined _NWPW_OPENCL_

/* can place opencl code from mac here */
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_SILENCE_DEPRECATION
#ifdef __APPLE__
#include <OpenCL/opencl.h>

#else
#include <CL/cl.h>
#endif


#define program1	"#pragma OPENCL EXTENSION cl_khr_fp64 : enable \n\n\
__kernel void NNmatmul(const int M, const int N, const int K,\n\
                     const __global double *A, \n\
                     const __global double *B, \n\
                     __global double *C) { \n\n\
    // Get the index of the current element  \n\
    int i = get_global_id(0); \n\
    int j = get_global_id(1); \n\n\
    // Do the operation  \n\
    double acc = 0.0;  \n\
    for (int l=0; l<K; l++) {  \n\
       acc += A[i + l*M]*B[l + j*K];  \n\
    }   \n\
    C[i+j*M] = acc; \n }\n"

#define program2	"#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n\
__kernel void TNmatmul(const int M, const int N, const int K,\n\
                     const __global double *A, \n\
                     const __global double *B, \n\
                     __global double *C) {\n\n\
    // Get the index of the current element\n\
    int i = get_global_id(0);\n\
    int j = get_global_id(1);\n\n\
    // Do the operation \n\
    double acc = 0.0; \n\
    for (int l=0; l<K; l++) { \n\
       acc += A[l + i*M]*B[l + j*K]; \n\
    } \n\
    C[i+j*M] = acc; \n}\n"



#else
/* standard host code */

#include        "blas.h"

class Gdevices {

public:

     void TN3_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_caa, double *host_cab, double *host_cbb)
     {
        DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_a,npack,host_a,npack,beta,host_caa,ne);
        DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_a,npack,host_b,npack,beta,host_cab,ne);
        DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_b,npack,host_b,npack,beta,host_cbb,ne);
     }

     void TN_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
        DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_a,npack,host_b,npack,beta,host_c,ne);
     }

     void NN_dgemm(int npack, int ne, double alpha, double *host_a, double *host_c, double beta, double *host_b) {
        DGEMM_PWDFT((char *) "N",(char *) "N",npack,ne,ne,alpha,host_a,npack,host_c,ne,beta,host_b,npack);
     }

};

#endif

#endif

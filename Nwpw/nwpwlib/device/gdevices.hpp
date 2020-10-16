#ifndef _GDEVICES_HPP_
#define _GDEVICES_HPP_

#ifdef NWPW_SYCL

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


#elif defined NWPW_OPENCL

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

#include        <iostream>
#include        "blas.h"

typedef struct {
   cl_uint num_platforms;
   cl_platform_id * platform_id;
   cl_uint * num_devices;
   cl_device_id ** device_id;
   cl_bool   ** avail;
   cl_bool   ** has_cl_khr_fp64;
   cl_uint   ** num_cores;
   cl_uint   ** freq;
   cl_uint   ** wdouble;
   cl_uint   ** wfloat;
   cl_ulong  ** mem;
   //Device_Type ** device;
} NWPW_GPU_Type;

class Gdevices {

   NWPW_GPU_Type gpu;

   int plat_indx,device_indx;
   cl_device_id     device_id_selected;
   cl_context       context;
   cl_command_queue command_queue;

public:
     Gdevices() {
        size_t size;
        char str[1000];

        // Get platforms 
        cl_int ret = clGetPlatformIDs(0, NULL, &(gpu.num_platforms));
        gpu.platform_id = (cl_platform_id *) malloc(sizeof(cl_platform_id)*gpu.num_platforms);
        ret = clGetPlatformIDs(gpu.num_platforms,gpu.platform_id,NULL);

        gpu.num_devices = (cl_uint *) malloc(sizeof(cl_uint)*gpu.num_platforms);
        gpu.device_id = (cl_device_id **) malloc(sizeof(cl_device_id *)*gpu.num_platforms);
        gpu.avail     = (cl_bool **) malloc(sizeof(cl_bool *)*gpu.num_platforms);
        gpu.has_cl_khr_fp64 = (cl_bool **) malloc(sizeof(cl_bool *)*gpu.num_platforms);
        gpu.num_cores = (cl_uint **) malloc(sizeof(cl_uint *)*gpu.num_platforms);
        gpu.freq      = (cl_uint **) malloc(sizeof(cl_uint *)*gpu.num_platforms);
        gpu.wdouble   = (cl_uint **) malloc(sizeof(cl_uint *)*gpu.num_platforms);
        gpu.wfloat    = (cl_uint **) malloc(sizeof(cl_uint *)*gpu.num_platforms);
        gpu.mem       = (cl_ulong **) malloc(sizeof(cl_ulong *)*gpu.num_platforms);
        for (cl_uint i=0; i<gpu.num_platforms; ++i)
        {
           ret = clGetDeviceIDs(gpu.platform_id[i], CL_DEVICE_TYPE_ALL, 0,NULL,&(gpu.num_devices[i]));

           gpu.device_id[i] = (cl_device_id *) malloc(sizeof(cl_device_id)*gpu.num_devices[i]);
           gpu.avail[i]     = (cl_bool *) malloc(sizeof(cl_bool)*gpu.num_devices[i]);
           gpu.has_cl_khr_fp64[i] = (cl_bool *) malloc(sizeof(cl_bool)*gpu.num_devices[i]);
           gpu.num_cores[i] = (cl_uint *) malloc(sizeof(cl_uint)*gpu.num_devices[i]);
           gpu.freq[i]      = (cl_uint *) malloc(sizeof(cl_uint)*gpu.num_devices[i]);
           gpu.wdouble[i]   = (cl_uint *) malloc(sizeof(cl_uint)*gpu.num_devices[i]);
           gpu.wfloat[i]    = (cl_uint *) malloc(sizeof(cl_uint)*gpu.num_devices[i]);
           gpu.mem[i]       = (cl_ulong *) malloc(sizeof(cl_ulong)*gpu.num_devices[i]);

           ret = clGetDeviceIDs(gpu.platform_id[i],CL_DEVICE_TYPE_ALL,gpu.num_devices[i],gpu.device_id[i],NULL);

           for (cl_uint j=0; j<gpu.num_devices[i]; ++j)
           {
              ret = clGetDeviceInfo(gpu.device_id[i][j],CL_DEVICE_AVAILABLE,sizeof(cl_bool),&(gpu.avail[i][j]),&size);
              ret = clGetDeviceInfo(gpu.device_id[i][j],CL_DEVICE_MAX_COMPUTE_UNITS,sizeof(cl_uint),&(gpu.num_cores[i][j]),&size);
              ret = clGetDeviceInfo(gpu.device_id[i][j],CL_DEVICE_MAX_CLOCK_FREQUENCY,sizeof(cl_uint),&(gpu.freq[i][j]),&size);
              ret = clGetDeviceInfo(gpu.device_id[i][j],CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE,sizeof(cl_uint),&(gpu.wdouble[i][j]),&size);
              ret = clGetDeviceInfo(gpu.device_id[i][j],CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT,sizeof(cl_uint),&(gpu.wfloat[i][j]),&size);
              ret = clGetDeviceInfo(gpu.device_id[i][j],CL_DEVICE_GLOBAL_MEM_SIZE,sizeof(cl_ulong),&(gpu.mem[i][j]),&size);
              ret = clGetDeviceInfo(gpu.device_id[i][j],CL_DEVICE_EXTENSIONS,1000*sizeof(char),str,&size);
              gpu.has_cl_khr_fp64[i][j] = (strstr(str, "cl_khr_fp64") != NULL);
           }
        }
        printf("Number of platforms = %d\n",gpu.num_platforms);
        plat_indx   = 0;
        device_indx = 0;
        for (int i=0; i<gpu.num_platforms; ++i)
        {
           printf(" - %d patform_id= %ld num_devices= %d\n",i, gpu.platform_id[i],gpu.num_devices[i]);
           for (int j=0; j<gpu.num_devices[i]; ++j)
           {
              printf("   -- %d device_id= %ld num_cores=%3d mem=%12ld  %4d MHz wfloat=%d wdouble=%d avail=%d has_cl_khr_fp64=%d\n",j,gpu.device_id[i][j],
                                                                gpu.num_cores[i][j],
                                                                (long) gpu.mem[i][j],
                                                                gpu.freq[i][j],
                                                                gpu.wfloat[i][j],
                                                                gpu.wdouble[i][j],
                                                                gpu.avail[i][j],
                                                                gpu.has_cl_khr_fp64[i][j]);
              if (gpu.avail[i][j] && gpu.wdouble[i][j] && gpu.has_cl_khr_fp64[i][j])
              {
                 plat_indx = i;
                 device_indx = j;
              }
           }
        }

        device_id_selected = gpu.device_id[plat_indx][device_indx];
        ret = clGetDeviceInfo(device_id_selected,CL_DEVICE_VENDOR,1000*sizeof(char),str,&size);
        printf("\n - Using platform_id=%ld device_id=%ld vendor=%s num_cores=%3d mem=%12ld  %4d MHz wfloat=%d wdouble=%d avail=%d has_cl_khr_fp64=%d\n",
               gpu.platform_id[plat_indx],gpu.device_id[plat_indx][device_indx],str,
                                                                gpu.num_cores[plat_indx][device_indx],
                                                                (long) gpu.mem[plat_indx][device_indx],
                                                                gpu.freq[plat_indx][device_indx],
                                                                gpu.wfloat[plat_indx][device_indx],
                                                                gpu.wdouble[plat_indx][device_indx],
                                                                gpu.avail[plat_indx][device_indx],
                                                                gpu.has_cl_khr_fp64[plat_indx][device_indx]);


        // Create an OpenCL context
        context = clCreateContext(NULL,1, &(device_id_selected), NULL, NULL, &ret);

        // Create a command queue
        command_queue = clCreateCommandQueue(context, device_id_selected, 0, &ret);


     }

     ~Gdevices() {
        cl_int ret = clReleaseCommandQueue(command_queue);
               ret = clReleaseContext(context);

        for (int i=0; i<gpu.num_platforms; ++i)
        {
           free(gpu.device_id[i]);
           free(gpu.avail[i]);
           free(gpu.has_cl_khr_fp64[i]);
           free(gpu.num_cores[i]);
           free(gpu.freq[i]);
           free(gpu.wdouble[i]);
           free(gpu.wfloat[i]);
           free(gpu.mem[i]);
        }
        free(gpu.platform_id);
        free(gpu.num_devices);
        free(gpu.device_id);
        free(gpu.avail);
        free(gpu.has_cl_khr_fp64);
        free(gpu.num_cores);
        free(gpu.freq);
        free(gpu.wdouble);
        free(gpu.wfloat);
        free(gpu.mem);

     }

     void TN3_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_caa, double *host_cab, double *host_cbb)
     {
        int one = 1;
        int shift1  = 0;
        int mshift1 = 0;

        for (auto k=1; k<=ne; ++k)
        {
           DGEMM_PWDFT((char *) "T",(char *) "N",k,one,npack,alpha,host_a,npack,&host_a[shift1],npack,beta,&host_caa[mshift1],k);
           DGEMM_PWDFT((char *) "T",(char *) "N",k,one,npack,alpha,host_a,npack,&host_b[shift1],npack,beta,&host_cab[mshift1],k);
           DGEMM_PWDFT((char *) "T",(char *) "N",k,one,npack,alpha,host_b,npack,&host_b[shift1],npack,beta,&host_cbb[mshift1],k);
           shift1  += npack;
           mshift1 += ne;
        }

        //DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_a,npack,host_a,npack,beta,host_caa,ne);
        //DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_a,npack,host_b,npack,beta,host_cab,ne);
        //DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_b,npack,host_b,npack,beta,host_cbb,ne);
        //std::cout << "In the OpenCL branch!!!" << std::endl;
        //std::cout << "program1=" << program1 << std::endl;

     }

     void TN_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
        DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_a,npack,host_b,npack,beta,host_c,ne);
     }

     void NN_dgemm(int npack, int ne, double alpha, double *host_a, double *host_c, double beta, double *host_b) {
        DGEMM_PWDFT((char *) "N",(char *) "N",npack,ne,ne,alpha,host_a,npack,host_c,ne,beta,host_b,npack);
     }

};



#else
/* standard host code */

#include        "blas.h"

class Gdevices {

public:

     void TN3_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_caa, double *host_cab, double *host_cbb)
     {
        int one = 1;
        int shift1  = 0;
        int mshift1 = 0;

        for (auto k=1; k<=ne; ++k)
        {
           DGEMM_PWDFT((char *) "T",(char *) "N",k,one,npack,alpha,host_a,npack,&host_a[shift1],npack,beta,&host_caa[mshift1],k);
           DGEMM_PWDFT((char *) "T",(char *) "N",k,one,npack,alpha,host_a,npack,&host_b[shift1],npack,beta,&host_cab[mshift1],k);
           DGEMM_PWDFT((char *) "T",(char *) "N",k,one,npack,alpha,host_b,npack,&host_b[shift1],npack,beta,&host_cbb[mshift1],k);
           shift1  += npack;
           mshift1 += ne;
        }

        //DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_a,npack,host_a,npack,beta,host_caa,ne);
        //DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_a,npack,host_b,npack,beta,host_cab,ne);
        //DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_b,npack,host_b,npack,beta,host_cbb,ne);
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

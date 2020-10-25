#ifndef _GDEVICES_HPP_
#define _GDEVICES_HPP_

#ifdef NWPW_SYCL
#pragma once

/* can place sycl mkl code here */
#include        <cstdio>
#include        <iostream>
#include        <limits>
#include        <CL/sycl.hpp>
#include        <oneapi/mkl.hpp>


class Gdevices {

   oneapi::mkl::transpose matT = oneapi::mkl::transpose::trans;
   oneapi::mkl::transpose matN = oneapi::mkl::transpose::nontrans;

    /* device memory */
    int    ndev_mem;
    bool   inuse[25];
    size_t ndsize_mem[25];
    double *dev_mem[25];

public:
     std::vector<cl::sycl::queue*> syclQueues;  // TODO: queues per device
     cl::sycl::queue* device_queue = nullptr; // default SYCL queue for now

     Gdevices();

    ~Gdevices() {
       for (auto i=0; i<ndev_mem; ++i)
          cl::sycl::free(dev_mem[i], *device_queue);
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
           dev_mem[ii]   = cl::sycl::malloc_device<double>(ndsize, *device_queue);
           ndev_mem += 1;
        }

        return ii;
     }

     void TN3_dgemm(const int npack, const int ne, double alpha, const
		    double *host_a, const double *host_b, const double beta,
                    double *host_caa, double *host_cab, double *host_cbb) {
        int ia = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne));
        int ib = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne));
        int icaa = fetch_dev_mem_indx(((size_t) ne) * ((size_t) ne));
        int icab = fetch_dev_mem_indx(((size_t) ne) * ((size_t) ne));
        int icbb = fetch_dev_mem_indx(((size_t) ne) * ((size_t) ne));
        try {
           device_queue->memcpy(dev_mem[ia], host_a, npack*ne*sizeof(double));
           device_queue->memcpy(dev_mem[ib], host_b, npack*ne*sizeof(double));

           oneapi::mkl::blas::gemm(*device_queue, matT, matN, ne, ne, npack, alpha, dev_mem[ia], npack, dev_mem[ia], npack, beta, dev_mem[icaa], ne);
           oneapi::mkl::blas::gemm(*device_queue, matT, matN, ne, ne, npack, alpha, dev_mem[ia], npack, dev_mem[ib], npack, beta, dev_mem[icab], ne);
           oneapi::mkl::blas::gemm(*device_queue, matT, matN, ne, ne, npack, alpha, dev_mem[ib], npack, dev_mem[ib], npack, beta, dev_mem[icbb], ne);

           device_queue->memcpy(host_caa, dev_mem[icaa], ne*ne*sizeof(double));
           device_queue->memcpy(host_cab, dev_mem[icab], ne*ne*sizeof(double));
           device_queue->memcpy(host_cbb, dev_mem[icbb], ne*ne*sizeof(double));
           device_queue->wait();
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

     void TN_dgemm(const int npack, const int ne, const int nprj, double alpha, const double *host_a, const double *host_b, const double beta, double *host_c) {
        int ia = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne));
        int ib = fetch_dev_mem_indx(((size_t) npack) * ((size_t) nprj));
        int ic = fetch_dev_mem_indx(((size_t) ne) * ((size_t) nprj));
        try {
           device_queue->memcpy(dev_mem[ia], host_a, npack*ne*sizeof(double));
           device_queue->memcpy(dev_mem[ib], host_b, npack*nprj*sizeof(double));

           oneapi::mkl::blas::gemm(*device_queue,matT,matN,ne,nprj,npack,alpha,dev_mem[ia],npack,dev_mem[ib],npack,beta,dev_mem[ic],ne);
           device_queue->memcpy(host_c, dev_mem[ic], ne*nprj*sizeof(double));
           device_queue->wait();
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
           device_queue->memcpy(dev_mem[ia], host_a, npack*ne*sizeof(double));
           device_queue->memcpy(dev_mem[ib], host_b, ne*ne*sizeof(double));

           oneapi::mkl::blas::gemm(*device_queue, matN, matN, npack,ne,ne, alpha, dev_mem[ia], npack, dev_mem[ib], ne, beta, dev_mem[ic], npack);
           device_queue->memcpy(host_c, dev_mem[ic], npack*ne*sizeof(double));
           device_queue->wait();
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

#define MAX_SOURCE_SIZE (0x100000)

#define NNmatmul_src	"#pragma OPENCL EXTENSION cl_khr_fp64 : enable \n\n\
__kernel void NNmatmul(const int M, const int N, const int K, \n\
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
    C[i+j*M] = acc; \n }\n \0\0"

#define TN3matmul_src	"#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n__kernel void TN3matmul(const int M, const int N,\n                     const __global double *A, \n                     const __global double *B, \n                     __global double *Caa) {\n    \n    // Get the index of the current element\n    int i = get_global_id(0);\n    int j = get_global_id(1);\n\n    // Do the operation\n    int NN = N*N;\n     double aa_acc = 0.0;\n    double ab_acc = 0.0;\n    double bb_acc = 0.0;\n    for (int l=0; l<M; ++l) {\n       aa_acc += A[l + i*M]*A[l + j*M];\n       ab_acc += A[l + i*M]*B[l + j*M];\n       bb_acc += B[l + i*M]*B[l + j*M];\n    }\n    Caa[i+j*N] = aa_acc;\n    Caa[i+j*N+NN] = ab_acc;\n    Caa[i+j*N+NN+NN] = bb_acc;\n}\n"

//#define TN3matmul_src	"#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n__kernel void TN3matmul(const int M, const int N,\n                     const __global double *A, \n                     const __global double *B, \n                     __global double *Caa) {\n    \n    // Get the index of the current element\n    int i = get_global_id(0);\n    int j = get_global_id(1);\n\n    // Do the operation\n    int NN = N*N;\n    double acc[3] = {0.0, 0.0, 0.0};\n    for (int l=0; l<M; ++l) {\n       acc[0] += A[l + i*M]*A[l + j*M];\n       acc[1] += A[l + i*M]*B[l + j*M];\n       acc[2] += B[l + i*M]*B[l + j*M];\n    }\n    Caa[i+j*N]       = acc[0];\n    Caa[i+j*N+NN]    = acc[1];\n    Caa[i+j*N+NN+NN] = acc[2];\n}"


#define NTmatmul_src	"#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n__kernel void NTmatmul(const int M, const int N, const int K,\n                     const __global double *A, \n                     const __global double *B, \n                     __global double *C) {\n    \n    // A(npack,nprj),  B(ne,nprj), C(npack,ne)\n\n    // Get the index of the current element, M=npack, N=ne, K=nprj\n    int i = get_global_id(0);\n    int j = get_global_id(1);\n\n    // Do the operation\n    //int la = 0; int lb = 0;\n    double acc = 0.0;\n    for (int l=0; l<K; ++l) {\n       acc += A[i+l*M]*B[j+l*N];\n    }\n    C[i+j*M] = acc;\n}"


#define TNmatmul_src	"#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n  \n__kernel void TNmatmul(const int M, const int N, const int K,\n                     const __global double *A,\n                     const __global double *B,\n                     __global double *C) {\n\n    //A(npack,ne), B(npack,nprj), C(ne,nprj), M=ne, N=nprj, K=npack\n    // Get the index of the current element\n    int i = get_global_id(0);\n    int j = get_global_id(1);\n\n    // Do the operation\n    double acc = 0.0;\n    for (int l=0; l<K; ++l)\n       acc += A[l + i*K]*B[l + j*K];\n    C[i+j*M] = acc;\n}"


#define Generate_projectors_src	"#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n__kernel void Generate_projectors(const int ii, const int ng0, const int nprj, const int nprjall, \n\t\t\t\t  const int nx, const int ny, const int nz,\n                                  const __global int    *indxi,\n                                  const __global int    *indxj,\n                                  const __global int    *indxk,\n                                  const __global double *phfacx,\n                                  const __global double *phfacy,\n                                  const __global double *phfacz,\n                                  const __global int    *sdfunction, \n                                  const __global double *vnl,\n                                  __global double *prj) {\n   int shftx = 2*ii*nx;\n   int shfty = 2*ii*ny;\n   int shftz = 2*ii*ny;\n   int ng    = 2*ng0;\n\n   int i = get_global_id(0); //ng\n   //int l = get_global_id(1); //nprj\n   \n   double ai = phfacx[shftx+2*indxi[i]]; double bi = phfacx[shftx+2*indxi[i]+1];\n   double aj = phfacy[shfty+2*indxj[i]]; double bj = phfacy[shfty+2*indxj[i]+1];\n   double ak = phfacz[shftz+2*indxk[i]]; double bk = phfacz[shftz+2*indxk[i]+1];\n   double c  = aj*ak - bj*bk;\n   double d  = aj*bk + ak*bj;\n   double rexi = (ai*c - bi*d);\n   double iexi = (ai*d + bi*c);\n\n   for (int l=0; l<nprj; ++l)\n   {\n      if (sdfunction[l])\n      {\n         prj[2*i   + (l+nprjall)*ng] = rexi * vnl[i+l*ng0];\n         prj[2*i+1 + (l+nprjall)*ng] = iexi * vnl[i+l*ng0];\n      } else { \n         prj[2*i   + (l+nprjall)*ng] = -iexi * vnl[i+l*ng0];\n         prj[2*i+1 + (l+nprjall)*ng] =  rexi * vnl[i+l*ng0];\n      }\n   }\n}"


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
   cl_program       NNmatmul_program,TN3matmul_program,NTmatmul_program,TNmatmul_program,Generate_projectors_program;
   cl_kernel        NNmatmul_kernel, TN3matmul_kernel, NTmatmul_kernel, TNmatmul_kernel, Generate_projectors_kernel;

   /* device memory */
   int    ndev_mem = 0;
   int    rw_mem[25]; //01 =1- read,10=2 - write, 11=3 - read/write
   bool   inuse[25];
   size_t ndsize_mem[25];
   cl_mem dev_mem[25];

   /* tmp memory */
   int    ntmp_mem=0;
   bool   tmpinuse[25];
   size_t tmpndsize_mem[25];
   double *tmp_mem[25];

   /* phafac memory */
   int nx_pf,ny_pf,nz_pf,npack_pf;
   int indxi_dev,indxj_dev,indxk_dev,exi_dev,exj_dev,exk_dev;

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
        context = clCreateContext(NULL,1, &(gpu.device_id[plat_indx][device_indx]), NULL, NULL, &ret); std::cout << " retcontex=" << ret << " context=" << context;


        // Create a command queue
        command_queue = clCreateCommandQueue(context, gpu.device_id[plat_indx][device_indx], 0, &ret); std::cout << " retcommand=" << ret << std::endl;

        // Create programs from the kernel source

        // Build the NNmatmul program
        char *source_str = (char*)malloc(MAX_SOURCE_SIZE);
        strcpy(source_str,NNmatmul_src);
        size_t source_size = strlen(source_str);

        NNmatmul_program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret); //std::cout << " retcreateprog=" << ret;
        ret = clBuildProgram(NNmatmul_program, 1, &(gpu.device_id[plat_indx][device_indx]), NULL, NULL, NULL);                    //std::cout << " retbuild=" << ret;

        size_t logSize;
        clGetProgramBuildInfo(NNmatmul_program, gpu.device_id[plat_indx][device_indx], CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
        char* messages = (char*)malloc((1+logSize)*sizeof(char));
        clGetProgramBuildInfo(NNmatmul_program, gpu.device_id[plat_indx][device_indx], CL_PROGRAM_BUILD_LOG, logSize, messages, NULL);
        messages[logSize] = '\0';
        if (logSize > 10) { printf(">>> NNmatmul Compiler message: %s\n", messages); }
        free(messages);

        NNmatmul_kernel = clCreateKernel(NNmatmul_program, "NNmatmul", &ret); //std::cout << " retkernel=" << ret << std::endl;

        // Build the TN3matmul program
        strcpy(source_str,TN3matmul_src);
        source_size = strlen(source_str);

        TN3matmul_program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret); //std::cout << " retcreateprog=" << ret;
        ret = clBuildProgram(TN3matmul_program, 1, &(gpu.device_id[plat_indx][device_indx]), NULL, NULL, NULL);                    //std::cout << " retbuild=" << ret;

        clGetProgramBuildInfo(TN3matmul_program, gpu.device_id[plat_indx][device_indx], CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
        messages = (char*)malloc((1+logSize)*sizeof(char));
        clGetProgramBuildInfo(TN3matmul_program, gpu.device_id[plat_indx][device_indx], CL_PROGRAM_BUILD_LOG, logSize, messages, NULL);
        messages[logSize] = '\0';
        if (logSize > 10) { printf(">>> TN3matmul Compiler message: %s\n", messages); }
        free(messages);

        TN3matmul_kernel = clCreateKernel(TN3matmul_program, "TN3matmul", &ret); //std::cout << " retTN3kernel=" << ret << std::endl;


        // Build the NTmatmul program
        strcpy(source_str,NTmatmul_src);
        source_size = strlen(source_str);

        NTmatmul_program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret); //std::cout << " retcreateprog=" << ret;
        ret = clBuildProgram(NTmatmul_program, 1, &(gpu.device_id[plat_indx][device_indx]), NULL, NULL, NULL);                    //std::cout << " retbuild=" << ret;

        clGetProgramBuildInfo(NTmatmul_program, gpu.device_id[plat_indx][device_indx], CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
        messages = (char*)malloc((1+logSize)*sizeof(char));
        clGetProgramBuildInfo(NTmatmul_program, gpu.device_id[plat_indx][device_indx], CL_PROGRAM_BUILD_LOG, logSize, messages, NULL);
        messages[logSize] = '\0';
        if (logSize > 10) { printf(">>> NTmatmul Compiler message: %s\n", messages); }
        free(messages);

        NTmatmul_kernel = clCreateKernel(NTmatmul_program, "NTmatmul", &ret); //std::cout << " retNTkernel=" << ret << std::endl;


        // Build the TNmatmul program
        strcpy(source_str,TNmatmul_src);
        source_size = strlen(source_str);

        TNmatmul_program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret); //std::cout << " retcreateprog=" << ret;
        ret = clBuildProgram(TNmatmul_program, 1, &(gpu.device_id[plat_indx][device_indx]), NULL, NULL, NULL);                    //std::cout << " retbuild=" << ret;

        clGetProgramBuildInfo(TNmatmul_program, gpu.device_id[plat_indx][device_indx], CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
        messages = (char*)malloc((1+logSize)*sizeof(char));
        clGetProgramBuildInfo(TNmatmul_program, gpu.device_id[plat_indx][device_indx], CL_PROGRAM_BUILD_LOG, logSize, messages, NULL);
        messages[logSize] = '\0';
        if (logSize > 10) { printf(">>> TNmatmul Compiler message: %s\n", messages); }
        free(messages);

        TNmatmul_kernel = clCreateKernel(TNmatmul_program, "TNmatmul", &ret); //std::cout << " retTNkernel=" << ret << std::endl;


        // Build the Generate_projectors program
        strcpy(source_str,Generate_projectors_src);
        source_size = strlen(source_str);

        Generate_projectors_program = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret); //std::cout << " retcreateprog=" << ret;
        ret = clBuildProgram(Generate_projectors_program, 1, &(gpu.device_id[plat_indx][device_indx]), NULL, NULL, NULL);                    //std::cout << " retbuild=" << ret;

        clGetProgramBuildInfo(Generate_projectors_program, gpu.device_id[plat_indx][device_indx], CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
        messages = (char*)malloc((1+logSize)*sizeof(char));
        clGetProgramBuildInfo(Generate_projectors_program, gpu.device_id[plat_indx][device_indx], CL_PROGRAM_BUILD_LOG, logSize, messages, NULL);
        messages[logSize] = '\0';
        if (logSize > 10) { printf(">>> Generate_projectors Compiler message: %s\n", messages); }
        free(messages);

        Generate_projectors_kernel = clCreateKernel(Generate_projectors_program, "Generate_projectors", &ret); //std::cout << " retGenerate_projectors_kernel=" << ret << std::endl;

        free(source_str);


        // Build the program2
/*
        strcpy(source_str,program2_src);
        source_size = strlen(source_str);
        program1 = clCreateProgramWithSource(context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
        ret = clBuildProgram(program2, 1, &device_id_selected, NULL, NULL, NULL);

        logSize;
        clGetProgramBuildInfo(program2, device_id_selected, CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);
        messages = (char*)malloc((1+logSize)*sizeof(char));
        clGetProgramBuildInfo(program2, device_id_selected, CL_PROGRAM_BUILD_LOG, logSize, messages, NULL);
        messages[logSize] = '\0';
        if (logSize > 10) { printf(">>> Compiler message: %s\n", messages); }
        free(messages);

        kernel2 = clCreateKernel(program2, "NN_matmul", &ret);
*/


     }

     ~Gdevices() {
        std::cout << "Deallocating Gdevices" << std::endl;

        for (auto i=0; i<ntmp_mem; ++i)
           free(tmp_mem[i]);

        cl_int ret;
        for (auto i=0; i<ndev_mem; ++i)
           ret = clReleaseMemObject(dev_mem[i]);

        ret = clReleaseKernel(Generate_projectors_kernel);
        ret = clReleaseProgram(Generate_projectors_program);

        ret = clReleaseKernel(TNmatmul_kernel);
        ret = clReleaseProgram(TNmatmul_program);

        ret = clReleaseKernel(NTmatmul_kernel);
        ret = clReleaseProgram(NTmatmul_program);

        ret = clReleaseKernel(TN3matmul_kernel);
        ret = clReleaseProgram(TN3matmul_program);

        ret = clReleaseKernel(NNmatmul_kernel);
        ret = clReleaseProgram(NNmatmul_program);
        //ret = clReleaseKernel(kernel2);
        //ret = clReleaseProgram(program2);
        ret = clReleaseCommandQueue(command_queue);
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

     int fetch_dev_mem_indx(const size_t ndsize, const int rw) {
        cl_int ret;
        int ii = 0;
        while (((ndsize!=ndsize_mem[ii]) || inuse[ii] || (rw!=rw_mem[ii])) && (ii<ndev_mem))
          ++ii;

        if (ii<ndev_mem) {
           inuse[ii] = true;
        } else {
           ii            = ndev_mem;
           inuse[ii]     = true;
           rw_mem[ii]    = rw;
           ndsize_mem[ii] = ndsize;
           if (rw==1)
              dev_mem[ii]   = clCreateBuffer(context, CL_MEM_READ_ONLY,  ndsize*sizeof(double), NULL, &ret);
           else if (rw==2)
              dev_mem[ii]   = clCreateBuffer(context, CL_MEM_WRITE_ONLY, ndsize*sizeof(double), NULL, &ret);
           else if (rw==3)
              dev_mem[ii]   = clCreateBuffer(context, CL_MEM_READ_WRITE, ndsize*sizeof(double), NULL, &ret);
           else if (rw==4)
              dev_mem[ii]   = clCreateBuffer(context, CL_MEM_READ_ONLY,  ndsize*sizeof(int), NULL, &ret);
           else if (rw==5)
              dev_mem[ii]   = clCreateBuffer(context, CL_MEM_WRITE_ONLY, ndsize*sizeof(int), NULL, &ret);
           else if (rw==6)
              dev_mem[ii]   = clCreateBuffer(context, CL_MEM_READ_WRITE, ndsize*sizeof(int), NULL, &ret);

           ndev_mem += 1;
        }

        return ii;
     }


     int fetch_tmp_mem_indx(const size_t tmpndsize) {
        int ii = 0;
        while (((tmpndsize!=tmpndsize_mem[ii]) || tmpinuse[ii] ) && (ii<ntmp_mem))
          ++ii;

        if (ii<ntmp_mem) {
           tmpinuse[ii] = true;
        } else {
           ii                = ntmp_mem;
           tmpinuse[ii]      = true;
           tmpndsize_mem[ii] = tmpndsize;
           tmp_mem[ii]       = (double *) malloc(tmpndsize*sizeof(double));
           ntmp_mem += 1;
        }

        return ii;
     }


     void TN3_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_caa, double *host_cab, double *host_cbb)
     {
#if 0
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
#else
        cl_int ret;
        cl_uint nevents = 2;
        cl_event events[5];
        cl_event wevent;
        int one = 1;
        int ia = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne),1);
        int ib = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne),1);
        int icaa = fetch_dev_mem_indx(((size_t) ne) * ((size_t) ne)*((size_t) 3),2);
        int iccaa = fetch_tmp_mem_indx( ((size_t) ne)*((size_t) ne)*((size_t) 3) );

        int nn = ne*ne;
        double *tmpc = tmp_mem[iccaa];

        int ifac1=1;
        for (int i=2; i<=16; ++i)
        {
           if ((ne%i)==0)    ifac1 = i;
        }
        ret = clEnqueueWriteBuffer(command_queue,dev_mem[ia],CL_FALSE,0,npack*ne*sizeof(double),host_a,0,NULL,&events[0]); //std::cout << " ret1=" << ret;
        ret = clEnqueueWriteBuffer(command_queue,dev_mem[ib],CL_FALSE,0,npack*ne*sizeof(double),host_b,0,NULL,&events[1]); //std::cout << " ret2=" << ret;

        const int MNK[2] = {npack,ne};

        ret = clSetKernelArg(TN3matmul_kernel,0,sizeof(int),(void*)&MNK[0]); //std::cout << " ret4=" << ret;
        ret = clSetKernelArg(TN3matmul_kernel,1,sizeof(int),(void*)&MNK[1]); //std::cout << " ret5=" << ret;
        ret = clSetKernelArg(TN3matmul_kernel,2,sizeof(cl_mem),(void *)&(dev_mem[ia])); //std::cout << " ret7=" << ret;
        ret = clSetKernelArg(TN3matmul_kernel,3,sizeof(cl_mem),(void *)&(dev_mem[ib])); //std::cout << " ret8=" << ret;
        ret = clSetKernelArg(TN3matmul_kernel,4,sizeof(cl_mem),(void *)&(dev_mem[icaa])); //std::cout << " ret9aa=" << ret;


        // Execute the OpenCL kernel on the list
        const size_t global_item_size[2] = {(size_t) ne, (size_t) ne};
        const size_t local_item_size[2]  = {(size_t) ifac1, (size_t) ifac1};
        ret = clEnqueueNDRangeKernel(command_queue,
                                     TN3matmul_kernel, 2, NULL,
                                     global_item_size,
                                     local_item_size, nevents,events, NULL);

        DSCAL_PWDFT(nn,beta,host_caa,one);
        DSCAL_PWDFT(nn,beta,host_cab,one);
        DSCAL_PWDFT(nn,beta,host_cbb,one);

        ret = clEnqueueReadBuffer(command_queue,dev_mem[icaa],CL_FALSE,0,3*ne*ne*sizeof(double),tmp_mem[iccaa],0,NULL,&wevent);
        ret = clWaitForEvents(1,&wevent);

        DAXPY_PWDFT(nn,alpha,tmpc,one,host_caa,one);
        DAXPY_PWDFT(nn,alpha,&(tmpc[nn]),one,host_cab,one);
        DAXPY_PWDFT(nn,alpha,&(tmpc[nn+nn]),one,host_cbb,one);

        inuse[ia] = false;
        inuse[ib] = false;
        inuse[icaa] = false;
        tmpinuse[iccaa] = false;
#endif
     }


     void TN_dgemm(int ne, int nprj, int npack, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
#if 0
        DGEMM_PWDFT((char *) "T",(char *) "N",ne,nprj,npack,alpha,host_a,npack,host_b,npack,beta,host_c,ne);
#else
        int one = 1;
        int ia = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne),1);
        int ib = fetch_dev_mem_indx(((size_t) npack) * ((size_t) nprj),1);
        int ic = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) nprj),2);
        int icc = fetch_tmp_mem_indx( ((size_t) ne)  * ((size_t) nprj) );

        cl_int ret;
        cl_uint nevents = 2;
        cl_event events[5];
        cl_event wevent;
        int ifac1=1;
        int ifac2=1;
        for (int i=2; i<=16; ++i)
        {
           if ((ne%i)==0)   ifac1 = i;
           if ((nprj%i)==0) ifac2 = i;
        }
        ret = clEnqueueWriteBuffer(command_queue,dev_mem[ia],CL_FALSE,0,npack*ne  *sizeof(double),host_a,0,NULL,&events[0]); //std::cout << " ret1=" << ret;
        ret = clEnqueueWriteBuffer(command_queue,dev_mem[ib],CL_FALSE,0,npack*nprj*sizeof(double),host_b,0,NULL,&events[1]); //std::cout << " ret2=" << ret;

        const int MNK[3] = {ne,nprj,npack};

        ret = clSetKernelArg(NNmatmul_kernel,0,sizeof(int),(void*)&MNK[0]); //std::cout << " ret4=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,1,sizeof(int),(void*)&MNK[1]); //std::cout << " ret5=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,2,sizeof(int),(void*)&MNK[2]); //std::cout << " ret6=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,3,sizeof(cl_mem),(void *)&(dev_mem[ia])); //std::cout << " ret7=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,4,sizeof(cl_mem),(void *)&(dev_mem[ib])); //std::cout << " ret8=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,5,sizeof(cl_mem),(void *)&(dev_mem[ic])); //std::cout << " ret9=" << ret << std::endl;;

        // Execute the OpenCL kernel on the list
        const size_t global_item_size[2] = {(size_t) ne,    (size_t) nprj};
        const size_t local_item_size[2]  = {(size_t) ifac1, (size_t) ifac2};
        ret = clEnqueueNDRangeKernel(command_queue,
                                     NNmatmul_kernel, 2, NULL,
                                     global_item_size,
                                     local_item_size, nevents,events, NULL);
        ret = clEnqueueReadBuffer(command_queue,dev_mem[ic],CL_FALSE,0,ne*nprj*sizeof(double),tmp_mem[icc],0,NULL,&wevent);
        ret = clWaitForEvents(1,&wevent);

        int nn = ne*nprj;
        DSCAL_PWDFT(nn,beta,host_c,one);
        DAXPY_PWDFT(nn,alpha,tmp_mem[icc],one,host_c,one);
        //std::cout << " host_c[0]=" << host_c[0] << " host_c=" << host_c[1] << std::endl << std::endl;;

        inuse[ia] = false;
        inuse[ib] = false;
        inuse[ic] = false;
        tmpinuse[icc] = false;
#endif
     }


     void NN_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
#if 0
        DGEMM_PWDFT((char *) "N",(char *) "N",npack,ne,ne,alpha,host_a,npack,host_b,ne,beta,host_c,npack);

        //DGEMM_PWDFT((char *) "N",(char *) "N",npack,ne,ne,alpha,host_a,npack,host_b,ne,beta,tmpc,npack);
        //std::cout << " ehost_c[0]=" << tmpc[0] << " ehost_c[1]=" << tmpc[1] << std::endl;

#else
        int one = 1;
        int ia = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne),1);
        int ib = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) ne),1);
        int ic = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne),2);
        int icc = fetch_tmp_mem_indx( ((size_t) npack)*((size_t) ne) );

        //std::cout << "ndev_mem=" << ndev_mem << " ia=" << ia << " ib=" << ib << " ic=" << ic << " alpha=" << alpha << " beta=" << beta;
        //std::cout << " npack=" << npack << " ne=" << ne;

        cl_int ret;
        cl_uint nevents = 2;
        cl_event events[5];
        cl_event wevent;
        int ifac1=1;
        int ifac2=1;
        for (int i=2; i<=16; ++i)
        {
           if ((npack%i)==0) ifac1 = i;
           if ((ne%i)==0)    ifac2 = i;
        }
        ret = clEnqueueWriteBuffer(command_queue,dev_mem[ia],CL_FALSE,0,npack*ne*sizeof(double),host_a,0,NULL,&events[0]); //std::cout << " ret1=" << ret;
        ret = clEnqueueWriteBuffer(command_queue,dev_mem[ib],CL_FALSE,0,ne*ne   *sizeof(double),host_b,0,NULL,&events[1]); //std::cout << " ret2=" << ret;
        //std::cout << "  ifac1=" << ifac1 << " ifac2=" << ifac2 << " events[0]=" << events[0] << " events[1]=" << events[1] << std::endl;

        const int MNK[3] = {npack,ne,ne};

        ret = clSetKernelArg(NNmatmul_kernel,0,sizeof(int),(void*)&MNK[0]); //std::cout << " ret4=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,1,sizeof(int),(void*)&MNK[1]); //std::cout << " ret5=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,2,sizeof(int),(void*)&MNK[2]); //std::cout << " ret6=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,3,sizeof(cl_mem),(void *)&(dev_mem[ia])); //std::cout << " ret7=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,4,sizeof(cl_mem),(void *)&(dev_mem[ib])); //std::cout << " ret8=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,5,sizeof(cl_mem),(void *)&(dev_mem[ic])); //std::cout << " ret9=" << ret << std::endl;;

        // Execute the OpenCL kernel on the list
        const size_t global_item_size[2] = {(size_t) npack, (size_t) ne};
        const size_t local_item_size[2]  = {(size_t) ifac1, (size_t) ifac2};
        ret = clEnqueueNDRangeKernel(command_queue,
                                     NNmatmul_kernel, 2, NULL,
                                     global_item_size,
                                     local_item_size, nevents,events, NULL);
        ret = clEnqueueReadBuffer(command_queue,dev_mem[ic],CL_FALSE,0,npack*ne*sizeof(double),tmp_mem[icc],0,NULL,&wevent);
        ret = clWaitForEvents(1,&wevent);

        int nn = npack*ne;
        DSCAL_PWDFT(nn,beta,host_c,one);
        DAXPY_PWDFT(nn,alpha,tmp_mem[icc],one,host_c,one);
        //std::cout << " host_c[0]=" << host_c[0] << " host_c=" << host_c[1] << std::endl << std::endl;;

        inuse[ia] = false;
        inuse[ib] = false;
        inuse[ic] = false;
        tmpinuse[icc] = false;

#endif

     }

     void NT_dgemm(int npack, int ne, int nprj, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
#if 0

        DGEMM_PWDFT((char *) "N",(char *) "T",npack,ne,nprj,alpha,host_a,npack,host_b,ne,beta,host_c,npack);
#else
        int one = 1;
        int ia = fetch_dev_mem_indx(((size_t) npack) * ((size_t) nprj),1);
        int ib = fetch_dev_mem_indx(((size_t) ne)    * ((size_t) nprj),1);
        int ic = fetch_dev_mem_indx(((size_t) npack) * ((size_t) ne),2);
        int icc = fetch_tmp_mem_indx( ((size_t) npack)*((size_t) ne) );

        //std::cout << "ndev_mem=" << ndev_mem << " ia=" << ia << " ib=" << ib << " ic=" << ic << " alpha=" << alpha << " beta=" << beta;
        //std::cout << " npack=" << npack << " ne=" << ne;
        cl_int ret;
        cl_uint nevents = 2;
        cl_event events[5];
        cl_event wevent;
        int ifac1=1;
        int ifac2=1;
        for (int i=2; i<=16; ++i)
        {
           if ((npack%i)==0) ifac1 = i;
           if ((ne%i)==0)    ifac2 = i;
        }
        ret = clEnqueueWriteBuffer(command_queue,dev_mem[ia],CL_FALSE,0,npack*nprj*sizeof(double),host_a,0,NULL,&events[0]); //std::cout << " ret1=" << ret;
        ret = clEnqueueWriteBuffer(command_queue,dev_mem[ib],CL_FALSE,0,ne*nprj   *sizeof(double),host_b,0,NULL,&events[1]); //std::cout << " ret2=" << ret;
        //std::cout << "  ifac1=" << ifac1 << " ifac2=" << ifac2 << " events[0]=" << events[0] << " events[1]=" << events[1] << std::endl;

        const int MNK[3] = {npack,ne,nprj};

        ret = clSetKernelArg(NNmatmul_kernel,0,sizeof(int),(void*)&MNK[0]); //std::cout << " ret4=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,1,sizeof(int),(void*)&MNK[1]); //std::cout << " ret5=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,2,sizeof(int),(void*)&MNK[2]); //std::cout << " ret6=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,3,sizeof(cl_mem),(void *)&(dev_mem[ia])); //std::cout << " ret7=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,4,sizeof(cl_mem),(void *)&(dev_mem[ib])); //std::cout << " ret8=" << ret;
        ret = clSetKernelArg(NNmatmul_kernel,5,sizeof(cl_mem),(void *)&(dev_mem[ic])); //std::cout << " ret9=" << ret << std::endl;;

        // Execute the OpenCL kernel on the list
        const size_t global_item_size[2] = {(size_t) npack, (size_t) ne};
        const size_t local_item_size[2]  = {(size_t) ifac1, (size_t) ifac2};
        ret = clEnqueueNDRangeKernel(command_queue,
                                     NTmatmul_kernel, 2, NULL,
                                     global_item_size,
                                     local_item_size, nevents,events, NULL);
        ret = clEnqueueReadBuffer(command_queue,dev_mem[ic],CL_FALSE,0,npack*ne*sizeof(double),tmp_mem[icc],0,NULL,&wevent);
        ret = clWaitForEvents(1,&wevent);

        int nn = npack*ne;
        DSCAL_PWDFT(nn,beta,host_c,one);
        DAXPY_PWDFT(nn,alpha,tmp_mem[icc],one,host_c,one);
        //std::cout

        inuse[ia] = false;
        inuse[ib] = false;
        inuse[ic] = false;
        tmpinuse[icc] = false;
#endif
     }


     void load_index_phafac(const int nion, const int npack, 
                            const int indxi[], 
                            const int indxj[], 
                            const int indxk[],
                            const int nx, const int ny, const int nz,
                            const double exi[], const double exj[], const double exk[])
     {
        cl_int ret;
        cl_uint nevents = 6;
        cl_event events[6];

        npack_pf = npack;
        nx_pf = nx; ny_pf = ny; nz_pf = nz;

        indxi_dev = fetch_dev_mem_indx(((size_t) npack),4);
        indxj_dev = fetch_dev_mem_indx(((size_t) npack),4);
        indxk_dev = fetch_dev_mem_indx(((size_t) npack),4);

        exi_dev = fetch_dev_mem_indx(((size_t) (2*nx*nion)),1);
        exj_dev = fetch_dev_mem_indx(((size_t) (2*ny*nion)),1);
        exk_dev = fetch_dev_mem_indx(((size_t) (2*nz*nion)),1);


        ret = clEnqueueWriteBuffer(command_queue,dev_mem[indxi_dev],CL_FALSE,0,npack*sizeof(int),indxi,0,NULL,&events[0]); 
        ret = clEnqueueWriteBuffer(command_queue,dev_mem[indxj_dev],CL_FALSE,0,npack*sizeof(int),indxi,0,NULL,&events[1]); 
        ret = clEnqueueWriteBuffer(command_queue,dev_mem[indxk_dev],CL_FALSE,0,npack*sizeof(int),indxi,0,NULL,&events[2]); 

        ret = clEnqueueWriteBuffer(command_queue,dev_mem[exi_dev],CL_FALSE,0,2*nx*nion*sizeof(double),exi,0,NULL,&events[3]); 
        ret = clEnqueueWriteBuffer(command_queue,dev_mem[exj_dev],CL_FALSE,0,2*nx*nion*sizeof(double),exj,0,NULL,&events[4]); 
        ret = clEnqueueWriteBuffer(command_queue,dev_mem[exk_dev],CL_FALSE,0,2*nx*nion*sizeof(double),exk,0,NULL,&events[5]); 


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

     void TN_dgemm(int ne, int nprj, int npack, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
        DGEMM_PWDFT((char *) "T",(char *) "N",ne,nprj,npack,alpha,host_a,npack,host_b,npack,beta,host_c,ne);
     }

     void NN_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
        DGEMM_PWDFT((char *) "N",(char *) "N",npack,ne,ne,alpha,host_a,npack,host_b,ne,beta,host_c,npack);
     }

     void NT_dgemm(int npack, int ne, int nprj, double alpha, double *host_a, double *host_b, double beta, double *host_c) {

        DGEMM_PWDFT((char *) "N",(char *) "T",npack,ne,nprj,alpha,host_a,npack,host_b,ne,beta,host_c,npack);
     }

};


#endif
#endif

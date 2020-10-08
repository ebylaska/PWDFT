#include <random>
#include <cstdio>
#include <iostream>
#include <limits>
#include <CL/sycl.hpp>
#include "mkl_blas_sycl.hpp"

#include <chrono>
using namespace std::chrono;
using namespace cl::sycl;


//////////////////////////////////////////////////////////////////////////////////////////
//                                  README
//
// Login to IRIS nodes and load the modules:
// qsub -I -n 1 -t 60 -q iris
// module load dpcpp/2020.07.30.001
// module load mkl/2020.07.30.002
//
//
// To build the code using:
// dpcpp -I$MKLROOT/include -Icommon -DMKL_ILP64 -std=c++11 -fsycl -O3 -o
// baseline_MKL_sycl baseline_MKL_sycl.cpp -L$MKLROOT/lib/intel64 -lmkl_sycl -lmkl_intel_ilp64
// -lmkl_sequential -lmkl_core -lOpenCL
//
// To execute: SYCL_BE=PI_OPENCL ./baseline_MKL_sycl
// 
// Change the Matrix sizes by appropriately setting the sizes with M, N and
// P values
//
//////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]) {
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(1.0, 2.0);

  // C = alpha * op(A) * op(B)  + beta * C
  oneapi::mkl::transpose matT = oneapi::mkl::transpose::trans;
  oneapi::mkl::transpose matN = oneapi::mkl::transpose::nontrans;

  // matrix data sizes
  int npack = 10000;
  int ne    = 100;
  if (argc>2)
  {
     /* npack*/
     sscanf(argv[1],"%d",&npack);
     /* ne*/
     sscanf(argv[2],"%d",&ne);
  }

  double nmultadds = (double) ne;
         nmultadds *= (double) ne;
         nmultadds *= (double) npack;

  // set scalar fp values
  double alpha = 1.0;
  double beta  = 0.0;

  // 1D arrays on host side
  double *host_a;
  double *host_b;
  double *host_c;
  host_a = new double[npack*ne]{};
  host_b = new double[npack*ne]{};
  host_c = new double[ne*ne]{};

  // prepare matrix data with ROW-major style
  // A(M, N)
  for (size_t j=0; j<ne; ++j)
  for (size_t i=0; i<npack; ++i)
  {
      host_a[i + j*npack] = dis(gen);
      host_b[i + j*npack] = dis(gen);
  }

  std::cout << "matrix sizes:" << std::endl;
  std::cout << "npack = " << npack << std::endl;
  std::cout << "ne    = " << ne << std::endl;
  std::cout << "c(" << ne << "," << ne << ") = a(" << npack << "," << ne << ")' * b(" << npack << "," << ne << ")" << std::endl;

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

  // Initializing the devices queue with the default selector
  // The device queue is used to enqueue the kernels and encapsulates
  // all the states needed for execution
  gpu_selector device_selector;
  queue device_queue(device_selector, asyncHandler);
  std::cout << "Device: " << device_queue.get_device().get_info<info::device::name>() << std::endl << std::endl;

  // Creating 1D buffers for matrices (double)
  double* dev_a = cl::sycl::malloc_device<double>(npack*ne, device_queue);
  double* dev_b = cl::sycl::malloc_device<double>(npack*ne, device_queue);
  double* dev_c = cl::sycl::malloc_device<double>(ne*ne,    device_queue);

  // Transfer info from CPU to GPU
  // copy host -> device
  //

  //device_queue.submit([&](cl::sycl::handler& cgh)
  //{
  //  cgh.memcpy(dev_a, host_a, npack*ne*sizeof(double));
  //});
  //device_queue.wait();
  //device_queue.submit([&](cl::sycl::handler& cgh)
  //{
  //  cgh.memcpy(dev_b, host_b, npack*ne*sizeof(double));
  //});
  //device_queue.wait();

  // device_queue.memcpy(dev_a, host_a, sizeof(double)*M*N);
  // device_queue.memcpy(dev_b, host_b, sizeof(double)*N*P);

  try {
//Warm up the kernel
      //for (int i=0; i<5; i++)
      oneapi::mkl::blas::gemm(device_queue, matT, matN, ne,ne,npack, alpha, dev_a, npack, dev_b, npack, beta, dev_c, ne);
      device_queue.wait();
      


      auto start = high_resolution_clock::now();
      device_queue.submit([&](cl::sycl::handler& cgh) { cgh.memcpy(dev_a, host_a, npack*ne*sizeof(double)); });
      device_queue.submit([&](cl::sycl::handler& cgh) { cgh.memcpy(dev_b, host_b, npack*ne*sizeof(double)); });
      device_queue.wait();

      oneapi::mkl::blas::gemm(device_queue, matT, matN, ne,ne,npack, alpha, dev_a, npack, dev_b, npack, beta, dev_c, ne);
      device_queue.memcpy(host_c, dev_c, ne*ne*sizeof(double));
      device_queue.wait();
      auto stop = high_resolution_clock::now();



      duration<double> deltatime = stop-start;
      double dt = (double) deltatime.count();

      std::cout << "GPU time (ms): " << dt*1000.0 << std::endl;
      std::cout << "GPU GFLOPs   : " << (nmultadds/dt)*1.0e-9 << std::endl;

      std::cout << std::endl;
      std::cout << "csv:" << npack << "," << ne << "," << dt*1000.0 << "," <<  (nmultadds/dt)*1.0e-9 << std::endl;
      std::cout << std::endl;
  }
  catch(cl::sycl::exception const& e) {
    std::cout << "\t\tSYCL exception during GEMM\n" << e.what() << std::endl << "OpenCL status: " << e.get_cl_code() << std::endl;
  }


  delete[] host_a;
  delete[] host_b;
  cl::sycl::free(dev_a, device_queue);
  cl::sycl::free(dev_b, device_queue);
  cl::sycl::free(dev_c, device_queue);

  return 0;
}

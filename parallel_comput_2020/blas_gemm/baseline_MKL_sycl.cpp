#include <random>
#include <iostream>
#include <limits>
#include <CL/sycl.hpp>
#include "mkl_blas_sycl.hpp"

#include <chrono>
using namespace std::chrono;
using namespace cl::sycl;

#define M 8192
#define N 8192
#define P 8192

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

int main() {
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(1.0, 2.0);

  // C = alpha * op(A) * op(B)  + beta * C
  oneapi::mkl::transpose transA = oneapi::mkl::transpose::nontrans;
  oneapi::mkl::transpose transB = oneapi::mkl::transpose::nontrans;

  // matrix data sizes
  int m = M;
  int n = P;
  int k = N;

  // leading dimensions of data
  int ldA = k;
  int ldB = n;
  int ldC = n;

  // set scalar fp values
  double alpha = 1.0;
  double beta  = 0.0;

  // 1D arrays on host side
  double *host_a;
  double *host_b;
  host_a = new double[M*N]{};
  host_b = new double[N*P]{};

  // prepare matrix data with ROW-major style
  // A(M, N)
  for (size_t i=0; i<M; i++)
    for (size_t j=0; j<N; j++)
      host_a[i*N + j] = dis(gen);
  // B(N, P)
  for (size_t i=0; i<N; i++)
    for (size_t j=0; j<P; j++)
      host_b[i*P + j] = dis(gen);

  std::cout << "Problem size: c(" << M << "," << P << ") = a(" << M << "," << N << ") * b(" << N << "," << P << ")" << std::endl;

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
  double* dev_a = cl::sycl::malloc_device<double>(M*N, device_queue);
  double* dev_b = cl::sycl::malloc_device<double>(N*P, device_queue);
  double* dev_c = cl::sycl::malloc_device<double>(M*P, device_queue);

  // Transfer info from CPU to GPU
  // copy host -> device
  device_queue.submit([&](cl::sycl::handler& cgh)
  {
    cgh.memcpy(dev_a, host_a, M*N*sizeof(double));
  });
  device_queue.wait();
  device_queue.submit([&](cl::sycl::handler& cgh)
  {
    cgh.memcpy(dev_b, host_b, N*P*sizeof(double));
  });
  device_queue.wait();

  // device_queue.memcpy(dev_a, host_a, sizeof(double)*M*N);
  // device_queue.memcpy(dev_b, host_b, sizeof(double)*N*P);

  try {
//Warm up the kernel
      for (int i=0; i<5; i++)
          oneapi::mkl::blas::gemm(device_queue, transB, transA, n, m, k, alpha, dev_b, ldB, dev_a, ldA, beta, dev_c, ldC);
      device_queue.wait();

      auto start = high_resolution_clock::now();
      oneapi::mkl::blas::gemm(device_queue, transB, transA, n, m, k, alpha, dev_b, ldB, dev_a, ldA, beta, dev_c, ldC);
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<microseconds>(stop - start);
      std::cout << "GPU time (ms): " << duration.count() << std::endl;
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

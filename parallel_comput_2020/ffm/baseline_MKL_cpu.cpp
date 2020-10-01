#include <random>
#include <iostream>
#include <limits>
#include <mkl.h>

#include <chrono>
using namespace std::chrono;

#define N_ROUNDS 10

//////////////////////////////////////////////////////////////////////////////////////////
//                                  README
//
// Login to IRIS nodes and load the modules:
// qsub -I -n 1 -t 60 -q iris
// module load oneapi
//
//
// To build the code using:
// icpx -o baseline_cpu_dbl baseline_MKL_cpu.cpp -mkl -lmkl_intel_ilp64 -lmkl_core -lpthread -ldl -lm
//
// To execute: ./baseline_cpu_db $Ne $Npack
//
//
//////////////////////////////////////////////////////////////////////////////////////////

int main(int arc, char ** argv) {

  int M = atoi(argv[1]); // Ne
  int K = atoi(argv[2]); // Npack
  int N = M;

  // std::cout << "
  std::random_device rd;  //Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(1.0, 2.0);

  // C = alpha * op(A) * op(B)  + beta * C

  // matrix data sizes
  int m = M;
  int n = N;
  int k = K;

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
  double *host_c;
  host_a = new double[M*K]{};
  host_b = new double[K*N]{};
  host_c = new double[M*N]{};

  // prepare matrix data with ROW-major style
  // A(M, K)
  for (size_t i=0; i<M; i++)
    for (size_t j=0; j<K; j++)
        host_a[i*K + j] = dis(gen);
  // B(K, N)
  for (size_t i=0; i<K; i++)
    for (size_t j=0; j<N; j++)
        host_b[i*N + j] = dis(gen);

  std::cout << "Problem size: c(" << M << "," << N << ") = a(" << M << "," << K << ") * b(" << K << "," << N << ")";

  //Warm up the kernel (dummy runs)
  for (int i=0; i<5; i++)
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, m, k,
                  alpha, host_b, ldB, host_a, ldA, beta, host_c, ldC);

  auto start = high_resolution_clock::now();
  for (int i=0; i<N_ROUNDS; i++)
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, m, k,
                  alpha, host_b, ldB, host_a, ldA, beta, host_c, ldC);
  auto stop = high_resolution_clock::now();
  double time = duration<double, std::micro>(stop - start).count();
  time /= N_ROUNDS;

  std::cout << "   avg CPU time (ms): " << time << std::endl;

  delete[] host_a;
  delete[] host_b;
  delete[] host_c;

  return 0;
}

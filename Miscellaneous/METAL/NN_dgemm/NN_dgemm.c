#include <Accelerate/Accelerate.h>
#include <stdio.h>
#include <stdlib.h>

// Function to perform NN_dgemm using Apple's Accelerate (BLAS)
void NN_dgemm(int M, int N, int K, 
              double alpha, 
              double *A, 
              double *B, 
              double beta, 
              double *C) 
{
    // Call BLAS `dgemm`
    cblas_dgemm(CblasRowMajor,  // Row-major storage
                CblasNoTrans,   // A is not transposed
                CblasNoTrans,   // B is not transposed
                M, N, K,        // Matrix dimensions
                alpha,          // Scalar alpha
                A, K,           // Matrix A and its leading dimension
                B, N,           // Matrix B and its leading dimension
                beta,           // Scalar beta
                C, N);          // Matrix C and its leading dimension
}

#include <mach/mach_time.h>  // High-precision timing for macOS

int main() {
    int M = 512, N = 512, K = 512;

    // Allocate matrices
    double *A = (double *)malloc(M * K * sizeof(double));
    double *B = (double *)malloc(K * N * sizeof(double));
    double *C = (double *)malloc(M * N * sizeof(double));

    // Initialize matrices with random values
    for (int i = 0; i < M * K; i++) A[i] = (double)(rand() % 10);
    for (int i = 0; i < K * N; i++) B[i] = (double)(rand() % 10);
    for (int i = 0; i < M * N; i++) C[i] = 0.0;  // Initialize C to zero

    // ✅ Start Timing
    uint64_t start = mach_absolute_time();

    // Call NN_dgemm (Accelerate BLAS)
    NN_dgemm(M, N, K, 1.0, A, B, 0.0, C);

    // ✅ End Timing
    uint64_t end = mach_absolute_time();

    // Convert to milliseconds
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);
    double elapsed_ms = (double)(end - start) * timebase.numer / timebase.denom / 1e6;

    printf("Matrix multiplication complete in %.3f ms\n", elapsed_ms);


    printf("Accelerate BLAS: Double precision matrix multiplication complete!\n");

    // Free allocated memory
    free(A);
    free(B);
    free(C);

    return 0;
}


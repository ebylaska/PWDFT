#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#import <stdio.h>
#include <stdlib.h>
#include <mach/mach_time.h>

// ✅ Global Metal Pipeline (Precompiled Once)
id<MTLDevice> globalDevice = nil;
id<MTLCommandQueue> globalCommandQueue = nil;
id<MTLComputePipelineState> globalPipelineState = nil;

// ✅ Metal Kernel (Embedded as a String)
const char *metalKernelSource = R"(
#include <metal_stdlib>
using namespace metal;

// Define FP64 emulation using float2 (hi = x, lo = y)
struct Double {
    float2 value;
};

// FP64 Addition (Emulated with FP32)
Double add_double(Double a, Double b) {
    Double result;
    float s1 = a.value.x + b.value.x;
    float v = s1 - a.value.x;
    float s2 = ((b.value.x - v) + (a.value.x - (s1 - v))) + a.value.y + b.value.y;
    result.value.x = s1 + s2;
    result.value.y = s2 - (result.value.x - s1);
    return result;
}

// FP64 Multiplication (Emulated with FP32)
Double mul_double(Double a, Double b) {
    Double result;
    float p1 = a.value.x * b.value.x;
    float p2 = a.value.x * b.value.y + a.value.y * b.value.x;
    result.value.x = p1 + p2;
    result.value.y = p2 - (result.value.x - p1);
    return result;
}

// ✅ Metal Kernel for FP64 Matrix Multiplication (Emulated)
kernel void matrix_multiply_fp64(
    const device Double *A,
    const device Double *B,
    device Double *C,
    constant uint &M,
    constant uint &N,
    constant uint &K,
    uint2 gid [[thread_position_in_grid]]) 
{
    if (gid.x >= M || gid.y >= N) return;

    Double sum;
    sum.value = float2(0.0, 0.0);

    for (uint k = 0; k < K; ++k) {
        Double a = A[gid.x * K + k];
        Double b = B[k * N + gid.y];
        sum = add_double(sum, mul_double(a, b));
    }

    C[gid.x * N + gid.y] = sum;
}
)";

// ✅ Function to Initialize Metal & Precompile Kernel
void initializeMetal() {
    NSArray<id<MTLDevice>> *devices = MTLCopyAllDevices();
    if (devices.count == 0) {
        fprintf(stderr, "Error: No Metal devices found!\n");
        exit(1);
    }
    globalDevice = devices[0];  // Select the first available GPU
    globalCommandQueue = [globalDevice newCommandQueue];

    NSError *error = nil;
    NSString *metalSource = [NSString stringWithUTF8String:metalKernelSource];
    id<MTLLibrary> library = [globalDevice newLibraryWithSource:metalSource options:nil error:&error];

    if (!library) {
        fprintf(stderr, "Error compiling Metal source: %s\n", error.localizedDescription.UTF8String);
        exit(1);
    }

    id<MTLFunction> function = [library newFunctionWithName:@"matrix_multiply_fp64"];
    globalPipelineState = [globalDevice newComputePipelineStateWithFunction:function error:&error];

    if (!globalPipelineState) {
        fprintf(stderr, "Error creating pipeline state: %s\n", error.localizedDescription.UTF8String);
        exit(1);
    }
    
    printf("✅ Metal Initialized & Kernel Precompiled on Device: %s\n", [[globalDevice name] UTF8String]);
}

// ✅ Struct for FP64 Emulation (Using Two FP32 Values)
struct Double {
    float hi;
    float lo;
};

// ✅ Function for FP64 Matrix Multiplication with FP32 Emulation
void NN_dgemm(int M, int N, int K, double *hostA, double *hostB, double *hostC) {
    if (!globalPipelineState) {
        fprintf(stderr, "Error: Metal pipeline not initialized!\n");
        return;
    }

    // Convert double arrays to Double (FP32 pairs)
    Double *A_fp32 = (Double *)malloc(M * K * sizeof(Double));
    Double *B_fp32 = (Double *)malloc(K * N * sizeof(Double));
    Double *C_fp32 = (Double *)malloc(M * N * sizeof(Double));

    for (int i = 0; i < M * K; i++) {
        A_fp32[i].hi = (float)hostA[i];
        A_fp32[i].lo = (float)(hostA[i] - A_fp32[i].hi);
    }
    for (int i = 0; i < K * N; i++) {
        B_fp32[i].hi = (float)hostB[i];
        B_fp32[i].lo = (float)(hostB[i] - B_fp32[i].hi);
    }

    // Create Metal Buffers
    id<MTLBuffer> bufferA = [globalDevice newBufferWithBytes:A_fp32 length:M*K*sizeof(Double) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufferB = [globalDevice newBufferWithBytes:B_fp32 length:K*N*sizeof(Double) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufferC = [globalDevice newBufferWithBytes:C_fp32 length:M*N*sizeof(Double) options:MTLResourceStorageModeShared];

    id<MTLCommandBuffer> commandBuffer = [globalCommandQueue commandBuffer];
    id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];

    [computeEncoder setComputePipelineState:globalPipelineState];
    [computeEncoder setBuffer:bufferA offset:0 atIndex:0];
    [computeEncoder setBuffer:bufferB offset:0 atIndex:1];
    [computeEncoder setBuffer:bufferC offset:0 atIndex:2];

    // Set Matrix Dimensions as Constants
    [computeEncoder setBytes:&M length:sizeof(uint) atIndex:3];
    [computeEncoder setBytes:&N length:sizeof(uint) atIndex:4];
    [computeEncoder setBytes:&K length:sizeof(uint) atIndex:5];

    // Dispatch Threads
    MTLSize gridSize = MTLSizeMake(M, N, 1);
    MTLSize threadGroupSize = MTLSizeMake(16, 16, 1);

    [computeEncoder dispatchThreads:gridSize threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];

    // Submit Command Buffer and Wait for Completion
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];

    // Convert Results Back to Double
    memcpy(C_fp32, [bufferC contents], M*N*sizeof(Double));
    for (int i = 0; i < M * N; i++) {
        hostC[i] = (double)C_fp32[i].hi + (double)C_fp32[i].lo;
    }

    free(A_fp32);
    free(B_fp32);
    free(C_fp32);
}

#include <Accelerate/Accelerate.h>  // Apple BLAS
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mach/mach_time.h>  // High-precision timing for macOS

// Global Metal Pipeline (Precompiled Once)
//extern void initializeMetal();  // Declare Metal initialization function
//extern void NN_dgemm(int M, int N, int K, double alpha, double *hostA, double *hostB, double beta, double *hostC);

int main() {
    //int M = 512, N = 512, K = 512;
    //int M = 1024, N = 1024, K = 1024;
    //int M = 2048, N = 2048, K = 2048;
    int M = 4096, N = 4096, K = 4096;
    //int M = 8192, N = 8192, K = 8192;
    //int M = 2048, N = 2048, K = 2048;  // Must use `int64_t` for ILP64

    double *A = (double *)malloc(M * K * sizeof(double));
    double *B = (double *)malloc(K * N * sizeof(double));
    double *C_ref = (double *)malloc(M * N * sizeof(double));  // Reference
    double *C_metal = (double *)malloc(M * N * sizeof(double));  // Metal Result

    printf("matrix dimensions M=%d N=%d K=%d\n",M,N,K);

    // ✅ Initialize Metal (Precompile Kernel Before Running)
    initializeMetal();

    // ✅ Initialize Matrices A and B
    //for (int i = 0; i < M * K; i++) A[i] = (double)(rand() % 10);
    //for (int i = 0; i < K * N; i++) B[i] = (double)(rand() % 10);
    for (int i = 0; i < M * K; i++) A[i] = ((double)(rand() % 1000))/1000.0;
    for (int i = 0; i < K * N; i++) B[i] = ((double)(rand() % 1000))/1000.0;

    //for (int i = 0; i < M * K; i++) A[i] = drand48() * 10.0;  // Scaled to [0,10)
    //for (int i = 0; i < K * N; i++) B[i] = drand48() * 10.0;
    
    // ✅ Compute C_ref using Accelerate BLAS (Reference)
    uint64_t start_ref = mach_absolute_time();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                M, N, K, 1.0, A, K, B, N, 0.0, C_ref, N);
    uint64_t end_ref = mach_absolute_time();
    double elapsed_ref = (double)(end_ref - start_ref) * 1e-6;
    printf("✅ Accelerate BLAS DGEMM completed in %.3f ms\n", elapsed_ref);

    // ✅ Compute C_metal using Metal FP64 Emulation
    uint64_t start_metal = mach_absolute_time();
    NN_dgemm(M, N, K, A, B, C_metal);
    uint64_t end_metal = mach_absolute_time();
    double elapsed_metal = (double)(end_metal - start_metal) * 1e-6;
    printf("✅ Metal FP64 Emulation DGEMM completed in %.3f ms\n", elapsed_metal);


    // ✅ Compute FLOP Rate
    double flops = 2.0 * M * N * K;
    double tflops_ref = (flops / (elapsed_ref   * 1e9));
    double tflops_metal = (flops / (elapsed_metal * 1e9));

    printf("✅ Accelerate BLAS DGEMM completed in %.3f ms (%.2f TFLOP/s %0.2f GFLOP/s)\n", elapsed_ref, tflops_ref, tflops_ref*1e3);
    printf("✅ Metal FP64 Emulation DGEMM completed in %.3f ms (%.2f TFLOP/s, %0.2f GFLOP/s)\n", elapsed_metal, tflops_metal, tflops_metal*1e3);



    // ✅ Compute Accuracy Metrics
    double max_error = 0.0;
    double sum_squared_error = 0.0;
    double sum_error = 0.0;
    
    for (int i = 0; i < M * N; i++) {
        double error = fabs(C_metal[i] - C_ref[i]) / fabs(C_ref[i] + 1e-12);
        max_error = fmax(max_error, error);
        sum_squared_error += error * error;
        sum_error += error;
    }
    
    double mean_error = sum_error / (M * N);
    double rms_error = sqrt(sum_squared_error / (M * N));
 
    printf("C_metal[100]=%20.16le C_ref[100]=%20.16le\n",C_metal[100], C_ref[100]);

    // ✅ Print Accuracy Metrics
    printf("\n🔍 **Accuracy Comparison (Metal FP64 Emulation vs. Accelerate BLAS)**\n");
    printf("🔹 Max Relative Error:   %.5le\n", max_error);
    printf("🔹 Mean Relative Error:  %.5le\n", mean_error);
    printf("🔹 RMS Error:            %.5le\n", rms_error);
    
    free(A);
    free(B);
    free(C_ref);
    free(C_metal);

    return 0;
}




#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#import <stdio.h>
#include <stdlib.h>
#include <mach/mach_time.h>
#include <Accelerate/Accelerate.h>  // Apple BLAS for comparison

// ✅ Global Metal Pipeline (Precompiled Once)
id<MTLDevice> globalDevice = nil;
id<MTLCommandQueue> globalCommandQueue = nil;
id<MTLComputePipelineState> globalPipelineState = nil;

// ✅ Metal Kernel (Embedded as a String)
const char *metalKernelSource = R"(
#include <metal_stdlib>
using namespace metal;

// ✅ Define FP96 emulation using three FP32 values (hi, mid, lo)
struct ExtendedDouble { float hi; float mid; float lo; };

// ✅ FP96 Addition (Emulated with FP32)
ExtendedDouble add_extended(ExtendedDouble a, ExtendedDouble b) {
    ExtendedDouble result;
    float s1 = a.hi + b.hi;
    float v = s1 - a.hi;
    float s2 = ((b.hi - v) + (a.hi - (s1 - v))) + a.mid + b.mid;
    result.hi = s1 + s2;
    result.mid = s2 - (result.hi - s1);
    result.lo = 0.0f; // Ignore `lo` for now
    return result;
}



// ✅ FP96 Multiplication (Emulated with FP32)
ExtendedDouble mul_extended(ExtendedDouble a, ExtendedDouble b) {
    ExtendedDouble result;
    float p1 = a.hi * b.hi;
    float p2 = a.hi * b.mid + a.mid * b.hi;
    result.hi = p1 + p2;
    result.mid = p2 - (result.hi - p1);
    result.lo = 0.0f; // Ignore `lo` for now
    return result;
}



// ✅ Metal Kernel for FP96 Matrix Multiplication (Emulated)
kernel void matrix_multiply_fp96(
    const device ExtendedDouble *A,
    const device ExtendedDouble *B,
    device ExtendedDouble *C,
    constant uint &M,
    constant uint &N,
    constant uint &K,
    uint2 gid [[thread_position_in_grid]]) 
{
    if (gid.x >= M || gid.y >= N) return;

    // ✅ Accumulate into all components
    ExtendedDouble sum = {0.0f, 0.0f, 0.0f};

    for (uint k = 0; k < K; k++) {
        ExtendedDouble a = A[gid.x * K + k];
        ExtendedDouble b = B[k * N + gid.y];
        sum = add_extended(sum, mul_extended(a, b));
    }


    // ✅ Explicitly store .mid separately for debugging
    ExtendedDouble result;
    result.hi = sum.hi;
    result.mid = sum.mid;  // 🔹 Ensure Metal stores this properly
    result.lo = 0.0f; // Ignore `.lo` for now


    // ✅ Store result in C (ensuring full FP96 precision)
    C[gid.x * N + gid.y] = result;
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

    id<MTLFunction> function = [library newFunctionWithName:@"matrix_multiply_fp96"];
    globalPipelineState = [globalDevice newComputePipelineStateWithFunction:function error:&error];

    if (!globalPipelineState) {
        fprintf(stderr, "Error creating pipeline state: %s\n", error.localizedDescription.UTF8String);
        exit(1);
    }
    
    printf("✅ Metal Initialized & Kernel Precompiled on Device: %s\n", [[globalDevice name] UTF8String]);
}

// ✅ Function for FP96 Matrix Multiplication with FP32 Emulation
void NN_dgemm(int M, int N, int K, double *hostA, double *hostB, double *hostC) {
    if (!globalPipelineState) {
        fprintf(stderr, "Error: Metal pipeline not initialized!\n");
        return;
    }

    // Convert double arrays to ExtendedDouble (FP96)
    struct ExtendedDouble { float hi; float mid; float lo; };
    ExtendedDouble *A_fp32 = (ExtendedDouble *)malloc(M * K * sizeof(ExtendedDouble));
    ExtendedDouble *B_fp32 = (ExtendedDouble *)malloc(K * N * sizeof(ExtendedDouble));
    ExtendedDouble *C_fp32 = (ExtendedDouble *)malloc(M * N * sizeof(ExtendedDouble));

    for (int i = 0; i < M * K; i++) {
        A_fp32[i].hi = (float)hostA[i];
        A_fp32[i].mid = (float)(hostA[i] - A_fp32[i].hi);
        A_fp32[i].lo = 0.0;
    }
    for (int i = 0; i < K * N; i++) {
        B_fp32[i].hi = (float)hostB[i];
        B_fp32[i].mid = (float)(hostB[i] - B_fp32[i].hi);
        B_fp32[i].lo = 0.0;
    }

    printf("\n🔍 **FP96 Representation of A (First 5 values)**\n");
    for (int i = 0; i < (5); i++) {
        printf("A_fp32[%d]: hi = %.8e, mid = %.8e, lo = %.8e\n", 
               i, A_fp32[i].hi, A_fp32[i].mid, A_fp32[i].lo);
    }
    for (int i = 0; i < (5); i++) {
        printf("B_fp32[%d]: hi = %.8e, mid = %.8e, lo = %.8e\n", 
               i, B_fp32[i].hi, B_fp32[i].mid, B_fp32[i].lo);
    }


    // Create Metal Buffers
    id<MTLBuffer> bufferA = [globalDevice newBufferWithBytes:A_fp32 length:M*K*sizeof(ExtendedDouble) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufferB = [globalDevice newBufferWithBytes:B_fp32 length:K*N*sizeof(ExtendedDouble) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufferC = [globalDevice newBufferWithBytes:C_fp32 length:M*N*sizeof(ExtendedDouble) options:MTLResourceStorageModeShared];

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


    // Create Metal Buffers
    MTLSize gridSize = MTLSizeMake(M, N, 1);
    MTLSize threadGroupSize = MTLSizeMake(16, 16, 1);
    //MTLSize threadGroupSize = MTLSizeMake(32, 32, 1);

    [computeEncoder dispatchThreads:gridSize threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];

    // Submit Command Buffer and Wait for Completion
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];

    ExtendedDouble *C_gpu = (ExtendedDouble *)[bufferC contents];

    for (int i = 0; i < (5); i++) {
       printf("C_gpu[%d] = (%.6e, %.6e, %.6e)\n", i, C_gpu[i].hi, C_gpu[i].mid, C_gpu[i].lo);
    }



    memcpy(C_fp32, [bufferC contents], M*N*sizeof(ExtendedDouble));


    for (int i = 0; i < M*N; i++) {
        hostC[i] = (double)C_fp32[i].hi + (double)C_fp32[i].mid + (double)C_fp32[i].lo;
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


// ✅ **Main Function**
int main() {
    int M = 2048, N = 2048, K = 2048;
    //int M = 4, N = 4, K = 4;

    double *A = (double *)malloc(M * K * sizeof(double));
    double *B = (double *)malloc(K * N * sizeof(double));
    double *C_ref   = (double *)malloc(M * N * sizeof(double));  // Reference
    double *C_metal = (double *)malloc(M * N * sizeof(double));

    initializeMetal();


    printf("matrix dimensions M=%d N=%d K=%d\n",M,N,K);

    //for (int i = 0; i < M * K; i++) A[i] = ((double)(rand() % 1000)) / 1000.0;
    //for (int i = 0; i < K * N; i++) B[i] = ((double)(rand() % 1000)) / 1000.0;
    for (int i = 0; i < M * K; i++) A[i] = drand48() * 2.0 - 1.0;  // Range: [-1, 1]
    for (int i = 0; i < K * N; i++) B[i] = drand48() * 2.0 - 1.0;

    //printf("\n🔍 **Initial Matrix A (First 5 values)**\n");
    //for (int i = 0; i < (M*K); i++) {
    //    printf("A[%d] = %.16e\n", i, A[i]);
    //}

    //printf("\n🔍 **Initial Matrix B (First 5 values)**\n");
    //for (int i = 0; i < (K*N); i++) {
    //    printf("B[%d] = %.16e\n", i, B[i]);
    //}



    // ✅ Compute C_ref using Accelerate BLAS (Reference)
    uint64_t start_ref = mach_absolute_time();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, M, N, K, 1.0, A, K, B, N, 0.0, C_ref, N);
    uint64_t end_ref = mach_absolute_time();
    double elapsed_ref = (double)(end_ref - start_ref) * 1e-6;
    printf("✅ Accelerate BLAS DGEMM completed in %.3f ms\n", elapsed_ref);


    // ✅ Compute C_metal using Metal FP64 Emulation
    uint64_t start_metal = mach_absolute_time();
    NN_dgemm(M, N, K, A, B, C_metal);
    uint64_t end_metal = mach_absolute_time();
    double elapsed_metal = (double)(end_metal - start_metal) * 1e-6;
    printf("✅ FP96 Emulation (float-based) matrix multiplication completed in %.3f ms\n", elapsed_metal);


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
 
  //  printf("C_metal[50]=%20.16lf C_ref[50]=%20.16lf\n",C_metal[50], C_ref[50]);

    printf("C_metal[100]=%20.16le C_ref[100]=%20.16le\n",C_metal[100], C_ref[100]);

    // ✅ Print Accuracy Metrics
    printf("\n🔍 **Accuracy Comparison (Metal FP64 Emulation vs. Accelerate BLAS)**\n");
    printf("🔹 Max Relative Error:   %.5le\n", max_error);
    printf("🔹 Mean Relative Error:  %.5le\n", mean_error);
    printf("🔹 RMS Error:            %.5le\n", rms_error);
    

    free(A);
    free(B);
    free(C_ref);
    //free(C_metal);

    return 0;
}


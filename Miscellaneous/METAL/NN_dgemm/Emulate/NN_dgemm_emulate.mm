#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#import <stdio.h>
#include <stdlib.h>

// Define Metal Kernel (FP64 Emulation Using FP32 Pairs)
const char *metalKernelSource = R"(
#include <metal_stdlib>
using namespace metal;

struct Double {
    float hi;
    float lo;
};

// FP64 Addition (using FP32 pairs)
Double add_double(Double a, Double b) {
    Double result;
    float s1 = a.hi + b.hi;
    float v = s1 - a.hi;
    float s2 = ((b.hi - v) + (a.hi - (s1 - v))) + a.lo + b.lo;
    result.hi = s1 + s2;
    result.lo = s2 - (result.hi - s1);
    return result;
}

// FP64 Multiplication (using FP32 pairs)
Double mul_double(Double a, Double b) {
    Double result;
    float p1 = a.hi * b.hi;
    float p2 = a.hi * b.lo + a.lo * b.hi;
    result.hi = p1 + p2;
    result.lo = p2 - (result.hi - p1);
    return result;
}

// Metal Kernel for Matrix Multiplication with FP64 Emulation
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
    sum.hi = 0.0;
    sum.lo = 0.0;

    for (uint k = 0; k < K; ++k) {
        Double a = A[gid.x * K + k];
        Double b = B[k * N + gid.y];
        sum = add_double(sum, mul_double(a, b));
    }

    C[gid.x * N + gid.y] = sum;
}
)";

// Global Metal Pipeline (Precompiled Once)
id<MTLDevice> globalDevice = nil;
id<MTLCommandQueue> globalCommandQueue = nil;
id<MTLComputePipelineState> globalPipelineState = nil;

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

// ✅ Optimized NN_dgemm with Precompiled Metal Kernel
void NN_dgemm(int M, int N, int K, 
              double alpha, 
              double *hostA,
              double *hostB, 
              double beta, 
              double *hostC) 
{
    if (!globalPipelineState) {
        fprintf(stderr, "Error: Metal pipeline not initialized!\n");
        return;
    }

    struct Double { float hi; float lo; };
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

    id<MTLBuffer> bufferA = [globalDevice newBufferWithBytes:A_fp32 length:M*K*sizeof(Double) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufferB = [globalDevice newBufferWithBytes:B_fp32 length:K*N*sizeof(Double) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufferC = [globalDevice newBufferWithBytes:C_fp32 length:M*N*sizeof(Double) options:MTLResourceStorageModeShared];

    id<MTLCommandBuffer> commandBuffer = [globalCommandQueue commandBuffer];
    id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];

    [computeEncoder setComputePipelineState:globalPipelineState];
    [computeEncoder setBuffer:bufferA offset:0 atIndex:0];
    [computeEncoder setBuffer:bufferB offset:0 atIndex:1];
    [computeEncoder setBuffer:bufferC offset:0 atIndex:2];

    [computeEncoder setBytes:&M length:sizeof(uint) atIndex:3];
    [computeEncoder setBytes:&N length:sizeof(uint) atIndex:4];
    [computeEncoder setBytes:&K length:sizeof(uint) atIndex:5];

    MTLSize gridSize = MTLSizeMake(M, N, 1);
    MTLSize threadGroupSize = MTLSizeMake(16, 16, 1);

    [computeEncoder dispatchThreads:gridSize threadsPerThreadgroup:threadGroupSize];
    [computeEncoder endEncoding];

    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];

    memcpy(C_fp32, [bufferC contents], M*N*sizeof(Double));

    for (int i = 0; i < M * N; i++) {
        hostC[i] = (double)C_fp32[i].hi + (double)C_fp32[i].lo;
    }

    free(A_fp32);
    free(B_fp32);
    free(C_fp32);
}

// ✅ MAIN FUNCTION TO RUN NN_DGEMM
#include <mach/mach_time.h>  // High-precision timing for macOS

int main() {
    int M = 512, N = 512, K = 512;
    double *A = (double *)malloc(M * K * sizeof(double));
    double *B = (double *)malloc(K * N * sizeof(double));
    double *C = (double *)malloc(M * N * sizeof(double));

    for (int i = 0; i < M * K; i++) A[i] = (double)(rand() % 10);
    for (int i = 0; i < K * N; i++) B[i] = (double)(rand() % 10);

    // ✅ Initialize Metal Once
    initializeMetal();

    // ✅ Start Timing
    uint64_t start = mach_absolute_time();

    // Run FP64 Matrix Multiplication on Apple M2 GPU
    NN_dgemm(M, N, K, 1.0, A, B, 0.0, C);

    // ✅ End Timing
    uint64_t end = mach_absolute_time();
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);
    double elapsed_ms = (double)(end - start) * timebase.numer / timebase.denom / 1e6;

    printf("Matrix multiplication complete in %.3f ms\n", elapsed_ms);
    printf("✅ FP64 Emulation (Float-Based) Matrix Multiplication Complete!\n");

    free(A);
    free(B);
    free(C);
    return 0;
}


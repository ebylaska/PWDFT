#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#import <stdio.h>

// ✅ Metal Kernel Source (Embedded in Code)
const char *metalKernelSource = R"(
#include <metal_stdlib>
using namespace metal;

kernel void matrix_multiply_float(
    const device float *A,
    const device float *B,
    device float *C,
    constant uint &M,
    constant uint &N,
    constant uint &K,
    uint2 gid [[thread_position_in_grid]]) 
{
    if (gid.x >= M || gid.y >= N) return;

    float sum = 0.0;
    for (uint k = 0; k < K; ++k) {
        sum += A[gid.x * K + k] * B[k * N + gid.y];
    }
    C[gid.x * N + gid.y] = sum;
}
)";

// ✅ Global Metal Objects for Precompiled Kernel
id<MTLDevice> device = nil;
id<MTLComputePipelineState> pipelineState = nil;
id<MTLCommandQueue> commandQueue = nil;

// ✅ Function to Compile Metal Kernel Before Calling `NN_dgemm`
void compileMetalKernel() {
    // Get Metal Device
    NSArray<id<MTLDevice>> *devices = MTLCopyAllDevices();
    if (devices.count == 0) {
        fprintf(stderr, "Error: No Metal devices found!\n");
        exit(1);
    }
    device = devices[0];  // Use the first available device
    commandQueue = [device newCommandQueue];

    // Compile Metal Source at Startup
    NSError *error = nil;
    NSString *metalSource = [NSString stringWithUTF8String:metalKernelSource];
    id<MTLLibrary> library = [device newLibraryWithSource:metalSource options:nil error:&error];

    if (!library) {
        fprintf(stderr, "Error compiling Metal source: %s\n", error.localizedDescription.UTF8String);
        exit(1);
    }

    id<MTLFunction> function = [library newFunctionWithName:@"matrix_multiply_float"];
    pipelineState = [device newComputePipelineStateWithFunction:function error:&error];
}

// ✅ Function to Perform Matrix Multiplication Using Precompiled Metal Kernel
void NN_dgemm(int M, int N, int K, 
              float alpha, 
              float *hostA,
              float *hostB, 
              float beta, 
              float *hostC) 
{
    // Ensure Kernel is Compiled
    if (!pipelineState) {
        fprintf(stderr, "Error: Metal kernel not compiled!\n");
        return;
    }

    // Create Buffers on GPU
    id<MTLBuffer> bufferA = [device newBufferWithBytes:hostA length:M*K*sizeof(float) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufferB = [device newBufferWithBytes:hostB length:K*N*sizeof(float) options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufferC = [device newBufferWithBytes:hostC length:M*N*sizeof(float) options:MTLResourceStorageModeShared];

    // Create Command Buffer and Encoder
    id<MTLCommandBuffer> commandBuffer = [commandQueue commandBuffer];
    id<MTLComputeCommandEncoder> computeEncoder = [commandBuffer computeCommandEncoder];

    [computeEncoder setComputePipelineState:pipelineState];
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

    // Copy Data Back to Host
    memcpy(hostC, [bufferC contents], M*N*sizeof(float));
}

// ✅ MAIN FUNCTION TO RUN NN_DGEMM
#include <mach/mach_time.h>  // High-precision timing for macOS

int main() {
    // ✅ Compile Metal Kernel Before Running NN_dgemm
    compileMetalKernel();

    int M = 512, N = 512, K = 512;
    float *A = (float *)malloc(M * K * sizeof(float));
    float *B = (float *)malloc(K * N * sizeof(float));
    float *C = (float *)malloc(M * N * sizeof(float));

    // Initialize Matrices A and B
    for (int i = 0; i < M * K; i++) A[i] = (float)(rand() % 10);
    for (int i = 0; i < K * N; i++) B[i] = (float)(rand() % 10);

    // ✅ Start Timing
    uint64_t start = mach_absolute_time();

    // Call the Precompiled Metal GEMM Function
    NN_dgemm(M, N, K, 1.0f, A, B, 0.0f, C);

    // ✅ End Timing
    uint64_t end = mach_absolute_time();

    // Convert to milliseconds
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);
    double elapsed_ms = (double)(end - start) * timebase.numer / timebase.denom / 1e6;

    printf("Matrix multiplication complete in %.3f ms\n", elapsed_ms);
    printf("✅ Single precision (float) matrix multiplication complete!\n");

    free(A);
    free(B);
    free(C);
    return 0;
}


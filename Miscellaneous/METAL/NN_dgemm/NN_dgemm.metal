#include <metal_stdlib>
using namespace metal;

#pragma clang extension metal::use_native_double : enable

// Metal kernel for double-precision matrix multiplication
kernel void matrix_multiply_double(const device double *A,
                                   const device double *B,
                                   device double *C,
                                   constant uint &M,
                                   constant uint &N,
                                   constant uint &K,
                                   uint2 gid [[thread_position_in_grid]]) 
{
    if (gid.x >= M || gid.y >= N) return;

    double sum = 0.0;
    for (uint k = 0; k < K; ++k) {
        sum += A[gid.x * K + k] * B[k * N + gid.y];
    }
    C[gid.x * N + gid.y] = sum;
}


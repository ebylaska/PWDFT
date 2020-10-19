#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void matmul(const int M, const int N, const int K,
                     const __global double *A, 
                     const __global double *B, 
                     __global double *C) {
    
    // Get the index of the current element
    int i = get_global_id(0);
    int j = get_global_id(1);

    // Do the operation
    double acc = 0.0;
    for (int l=0; l<K; l++) {
       acc += A[i + l*M]*B[l + j*K];
    }
    C[i+j*M] = acc;
}

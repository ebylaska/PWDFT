#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void TN3matmul(const int M, const int N,
                     const __global float *A, 
                     const __global float *B, 
                     __global float *Caa) {
    
    // Get the index of the current element
    int i = get_global_id(0);
    int j = get_global_id(1);

    // Do the operation
    int NN = N*N;
     float aa_acc = 0.0;
    float ab_acc = 0.0;
    float bb_acc = 0.0;
    for (int l=0; l<M; ++l) {
       aa_acc += A[l + i*M]*A[l + j*M];
       ab_acc += A[l + i*M]*B[l + j*M];
       bb_acc += B[l + i*M]*B[l + j*M];
    }
    Caa[i+j*N] = aa_acc;
    Caa[i+j*N+NN] = ab_acc;
    Caa[i+j*N+NN+NN] = bb_acc;
}

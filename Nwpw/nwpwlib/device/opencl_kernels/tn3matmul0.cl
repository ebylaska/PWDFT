#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void TN3matmul(const int M, const int N,
                     const __global double *A, 
                     const __global double *B, 
                     __global double *Caa) {
    
    // Get the index of the current element
    int i = get_global_id(0);
    int j = get_global_id(1);

    // Do the operation
    int NN = N*N;
     double aa_acc = 0.0;
    double ab_acc = 0.0;
    double bb_acc = 0.0;
    for (int l=0; l<M; ++l) {
       aa_acc += A[l + i*M]*A[l + j*M];
       ab_acc += A[l + i*M]*B[l + j*M];
       bb_acc += B[l + i*M]*B[l + j*M];
    }
    Caa[i+j*N] = aa_acc;
    Caa[i+j*N+NN] = ab_acc;
    Caa[i+j*N+NN+NN] = bb_acc;
}

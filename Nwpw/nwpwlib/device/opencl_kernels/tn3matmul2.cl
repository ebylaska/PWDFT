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
    double acc[3] = {0.0, 0.0, 0.0};
    for (int l=0; l<M; ++l) {
       acc[0] += A[l + i*M]*A[l + j*M];
       acc[1] += A[l + i*M]*B[l + j*M];
       acc[2] += B[l + i*M]*B[l + j*M];
    }
    Caa[i+j*N]       = acc[0];
    Caa[i+j*N+NN]    = acc[1];
    Caa[i+j*N+NN+NN] = acc[2];
}

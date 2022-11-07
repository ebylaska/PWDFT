#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void vector_add(__global double *A, __global double *B, __global double *C) {
    
    // Get the index of the current element
    int i = get_global_id(0);

    // Do the operation
    //C[i] = A[i] + B[i];
    //C[i] = cos(3.4934*A[i]) * sin(B[i]*B[i]) + A[i];
     C[i] = sin(A[i]) + cos(B[i]) + 3.141949;
}

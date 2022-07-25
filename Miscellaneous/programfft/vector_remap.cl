#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void butter_xy(int nf, __global int m, int s,  __global double *x, __global double *y) {
    
    // Get the index of the current element
    int p = get_global_id(0);
    int q = get_global_id(1);
    int f = get_global_id(2);

    ar = x[f*nf + 2*(q + s*(p + 0))]; 
    ai = x[f*nf + 2*(q + s*(p + 0))+1];
    br = x[f*nf + 2*(q + s*(p + m))]; b
    bi = x[f*nf + 2*(q + s*(p + m))+1];

    y[f*nf + 2*(q + s*(2*p + 0))]   = ar+br; 
    y[f*nf + 2*(q + s*(2*p + 0))+1] = ai+bi;
    y[f*nf + 2*(q + s*(2*p + 1))]   = (ar-br)*wpr - (ai-bi)*wpi; 
    y[f*nf + 2*(q + s*(2*p + 1))+1] = (ar-br)*wpi + (ai-bi)*wpr;

}

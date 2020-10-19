#include "clFFT.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/////////////////////////////////////
// OpenCL FFT 1D function ///////////
/////////////////////////////////////

int FFT_1D_OpenCL(float *tab[], const char* direction, int size) {

 // Index
 int i;

 // OpenCL variables
 cl_int err;
 cl_platform_id platform = 0;
 cl_device_id device = 0;
 cl_context ctx = 0;
 cl_command_queue queue = 0;

 // Input and Output buffer
 cl_mem buffersIn[2]  = {0, 0};
 cl_mem buffersOut[2] = {0, 0};

 // Temporary buffer
 cl_mem tmpBuffer = 0;

 // Size of temp buffer
 size_t tmpBufferSize = 0;
 int status = 0;
 int ret = 0;

 // Size of FFT
 size_t N = size;

 // FFT library realted declarations
 clfftPlanHandle planHandle;
 clfftDim dim = CLFFT_1D;
 size_t clLengths[1] = {N};

 // Setup OpenCL environment
 err = clGetPlatformIDs(1, &platform, NULL);
 err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);

 // Create an OpenCL context
 ctx = clCreateContext(NULL, 1, &device, NULL, NULL, &err);

 // Create a command queue
 queue = clCreateCommandQueueWithProperties(ctx, device, 0, &err);

 // Setup clFFT
 clfftSetupData fftSetup;
 err = clfftInitSetupData(&fftSetup);
 err = clfftSetup(&fftSetup);

 // Create a default plan for a complex FFT
 err = clfftCreateDefaultPlan(&planHandle, ctx, dim, clLengths);

 // Set plan parameters
 err = clfftSetPlanPrecision(planHandle, CLFFT_SINGLE);
 err = clfftSetLayout(planHandle, CLFFT_COMPLEX_PLANAR, CLFFT_COMPLEX_PLANAR);
 err = clfftSetResultLocation(planHandle, CLFFT_OUTOFPLACE);

 // Bake the plan
 err = clfftBakePlan(planHandle, 1, &queue, NULL, NULL);

 // Real and Imaginary arrays
 cl_float* inReal  = (cl_float*) malloc (N * sizeof (cl_float));
 cl_float* inImag  = (cl_float*) malloc (N * sizeof (cl_float));
 cl_float* outReal = (cl_float*) malloc (N * sizeof (cl_float));
 cl_float* outImag = (cl_float*) malloc (N * sizeof (cl_float));

 // Initialization of inReal, inImag, outReal and outImag
 for(i=0; i<N; i++)
 {
  inReal[i]  = tab[0][i];
  inImag[i]  = 0.0f;
  outReal[i] = 0.0f;
  outImag[i] = 0.0f;
 }

 // Create temporary buffer
 status = clfftGetTmpBufSize(planHandle, &tmpBufferSize);

 if ((status == 0) && (tmpBufferSize > 0)) {
  tmpBuffer = clCreateBuffer(ctx, CL_MEM_READ_WRITE, tmpBufferSize, 0, &err);
  if (err != CL_SUCCESS)
   printf("Error with tmpBuffer clCreateBuffer\n");
 }

 // Prepare OpenCL memory objects : create buffer for input
 buffersIn[0] = clCreateBuffer(ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
   N * sizeof(cl_float), inReal, &err);
 if (err != CL_SUCCESS)
  printf("Error with buffersIn[0] clCreateBuffer\n");

 // Enqueue write tab array into buffersIn[0]
 err = clEnqueueWriteBuffer(queue, buffersIn[0], CL_TRUE, 0, N *
   sizeof(float),
   inReal, 0, NULL, NULL);
 if (err != CL_SUCCESS)
  printf("Error with buffersIn[0] clEnqueueWriteBuffer\n");

 // Prepare OpenCL memory objects : create buffer for input
 buffersIn[1] = clCreateBuffer(ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
   N * sizeof(cl_float), inImag, &err);
 if (err != CL_SUCCESS)
  printf("Error with buffersIn[1] clCreateBuffer\n");

 // Enqueue write tab array into buffersIn[1]
 err = clEnqueueWriteBuffer(queue, buffersIn[1], CL_TRUE, 0, N * sizeof(float),
   inImag, 0, NULL, NULL);
 if (err != CL_SUCCESS)
  printf("Error with buffersIn[1] clEnqueueWriteBuffer\n");

 // Prepare OpenCL memory objects : create buffer for output
 buffersOut[0] = clCreateBuffer(ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, N *
   sizeof(cl_float), outReal, &err);
 if (err != CL_SUCCESS)
  printf("Error with buffersOut[0] clCreateBuffer\n");

 buffersOut[1] = clCreateBuffer(ctx, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, N *
   sizeof(cl_float), outImag, &err);
 if (err != CL_SUCCESS)
  printf("Error with buffersOut[1] clCreateBuffer\n");

 // Execute Forward or Backward FFT
 if(strcmp(direction,"forward") == 0)
 {
  // Execute the plan
  err = clfftEnqueueTransform(planHandle, CLFFT_FORWARD, 1, &queue, 0, NULL, NULL,
    buffersIn, buffersOut, tmpBuffer);
 }
 else if(strcmp(direction,"backward") == 0)
 {
  // Execute the plan
  err = clfftEnqueueTransform(planHandle, CLFFT_BACKWARD, 1, &queue, 0, NULL, NULL,
    buffersIn, buffersOut, tmpBuffer);
 }

 // Wait for calculations to be finished
 err = clFinish(queue);

 // Fetch results of calculations : Real and Imaginary
 err = clEnqueueReadBuffer(queue, buffersOut[0], CL_TRUE, 0, N * sizeof(float), tab[0],
   0, NULL, NULL);
 err = clEnqueueReadBuffer(queue, buffersOut[1], CL_TRUE, 0, N * sizeof(float), tab[1],
   0, NULL, NULL);

 // Release OpenCL memory objects
 clReleaseMemObject(buffersIn[0]);
 clReleaseMemObject(buffersIn[1]);
 clReleaseMemObject(buffersOut[0]);
 clReleaseMemObject(buffersOut[1]);
 clReleaseMemObject(tmpBuffer);

 // Release the plan
 err = clfftDestroyPlan(&planHandle);

 // Release clFFT library
 clfftTeardown();

 // Release OpenCL working objects
 clReleaseCommandQueue(queue);
 clReleaseContext(ctx);

 return ret;
}

int main(void) {

 // Index
 int i;

 // Signal array and FFT output array
 float *Array[2];

 // Number of sampling points
 int size = 100;

 // Cumulative time
 float h = 0;

 // Signal frequency
 float frequency_signal = 10;

 // Sampling frequency : points between 0 and T_signal
 float frequency_sampling = size*frequency_signal;

 // Step = T_sampling
 float step = 1.0/frequency_sampling;

 // File for saving outputs
 FILE *FFT_1D;

 // Allocation of Array
 Array[0] = (float*) malloc(size*sizeof(float));
 Array[1] = (float*) malloc(size*sizeof(float));

 // Initialization of 1D ArrayInput
 FFT_1D = fopen("FFT_1D_OpenCL_Input.dat","w");
 for(i=0; i<size; i++)
 {
  Array[0][i] = cos(2*M_PI*frequency_signal*h);
  Array[1][i] = 0.0f;
  fprintf(FFT_1D,"%f %e\n",i/(frequency_signal*size),Array[0][i]);
  h = h + step;
 }
 fclose(FFT_1D);

 // Perform Forward FFT
 if (FFT_1D_OpenCL(Array,"forward", size) == 0)
  printf("FFT passed !\n");

 // Save Output Array
 FFT_1D = fopen("FFT_1D_OpenCL_Forward.dat","w");
 for (i=0; i<size; i++)
  fprintf(FFT_1D,"%f %e\n", i*frequency_sampling/size, Array[0][i]);
 fclose(FFT_1D);

 // Perform Backward FFT
 if (FFT_1D_OpenCL(Array,"backward", size) == 0)
  printf("IFFT passed !\n");

 // Save Output Array
 FFT_1D = fopen("FFT_1D_OpenCL_Backward.dat","w");
 for (i=0; i<size; i++)
  fprintf(FFT_1D,"%f %e\n", i/(size*frequency_signal), Array[0][i]);
 fclose(FFT_1D);

 return 0;
}

#include "clFFT.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/////////////////////////////////////
// OpenCL FFT 3D function ///////////
/////////////////////////////////////

int FFT_3D_OpenCL(float *tab[], const char* direction, int sizex, int sizey, int sizez) {

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

 // Total size of FFT
 size_t N = sizex*sizey*sizez;

 // FFT library realted declarations
 clfftPlanHandle planHandle;
 clfftDim dim = CLFFT_3D;
 size_t clLengths[3] = {sizex, sizey, sizez};

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

 // Indices
 int i,j,k;

 // Signal array and FFT output array
 float *Array[2];

 // Number of sampling points
 int sizex = 10;
 int sizey = 20;
 int sizez = 40;

 // Total size of FFT
 int N = sizex*sizey*sizez;

 // Cumulative time
 float hx = 0;
 float hy = 0;
 float hz = 0;

 // Signal frequency
 float frequency_signalx = 1;
 float frequency_signaly = 2;
 float frequency_signalz = 4;

 // Sampling frequency : points between 0 and T_signal
 float frequency_samplingx = sizex*frequency_signalx;
 float frequency_samplingy = sizey*frequency_signaly;
 float frequency_samplingz = sizez*frequency_signalz;

 // Step = T_sampling
 float stepx = 1.0/frequency_samplingx;
 float stepy = 1.0/frequency_samplingy;
 float stepz = 1.0/frequency_samplingz;

 // File for saving outputs
 FILE *FFT_3D;

 // Allocation of Array
 Array[0] = (float*) malloc(N*sizeof(float));
 Array[1] = (float*) malloc(N*sizeof(float));

 // Initialization of 2D (real + imaginary) ArrayInput
 FFT_3D = fopen("FFT_3D_OpenCL_Input.dat","w");
 for(k=0; k<sizez; k++)
 {
  for(j=0; j<sizey; j++)
  {
   for(i=0; i<sizex; i++)
   {
    Array[0][k*sizex*sizey+j*sizex+i] = cos(2*M_PI*frequency_signalx*hx)*
     cos(2*M_PI*frequency_signaly*hy)*
     cos(2*M_PI*frequency_signalz*hz);
    Array[1][k*sizex*sizey+j*sizex+i] = 0.0f;
    fprintf(FFT_3D,"%f %f %f %e\n", i/(frequency_signalx*sizex),
      j/(frequency_signaly*sizey),
      k/(frequency_signalz*sizez),
      Array[0][k*sizex*sizey+j*sizex+i]);
    hx = hx + stepx;
   }
   hx = 0.0f;
   hy = hy + stepy;
  }
  hy=0.0f;
  hz = hz + stepz;
 }
 fclose(FFT_3D);

 // Perform Forward FFT
 if (FFT_3D_OpenCL(Array,"forward", sizex, sizey, sizez) == 0)
  printf("FFT passed !\n");

 // Save Output Array
 FFT_3D = fopen("FFT_3D_OpenCL_Forward.dat","w");
 for(k=0; k<sizez; k++)
  for(j=0; j<sizey; j++)
   for(i=0; i<sizex; i++)
    fprintf(FFT_3D,"%f %f %f %e\n", i*frequency_samplingx/sizex,
      j*frequency_samplingy/sizey,
      k*frequency_samplingz/sizez,
      Array[0][k*sizex*sizey+j*sizex+i]);
 fclose(FFT_3D);

 // Perform Backward FFT
 if (FFT_3D_OpenCL(Array,"backward", sizex, sizey, sizez) == 0)
  printf("IFFT passed !\n");

 // Save Output Array
 FFT_3D = fopen("FFT_3D_OpenCL_Backward.dat","w");
 for(k=0; k<sizez; k++)
  for(j=0; j<sizey; j++)
   for(i=0; i<sizex; i++)
    fprintf(FFT_3D,"%f %f %f %e\n", i/(frequency_signalx*sizex),
      j/(frequency_signaly*sizey),
      k/(frequency_signalz*sizez),
      Array[0][k*sizex*sizey+j*sizex+i]);
 fclose(FFT_3D);

 return 0;
}

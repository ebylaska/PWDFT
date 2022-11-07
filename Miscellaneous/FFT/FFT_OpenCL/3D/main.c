#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

/////////////////////////////////////
// fftw3 FFT 3D function ////////////
/////////////////////////////////////

int main(void) {

 // Indices
 int i,j,k;

 // Signal 1D signal array and FFT 1D output array
 fftw_complex *Array1DIn;
 fftw_complex *Array1DOut;

 // Number of sampling points
 int sizex = 10;
 int sizey = 20;
 int sizez = 40;

 // Total size of FFT
 int N = sizex*sizey*sizez;

 // Cumulative time
 double hx = 0.0f;
 double hy = 0.0f;
 double hz = 0.0f;

 // Signal frequency
 double frequency_signalx = 1;
 double frequency_signaly = 2;
 double frequency_signalz = 4;

 // Sampling frequency
 double frequency_samplingx = sizex*frequency_signalx;
 double frequency_samplingy = sizey*frequency_signaly;
 double frequency_samplingz = sizez*frequency_signalz;

 // Step = T_sampling
 double stepx = 1.0/frequency_samplingx;
 double stepy = 1.0/frequency_samplingy;
 double stepz = 1.0/frequency_samplingz;

 // File for saving outputs
 FILE *FFT_3D;

 // FFTW plan signal
 fftw_plan fft_signal;

 // Allocation of Arrays
 Array1DIn = (fftw_complex*) fftw_malloc(N*sizeof(fftw_complex));
 Array1DOut = (fftw_complex*) fftw_malloc(N*sizeof(fftw_complex));

 // Initialization of 2D (real + imaginary) ArrayInput
 FFT_3D = fopen("FFT_3D_Input.dat","w");
 for(k=0; k<sizez; k++)
 {
  for(j=0; j<sizey; j++)
  {
   for(i=0; i<sizex; i++)
   {
    Array1DIn[k*sizex*sizey+j*sizex+i][0] = cos(2*M_PI*frequency_signalx*hx)*
     cos(2*M_PI*frequency_signaly*hy)*
     cos(2*M_PI*frequency_signalz*hz);
    Array1DIn[k*sizex*sizey+j*sizex+i][1] = 0.0f;
    fprintf(FFT_3D,"%f %f %f %e\n", i/(frequency_signalx*sizex),
      j/(frequency_signaly*sizey),
      k/(frequency_signalz*sizez),
      Array1DIn[k*sizex*sizey+j*sizex+i][0]);
    hx = hx + stepx;
   }
   hx = 0.0f;
   hy = hy + stepy;
  }
  hy = 0.0f;
  hz = hz + stepz;
 }
 fclose(FFT_3D);

 // FFTW plan Forward 3D
 fft_signal = fftw_plan_dft_3d(sizez, sizey, sizex, Array1DIn, Array1DOut, FFTW_FORWARD, FFTW_ESTIMATE);
 // Perform Forward FFT
 fftw_execute(fft_signal);
 printf("FFT passed !\n");

 // Save Output Array
 FFT_3D = fopen("FFT_3D_Forward.dat","w");
 for(k=0; k<sizez; k++)
  for(j=0; j<sizey; j++)
   for(i=0; i<sizex; i++)
    fprintf(FFT_3D,"%f %f %f %e\n", i*frequency_samplingx/sizex,
      j*frequency_samplingy/sizey,
      k*frequency_samplingz/sizez,
      Array1DOut[k*sizex*sizey+j*sizex+i][0]);
 fclose(FFT_3D);

 // FFTW plan Backward 3D
 fft_signal = fftw_plan_dft_3d(sizez, sizey, sizex, Array1DOut, Array1DIn, FFTW_BACKWARD, FFTW_ESTIMATE);
 // Perform Backward FFT
 fftw_execute(fft_signal);
 printf("IFFT passed !\n");

 // Destroy plan signal
 fftw_destroy_plan(fft_signal);

 // Normalization
 for(k=0; k<sizez; k++)
  for(j=0; j<sizey; j++)
   for(i=0; i<sizex; i++)
    Array1DIn[k*sizex*sizey+j*sizex+i][0] = Array1DIn[k*sizex*sizey+j*sizex+i][0]/(sizex*sizey*sizez);

 // Save Output Array
 FFT_3D = fopen("FFT_3D_Backward.dat","w");
 for(k=0; k<sizez; k++)
  for(j=0; j<sizey; j++)
   for(i=0; i<sizex; i++)
    fprintf(FFT_3D,"%f %f %f %e\n", i/(sizex*frequency_signalx),
      j/(sizey*frequency_signaly),
      k/(sizez*frequency_signalz),
      Array1DIn[k*sizex*sizey+j*sizex+i][0]);
 fclose(FFT_3D);

 return 0;
}

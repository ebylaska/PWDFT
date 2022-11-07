#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>

/////////////////////////////////////
// fftw3 FFT 1D function ////////////
/////////////////////////////////////

int main(void) {

 // Index
 int i;

 // Signal 1D signal array and FFT 1D output array
 fftw_complex *Array1DIn;
 fftw_complex *Array1DOut;

 // Number of sampling points
 int size = 100;

 // Cumulative time
 double h = 0;

 // Signal frequency
 double frequency_signal = 10;

 // Sampling frequency
 double frequency_sampling = size*frequency_signal;

 // step = T_sampling
 double step = 1.0/frequency_sampling;

 // File for saving outputs
 FILE *FFT_1D;

 // FFTW plan signal
 fftw_plan fft_signal;

 // Allocation of Arrays
 Array1DIn = (fftw_complex*) fftw_malloc(size*sizeof(fftw_complex));
 Array1DOut = (fftw_complex*) fftw_malloc(size*sizeof(fftw_complex));

 // FFTW plan Forward dft 1D
 fft_signal = fftw_plan_dft_1d(size, Array1DIn, Array1DOut, FFTW_FORWARD, FFTW_ESTIMATE);

 // Initialization of 2D (real + imaginary) ArrayInput
 FFT_1D = fopen("FFT_1D_Input.dat","w");
 for(i=0; i<size; i++)
 {
  Array1DIn[i][0] = cos(2*M_PI*frequency_signal*h);
  Array1DIn[i][1] = 0.0f;
  fprintf(FFT_1D,"%f %e\n", i/(frequency_signal*size), Array1DIn[i][0]);
  h = h + step;
 }
 fclose(FFT_1D);

 // Perform Forward FFT
 fftw_execute(fft_signal);
 printf("FFT passed !\n");

 // Save Output Array
 FFT_1D = fopen("FFT_1D_Forward.dat","w");
 for(i=0; i<size; i++)
  fprintf(FFT_1D,"%f %e\n", i*frequency_sampling/size, Array1DOut[i][0]);
 fclose(FFT_1D);

 // FFTW plan Backward dft 1D
 fft_signal = fftw_plan_dft_1d(size, Array1DOut, Array1DIn, FFTW_BACKWARD, FFTW_ESTIMATE);

 // Perform Forward FFT
 fftw_execute(fft_signal);
 printf("IFFT passed !\n");

 // Destroy plan signal
 fftw_destroy_plan(fft_signal);

 // Normalization
 for(i=0; i<size; i++)
  Array1DIn[i][0] = Array1DIn[i][0]/size;

 // Save Output Array
 FFT_1D = fopen("FFT_1D_Backward.dat","w");
 for(i=0; i<size; i++)
  fprintf(FFT_1D,"%f %e\n", i/(size*frequency_signal), Array1DIn[i][0]);
 fclose(FFT_1D);

 return 0;
}

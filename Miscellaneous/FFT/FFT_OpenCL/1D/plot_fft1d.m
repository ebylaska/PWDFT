%%%% Plot Input Signal, Forward FFT and Backward FFT (OpenCL version) %%%%

%%%% Parameters must be identical to main_OpenCL.c %%%%
frequency_signal=10;
number_sampling=100;
frequency_sampling=frequency_signal*number_sampling;
x_temp1=(1:number_sampling)/number_sampling*1/frequency_signal;
y_signal1=cos(2*pi*frequency_signal*x_temp1);

%%%% Plot Input Signal %%%%
subplot(131);
stem(x_temp1,y_signal1);
title('Input Signal');

%%%% Plot Forward FFT %%%%
subplot(132);
data=load('FFT_1D_OpenCL_Forward.dat');
x_frequency=data(:,1);
y_fft=data(:,2);
stem(x_frequency,y_fft);
title('Forward FFT');

%%%% Plot Backward FFT %%%%
subplot(133);
data=load('FFT_1D_OpenCL_Backward.dat');
x_temp2=data(:,1);
y_signal2=data(:,2);
stem(x_temp2,y_signal2);
title('Backward FFT');

%%%% Plot Input Signal, Forward FFT and Backward FFT (OpenCL version) %%%%

%%%% Parameters must be identical to main_OpenCL.c %%%%
frequency_signalx=10;
frequency_signaly=20;
number_samplingx=100;
number_samplingy=200;
frequency_samplingx=frequency_signalx*number_samplingx;
frequency_samplingy=frequency_signaly*number_samplingy;

%%%% Plot Input Signal %%%%
subplot(131);
x11=(1:number_samplingx)/number_samplingx*1/frequency_signalx;
x12=(1:number_samplingy)/number_samplingy*1/frequency_signaly;
[x y]=meshgrid(x11,x12);
y_signal1=cos(2*pi*frequency_signalx*x).*cos(2*pi*frequency_signaly*y);
z=reshape(y_signal1,number_samplingy,number_samplingx);
surf(x,y,z,'Facecolor','interp');
shading interp;
title('Input Signal');

%%%% Plot OpenCL Forward FFT %%%%
subplot(132);
data1=load('FFT_2D_OpenCL_Forward.dat');
x_frequency=reshape(data1(:,1),number_samplingx,number_samplingy);
y_frequency=reshape(data1(:,2),number_samplingx,number_samplingy);
z_fft=reshape(data1(:,3),number_samplingx,number_samplingy);
surf(x_frequency,y_frequency,z_fft);
shading interp;
title('OpenCL Forward FFT');

%%%% Plot OpenCL Forward + Backward FFT %%%%
subplot(133)
data2=load('FFT_2D_OpenCL_Backward.dat');
x21=reshape(data2(:,1),number_samplingx,number_samplingy);
x22=reshape(data2(:,2),number_samplingx,number_samplingy);
y_signal2=reshape(data2(:,3),number_samplingx,number_samplingy);
surf(x21,x22,y_signal2);
shading interp;
title('OpenCL Forward + Backward FFT');


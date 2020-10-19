%%%% Plot Input Signal, Forward FFT and Backward FFT %%%%

%%%% Parameters must be identical to main_OpenCL.c %%%%
frequency_signalx=1;
frequency_signaly=2;
frequency_signalz=4;
tperiod_signalx=1/frequency_signalx;
tperiod_signaly=1/frequency_signaly;
tperiod_signalz=1/frequency_signalz;

number_samplingx=10;
number_samplingy=20;
number_samplingz=40;
frequency_samplingx=frequency_signalx*number_samplingx;
frequency_samplingy=frequency_signaly*number_samplingy;
frequency_samplingz=frequency_signalz*number_samplingz;

%%%% Plot Input Signal from Matlab %%%%
x11=(0:number_samplingx-1)/number_samplingx*tperiod_signalx;
x12=(0:number_samplingy-1)/number_samplingy*tperiod_signaly;
x13=(0:number_samplingz-1)/number_samplingz*tperiod_signalz;
[ym1 xm1 zm1]=ndgrid(x12,x11,x13);
y_signal1=cos(2*pi*frequency_signalx*xm1).*cos(2*pi*frequency_signaly*ym1).*cos(2*pi*frequency_signalz*zm1);
hFig = figure(1);
set(hFig, 'Position', [400 400 750 600]);
slice(xm1,ym1,zm1,y_signal1,[tperiod_signalx/2 tperiod_signalx-tperiod_signalx/number_samplingx],...
                            [tperiod_signaly/2 tperiod_signaly-tperiod_signaly/number_samplingy],...
                            [0 tperiod_signalz/2]);
shading faceted;
view([-42,22]);
hc=colorbar;
set(hc,'position',[0.932 0.3 0.02 0.6]);
title('Slice of Input Signal From Matlab');
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');

%%%% Plot Input Signal from FFT_3D_Input.dat %%%%
data=load('FFT_3D_Input.dat');
xm2=reshape(data(:,1),number_samplingx,number_samplingy,number_samplingz);
ym2=reshape(data(:,2),number_samplingx,number_samplingy,number_samplingz);
zm2=reshape(data(:,3),number_samplingx,number_samplingy,number_samplingz);
y_signal2=data(:,4);
z_fft1=reshape(y_signal2,number_samplingx,number_samplingy,number_samplingz);
P = [2 1 3];
xm2= permute(xm2,P);
ym2= permute(ym2,P);
zm2= permute(zm2,P);
z_fft1= permute(z_fft1,P);
hFig = figure(2);
set(hFig, 'Position', [400 400 750 600]);
slice(xm2,ym2,zm2,z_fft1,[tperiod_signalx/2 tperiod_signalx-tperiod_signalx/number_samplingx],...
                         [tperiod_signaly/2  tperiod_signaly-tperiod_signaly/number_samplingy],...
                         [0 tperiod_signalz/2]);
shading faceted;
view([-42,22]);
hc=colorbar;
set(hc,'position',[0.932 0.3 0.02 0.6]);
title('Slice of Input Signal From FFT3D Input.dat');
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');

%%%% Plot OpenCL Forward FFT %%%%
data=load('FFT_3D_OpenCL_Forward.dat');
x=data(:,1);
y=data(:,2);
z=data(:,3);
z_fft=data(:,4);
z_fft_max=max(z_fft);
hFig = figure(3);
set(hFig, 'Position', [400 400 750 600]);
view([-42,22]);
scatter3(x,y,z,15,z_fft,'filled');
colormap(jet);
hc=colorbar;
set(hc,'position',[0.932 0.3 0.02 0.6]);
x(z_fft<0.9*z_fft_max)=NaN;
y(z_fft<0.9*z_fft_max)=NaN;
z(z_fft<0.9*z_fft_max)=NaN;
z_fft(z_fft<0.9*z_fft_max)=NaN;
hold on;
scatter3(x,y,z,250,z_fft,'filled')
title('OpenCL Forward FFT');
xlabel('x frequency');
ylabel('y frequency');
zlabel('z frequency');

%%%% Plot OpenCL Forward + Backward FFT %%%%
data=load('FFT_3D_OpenCL_Backward.dat');
x21=(0:number_samplingx-1)*tperiod_signalx/number_samplingx;
x22=(0:number_samplingy-1)*tperiod_signaly/number_samplingy;
x23=(0:number_samplingz-1)*tperiod_signalz/number_samplingz;
[xm2 ym2 zm2]=meshgrid(x21,x22,x23);
y_signal2=data(:,4);
z_fft2=reshape(y_signal2,number_samplingx,number_samplingy,number_samplingz);
P = [2 1 3];
z_fft2= permute(z_fft2,P);
hFig = figure(4);
set(hFig, 'Position', [400 400 750 600]);
slice(xm2,ym2,zm2,z_fft2,[tperiod_signalx/2 tperiod_signalx-tperiod_signalx/number_samplingx],...
                         [tperiod_signaly/2 tperiod_signaly-tperiod_signaly/number_samplingy],...
                         [0 tperiod_signalz/2]);
shading faceted;
view([-42,22]);
hc=colorbar;
set(hc,'position',[0.932 0.3 0.02 0.6]);
title('OpenCL Forward + Backward FFT');
xlabel('x axis');
ylabel('y axis');
zlabel('z axis');

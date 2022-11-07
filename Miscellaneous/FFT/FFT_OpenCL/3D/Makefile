####### Compiler, tools and options

CC            = gcc
CXX           = g++-4.4
CFLAGS        = -m64 -pipe -O2 -Wno-unused-but-set-parameter -Wall
CXXFLAGS      = -m64 -pipe -O2 -Wno-unused-but-set-parameter -Wall
OPENCL        = /opt/AMDAPP
CLFFT         = /opt/clFFT
INC_OPENCL    = -I$(OPENCL)/include
INC_CLFFT     = -I$(CLFFT)/src/include
INCPATH       = $(INC_OPENCL) $(INC_CLFFT)
LIB_FFTW      = -lfftw3
LIB_CLFFT     = -L$(CLFFT)/src/library -lclFFT -L$(OPENCL)/lib/x86_64/sdk -lOpenCL
LIB_MATH      = -lm
DEL_FILE      = rm -f

####### Build rules

all: FFT_3D FFT_3D_OpenCL

# Classic FFT with fftw
FFT_3D: main.c
	$(CC) $(CFLAGS) $(LIB_FFTW) $(LIB_MATH) -o FFT_3D main.c

# OpenCL FFT
FFT_3D_OpenCL: main_OpenCL.c
	$(CC) $(CFLAGS) $(INCPATH) $(LIB_CLFFT) $(LIB_MATH) -o FFT_3D_OpenCL main_OpenCL.c

clean:
	$(DEL_FILE) FFT_3D FFT_3D_OpenCL
	$(DEL_FILE) *.o

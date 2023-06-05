# PWDFT
PW-DFT development for NWChemEx

# CMAKE - Generate a Project Buildsystem
```
cmake -S Nwpw/ -B build
cd build
make

Alternatively you can build :
mkdir build
cd build
cmake ../Nwpw
make

```

Standard cmake build commands
```
cmake [<options>] <path-to-source>
$ mkdir build ; cd build
$ cmake ../src
 
cmake [<options>] -S <path-to-source> -B <path-to-build>
$ cmake -S src -B build 

cmake [<options>] <path-to-existing-build>
$ cd build
$ cmake .
```

# Timings 
## laptop timings
| machine     | ncpus  |   cputime | non-local |       ffm |       fmf |       fft | diagonalize |  mmm_mult |
| :----:      | :----: |       ---:|        --:|        --:|        --:|      ---: |          --:|        --:|
| **QA/CCO-Cu_surface**
| mac-m1      | 1      | 4.763e+01 | 1.413e+01 | 1.794e+01 | 1.217e+01 | 1.353e+00 | 3.725e-02 | 1.642e-01 |
| mac-m1      | 2      | 2.523e+01 | 7.349e+00 | 9.440e+00 | 6.547e+00 | 7.941e-01 | 4.310e-02 | 1.684e-01 |
| mac-m1      | 4      | 1.405e+01 | 3.887e+00 | 4.742e+00 | 4.125e+00 | 4.404e-01 | 4.183e-02 | 1.715e-01 |
| mac-m1      | 6      | 1.018e+01 | 2.415e+00 | 3.311e+00 | 3.449e+00 | 3.274e-01 | 4.599e-02 | 1.785e-01 | 
| mac-m1      | 8      | 8.539e+00 | 1.848e+00 | 2.417e+00 | 3.187e+00 | **2.792e-01** | 4.307e-02 | 1.605e-01 |
|      |      |  |  |  |  |  |  | |
| WE45090     | 1      | 5.029e+01 | 1.198e+01 | 1.262e+01 | 2.359e+01 | 1.780e+00 | 3.892e-02 | 7.100e-02 |
| WE45090     | 2      | 2.701e+01 | 6.422e+00 | 6.663e+00 | 1.270e+01 | 9.692e-01 | 4.015e-02 | 7.448e-02 |
| WE45090     | 4      | 1.900e+01 | 4.062e+00 | 5.126e+00 | 8.980e+00 | 6.095e-01 | 4.154e-02 | 9.022e-02 |
| WE45090     | 6      | 1.795e+01 | 3.499e+00 | 5.071e+00 | 8.675e+00 | 4.499e-01 | 4.439e-02 | 8.899e-02 |
| WE45090     | 8      | 1.756e+01 | 3.279e+00 | 4.983e+00 | 8.534e+00 |**4.010e-01**| 4.633e-02 | 9.845e-02 |
|      |      |  |  |  |  |  |  | |
| **WE45090-GPU** | **1**  | **5.878e+00** | **1.050e+00** |**7.224e-01**  |**1.180e+00**| 2.614e+00 | 3.999e-02 | 7.225e-02 |
| WE45090-GPU | 2      | 5.928e+00 | 1.053e+00 | 1.062e+00 | 1.346e+00 | 2.217e+00 | 4.219e-02 | 7.679e-02 |
| WE45090-GPU | 4      | 6.895e+00 | 1.078e+00 | 1.193e+00 | 1.459e+00 | 2.922e+00 | 4.405e-02 | 8.771e-02 |
| WE45090-GPU | 6      | 7.884e+00 | 1.037e+00 | 9.116e-01 | 1.486e+00 | 4.183e+00 | 4.977e-02 | 1.044e-01 |
| WE45090-GPU | 8      |  |  |  |  |  |    |  |
|      |      |  |  |  |  |  |  | |
| perlmutter-GPU | 1   | 2.924e+00 | 3.580e-01 | 3.084e-02 | 3.595e-02 | 2.188e+00 | 6.219e-03  |  |
| perlmutter-GPU | 2   | 1.639e+00 | 2.246e-01 | 2.258e-02 | 2.222e-02 | 1.219e+00 | 5.878e-03  |  |
| perlmutter-GPU | 3   | 1.246e+00 | 1.575e-01 | 2.566e-02 | 1.520e-02 | 9.523e-01 | 4.568e-03  |  |
| perlmutter-GPU | 4   | 1.080e+00 | 1.657e-01 | 1.588e-02 | 1.874e-02 | 7.948e-01 | 6.037e-03  |  |
|      |      |  |  |  |  |  |  | |
| perlmutter-CPU | 1   | 1.750e+01 | 7.615e+00 | 4.218e+00 | 2.856e+00 | 2.000e+00 | 5.851e-03 |  |
| perlmutter-CPU | 2   | 1.018e+01 | 3.822e+00 | 2.135e+00 | 1.371e+00 | 2.127e+00 | 6.039e-03 |  |
| perlmutter-CPU | 3   | 6.945e+00 | 2.520e+00 | 1.351e+00 | 9.250e-01 | 1.499e+00 | 6.011e-03 |  |
| perlmutter-CPU | 4   | 5.765e+00 | 2.103e+00 | 1.067e+00 | 7.114e-01 | 1.254e+00 | 6.465e-03 |  |
| perlmutter-CPU | 6   | 3.753e+00 | 1.138e+00 | 6.549e-01 | 3.859e-01 | 8.889e-01 | 6.006e-03 |  |
| perlmutter-CPU | 8   | 2.739e+00 | 8.779e-01 | 5.078e-01 | 2.521e-01 | 4.427e-01 | 6.179e-03 |  |
| perlmutter-CPU | 16  | 1.789e+00 | 4.439e-01 | 2.408e-01 | 1.404e-01 | 3.154e-01 | 7.301e-03 |  |
| perlmutter-CPU | 32  | 1.343e+00 | 2.417e-01 | 1.227e-01 | 7.259e-02 | 2.748e-01 | 7.184e-03 |  |
| perlmutter-CPU | 64  | 1.126e+00 | 1.514e-01 | 9.311e-02 | 3.653e-02 | 1.691e-01 | 7.429e-03 |  |
| perlmutter-CPU | 128 | 1.335e+00 | 1.298e-01 | 8.846e-02 | 2.510e-02 | 1.530e-01 | 9.429e-03 |  |

These timings suggest that parallel FFTs should be implemented using hybrid MPI-OpenMP code, and the large DGEMMs should use GPUs.  This is somewhat justified, since the cost of parallel FFTs is mostly due to data movement, i.e. FFTs are memory bound rather than computationlly bound.  However, we need to test the competiveness of pipelining FFT data to GPUs, and using Stockholm FFT kernels (https://github.com/ebylaska/PWDFT/tree/master/Miscellaneous/programfft), versus an MPI-only algorithm.



# Compiling and Running Instructions 

## Build instructions on Sunspot

### Required Modules

```
module add gcc/11.2.0
module add cray-libpals/1.2.3
module add intel_compute_runtime/release/pvc-prq-66
module add oneapi-prgenv/2022.10.15.006.001
module add mpich/50.1/icc-all-pmix-gpu
module add spack/linux-sles15-x86_64-ldpath
module add oneapi/eng-compiler/2022.10.15.006
module add ncurses/6.1.20180317-gcc-11.2.0-zedoshf
module add libfabric/1.15.2.0
module add openssl/1.1.1d-gcc-11.2.0-amlvxob
module add cray-pals/1.2.3
module add cmake/3.24.2-gcc-11.2.0-pcasswq
```

### Getting the code and building instrunctions

```
export HTTP_PROXY=http://proxy.alcf.anl.gov:3128
export HTTPS_PROXY=http://proxy.alcf.anl.gov:3128
export http_proxy=http://proxy.alcf.anl.gov:3128
export https_proxy=http://proxy.alcf.anl.gov:3128
git config --global http.proxy http://proxy.alcf.anl.gov:3128
git clone https://github.com/alvarovm/PWDFT.git

cd PWDFT
cmake -H. -Bbuild_sycl -DNWPW_SYCL=On -DCMAKE_CXX_COMPILER=dpcpp ./Nwpw
cd build_sycl
make 
```
### Running
```
qsub -l select=1 -l walltime=30:00 -A Aurora_deployment -q workq -I
qsub -l select=1 -l walltime=30:00 -A catalysis_aesp_CNDA -q workq -I
```

## Examples on JSLE -  `SYCL` backend
### Required Modules
```
export MODULEPATH=$MODULEPATH:/soft/modulefiles:/soft/restricted/CNDA/modules
module load oneapi
module load cmake
```

### Build Instructions (for `SYCL` backend)
```
cd PWDFT
cmake -H. -Bbuild_sycl -DNWPW_SYCL=On -DCMAKE_CXX_COMPILER=dpcpp ./Nwpw
make -j4
```

### Running on Articus
```
qsub -I -n 1 -t 60 -q arcticus
```

## CUDA backend
### Required Modules
```
export MODULEPATH=/soft/modulefiles:/usr/share/Modules/modulefiles:/etc/modulefiles:/usr/share/modulefiles
module add cmake/3.20.3   cuda/11.6.2    gcc/9.5.0
module add openmpi/4.1.1-gcc
```

### Build Instructions (for CUDA backend)
```
cd PWDFT
mkdir build_cuda
cd build_cuda
cmake -DNWPW_CUDA=ON ../Nwpw/
make -j4
```

### Running on JSLE in V100
```
qsub  -I -t 30  -n 1 -q gpu_v100_smx2_debug
```


##  Build Instructions for NERSC Perlmutter
<details>
<summary>Toggle for Details</summary>

CMake Build/Install
```
cmake -H. -Bbuild_cuda -DNWPW_CUDA=ON ./Nwpw -DCMAKE_CUDA_ARCHITECTURES=80 -DCUDA_cublas_LIBRARY=${CRAY_CUDATOOLKIT_DIR}/../../math_libs/${CRAY_CUDATOOLKIT_VERSION#*_}/lib64/libcublas.so -DCUDA_cufft_LIBRARY=${CRAY_CUDATOOLKIT_DIR}/../../math_libs/${CRAY_CUDATOOLKIT_VERSION#*_}/lib64/libcufft.so -DCUDA_cusolver_LIBRARY=${CRAY_CUDATOOLKIT_DIR}/../../math_libs/${CRAY_CUDATOOLKIT_VERSION#*_}/lib64/libcusolver.so
cd build_cuda
make
```
</details>

##  Build Instructions on NERSC Cori-Haswell

### required modules on Cori-Haswell
```
module purge
module load cmake/3.20.2
module load PrgEnv-intel
module load craype-haswell
module load openmpi
```

module purge;
module load cmake/3.20.2;
module load PrgEnv-intel;
module load craype-haswell;
module load openmpi

### Build Instructions on Cori-Haswell (starting from PWDFT directory)
```
mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=CC ../Nwpw/
```

### Running on Cori-Haswell
```
salloc --nodes 1 --qos interactive --time 01:00:00 --constraint haswell
srun -n <num_mpi_processes> -c <cpus_per_task> a.out
srun -n <num_mpi_processes> -c <cpus_per_task> pwdft
```

### runnin on haswell??
```
salloc --nodes 1 --qos interactive --time 01:00:00 --constraint haswell
cd PWDFT/QA/C2_steepest_descent
srun -n 24 ../../build/pwdft c2-sd.nw
```

## Build instructions on Polaris-CUDA

### Required Modules on Polaris-CUDA
```
module purge
module load cmake cpe-cuda aocl
module load PrgEnv-nvhpc
module unload craype-accel-nvidia80

```

### Build Instructions on Polaris-CUDA (starting from PWDFT directory)
```
mkdir build_cuda
cd build-cuda
cmake -DNWPW_CUDA=ON  -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_Fortran_COMPILER=ftn ../Nwpw/
```
### Running on Polaris-CUDA
```
qsub -q  debug -I -l walltime=01:00:00 -lselect=1 -A myproject -l filesystems=home:eagle:grand
module purge
module load cmake cpe-cuda aocl
module load PrgEnv-nvhpc
module unload craype-accel-nvidia80

mpiexec -n 2 --ppn 2 --cpu-bind=verbose --cpu-bind depth--env CUDA_VISIBLE_DEVICES=2 ./pwdft nwinput.nw

```

## Build instructions on NERSC Cori-CUDA

### Required Modules on Cori-CUDA
```
module unload impi
module load PrgEnv-intel
module load cmake
module load cudatoolkit
```

### Build Instructions on Cori-CUDA (starting from PWDFT directory)
```
mkdir build_cuda
cd build-cuda
cmake -DNWPW_CUDA=ON ../Nwpw/
```

### Running on Cori-CUDA
```
module load cgpu
salloc -C gpu -t 60 -c 10 -G 1 -q interactive -A <account>
salloc -C gpu -t 60 -c 10 -G 1 -q interactive -A mp119
```

# Making shared library
To generate a library clean the build directory and then regenerate cmake with
```
cmake ../Nwpw -DMAKE_LIBRARY=true
```
on Linux use: 
```
cmake ../Nwpw -DMAKE_LIBRARY=true -DCMAKE_POSITION_INDEPENDENT_CODE=ON
```
Compile and generate the shared library in the build directory.
```
make
```

The shared library, libpwdft.dylib, should be generated and present in the build directory.
```
prompt% ls
CMakeCache.txt       Makefile             cmake_install.cmake  nwpwlib/
CMakeFiles/          NwpwConfig.h         libpwdft.dylib*      pspw/
```


Example header to make function calls
```
#include <string>
#include "mpi.h"

namespace pwdft {
using namespace pwdft;

extern char *util_date();
extern void seconds(double *);
extern int cpsd(MPI_Comm, std::string&);
extern int cpmd(MPI_Comm, std::string&);
extern int pspw_minimizer(MPI_Comm, std::string&);
extern int pspw_geovib(MPI_Comm, std::string&);
}
```

Example function call
```
ierr += pwdft::pspw_geovib(MPI_COMM_WORLD,nwinput);
```

## Using shared library on a MAC
```
mpic++ test.cpp ../build-shared/libpwdft.dylib
setenv DYLD_LIBRARY_PATH /Users/bylaska/Codes/PWDFT/build-shared
a.out
```

## Using shared library on LINUX
```
mpic++ test.cpp ../build-shared/libpwdft.so
setenv DYLD_LIBRARY_PATH /Users/bylaska/Codes/PWDFT/build-shared
a.out
```

 [2:17 PM] Bagusetty, Abhishek
 add_library(pwdft SHARED nwpw.cpp)

 [2:18 PM] Bagusetty, Abhishek
 CMakeLists.txt (right after this line add_executable(pwdft nwpw.cpp))

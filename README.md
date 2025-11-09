# PWDFT
PW-DFT development for NWChemEx

Web location:
https://ebylaska.github.io/PWDFT/


<img width="1044" height="1814" alt="image" src="https://github.com/user-attachments/assets/623b3dc8-5882-4036-a824-b6ae3ac484ad" />


# ![PWDFT-QR Code](https://raw.githubusercontent.com/ebylaska/PWDFT/gh-pages/qr-code-pwdft.png)



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
| WE45090     | 1      | 2.786e+01 | 1.298e+01 | 6.825e+00 | 4.932e+00 | 2.154e+00 | 6.489e-03 | 7.100e-02 |
| WE45090     | 2      | 1.514e+01 | 6.859e+00 | 3.674e+00 | 2.543e+00 | 1.136e+00 | 6.998e-03 | 7.448e-02 |
| WE45090     | 4      | 1.055e+01 | 4.408e+00 | 2.725e+00 | 1.775e+00 | 6.635e-01 | 6.676e-03 | 9.022e-02 |
| WE45090     | 6      | 9.858e+00 | 3.886e+00 | 2.660e+00 | 1.709e+00 | 4.967e-01 | 7.285e-03 | 8.899e-02 |
| WE45090     | 8      | 9.450e+00 | 3.542e+00 | 2.594e+00 | 1.669e+00 | 4.279e-01 | 7.650e-03 | 9.845e-02 |
|      |      |  |  |  |  |  |  | |
| WE45090-GPU | 1      | 4.773e+00 | 1.140e+00 | 7.054e-01 | 2.313e-01 | 2.481e+00 | 6.870e-03 | 7.225e-02 |
| **WE45090-GPU** | 2  | 4.387e+00 | 1.037e+00 | 9.929e-01 | 2.427e-01 | 1.911e+00 | 7.432e-03 | 7.679e-02 |
| WE45090-GPU | 4      | 4.164e+00 | 9.757e-01 | 9.960e-01 | 2.452e-01 | 1.667e+00 | 8.277e-03 | 8.771e-02 |
| WE45090-GPU | 6      | 4.765e+00 | 1.076e+00 | 8.142e-01 | 2.564e-01 | 2.217e+00  | 7.819e-03 | 1.044e-01 |
| WE45090-GPU | 8      |  |  |  |  |  |    |  |
|      |      |  |  |  |  |  |  | |
| WE45090-nwchem | 8   | 5.194e+00 |  |  |  |  |    |  |
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

## Build instructions on ALCF Aurora/Sunspot

### Required Modules

```
module restore
module load cmake
```
### Getting the code and building instructions

```
export HTTP_PROXY=http://proxy.alcf.anl.gov:3128
export HTTPS_PROXY=http://proxy.alcf.anl.gov:3128
export http_proxy=http://proxy.alcf.anl.gov:3128
export https_proxy=http://proxy.alcf.anl.gov:3128
git config --global http.proxy http://proxy.alcf.anl.gov:3128
git clone https://github.com/ebylaska/PWDFT.git

cd PWDFT
cmake -H. -Bbuild_sycl -DNWPW_SYCL=On -DCMAKE_CXX_COMPILER=icpx -DCMAKE_C_COMPILER=icx -DCMAKE_Fortran_COMPILER=ifx ./Nwpw
```
### Running
```
qsub -l select=1 -l walltime=30:00 -A catalysis_aesp_CNDA -q lustre_scaling -I
qsub -l select=4 -l walltime=30:00 -l filesystems=flare -A ExaCatChem -q debug-scaling -I  
```
```
export MPIR_CVAR_ENABLE_GPU=0
export OMP_NUM_THREADS=1
mpiexec -n 12 --ppn 12 --cpu-bind list:0-7:8-15:16-23:24-31:32-39:40-47:52-59:60-67:68-75:76-83:84-91:92-99 --mem-bind list:0:0:0:0:0:0:1:1:1:1:1:1 --env OMP_NUM_THREADS=1 gpu_tile_compact.sh ../../build_sycl/pwdft cco-cu_surf30.nw

mpiexec -n 12 --ppn 6 --cpu-bind list:0-15:16-31:32-47:52:67:68-83:84-99 --mem-bind list:0:0:0:1:1:1 --env OMP_NUM_THREADS=1 gpu_tile_compact.sh ../../build_sycl/pwdft cco-cu_surf30.nw

```

##  Instructions for OLCF Frontier
<details>
<summary>Toggle for Details</summary>

1. Login to Frontier:
```
ssh frontier.olcf.ornl.gov
```
2. Modules to load:
```
module load amd-mixed
```
3. CMake Build/Install
```
cmake -Bbuild_hip -DNWPW_HIP=ON ./Nwpw -DGPU_TARGETS=gfx90a
cd build_hip
make -j
```
4. Job submission script via `sbatch job_submit.sbatch`
```
#!/bin/bash

#SBATCH -A
#SBATCH -J
#SBATCH -o %x-%j.out
#SBATCH -t 00:15:00
#SBATCH -N 1
#SBATCH -C nvme
#SBATCH --mail-user=
#SBATCH --mail-type=END

module load amd-mixed
module list

export MPICH_GPU_SUPPORT_ENABLED=0
export OMP_NUM_THREADS=1
export CRAYPE_LINK_TYPE=dynamic

date

NNODES=1
NRANKS_PER_NODE=8
NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))

PWDFT_EXE=
PWDFT_INPUT=

srun -N${NNODES} -n${NTOTRANKS} -c1 --ntasks-per-gpu=1 --gpus-per-node=8 --gpu-bind=closest ${PWDFT_EXE} ${PWDFT_INPUT}
```
</details>

##  Instructions for NERSC Perlmutter
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

##  Instructions for ALCF Polaris
<details>
<summary>Toggle for Details</summary>

1. Login to Frontier:
```
ssh polaris.alcf.anl.gov
```
2. Modules to load:
```
module load PrgEnv-gnu cudatoolkit-standalone cmake
```
3. CMake Build/Install
```
cmake -Bbuild_cuda -DNWPW_CUDA=ON ./Nwpw -DCMAKE_CUDA_ARCHITECTURES=80
cd build_cuda
make -j
```
4. Job submission script via `qsub polaris_submit.pbs`. Please note that it needs an additional bash-script to bind the GPU-IDs to MPI-ranks.
The script below is the MPI-ranks with GPUs affinity script:
```
#!/bin/bash
num_gpus=4
# need to assign GPUs in reverse order due to topology
# See Polaris Device Affinity Information https://www.alcf.anl.gov/support/user-guides/polaris/hardware-overview/machine-overview/index.html
gpu=$((${num_gpus} - 1 - ${PMI_LOCAL_RANK} % ${num_gpus}))

unset CUDA_VISIBLE_DEVICES
if [ ${PMI_LOCAL_RANK} -ne 4 ]; then
   export CUDA_VISIBLE_DEVICES=$gpu
fi
#echo "RANK= ${PMI_RANK} LOCAL_RANK= ${PMI_LOCAL_RANK} gpu= ${CUDA_VISIBLE_DEVICES}"
exec "$@"
```
The actual job-script is below:
```
#!/bin/bash

#PBS -N develop
#PBS -l select=8:system=polaris
#PBS -l place=scatter
#PBS -l walltime=01:00:00
#PBS -l filesystems=home:eagle
#PBS -A
#PBS -q workq

module load PrgEnv-gnu cudatoolkit-standalone
module list

export MPICH_GPU_SUPPORT_ENABLED=0
export CRAYPE_LINK_TYPE=dynamic
env
nvidia-smi topo -m

cd ${PBS_O_WORKDIR}

NNODES=`wc -l < $PBS_NODEFILE`
NRANKS_PER_NODE=4
NTHREADS=1

NTOTRANKS=$(( NNODES * NRANKS_PER_NODE ))
echo "NUM_OF_NODES= ${NNODES} TOTAL_NUM_RANKS= ${NTOTRANKS} RANKS_PER_NODE= ${NRANKS_PER_NODE} THREADS_PER_RANK= ${NTHREADS}"

PWDFT_EXE=
PWDFT_INPUT=

mpiexec -n ${NTOTRANKS} --ppn ${NRANKS_PER_NODE} --mem-bind list:0:1:2:3 --cpu-bind list:0-7:8-15:16-23:24-31 --env OMP_NUM_THREADS=${NTHREADS} ./gpu_bind_affinity.sh ${PWDFT_EXE} ${PWDFT_INPUT}
```
</details>

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

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

# Build instructions on JLSE

## `SYCL` backend
### Required Modules
```
export MODULEPATH=$MODULEPATH:/soft/modulefiles:/soft/restricted/CNDA/modules
module load oneapi
module load cmake
```

### Build Instructions (for `SYCL` backend)
```
cd PWDFT
```
```
cmake -H. -Bbuild_sycl -DNWPW_SYCL=On -DCMAKE_CXX_COMPILER=dpcpp ./Nwpw
```
```
make -j4
```

### Running on JSLE
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


# Build Instructions on NERSC Cori-Haswell

## required modules on Cori-Haswell
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

## Build Instructions on Cori-Haswell (starting from PWDFT directory)
```
mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=CC ../Nwpw/
```

## Running on Cori-Haswell
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

# Build instructions on Polaris-CUDA

## Required Modules on Polaris-CUDA
```
module purge
module load cmake cpe-cuda aocl
module load PrgEnv-nvhpc
module unload craype-accel-nvidia80

```

## Build Instructions on Polaris-CUDA (starting from PWDFT directory)
```
mkdir build_cuda
cd build-cuda
cmake -DNWPW_CUDA=ON  -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_Fortran_COMPILER=ftn ../Nwpw/
```
## Running on Polaris-CUDA
```
qsub -q  debug -I -l walltime=01:00:00 -lselect=1 -A myproject -l filesystems=home:eagle:grand
module purge
module load cmake cpe-cuda aocl
module load PrgEnv-nvhpc
module unload craype-accel-nvidia80

mpiexec -n 2 --ppn 2 --cpu-bind=verbose --cpu-bind depth--env CUDA_VISIBLE_DEVICES=2 ./pwdft nwinput.nw

```

# Build instructions on NERSC Cori-CUDA

## Required Modules on Cori-CUDA
```
module unload impi
module load PrgEnv-intel
module load cmake
module load cudatoolkit
```

## Build Instructions on Cori-CUDA (starting from PWDFT directory)
```
mkdir build_cuda
cd build-cuda
cmake -DNWPW_CUDA=ON ../Nwpw/
```

## Running on Cori-CUDA
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

## Running on a MAC
```
mpic++ test.cpp ../build-shared/libpwdft.dylib
setenv DYLD_LIBRARY_PATH /Users/bylaska/Codes/PWDFT/build-shared
a.out
```

## Running on a LINUX
```
mpic++ test.cpp ../build-shared/libpwdft.so
setenv DYLD_LIBRARY_PATH /Users/bylaska/Codes/PWDFT/build-shared
a.out
```

 [2:17 PM] Bagusetty, Abhishek
 add_library(pwdft SHARED nwpw.cpp)

 [2:18 PM] Bagusetty, Abhishek
 CMakeLists.txt (right after this line add_executable(pwdft nwpw.cpp))

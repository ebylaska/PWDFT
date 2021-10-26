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

## Required Modules
```
export MODULEPATH=$MODULEPATH:/soft/modulefiles:/soft/restricted/CNDA/modules
module load oneapi
module load mpi/aurora_mpich
module load cmake
```

## Build Instructions (for `SYCL` backend)
```
cd PWDFT
cmake -H. -Bbuild_sycl -DNWPW_SYCL=On -DCMAKE_CXX_COMPILER=dpcpp ./Nwpw
make -j4
```

## Running on JSLE
```
qsub -I -n 1 -t 60 -q arcticus
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

 [2:17 PM] Bagusetty, Abhishek
 add_library(pwdft SHARED nwpw.cpp)

 [2:18 PM] Bagusetty, Abhishek
 CMakeLists.txt (right after this line add_executable(pwdft nwpw.cpp))

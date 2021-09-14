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
module load cmake/3.18.2
module load PrgEnv-intel
module load craype-haswell
module load openmpi
```

## Build Instructions on Cori-Haswell (starting from PWDFT directory)
```
mkdir build
cd build
cmake ../Nwpw/
```

## Running on Cori-Haswell



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
```
To generate a library clean the build directory and then regenerate cmake with

cmake ../Nwpw -DMAKE_LIBRARY=true
```
 [2:17 PM] Bagusetty, Abhishek
 add_library(pwdft SHARED nwpw.cpp)

 [2:18 PM] Bagusetty, Abhishek
 CMakeLists.txt (right after this line add_executable(pwdft nwpw.cpp))

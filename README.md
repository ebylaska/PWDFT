# PWDFT
PW-DFT development for NWChemEx

# Build instructions on JLSE

## Required Modules
```
module load oneapi
module load mpi/aurora
module load cmake
```

## Build Instructions (for `SYCL` backend)
```
cd PWDFT
cmake -H. -Bbuild_sycl -DNWPW_SYCL=On -DCMAKE_CXX_COMPILER=dpcpp ./Nwpw
make -j4
```
# PWDFT
PW-DFT development for NWChemEx

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
qsub -I -n 1 -t 60 -q iris
qsub -I -n 1 -t 60 -q R.coe_dungeon_iris
```

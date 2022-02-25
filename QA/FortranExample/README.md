```
cd /Users/bylaska/Codes/PWDFT
```

then excecute the following steps:

 1) mkdir build_library
 2) cd build_library
 3) cmake ../Nwpw -DMAKE_LIBRARY=true
 4) On linux: cmake ../Nwpw -DMAKE_LIBRARY=true -DCMAKE_POSITION_INDEPENDENT_CODE=ON
 5) make

Go back to the QA directory, e.g.,

```
cd /Users/bylaska/Codes/PWDFT/QA/FortranExample
```

## Using macOS ##
Set location of the DYLD_LIBRARY_PATH

#### with csh ####
1) setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:/Users/bylaska/Codes/PWDFT/build-shared

#### with bash ####
1) export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:/Users/bylaska/Codes/PWDFT/build-shared

#### Compile fortran_test.f using mpif90 ####
2) mpif90 fortran_test.f /Users/bylaska/Codes/PWDFT/build_library/libpwdft.dylib 


## Compiling on LINUX ##
Add location of shared library location to LD_LIBRARY_PATH, e.g., 

#### with csh ####
1) setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/qfs/people/bylaska/lib

#### with bash ####
1) export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/qfs/people/bylaska/lib

#### Compile test.cpp using mpic++ ####
2) mpif90 fortran_test.f /Users/bylaska/Codes/PWDFT/build_library/libpwdft.so

## How to Run ##
mpirun -np 8 a.out



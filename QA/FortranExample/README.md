```
cd /Users/bylaska/Codes/PWDFT   (mac laptop location)
cd /home/bylaska/Codes/PWDFT    (ubuntu laptop location)
```

then excecute the following steps:

 1) mkdir build_library
 2) cd build_library
 3) cmake ../Nwpw -DMAKE_LIBRARY=true
 4) On linux: cmake ../Nwpw -DMAKE_LIBRARY=true -DCMAKE_POSITION_INDEPENDENT_CODE=ON
 5) make

Next, go back to the QA directory, e.g.,

```
cd /Users/bylaska/Codes/PWDFT/QA/FortranExample   (mac laptop location)
cd /home/bylaska/Codes/PWDFT/QA/FortranExample    (ubuntu laptop location)
```

## Using macOS ##
Set location of the DYLD_LIBRARY_PATH
#### with csh ####
```
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:/Users/bylaska/Codes/PWDFT/build-shared
```
#### with bash ####
```
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:/Users/bylaska/Codes/PWDFT/build-shared
```
then compile fortran_test.f using mpif90
```
mpif90 fortran_test.f /Users/bylaska/Codes/PWDFT/build_library/libpwdft.dylib 
```
## Using LINUX ##
Add location of shared library location to LD_LIBRARY_PATH, e.g., 
#### with csh ####
```
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/qfs/people/bylaska/lib
```
#### with bash ####
```
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/qfs/people/bylaska/lib
```
then compile test.cpp using mpic++
```
mpif90 fortran_test.f /Users/bylaska/Codes/PWDFT/build_library/libpwdft.so

or

mpif90 fortran_test.f /home/bylaska/Codes/PWDFT/build_library/libpwdft.so
```

## Running Example ##
To run the program just type 
```
a.out
```
or run using mpirun, e,g,
```
 mpirun -np 4 a.out
```



```
cd /Users/bylaska/Codes/PWDFT
```

then excecute the following steps:

 1) mkdir build_library
 2) cd build_library
 3) cmake ../Nwpw -DMAKE_LIBRARY=true
 4) On linux: cmake ../Nwpw -DMAKE_LIBRARY=true -DCMAKE_POSITION_INDEPENDENT_CODE=ON
 5) make


## Compiling on macOS ##

Go back to the QA directory, e.g.,

```
cd /Users/bylaska/Codes/PWDFT/QA/LAMMPSExample
```

Set location of the DYLD_LIBRARY_PATH

 1) setenv DYLD_LIBRARY_PATH /${DYLD_LIBRARY_PATH}:/Users/bylaska/Codes/PWDFT/build-shared

Compile test.cpp using mpic++

 2) mpic++ test.cpp /Users/bylaska/Codes/PWDFT/build_library/libpwdft.dylib 


## Compiling on LINUX ##
Add location of shared library location to LD_LIBRARY_PATH, e.g., 

### csh ###
1) setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/qfs/people/bylaska/lib

### bash ###

Compile test.cpp using mpic++

 2) mpic++ test-lammps.cpp /Users/bylaska/Codes/PWDFT/build_library/libpwdft.so


## Running Example ##

To run the program just type 

```
a.out
```

or run using mpirun, e,g,

```
 mpirun -np 4 a.out
```

## Example Output ##

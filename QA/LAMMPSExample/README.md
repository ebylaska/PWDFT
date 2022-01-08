```
cd /Users/bylaska/Codes/PWDFT
```

then excecute the following steps:

 1) mkdir build_library
 2) cd build_library
 3) cmake ../Nwpw -DMAKE_LIBRARY=true
 4) make


## Compiling on macOS ##

Go back to the QA directory, e.g.,

```
cd /Users/bylaska/Codes/PWDFT/QA/LibExample
```


Set location of the DYLD_LIBRARY_PATH

 1) setenv DYLD_LIBRARY_PATH /Users/bylaska/Codes/PWDFT/build-shared

Compile test.cpp using mpic++

 2) mpic++ test.cpp /Users/bylaska/Codes/PWDFT/build_library/libpwdft.dylib 


## Compiling on LINUX ##
Add location of shared library location to LD_LIBRARY_PATH, e.g., 

1) setenv LD_LIBRARY_PATH setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/qfs/people/bylaska/lib

Compile test.cpp using mpic++

 2) mpic++ test-lammps.cpp /Users/bylaska/Codes/PWDFT/build_library/libpwdft.dylib 


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

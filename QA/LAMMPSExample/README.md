```
cd /Users/bylaska/Codes/PWDFT    (mac laptop location)
cd /home/bylaska/Codes/PWDFT     (ubuntu laptop location)
```

then excecute the following steps:

 1) mkdir build_library
 2) cd build_library
 3) cmake ../Nwpw -DMAKE_LIBRARY=true
 4) On linux: cmake ../Nwpw -DMAKE_LIBRARY=true -DCMAKE_POSITION_INDEPENDENT_CODE=ON
 5) make

Next, go back to the QA directory, e.g.,

```
cd /Users/bylaska/Codes/PWDFT/QA/LAMMPSExample   (mac laptop location)
cd /home/bylaska/Codes/PWDFT/QA/LAMMPSExample    (ubuntun laptop location)
```


## Using macOS ##
Set location of the DYLD_LIBRARY_PATH

#### with csh ####
```
setenv DYLD_LIBRARY_PATH ${DYLD_LIBRARY_PATH}:/Users/bylaska/Codes/PWDFT/build-shared   (mac laptop example)
```

#### with bash ####
```
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:/Users/bylaska/Codes/PWDFT/build-shared   (mac laptop example)
```

then compile test.cpp using mpic++
```
mpic++ test.cpp /Users/bylaska/Codes/PWDFT/build_library/libpwdft.dylib 
```

## Using LINUX ##
Add location of shared library location to LD_LIBRARY_PATH, e.g., 

#### with csh ####
```
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/qfs/people/bylaska/lib                  (constance example)
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/home/bylaska/Codes/PWDFT/build_library  (ubuntu laptop example)
```

#### with bash ####
```
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/qfs/people/bylaska/lib                  (constance example)
export LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/home/bylaska/Codes/PWDFT/build_library  (ubuntu laptop example)
```

then compile test.cpp using mpic++
```
mpic++ test-lammps.cpp /Users/bylaska/Codes/PWDFT/build_library/libpwdft.so
```

## Running Example ##

To run the program just type 

```
a.out
a.out w2.nw
```

or run using mpirun, e,g,

```
 mpirun -np 4 a.out
 mpirun -np 4 a.out w2.nw
```


## Example Output ##

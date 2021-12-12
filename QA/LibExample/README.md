## Generate pwdft shared library ##

First go to the PWDFT home directory, e.g.,

cd /Users/bylaska/Codes/PWDFT

then excecute the following steps:

1) mkdir build_library
2) cd build_library
3) cmake ../Nwpw -DMAKE_LIBRARY=true
4) make


## compiling ##
Go back to the QA directory, e.g.,

cd /Users/bylaska/Codes/PWDFT/QA/LibExample


### On a mac ###

Set location of the DYLD_LIBRARY_PATY
1) setenv DYLD_LIBRARY_PATH /Users/bylaska/Codes/PWDFT/build-shared

Compile test.cpp using mpic++
2) mpic++ test.cpp /Users/bylaska/Codes/PWDFT/build_library/libpwdft.dylib 


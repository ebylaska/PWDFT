## Generate pwdft shared library ##

First go to the PWDFT home directory, e.g.,

cd /Users/bylaska/Codes/PWDFT

then excecute the following steps:

1) mkdir build_library
2) cd build_library
3) cmake ../Nwpw -DMAKE_LIBRARY=true
4) make


## compiling ##
setenv DYLD_LIBRARY_PATH /Users/bylaska/Codes/PWDFT/build-shared
mpic++ test.cpp /Users/bylaska/Codes/PWDFT/build_library/libpwdft.dylib 


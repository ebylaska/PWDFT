
cd $PWDFTHOME

mkdir build_library
cd build_library

cmake ../Nwpw -DMAKE_LIBRARY=true
(on linux use: cmake ../Nwpw -DMAKE_LIBRARY=true -DCMAKE_POSITION_INDEPENDENT_CODE=ON)

make


## compiling ##
setenv DYLD_LIBRARY_PATH /Users/bylaska/Codes/PWDFT/build-shared
mpic++ test.cpp /Users/bylaska/Codes/PWDFT/build_library/libpwdft.dylib 


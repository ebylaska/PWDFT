# PWDFT
PW DFT development for NWChemEx

# Build instructions on JLSE
module load dpcpp
module load mkl
module load mpi

cd PWDFT
mkdir build && cd build
cmake -DCMAKE_CXX_COMPILER=dpcpp ../Nwpw
make -j4

# PWDFT
PW DFT development for NWChemEx

# Build instructions on JLSE
module load dpcpp
module load mkl
module load cmake/3.16.3
module load mpi/2019.4.243

cd pwdft-dga/Nwpw
mkdir build
cd build
cmake -DCMAKE_CXX_COMPILER=dpcpp ../
make -j4

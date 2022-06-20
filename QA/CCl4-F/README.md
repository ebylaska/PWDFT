

mpic++ bomd-qmmm3.cpp ../../build_library/libpwdft.so -o bomd_qmmm3.x
mpirun -np 4 ./bomd_qmmm3.x ccl4.nw  | tee bomd_qmmm3.out


mpic++ bomd-qmmm4.cpp ../../build_library/libpwdft.so -o bomd_qmmm4.x
mpirun -np 4 ./bomd_qmmm4.x ccl4-f.nw  | tee bomd_qmmm4.out

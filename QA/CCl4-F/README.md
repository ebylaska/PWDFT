
## CCl4-F QM/MM example
```
mpic++ -std=c++17 bomd-qmmm3.cpp ../../build_library/libpwdft.so -o bomd_qmmm3.x
mpirun -np 4 ./bomd_qmmm3.x ccl4.nw  | tee bomd_qmmm3.out
fgrep @@ bomd_qmmm3.out00 | tee bomd_qmmm3.dat
```

## CCl4-F QM/MM example - read in MM atoms
```
 mpic++ -std=c++17 bomd-qmmm4.cpp ../../build_library/libpwdft.so -o bomd_qmmm4.x
 mpirun -np 4 ./bomd_qmmm4.x ccl4-f.nw  | tee bomd_qmmm4.out
 fgrep @@ bomd_qmmm4.out00 | tee bomd_qmmm4.dat
```

## CCl4-F QM/MM example - read in MM atoms - qm/mm coulomb interactions removed and added back in
```
 mpic++ -std=c++17 bomd-qmmm5.cpp ../../build_library/libpwdft.so -o bomd_qmmm5.x
 mpirun -np 4 ./bomd_qmmm5.x ccl4-f.nw  | tee bomd_qmmm5.out
 fgrep @@ bomd_qmmm5.out00 | tee bomd_qmmm5.dat
```

## CCl4-F QM/MM example - read in MM atoms - qm/mm and qm/qm coulomb interactions removed and added back in incorrectly like LJ pairs
```
 mpic++ -std=c++17 bomd-qmmm6.cpp ../../build_library/libpwdft.so -o bomd_qmmm6.x
 mpirun -np 4 ./bomd_qmmm6.x ccl4-f.nw  | tee bomd_qmmm6.out
 fgrep @@ bomd_qmmm6.out00 | tee bomd_qmmm6.dat
```

## CCl4-F QM/MM example - read in MM atoms - qm/mm and qm/qm coulomb interactions removed and added back in correctly
```
 mpic++ -std=c++17 bomd-qmmm7.cpp ../../build_library/libpwdft.so -o bomd_qmmm7.x
 mpirun -np 4 ./bomd_qmmm7.x ccl4-f.nw  | tee bomd_qmmm7.out
 fgrep @@ bomd_qmmm7.out00 | tee bomd_qmmm7.dat
```


## CCl4-F QM/MM example - read in MM atoms - qm/mm and qm/qm coulomb interactions removed and added back in correctly - io turned off
interface definitions:

- extern int  lammps_pspw_qmmm_minimizer_filename(MPI_Comm, double*, double*, double*, double*, double*, bool, bool, std::string&);
- extern void lammps_pspw_input_filename(MPI_Comm, std::string&, std::string&);

```
 mpic++ -std=c++17 bomd-qmmm8.cpp ../../build_library/libpwdft.so -o bomd_qmmm8.x
 mpirun -np 4 ./bomd_qmmm8.x ccl4-f.nw  | tee bomd_qmmm8.out
 fgrep @@ bomd_qmmm8.out00 | tee bomd_qmmm8.dat
```

## CCl4-F QM/MM example - read in MM atoms - qm/mm and qm/qm coulomb interactions removed and added back in correctly - io redirected
interface definitions:

- extern int  lammps_pspw_qmmm_minimizer_filename(MPI_Comm, double*, double*, double*, double*, double*, bool, bool, std::string&);
- extern void lammps_pspw_input_filename(MPI_Comm, std::string&, std::string&);
- io redirected to qmmm9.nwout

```
 mpic++ -std=c++17 bomd-qmmm9.cpp ../../build_library/libpwdft.so -o bomd_qmmm9.x
 mpirun -np 4 ./bomd_qmmm9.x ccl4-f.nw  | tee bomd_qmmm9.out
```



## Plot energies
```
gnuplot

gnuplot> plot "bomd_qmmm5.dat" using 3:4 w l,"bomd_qmmm5.dat" using 3:5 w l,"bomd_qmmm4.dat" using 3:4 w l,"bomd_qmmm4.dat" using 3:5 w l,"bomd_qmmm3.dat"
using 3:4 w l,"bomd_qmmm3.dat" using 3:5 w l,"bomd_qmmm6.dat" using 3:4,"bomd_qmmm6.dat" using 3:5, "bomd_qmmm7.dat" using 3:4,"bomd_qmmm7.dat" using 3:5

```

<p align="center">
  <img src="./qmmm.svg" width="600" title="hover text">
</p>

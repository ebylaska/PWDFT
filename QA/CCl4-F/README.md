
## CCl4-F QM/MM example
```
mpic++ bomd-qmmm3.cpp ../../build_library/libpwdft.so -o bomd_qmmm3.x
mpirun -np 4 ./bomd_qmmm3.x ccl4.nw  | tee bomd_qmmm3.out
fgrep @@ bomd_qmmm3.out00 | tee bomd_qmmm3.dat
```

## CCl4-F QM/MM example - read in MM atoms
```
 mpic++ bomd-qmmm4.cpp ../../build_library/libpwdft.so -o bomd_qmmm4.x
 mpirun -np 4 ./bomd_qmmm4.x ccl4-f.nw  | tee bomd_qmmm4.out
 fgrep @@ bomd_qmmm4.out00 | tee bomd_qmmm4.dat
```



## Plot energies
```
 gnuplot
 plot "bomd_qmmm4.dat" using 3:4,"bomd_qmmm4.dat" using 3:5,"bomd_qmmm3.dat" using 3:4 w l,"bomd_qmmm3.dat" using 3:5 w l
```

<p align="center">
  <img src=".qmmm.svg" width="350" title="hover text">
</p>

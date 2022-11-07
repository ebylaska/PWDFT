# QM/MM Examples #

## QM/MM Examples - 1 QM water and 1 MM water  ##
The following examples run a QM/MM two water example.

NWChem input deck used: https://github.com/ebylaska/PWDFT/blob/master/QA/QMMM2/water2.nw

### Example 5 - Baseline implementation - Not ideal for LAMMPs ###

Program filename: https://github.com/ebylaska/PWDFT/blob/master/QA/QMMM2/qmmm-example05.cpp

Program library files: https://github.com/ebylaska/PWDFT/blob/master/QA/QMMM2/qmmm.cpp
                       https://github.com/ebylaska/PWDFT/blob/master/QA/QMMM2/qmmm.hpp
                       https://github.com/ebylaska/PWDFT/blob/master/QA/QMMM2/parsestring.hpp

How to compile: mpic++ -O3 qmmm-example05.cpp qmmm.cpp ../../build_library/libpwdft.so 

<img src="./xampl5.svg" width=600>

This example runs a QM/MM two water example in which the call to `c_lammps_pspw_qmmm_minimizer_filename` contains all the QM/QM Coulomb energies and QM/MM Coulomb energies, and the energies and forces between the QM/QM atoms.  User qmmm codes needs to include functions that calculate the electrostatic potentials on the QM atoms from the MM atoms, QM/MM forces, and the MM/MM energies and forces. 

This example calls:
- `extern void c_lammps_pspw_input_filename(MPI_Comm comm_world, const char *cnwfilename, const char *cfilename);`
- `extern int c_lammps_pspw_qmmm_minimizer_filename(MPI_Comm comm_world, double *rion, double *uion, double *fion, double *qion, double *E,
                                                   bool removeqmmmcoulomb, bool removeqmqmcoulomb, const char *cfilename);`
    - where removeqmmmcoulomb = false, and removeqmqmcoulomb = false.
    
- Functions used that are contained in qmmm.cpp:
    - `qmmm.QMMM_electrostatic_potential(qion,rion1,uion);`
    - `qmmm.QMMM_electrostatic_force(qion,rion1,fion)`
    - `ELJ=qmmm.QMMM_LJ_Energy(rion1);`
    - `qmmm.QMMM_LJ_Force(rion1,fion);`
    - `Espring=qmmm.spring_Energy(rion1);`
    - `qmmm.spring_Force(rion1,fion);`


To use this formulation of QM/MM requires the MD code to calculate
- The forces on the QM atoms are included in the call to `c_lammps_pspw_input_filename`
- The electrostatic potential on the QM atoms from the MM atoms
    - $U_I = \sum_{j=MM} {q_{j} \over |R_I - r_j| }$  
    - function: `qmmm.QMMM_electrostatic_potential(qion,rion1,uion)`
- Forces on the QM and MM atoms from the electrostatic (Coulomb) forces between the QM and MM atoms
    - ${\vec F_I} = - Q_{I} \sum_{j=MM} q_{j} {{\vec R_I} - {\vec r_j} \over |R_I - r_j|^3 }$
    - ${\vec F_{j}} = - q_{j} \sum_{I=QM} Q_{I} {{\vec r_j} - {\vec R_I} \over |r_j - R_I|^3}$
    - function: `qmmm.QMMM_electrostatic_force(qion,rion1,fion)`
    - Note the electostatic energy between the QM and MM atoms is already included in the call to `c_lammps_pspw_qmmm_minimizer_filename`
- Forces on MM atoms from the electrostatic (Coulomb) forces between the MM atoms
    - This interaction not needed for the two water example.
   
- The electrostatic (Coulomb) energy and forces between the MM atoms.
    - This interaction is not needed for the two water example. 
  
- The LJ energy and forces between the QM and MM atoms
    - functions: `ELJ=qmmm.QMMM_LJ_Energy(rion1)` and `qmmm.QMMM_LJ_Force(rion1,fion)`

- The LJ energy and forces between the MM atoms
    - This interaction is not needed for the two water example.

- The spring energy and forces between the MM atoms
    - functions: `Espring=qmmm.spring_Energy(rion1)` and `qmmm.spring_Force(rion1,fion)`



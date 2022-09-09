# QM/MM Examples #

## 2 water QM/MM Examples - 1 QM water and 1 MM water ##
The following examples run a QM/MM two water example.


### Example 1 - Baseline implementation - Not ideal for LAMMPs ###

This example runs a QM/MM two water example in which the the call to c_lammps_pspw_qmmm_minimizer_filename contains all the QM/QM energies and QM/MM energies, and the forces on between the QM/QM atoms.  

This example calls 

extern int c_lammps_pspw_qmmm_minimizer_filename(MPI_Comm comm_world, double *rion, double *uion, double *fion, double *qion, double *E,
                                                 bool removeqmmmcoulomb, bool removeqmqmcoulomb, const char *cfilename)

where removeqmmmcoulomb = false, and removeqmqmcoulomb = false.

To use this formulation of QM/MM requires the MD code to calculate
- The electrostatic potential on the QM atoms from the MM atoms
- The electrostatic (Coulomb), LJ and spring energy and forces between the MM atoms
- The electrostatic (Coulomb) forces between the QM and MM atoms
  - Note the electostatic energy between the QM and MM atoms is already included in the call to c_lammps_pspw_qmmm_minimizer_filename

Filename: https://github.com/ebylaska/PWDFT/blob/master/QA/QMMM/qmmm-example01.cpp

Library file used: https://github.com/ebylaska/PWDFT/blob/master/QA/QMMM/qmmm.cpp
                   https://github.com/ebylaska/PWDFT/blob/master/QA/QMMM/qmmm.hpp



### Example 2 - ###


## Formulas for Energies and forces ##

The Coulomb energy between the QM and MM atoms is written as

$$ E_{QMMM} = {\sum_{I=QM} \sum_{j=MM} {Q_{I}  q_{j} \over |R_I - r_j| } }  = {\sum_{I=QM} Q_I \left(\sum_{j=MM} {q_{j} \over |R_I - r_j| }\right) } ={\sum_{I=QM} Q_I U_I} , $$

where $U_I$ is the electrostatic potential located at the QM atom $I$ from the MM atom charges, $q_j$.  Note, the functional dependence of the the two factors in the ${\sum_I Q_I U_I}$ sum are connected to just the QM atom positions, and the MM atom positions positions are embedded in the electrostatic potentials, $U_I$.

, the Blochl charge, $Q_I$, at QM atom is a function of the QM position $R_I$, and the electrostatic potential, $U_I$, is a function of the QM positions, $R_I$, and MM atom positions, $r_j$.  

The forces on the QM atoms are

$$ \frac{\partial E_{QM-mm}}{\partial R} 
= \frac{\partial}{\partial R} \left({\sum_I Q_I U_I}\right)
= {\sum_I \left( \frac{\partial Q_I}{\partial R} U_I + Q_I \frac{\partial U_I}{\partial R} \right)} $$

and the forces and MM atoms are

$$ \frac{\partial E_{QM-mm}}{\partial r} 
$$


The overall Coulomb energy can be written in the following form.

$$ E_{Coulomb} = { {\frac{1}{2}} \sum_{{I,J=QM} \atop {I \ne J}} \frac{Q_I Q_J}{|R_I - R_J|} } + { {\frac{1}{2}} \sum_{{i,j=MM} \atop {i \ne j}}  \frac{q_i q_j}{|r_i - r_j|} } + {\sum_{I=QM} \sum_{j=MM} {Q_{I}  q_{j} \over |R_I - r_j| } } 
$$

$$ E_{Coulomb} = { {\frac{1}{2}} \sum_{{a,b=Atoms} \atop {a \ne b}} \frac{\eta_a \eta_b}{|x_a - x_b|} } 
$$






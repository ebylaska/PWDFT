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

https://jsfiddle.net/8ndx694g/

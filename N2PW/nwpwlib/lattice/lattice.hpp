#ifndef _LATTICE_H_
#define _LATTICE_H_
/* lattice.hpp
   Author - Eric Bylaska

*/

extern void   lattice_init();
extern float lattice_unita(const int, const int);
extern float lattice_unitg(const int, const int);
extern float lattice_ecut();
extern float lattice_wcut();
extern float lattice_omega();
extern float lattice_eggcut();
extern float lattice_wggcut();

#endif

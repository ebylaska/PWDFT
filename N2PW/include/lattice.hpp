#ifndef _LATTICE_H_
#define _LATTICE_H_
/* lattice.hpp
   Author - Eric Bylaska

*/

extern void   lattice_init();
extern double lattice_unita(const int, const int);
extern double lattice_unitg(const int, const int);
extern double lattice_ecut();
extern double lattice_wcut();
extern double lattice_omega();
extern double lattice_eggcut();
extern double lattice_wggcut();

#endif

#ifndef _UTIL_CELL_HPP_
#define _UTIL_CELL_HPP_


namespace pwdft {

extern void util_cell_unita_abc_abg(const double unita_f[9], double lattice[6]);
extern void util_cell_abc_abg_unita(const double lattice[6], double lattice_a_f[9]);
extern void util_cell_lattice_gradient(const double stressin_f[9], const double unita_f[9], double lstress[6]);

}

#endif


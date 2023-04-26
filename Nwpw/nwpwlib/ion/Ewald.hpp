#ifndef _EWALD_HPP_
#define _EWALD_HPP_

#include "Control2.hpp"
#include "Ion.hpp"
#include "Lattice.hpp"
#include "Parallel.hpp"
//#include	"Pseudopotential.hpp"

namespace pwdft {

class Ewald {
  int encut, enx, eny, enz, enshl3d, enpack, enpack_all, enida;
  int *i_indx, *j_indx, *k_indx;
  double *vg, *eG, *vcx, *zv, *ss, *exi, *tmp3, *ftmp;
  double *ewx1, *ewy1, *ewz1;
  double unita[9], unitg[9], ercut, cewald, alpha;
  double eecut;

public:
  Parallel *ewaldparall;
  Ion *ewaldion;
  Lattice *ewaldlattice;

  double *rcell;

  /* Constructors */
  // Ewald(Parallel *, Ion *, Pseudopotential *);
  Ewald(Parallel *, Ion *, Lattice *, Control2 &, double *);

  /* destructor */
  ~Ewald() {
    delete[] i_indx;
    delete[] j_indx;
    delete[] k_indx;
    delete[] vg;
    delete[] vcx;
    delete[] ss;
    delete[] exi;
    delete[] tmp3;
    delete[] ftmp;
    delete[] rcell;
    delete[] eG;
    delete[] zv;
    delete[] ewx1;
    delete[] ewy1;
    delete[] ewz1;
  }

  void phafac();
  int ncut() { return encut; }
  int nida() { return enida; }
  int npack() { return enpack; }
  int npack_all() { return enpack_all; }
  int nshl3d() { return enshl3d; }
  int nx() { return enx; }
  int ny() { return eny; }
  int nz() { return enz; }
  double ecut() { return eecut; }
  double rcut() { return ercut; }
  double mandelung() { return alpha; }
  double energy();
  void force(double *);

  double rs() {
    return pow(3.0 * ewaldlattice->omega() / (16 * atan(1.0)), 1.0 / 3.0);
  }
  double rsalpha() {
    return (alpha * pow((16.0 * atan(1.0) / 3.0), 1.0 / 3.0));
  }
};
} // namespace pwdft

#endif

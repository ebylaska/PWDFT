#ifndef _PSP1D_HAMANN_HPP_
#define _PSP1D_HAMANN_HPP_

#include "PGrid.hpp"
#include "CGrid.hpp"
#include "Parallel.hpp"

namespace pwdft {

class Psp1d_Hamann {

public:
  int version, nrho, nmax, lmax0, lmax, locp, psp_type;
  int nprj, n_extra, n_expansion[10];
  double drho, rlocal, amass, zv;

  double rc[10];
  char atom[2];
  char comment[80];

  double *rho;
  double *vp;
  double *wp;
  double *vnlnrm;
  int *n_prj, *l_prj, *m_prj, *b_prj;

  double *up;
  double *r3_matrix;

  bool semicore;
  int isemicore;
  double rcore;
  double *rho_sc_r;

  /* Constructors */
  Psp1d_Hamann(Parallel *, const char *, const int);

  /* destructor */
  ~Psp1d_Hamann() {
    delete[] rho;
    delete[] vp;
    delete[] wp;
    delete[] vnlnrm;
    if (nprj > 0) {
      delete[] n_prj;
      delete[] l_prj;
      delete[] m_prj;
      delete[] b_prj;
    }
    if (psp_type == 9) {
      delete[] up;
      delete[] r3_matrix;
    }

    if (semicore)
      delete[] rho_sc_r;
  }

  /* G integration routines */
  void vpp_generate_ray(Parallel *, int, double *, double *, double *, double *);
  void vpp_generate_spline(PGrid *, int, double *, double *, double *, double *, double *, double *, double *);

  void cpp_generate_ray(Parallel *, int, double *, double *, double *, double *);
  void cpp_generate_spline(CGrid *, int, double *, double *, double *, double *, double *, double *, double *);
};

} // namespace pwdft

#endif

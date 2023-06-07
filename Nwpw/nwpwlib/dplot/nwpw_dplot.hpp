#ifndef _NWPW_DPLOT_HPP_
#define _NWPW_DPLOT_HPP_

#pragma once

// ********************************************************************
// *                                                                  *
// *       nwpw_dplot: used to generate Gaussian cube files from a    *
// *                   plane-wave density.                            *
// *                                                                  *
// *   The algorithms used in this module are based on the work of    *
// *                                                                  *
// *                                                                  *
// ********************************************************************
#include "blas.h"
#include <iostream>

#include "Control2.hpp"
#include "Ion.hpp"
#include "Pneb.hpp"
#include "Strfac.hpp"
//#include "gdevice.hpp"

namespace pwdft {

class nwpw_dplot {

  Pneb *mypneb;
  Ion *myion;
  std::string permanent_dir_str;

  double dv;
  int n2ft3d, ispin, *ne;

  int num_cubefiles = 0;
  std::vector<std::string> filename;
  std::vector<int> cubetype;

public:
  int ncell[3] = {0, 0, 0};
  double position_tolerance = 0.0;
  double origin[3] = {0.0, 0.0, 0.0};

  /* constructor */
  nwpw_dplot(Ion *, Pneb *, Control2 &);

  /* destructor */
  ~nwpw_dplot() {}

  void gcube_write(std::string, const int, std::string, double *);
  // std::string gcube_write3d(number,cube_comment,rho);
  // std::string gcube_write1d(number,cube_comment,rho);

  void set_position_tolerance(double tol) { position_tolerance = tol; };
  void set_ncell(const int n[]) {
    ncell[0] = n[0];
    ncell[1] = n[1];
    ncell[2] = n[2];
  };

  int number_cubefiles() { return num_cubefiles; }
  int cubetype_cubefile(const int i) { return cubetype[i]; }
  std::string filename_cubefile(const int i) { return filename[i]; }
};

} // namespace pwdft

#endif

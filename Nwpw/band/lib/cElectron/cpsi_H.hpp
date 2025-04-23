#ifndef _CPSI_HPP_
#define _CPSI_HPP_

#pragma once

#include "cKinetic.hpp"
#include "Cneb.hpp"
#include "CPseudopotential.hpp"

namespace pwdft {

//extern void cpsi_H(Cneb *, cKinetic_Operator *, CPseudopotential *, double *,
//                  double *, double *, double *, double *, double *, bool,
//                  double *);
extern void cpsi_H(Cneb *mygrid, cKinetic_Operator *myke, CPseudopotential *mypsp,
                   double *psi, double *psi_r, double *vl, double *vc, double *xcp,
                   double *Hpsi, bool move, double *fion, double *occ = nullptr);



void cpsi_H_orb(const int, Cneb *, cKinetic_Operator *, CPseudopotential *,
                double *, double *, double *, double *);

} // namespace pwdft
#endif

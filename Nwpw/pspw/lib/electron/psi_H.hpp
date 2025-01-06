#ifndef _PSI_HPP_
#define _PSI_HPP_

#pragma once

#include "Kinetic.hpp"
#include "Pneb.hpp"
#include "Pseudopotential.hpp"

namespace pwdft {

//extern void psi_H(Pneb *, Kinetic_Operator *, Pseudopotential *, double *,
//                  double *, double *, double *, double *, double *, bool,
//                  double *);
//
//extern void psi_H(Pneb *, Kinetic_Operator *, Pseudopotential *, double *,
//                  double *, double *, double *, double *, double *, bool,
//                  double *, bool, double *);

extern void psi_H(Pneb *mygrid, Kinetic_Operator *myke, Pseudopotential *mypsp,
                  double *psi, double *psi_r, double *vl, double *vc, double *xcp,
                  double *Hpsi, bool move, double *fion, double *occ = nullptr);

extern void psi_Hv4(Pneb *mygrid, Kinetic_Operator *myke, Pseudopotential *mypsp,
                    double *psi, double *psi_r, double *vsr_l, double *vlr_l,
                    double *vc, double *xcp, double *Hpsi, bool move, double *fion,
                    double *occ = nullptr);

extern void psi_H_orb(Pneb *, Kinetic_Operator *, Pseudopotential *, double *, double *, double *, double *);

} // namespace pwdft
#endif

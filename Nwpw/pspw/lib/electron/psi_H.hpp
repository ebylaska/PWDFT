#ifndef _PSI_HPP_
#define _PSI_HPP_

#pragma once

#include "Kinetic.hpp"
#include "Pneb.hpp"
#include "Pseudopotential.hpp"

namespace pwdft {

extern void psi_H(Pneb *, Kinetic_Operator *, Pseudopotential *, double *,
                  double *, double *, double *, double *, double *, bool,
                  double *);

extern void psi_Hv4(Pneb *, Kinetic_Operator *, Pseudopotential *, double *,
                    double *, double *, double *, double *, double *, double *,
                    bool, double *);

} // namespace pwdft
#endif

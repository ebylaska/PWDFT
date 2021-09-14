#ifndef _PSI_H_
#define _PSI_H_

#include        "Pneb.hpp"
#include        "Kinetic.hpp"
#include        "Pseudopotential.hpp"

namespace pwdft {
using namespace pwdft;

extern void psi_H(Pneb *,
                  Kinetic_Operator *,
                  Pseudopotential *,
                  double *, double *, 
                  double *, double *, double *,
                  double *, 
                  bool, double *);
}
#endif

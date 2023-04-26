#ifndef _INNERLOOP_MD_HPP_
#define _INNERLOOP_MD_HPP_

#include "Control2.hpp"
#include "Coulomb12.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "Kinetic.hpp"
#include "PGrid.hpp"
#include "Pneb.hpp"
#include "Pseudopotential.hpp"
#include "Strfac.hpp"
#include "exchange_correlation.hpp"
#include "nwpw_Nose_Hoover.hpp"

namespace pwdft {
extern void inner_loop_md(const bool, double *, Control2 &, Pneb *, Ion *,
                          nwpw_Nose_Hoover *, Kinetic_Operator *,
                          Coulomb12_Operator *, XC_Operator *,
                          Pseudopotential *, Strfac *, Ewald *, double *,
                          double *, double *, double *, double *, double *,
                          double *, double *, const int, double *);
}
#endif

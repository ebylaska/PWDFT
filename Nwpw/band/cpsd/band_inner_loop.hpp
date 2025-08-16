#ifndef _BANDINNERLOOP_HPP_
#define _BANDINNERLOOP_HPP_

#include "Control2.hpp"
#include "cCoulomb.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "cKinetic.hpp"
#include "CGrid.hpp"
#include "Cneb.hpp"
#include "CPseudopotential.hpp"
#include "CStrfac.hpp"
#include "cExchange_Correlation.hpp"

namespace pwdft {

extern void band_inner_loop(Control2 &, Cneb *, Ion *, cKinetic_Operator *,
                            cCoulomb_Operator *, cXC_Operator *, CPseudopotential *,
                            CStrfac *, Ewald *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *,
                            double *, double *, double *);
}
#endif

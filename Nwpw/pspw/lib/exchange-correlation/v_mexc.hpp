#ifndef _VMEXC_HPP_
#define _VMEXC_HPP_

#include "Pneb.hpp"
#include "vdw_DF.hpp"

namespace pwdft {

extern void v_mexc(const int, Pneb *, vdw_DF *, const double *, double *,
                   const double, double *, double *, double *, double *,
                   double *, double *, double *, double *, double *);

}
#endif

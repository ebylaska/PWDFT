#ifndef _VBWEXC_HPP_
#define _VBWEXC_HPP_

#include "Pneb.hpp"
#include "vdw_DF.hpp"

namespace pwdft {

extern void v_bwexc(const int, Pneb *, vdw_DF *, const double *, const double,
                    const double, double *, double *, double *, double *,
                    double *, double *, double *, double *, double *);
}
#endif

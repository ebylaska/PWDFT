#ifndef _VCWEXC_HPP_
#define _VCWEXC_HPP_

#include "Cneb.hpp"
#include "cvdw_DF.hpp"

namespace pwdft {

extern void v_cwexc(const int, Cneb *, cvdw_DF *, const double *, const double,
                    const double, double *, double *, double *, double *,
                    double *, double *, double *, double *, double *);

void v_cwexc_euv(const int, Cneb *, cvdw_DF *,  const double *,
                 const double, const double, double *);

}

#endif

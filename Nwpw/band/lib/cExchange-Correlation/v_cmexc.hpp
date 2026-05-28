#ifndef _VCMEXC_HPP_
#define _VCMEXC_HPP_

#include "Cneb.hpp"
#include "cvdw_DF.hpp"

namespace pwdft {

extern void v_cmexc(const int, Cneb *, cvdw_DF *, const double *, double *,
                    const double, const double, double *, double *, double *, 
                    double *, double *, double *, double *, 
                    double *, double *, double *);

}
#endif


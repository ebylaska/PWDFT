#ifndef _VBWEXC_HPP_
#define _VBWEXC_HPP_

#include        "Pneb.hpp"

namespace pwdft {
using namespace pwdft;

extern void v_bwexc(const int, Pneb *, const double *, const double, const double,
                    double *, double *, 
                    double *, double *, double *, double *,
                    double *, double *, double *);
}
#endif

#ifndef _BEEF_HPP_
#define _BEEF_HPP_

namespace pwdft {
using namespace pwdft;

extern void gen_BEEF_BW_unrestricted(const int, 
                                      double *, double *,
                                      const double, const double, const double,
                                      double *, double *, double *);

extern void gen_BEEF_BW_restricted(const int,
                                    double *, double *,
                                    const double, const double, const double,
                                    double *, double *, double *);

}
#endif

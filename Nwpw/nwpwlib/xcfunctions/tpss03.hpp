#ifndef _TPSS03_HPP_
#define _TPSS03_HPP_

namespace pwdft {

extern void gen_TPSS03_unrestricted(const int, double *, double *, double *,
                                  const double, const double, 
                                  double *, double *, double *, double *);

extern void gen_TPSS03_restricted(const int, double *, double *, double *, 
                                  const double, const double, 
                                  double *, double *, double *, double *);

} // namespace pwdft
#endif

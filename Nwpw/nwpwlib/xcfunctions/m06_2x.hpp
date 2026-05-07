#ifndef _M06_2X_HPP_
#define _M06_2X_HPP_

namespace pwdft {

extern void gen_M06_2x_unrestricted(const int, double *, double *, double *,
                                  const double, const double, 
                                  double *, double *, double *, double *);

extern void gen_M06_2x_restricted(const int, double *, double *, double *, 
                               const double, const double, 
                               double *, double *, double *, double *);

} // namespace pwdft
#endif

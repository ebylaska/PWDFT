#ifndef _M06L_HPP_
#define _M06L_HPP_

namespace pwdft {

extern void gen_M06L_unrestricted(const int, double *, double *, double *,
                                  const double, const double, 
                                  double *, double *, double *, double *);

extern void gen_M06L_restricted(const int, double *, double *, double *, 
                                const double, const double, 
                                double *, double *, double *, double *);

} // namespace pwdft
#endif

#ifndef _PKZB_HPP_
#define _PKZB_HPP_

namespace pwdft {

extern void gen_PKZB_unrestricted(const int, double *, double *, double *,
                                  const double, const double, 
                                  double *, double *, double *, double *);

extern void gen_PKZB_restricted(const int, double *, double *, double *, 
                                const double, const double, 
                                double *, double *, double *, double *);

} // namespace pwdft
#endif

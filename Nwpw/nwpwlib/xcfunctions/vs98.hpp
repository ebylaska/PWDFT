#ifndef _VS98_HPP_
#define _VS98_HPP_

namespace pwdft {

extern void gen_VS98_unrestricted(const int, double *, double *, double *,
                                  const double, const double, 
                                  double *, double *, double *, double *);

extern void gen_VS98_restricted(const int, double *, double *, double *, 
                                const double, const double, 
                                double *, double *, double *, double *);

} // namespace pwdft
#endif

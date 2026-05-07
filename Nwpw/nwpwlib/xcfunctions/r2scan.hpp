#ifndef _R2SCAN_HPP_
#define _R2SCAN_HPP_

namespace pwdft {

extern void gen_r2SCAN_unrestricted(const int, double *, double *, double *,
                                  const double, const double, 
                                  double *, double *, double *, double *);

extern void gen_r2SCAN_restricted(const int, double *, double *, double *, 
                                const double, const double, 
                                double *, double *, double *, double *);

} // namespace pwdft
#endif

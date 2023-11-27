#ifndef _B3LYP_HPP_
#define _B3LYP_HPP_

namespace pwdft {

extern void gen_B3LYP_BW_unrestricted(const int, double *, double *,
                                      const double, const double, double *,
                                      double *, double *);

extern void gen_B3LYP_BW_restricted(const int, double *, double *, const double,
                                    const double, double *, double *, double *);
} // namespace pwdft

#endif

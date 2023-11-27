#ifndef _HSEPBE_HPP_
#define _HSEPBE_HPP_

namespace pwdft {

extern void gen_HSE_BW_unrestricted(const int, double *, double *, const double,
                                    const double, double *, double *, double *);

extern void gen_HSE_BW_restricted(const int, double *, double *, const double,
                                  const double, double *, double *, double *);
} // namespace pwdft

#endif

#ifndef _revPBE_HPP_
#define _revPBE_HPP_

namespace pwdft {

extern void gen_revPBE_BW_unrestricted(const int, double *, double *,
                                       const double, const double, double *,
                                       double *, double *);

extern void gen_revPBE_BW_restricted(const int, double *, double *,
                                     const double, const double, double *,
                                     double *, double *);
} // namespace pwdft

#endif

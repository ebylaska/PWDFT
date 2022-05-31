#ifndef _BLYP_HPP_
#define _BLYP_HPP_

namespace pwdft {


extern void gen_BLYP_BW_unrestricted(const int, 
                                     double *, double *,
                                     const  double, const double,
                                     double *, double *, double *);

extern void gen_BLYP_BW_restricted(const int,
                                   double *, double *,
                                   const double, const double,
                                   double *, double *, double *);
}

#endif

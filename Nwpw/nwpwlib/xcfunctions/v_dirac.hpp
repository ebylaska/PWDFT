// v_dirac.hpp
#ifndef _VDIRAC_HPP_
#define _VDIRAC_HPP_

namespace pwdft {

void v_dirac(const int ispin,
             const int n2ft3d,
             double *dn,
             double *xcp,
             double *xce,
             double *x);

} // namespace pwdft

#endif


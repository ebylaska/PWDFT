#ifndef _ION_ION_HPP_
#define _ION_ION_HPP_

namespace pwdft {
using namespace pwdft;

extern double ion_ion_e(const int, const double *, const double *);
extern void ion_ion_f(const int, const double *, const double *, double *);
extern void ion_ion_m_f(const int, const double *, const double *, double *);

} // namespace pwdft

#endif

#ifndef _SYMMETRY_ELEMENTS_HPP_
#define _SYMMETRY_ELEMENTS_HPP_

#include <string>

namespace pwdft {
using namespace pwdft;

extern void determine_point_group(const double *, const double *, const int, const double,
                                  std::string&, int&, std::string&, double *, double *, double *, double *);

} // namespace pwdft

#endif


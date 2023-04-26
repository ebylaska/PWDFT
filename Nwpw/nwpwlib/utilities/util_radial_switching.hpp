#ifndef _util_RADIAL_SWITCHING_HPP_
#define _util_RADIAL_SWITCHING_HPP_

#pragma once

namespace pwdft {
using namespace pwdft;

extern double util_radial_switching(const double, const double, const double);
extern void util_radial_dswitching(const double, const double, const double,
                                   double *, double *);

} // namespace pwdft

#endif

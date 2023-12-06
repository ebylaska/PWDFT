#ifndef _CPSI_HPP_
#define _CPSI_HPP_

#pragma once

#include "cKinetic.hpp"
#include "Cneb.hpp"
#include "CPseudopotential.hpp"

namespace pwdft {

extern void cpsi_H(Cneb *, cKinetic_Operator *, CPseudopotential *, double *,
                  double *, double *, double *, double *, double *, bool,
                  double *);


} // namespace pwdft
#endif

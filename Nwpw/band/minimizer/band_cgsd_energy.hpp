#ifndef _BAND_CGSD_ENERGY_HPP_
#define _BAND_CGSD_ENERGY_HPP_

#pragma once

namespace pwdft {

#include "Solid.hpp"

extern double band_cgsd_noit_energy(Solid &, bool, std::ostream &);
extern double band_cgsd_energy(Control2 &, Solid &, bool, std::ostream &);
extern void band_cgsd_energy_gradient(Solid &, double *);

} // namespace pwdft
#endif

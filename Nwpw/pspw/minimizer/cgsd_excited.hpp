#ifndef _CGSD_EXCITED_HPP_
#define _CGSD_EXCITED_HPP_

#pragma once

namespace pwdft {

#include "Molecule.hpp"

extern void cgsd_excited(Control2 &, Molecule &, bool, std::ostream &);

} // namespace pwdft
#endif

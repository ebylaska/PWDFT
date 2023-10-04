#ifndef _FILON_FILTER_HPP_
#define _FILON_FILTER_HPP_

#pragma once

#include <cmath>
#include <cstring> //memset
#include <iostream>

#include "Pneb.hpp"

namespace pwdft {

extern void coulomb_filter(Pneb *, double *, const std::string);
}

#endif

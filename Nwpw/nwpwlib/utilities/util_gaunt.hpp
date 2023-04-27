#ifndef _util_GAUNT_HPP_
#define _util_GAUNT_HPP_

#pragma once

namespace pwdft {

extern void util_gaunt_init(const bool, const int);
extern void util_gaunt_end();

extern double util_gaunt(const bool, const int, const int, const int, const int,
                         const int, const int);
extern double util_gaunt2(const bool, const int, const int, const int,
                          const int, const int, const int);
extern double util_gaunt3(const bool, const int, const int, const int,
                          const int, const int, const int);

} // namespace pwdft

#endif

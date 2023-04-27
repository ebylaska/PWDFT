/* psp_library.h -
   author - Eric Bylaska

*/

#pragma once

// extern "C" {
//#include        "compressed_io.h"
// }
#include "Control2.hpp"
#include "compressed_io.hpp"
#include <map>
#include <string>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace pwdft {

class psp_library {

public:
  std::string nwpw_permanent_dir;
  std::string nwpw_libraryps_dir;
  std::string default_library;
  std::map<std::string, std::string> libraries;

  /* Constructors */
  psp_library(const std::string);
  psp_library(Control2 &);

  /* other */
  void psp_check(const char *, Control2 &, double *, std::ostream &);

  std::string psp_libname(const char *);
  int psp_type(const char *);
  int psp_lmax(const char *);
  int psp_locp(const char *);
  double psp_rlocal(const char *);

  void psp_generator_auto(const char *, Control2 &, std::ostream &);

  // void psp_generate(const char *);

  // void geninput(const char *, char *);

  // void genpspfile(const char *);
};

} // namespace pwdft

#ifndef _PSP_FILE_CHECK_HPP_
#define _PSP_FILE_CHECK_HPP_

#include "Control2.hpp"
#include "Ion.hpp"
#include "Parallel.hpp"
#include <iostream>

namespace pwdft {

extern void psp_file_check(Parallel *, Ion *, Control2 &, std::ostream &);

}

#endif

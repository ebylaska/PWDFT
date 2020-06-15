#ifndef _PSP_FORMATTER_H_
#define _PSP_FORMATTER_H_

#include        "Parallel.hpp"
#include        "Control2.hpp"
#include        "Ion.hpp"

extern void   psp_formatter_check(Parallel *, Ion *, Control2&);
extern double psp_formatter_auto(Parallel *, Control2&, char *a);

#endif

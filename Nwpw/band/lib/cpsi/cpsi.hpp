#ifndef _CPSIGETHEADER_HPP_
#define _CPSIGETHEADER_HPP_

#pragma once

#include "Parallel.hpp"
#include "Cneb.hpp"

namespace pwdft {

extern void cpsi_get_header(Parallel *, int *, int *, double *, int *, int *, int *, char *);

extern void cpsi_read0(Cneb *, int *, int *, double *, int *, int *, int *, double *, char *);
extern bool cpsi_read(Cneb *, char *, bool, double *, std::ostream &);

extern void cpsi_write(Cneb *, int *, int *, double *, int *, int *, int *, double *, char *, std::ostream &);
extern bool cpsi_filefind(Cneb *, char *);

// extern void v_psi_read(Pneb *, int *, int *, double *, int *, int *,double
// *); extern void v_psi_write(Pneb *, int *, int *, double *, int *, int
// *,double *);

} // namespace pwdft
#endif

#ifndef _PSIGETHEADER_HPP_
#define _PSIGETHEADER_HPP_

#pragma once

#include "Parallel.hpp"
#include "Pneb.hpp"

namespace pwdft {

extern void psi_get_header(Parallel *, int *, int *, double *, int *, int *, char *);

extern void psi_get_header_occupation(Parallel *, int *, int *, double *, int *, int *, int *, char *);

extern void psi_read0(Pneb *, int *, int *, double *, int *, int *, double *, int *, double *,char *, bool);
extern void psi_read0(Pneb *, int *, int *, double *, int *, int *, double *, char *);
extern bool psi_read(Pneb *, char *, bool, double *, int *, double *, std::ostream &);
extern bool psi_read(Pneb *, char *, bool, double *, std::ostream &);

extern void psi_write(Pneb *, int *, int *, double *, int *, int *, double *, int *, double *, char *, std::ostream &);
extern void psi_write(Pneb *, int *, int *, double *, int *, int *, double *, char *, std::ostream &);
extern bool psi_filefind(Pneb *, char *);

extern bool epsi_read(Pneb *, char *, bool, const int *, double *, std::ostream &);
extern void epsi_write(Pneb *, int *, int *, double *, int *, int *, double *, char *, std::ostream &);

// extern void v_psi_read(Pneb *, int *, int *, double *, int *, int *,double
// *); extern void v_psi_write(Pneb *, int *, int *, double *, int *, int
// *,double *);

} // namespace pwdft
#endif

#ifndef _PSIGETHEADER_H_
#define _PSIGETHEADER_H_

#include	"Parallel.hpp"
#include	"Pneb.hpp"

extern void psi_get_header(Parallel *, int *, int *, float *, int *, int *, char *);

extern void psi_read0(Pneb *, int *, int *, float *, int *, int *,float *, char *);
extern void psi_read(Pneb *, char *, float *);

extern void psi_write(Pneb *, int *, int *, float *, int *, int *,float *, char *);
extern bool psi_filefind(Pneb *, char *);

//extern void v_psi_read(Pneb *, int *, int *, float *, int *, int *,float *);
//extern void v_psi_write(Pneb *, int *, int *, float *, int *, int *,float *);

#endif

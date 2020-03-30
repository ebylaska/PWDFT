#ifndef _INNERLOOP_MD_H_
#define _INNERLOOP_MD_H_

#include	"Control2.hpp"
#include        "Pneb.hpp"
#include        "PGrid.hpp"
#include        "Ion.hpp"
#include        "Ewald.hpp"
#include	"Kinetic.hpp"
#include	"Coulomb.hpp"
#include	"Strfac.hpp"
#include        "Pseudopotential.hpp"

extern void inner_loop_md(bool, Control2&, Pneb *, Ion *, 
                       Kinetic_Operator *, Coulomb_Operator *,
                       Pseudopotential *, Strfac *, Ewald *,
                       double *, double *, double *, double *,
                       double *, double *, double *,
                       double *, double *, double *, double *);
#endif
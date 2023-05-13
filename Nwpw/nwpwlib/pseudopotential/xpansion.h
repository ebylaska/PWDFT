/*
 $Id$
*/
#ifndef _GET_CS_H_
#define _GET_CS_H_
/* get_cs.h
 */

extern void p_xpansion(double []);
extern void dp_xpansion(double []);
extern void ddp_xpansion(double []);
extern void dddp_xpansion(double []);
extern void ddddp_xpansion(double []);
extern void psp_xpansion(void);
extern void psi_xpansion(void);
extern void dpsi_xpansion(void);
extern void init_xpansion(int, int, double, double, double *, double *,
                          double *, double *, double *, double *);

#endif

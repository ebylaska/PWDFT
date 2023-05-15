/*
 $Id$
*/
#ifndef _SPLINE_H_
#define _SPLINE_H_
/* spline.h -
 */

extern void init_Linear(char *);
extern void end_Linear(void);
extern int nrl_Linear(void);
extern double drl_Linear(void);
extern void Log_to_Linear(double [], double [], double []);
extern void Log_to_Linear_zero(double [], double [], double []);
extern void normalize_Linear(double *);
extern void normalize_Linear2(double *, double *);

#endif

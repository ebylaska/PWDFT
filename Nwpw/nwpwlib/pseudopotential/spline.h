/*
 $Id$
*/
#ifndef _SPLINE_H_
#define _SPLINE_H_
/* spline.h -
*/

extern void	init_Linear(char*);
extern void	end_Linear();
extern	int	nrl_Linear();
extern	float	drl_Linear();
extern void	Log_to_Linear();
extern void	Log_to_Linear_zero();
extern void	normalize_Linear(float*);
extern void	normalize_Linear2(float*,float*);

#endif


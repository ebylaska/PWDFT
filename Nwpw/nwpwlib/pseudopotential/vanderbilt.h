#ifndef _VANDERBILT_H_
#define _VANDERBILT_H_
/* vanderbilt.h -
   Author - Eric Bylaska

*/

void Suggested_Param_Vanderbilt(int *, int [], int [], double [], double [],
                                double [], double *, double *);
void solve_Vanderbilt(int, int [], int [], double [], double [], double [],
                      double, double, int ns[10], int indx_il[4][10], int indx_ijl[4][4][10],
                      double **, double **, double **, double *, double *,
                      double **, double *, double *, double *, double *,
                      double *, double *, double *, double *, double *, double *);

#endif
/* $Id$ */

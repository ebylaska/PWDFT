#ifndef _DIRAC_H_
#define _DIRAC_H_
/* dirac.c */
extern void R_Dirac (int n, int l, int s2, float Z, const float *v,
                         int *mch, float *Eig, float *u, float *uprime);
extern void R_Dirac_Fixed_E (int n, int l, int s2, float Z, const float *v,
                                 int match, float E, float *u, float *uprime);
extern void R_Dirac_FixedLogDeriv(int n, int l, int s2, float Z, const float *v, 
	int match, float u_logderiv, float *Eig, float *u, float *uprime);

#endif
/* $Id$ */

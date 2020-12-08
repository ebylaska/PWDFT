#ifndef _ZORA_H_
#define _ZORA_H_
/* dirac.c */
extern void R_ZORA (int n, int l, int s2, float Z, const float *v,
                         int *mch, float *Eig, float *u, float *uprime);
extern void R_ZORA_Fixed_E (int n, int l, int s2, float Z, const float *v,
                                 int match, float E, float *u, float *uprime);
#endif
/* $Id$ */

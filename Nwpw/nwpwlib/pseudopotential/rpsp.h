#ifndef _RPSP_H_
#define _RPSP_H_
#include <stdio.h>

/* rpsp.c */
extern void init_RelPsp(char *filename);
extern void solve_RelPsp(void);
extern char *solver_Name_RelPsp(void);
extern void print_RelPsp(FILE *fp);
extern void set_Solver_RelPsp(int solver);
extern int Vanderbilt_RelPsp(void);
extern int NormConserving_RelPsp(void);
extern float E_RelPsp(void);
extern float eigenvalue_RelPsp(int i);
extern float *rho_RelPsp(void);
extern float *rho_semicore_RelPsp(void);
extern float *drho_semicore_RelPsp(void);
extern float r_semicore_RelPsp(void);
extern float *Beta_RelPsp(int i, int l);
extern float *r_psi_il_RelPsp(int i, int l);
extern float *r_hard_psi_il_RelPsp(int i, int l);
extern int ns_RelPsp(int l);
extern float D0_RelPsp(int i, int j, int l);
extern float q_RelPsp(int i, int j, int l);
extern float *Vlocal_RelPsp(void);
extern float *V_RelPsp(int i);
extern float *r_psi_RelPsp(int i);
extern int n_RelPsp(int i);
extern int l_RelPsp(int i);
extern int lmax_RelPsp(void);
extern float fill_RelPsp(int i);
extern int Nvalence_RelPsp(void);
extern float peak_RelPsp(int i);
extern float rcut_RelPsp(int i);
extern float rcut_il_RelPsp(int i, int l);
extern float Zion_RelPsp(void);
extern int state_RelPsp(int nt, int lt, int st);
extern char *comment_RelPsp(void);
#endif
/* $Id$ */

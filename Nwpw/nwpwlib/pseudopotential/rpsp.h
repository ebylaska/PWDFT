#ifndef _RPSP_H_
#define _RPSP_H_
#include <stdio.h>

/* rpsp.c */
extern void init_RelPsp(char *);
extern void solve_RelPsp(void);
extern char *solver_Name_RelPsp(void);
extern void print_RelPsp(FILE *);
extern void set_Solver_RelPsp(int);
extern int Vanderbilt_RelPsp(void);
extern int NormConserving_RelPsp(void);
extern double E_RelPsp(void);
extern double eigenvalue_RelPsp(int);
extern double *rho_RelPsp(void);
extern double *rho_semicore_RelPsp(void);
extern double *drho_semicore_RelPsp(void);
extern double r_semicore_RelPsp(void);
extern double *Beta_RelPsp(int, int);
extern double *r_psi_il_RelPsp(int, int);
extern double *r_hard_psi_il_RelPsp(int, int);
extern int ns_RelPsp(int);
extern double D0_RelPsp(int, int, int);
extern double q_RelPsp(int, int, int);
extern double *Vlocal_RelPsp(void);
extern double *V_RelPsp(int);
extern double *r_psi_RelPsp(int);
extern int n_RelPsp(int);
extern int l_RelPsp(int);
extern int lmax_RelPsp(void);
extern double fill_RelPsp(int);
extern int Nvalence_RelPsp(void);
extern double peak_RelPsp(int);
extern double rcut_RelPsp(int);
extern double rcut_il_RelPsp(int, int);
extern double Zion_RelPsp(void);
extern int state_RelPsp(int, int, int);
extern char *comment_RelPsp(void);
#endif
/* $Id$ */

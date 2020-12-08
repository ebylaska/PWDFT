#ifndef _RHAMANN_H_
#define _RHAMANN_H_
/* rhamann.c */
extern void Suggested_Param_RelHamann (int *num_states_psp,
				       int n_psp[], int l_psp[], int s_psp[],
				       float e_psp[], float fill_psp[],
				       float rcut_psp[]);

extern void solve_RelHamann (int num_psp, int *n_psp, int *l_psp,
			     int *s_psp, float *e_psp, float *fill_psp,
			     float *rcut_psp, float **r_psi_psp,
			     float **r_psi_prime_psp, float *rho_psp,
			     float *rho_semicore, float **V_psp,
			     float *eall_psp, float *eh_psp, float *ph_psp,
			     float *ex_psp, float *px_psp, float *ec_psp,
			     float *pc_psp);

#endif
/* $Id$ */

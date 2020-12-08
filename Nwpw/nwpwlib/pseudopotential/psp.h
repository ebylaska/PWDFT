#ifndef	_PSP_H_
#define _PSP_H_
/* psp.h -
   author - Eric Bylaska

*/


/* Solver type: Solve_Type */
#define	Hamann		-8401
#define	Troullier	-8402
#define Vanderbilt      -8403

extern void	init_Psp(char*);
extern void	end_Psp();
extern void	solve_Psp();
extern void	print_Psp();
extern float	E_Psp();
extern float	eigenvalue_Psp(int);
extern float	*V_Psp();
extern float	*rho_Psp();
extern float	*rho_semicore_Psp();
extern float	*drho_semicore_Psp();
extern float	r_semicore_Psp();
extern float	*r_psi_Psp(int );
extern int	n_Psp(int);
extern int	ns_Psp(int);
extern int	l_Psp(int);
extern int	lmax_Psp();
extern float	fill_Psp(int);
extern int	Nvalence_Psp();
extern float	peak_Psp(int);
extern float	rcut_Psp(int);
extern float	Zion_Psp();
extern int	state_Psp(int, int);
extern int      Vanderbilt_Psp();
extern float   *Vlocal_Psp();
extern float   *Beta_Psp();
extern float   *r_hard_psi_il_Psp();
extern float   *r_psi_il_Psp();
extern float   D0_Psp();
extern float   q_Psp();
extern float   rcut_il_Psp(int,int);
extern int      Vanderbilt_Psp();
extern int      NormConserving_Psp();
extern char     *comment_Psp();

extern int      kb_extra_Psp();
extern int      kb_expansion_Psp(int);
extern float   *r_psi_extra_Psp(int );


/* used for setting solver parameters */
extern void	set_Solver_Psp();

#endif
/* $Id$ */

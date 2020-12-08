#ifndef _PAW_LOG_GRID_H_
#define _PAW_LOG_GRID_H_
/*
   $Id$
*/

#include        <stdio.h>


extern void     paw_init_LogGrid_from_file( float Z, FILE *fp);
extern float   paw_r0_LogGrid();
extern float   *paw_alloc_LogGrid();
extern void	    paw_dealloc_LogGrid(float *grid);
extern float   *paw_r_LogGrid();
extern float   *paw_r2_LogGrid();
extern float   *paw_r3_LogGrid();
extern float   *paw_scratch_LogGrid();
extern int      paw_N_LogGrid();
extern float   paw_log_amesh_LogGrid();
extern float   paw_amesh_LogGrid();
extern float	  paw_Integrate_LogGrid();
extern float      paw_Def_Integr(float fpow,float *f,float rpow,int Nrange);
extern void	    paw_Zero_LogGrid();
extern void	    paw_Copy_LogGrid();
extern void     paw_Copy_spin_LogGrid(float **gridnew, float **gridold);
extern float paw_Norm_LogGrid(int M, float gamma, float *u);
extern void    paw_Derivative_LogGrid(float *f,float *df);
extern float*	paw_Indef_Integr();
extern float   paw_dot_product(float*, float*);
extern float   paw_dot_product1(int ,float*, float*);
extern int      paw_get_grid_index(float r);
extern int      paw_get_grid_index(float r);
extern void     paw_print_loggrid_information();
extern float paw_dot_product1(int n, float *f, float *g);

/* this one was missing before */
extern void  paw_end_LogGrid();

#endif



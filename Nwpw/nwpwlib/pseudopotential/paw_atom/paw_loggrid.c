/*
   $Id$
*/

//#include        <stdio.h>
#include        <stdlib.h>
#include        <string.h>
#include        <math.h>

#include        "paw_pred_cor.h"
#include        "paw_loggrid.h"
#include        "paw_utilities.h"
#include        "paw_my_memory.h"
#include        "paw_sdir.h"

/***********************/
/* LogGrid data structure */
/***********************/

/* Hamman definitions */
/*static float amesh = 1.0247*/
/*static float r0Z   =  0.00625 */
/*static float Lmax  = 45.0  */
/*static float Lmax  = 45.0; */


int     Ngrid;
float  log_amesh;
float  amesh;
float  *rgrid;
float  *rgrid2;
float  *rgrid3;
float  *scratch;


float Lmax;
static float r0Z   = 0.00025;
float  r0;

void  paw_end_LogGrid()
{
    paw_dealloc_LogGrid(rgrid);
    paw_dealloc_LogGrid(rgrid2);
    paw_dealloc_LogGrid(rgrid3);
    paw_dealloc_LogGrid(scratch);
}


/****************************************
Function name	  : paw_init_LogGrid
Description	    :
Return type		  : void
Argument        : float Z -> ion charge
Argument        : FILE *fp
Author     		  : Marat Valiev
Date & Time		  : 1/7/99 4:26:57 PM
****************************************/
void  paw_init_LogGrid_from_file( float Z, FILE *fp)
{
    int  i;
    char input[30];

    strcpy(input,"<grid>");
    if (paw_find_word(input,fp) != 0)
    {
        if (paw_debug()) printf("Using default parameters\n");
        Lmax=25.0;
        amesh=1.005;
        log_amesh = log(amesh);

        r0 = r0Z/Z;
        Ngrid = (int) floor(log(Lmax/r0)/log_amesh)+1;

        /* make sure Ngrid is odd */
        if ((Ngrid%2)==0) Ngrid += 1;

    }
    else
    {
        fscanf(fp,"%e",&Lmax);
        fscanf(fp,"%d", &Ngrid);
        fscanf(fp,"%e",&r0);

        log_amesh = log(Lmax/r0)/(Ngrid-1);
        amesh     = exp(log_amesh);
    }

    Lmax = r0*pow(amesh,Ngrid-1);
    rgrid   = paw_alloc_LogGrid();
    rgrid2  = paw_alloc_LogGrid();
    rgrid3  = paw_alloc_LogGrid();
    scratch = paw_alloc_LogGrid();

    /* define rgrid */
    rgrid[0] = r0;
    for (i=1; i <= Ngrid-1; ++i)
        rgrid[i] = amesh*rgrid[i-1];

    for (i=0; i <= Ngrid-1; ++i)
    {
        rgrid2[i] = rgrid[i]*rgrid[i];
        rgrid3[i] = rgrid[i]*rgrid[i]*rgrid[i];
    }

}



/****************************************
  Function name	  : paw_alloc_LogGrid
  Description	    : creates a loggrid array
  Return type		  : float*
  Author     		  : Eric Bylaska & Marat Valiev
  Date & Time		  : 1/7/99 4:31:15 PM
****************************************/
float* paw_alloc_LogGrid()
{
    float *tt;

    tt = paw_alloc_1d_array(Ngrid);

    return tt;

} /* paw_alloc_LogGrid */




/****************************************
  Function name	  : paw_dealloc_LogGrid
  Description	    : deallocates the LogGrid array
  Return type		  : void
  Argument         : grid
  Author     		  : Eric Bylaska & Marat Valiev
  Date & Time		  : 1/7/99 4:31:39 PM
****************************************/
void  paw_dealloc_LogGrid(float *grid)
{

    paw_dealloc_1d_array(grid);

} /* dealloc_LogGrid */





/****************************************
Function name	  : paw_r_LogGrid
Description	    : returns the pointer to the rgrid array
Return type		  : float*
Author     		  : Eric Bylaska & Marat Valiev
Date & Time		  : 1/7/99 4:35:15 PM
****************************************/
float* paw_r_LogGrid()
{
    return rgrid;

} /* r_LogGrid */


float* paw_r2_LogGrid()
{
    return rgrid2;

}

float* paw_r3_LogGrid()
{
    return rgrid3;

}
/****************************************
  Function name	  : paw_N_LogGrid
  Description	    :
  Return type		  : int
  Author     		  : Eric Bylaska & Marat Valiev
  Date & Time		  : 1/7/99 4:35:40 PM
****************************************/
int  paw_N_LogGrid()
{
    return Ngrid;

} /* paw_N_LogGrid */




/****************************************
Function name	  : paw_r0_LogGrid
Description	    : returns the first nonzero coordinate
of a log grid
Return type		  : float 
Author     		  : Eric Bylaska & Marat Valiev
Date & Time		  : 1/7/99 4:36:10 PM
****************************************/
float paw_r0_LogGrid()
{
    return r0;
}



/****************************************
Function name	  : paw_log_amesh_LogGrid
Description	    : returns the value log(amesh)
Return type		  : float
Author     		  : Eric Bylaska & Marat Valiev
Date & Time		  : 1/7/99 4:37:11 PM
****************************************/
float paw_log_amesh_LogGrid()
{
    return log_amesh;

} /* paw_log_amesh_LogGrid */




/****************************************
  Function name	  : paw_amesh_LogGrid
  Description	    : returns the value (amesh)
  Return type		  : float
  Author     		  : Eric Bylaska & Marat Valiev
  Date & Time		  : 1/7/99 4:37:43 PM
****************************************/
float paw_amesh_LogGrid()
{
    return(amesh);

} /* paw_amesh_LogGrid */





/****************************************
  Function name	  : paw_Def_Integr
  Description	    : calculates definite integral of
  the function f with the weight
  r**(rpow) from 0 to Nrange.
  Function f is assumed to behave
  as r**(fpow) near 0
  Return type		  : float
  Argument         : float fpow
  Argument         : float *f
  Argument         : float rpow
  Argument         : int Nrange
  Author     		  : Marat Valiev
  Date & Time		  : 1/7/99 4:38:43 PM
****************************************/
float  paw_Def_Integr(float fpow,float *f,float rpow,int Nrange)

{
    int i;
    float sum;

    sum = (   9.0*f[0]*pow(rgrid[0],rpow+1)
              + 23.0*f[1]*pow(rgrid[1],rpow+1)
              + 28.0*f[2]*pow(rgrid[2],rpow+1)
          )/28.0;

    for (i=3; i<Nrange; ++i)
    {
        sum += f[i]*pow(rgrid[i],rpow+1);
    }

    sum = log_amesh*sum + f[0]*pow(rgrid[0],rpow+1)/(rpow+fpow+1);

    return sum;

}


/****************************************
Function name	  : paw_Integrate_LogGrid
Description	    :  returns a definite integral of
of the given function f times
r squared
Return type		  : float
Argument         : float f[]
Author     		  : Eric Bylaska & Marat Valiev
Date & Time		  : 1/7/99 4:40:59 PM
****************************************/

float  paw_Integrate_LogGrid(float *f)
{
    int    i;
    float sum;

    sum = (   9.0*f[0]*(rgrid[0]*rgrid[0]*rgrid[0])
              + 23.0*f[1]*(rgrid[1]*rgrid[1]*rgrid[1])
              + 28.0*f[2]*(rgrid[2]*rgrid[2]*rgrid[2])
          )/28.0;

    for (i=3; i<=Ngrid-1; i++)
    {
        sum += f[i]*(rgrid[i]*rgrid[i]*rgrid[i]);
    }

    sum = log_amesh*sum + f[0]*(rgrid[0]*rgrid[0]*rgrid[0])/3.0;

    return sum;

}




/****************************************
Function name	  :   paw_Zero_LogGrid
Description	    :
Return type		  : void
Argument         : float* grid
Author     		  : Eric Bylaska & Marat Valiev
Date & Time		  : 1/7/99 4:42:43 PM
****************************************/
void    paw_Zero_LogGrid(float *grid)

{
    int i;

    for (i=0; i<Ngrid; ++i)
        grid[i] = 0.0;

} /* paw_Zero_LogGrid */


/****************************************
  Function name	  : paw_Copy_LogGrid
  Description	    :
  Return type		  : void
  Argument         : float *gridnew
  Argument         : float *gridold
  Author     		  : Eric Bylaska & Marat Valiev
  Date & Time		  : 1/7/99 4:44:22 PM
****************************************/
void paw_Copy_LogGrid(float *gridnew, float *gridold)

{
    int i;
    for (i=0; i<Ngrid; ++i)
        gridnew[i] = gridold[i];
}



/****************************************
Function name	  : paw_Norm_LogGrid
Description	    : This routine calculates the Norm
of a wavefunction assuming that
the wavefunction decays like an
exponential as r goes to  infinity
Return type		  : float
Argument         : int M -> endpoint
Argument         : float gamma -> power of u near the 0
Argument         : float *u
Author     		  : Eric Bylaska & Marat Valiev
Date & Time		  : 1/7/99 4:46:08 PM
****************************************/
float paw_Norm_LogGrid(int M, float gamma, float *u)
{
    int   i;
    float sum;


    sum = (   9.0*u[0]*u[0]*rgrid[0]
              + 23.0*u[1]*u[1]*rgrid[1]
              + 28.0*u[2]*u[2]*rgrid[2]
          )/28.0;

    for (i=3; i<=M; ++i)
    {
        sum += u[i]*u[i]*rgrid[i];
    }

    sum = log_amesh*sum + u[0]*u[0]*rgrid[0]/(2.0*gamma+1.0);

    return sum;


} /* paw_Norm_LogGrid */



/****************************************
  Function name	  :   paw_Derivative_LogGrid
  Description	    : calculates the derivative
  of the function defined on the
  loggrid array
  Return type		  : void
  Argument         : float *f  -> original function
  Argument         : float *df -> derivative
  Author     		  : Eric Bylaska & Marat Valiev
  Date & Time		  : 1/7/99 4:47:29 PM
****************************************/
void    paw_Derivative_LogGrid(float *f,float *df)
{
    int i;

    df[0] = paw_Derivative5_1(0,f)/(log_amesh*rgrid[0]);
    df[1] = paw_Derivative5_2(1,f)/(log_amesh*rgrid[1]);

    for (i=2; i<Ngrid-2; ++i)
        df[i] = paw_Derivative5_3(i,f)/(log_amesh*rgrid[i]);


    df[Ngrid-2] = paw_Derivative5_4(Ngrid-2,f)/(log_amesh*rgrid[Ngrid-2]);
    df[Ngrid-1] = paw_Derivative5_5(Ngrid-1,f)/(log_amesh*rgrid[Ngrid-1]);

} /* paw_Derivative_LogGrid */


/****************************************
  Function name	  : paw_dot_product
  Description	    :
  Return type		  : float
  Argument         : float *f
  Argument         : float *g
  Author     		  : Marat Valiev
  Date & Time		  : 2/7/99 4:53:24 PM
****************************************/
float paw_dot_product(float *f, float *g)
{

    int k;
    float norm;

    norm =0.0;

    /* Integrate from 0 to r0 */
    norm = 0.5 * f[0] * g[0] * rgrid[0];

    for (k = 0; k < Ngrid; ++k)
        norm += f[k] * g[k] * rgrid[k] * log_amesh;

    return norm;

}

float paw_dot_product1(int n, float *f, float *g)
{

    int k;
    float norm;

    norm =0.0;

    /* Integrate from 0 to r0 */
    norm = 0.5 * f[0] * g[0] * rgrid[0];

    for (k = 0; k < n; ++k)
        norm += f[k] * g[k] * rgrid[k] * log_amesh;

    return norm;

}
/****************************************
 Function name	  : paw_print_loggrid_information
 Description	    :
 Return type		  : void
 Argument         : FILE *fp
 Author     		  : Marat Valiev
 Date & Time		  : 3/31/99 3:00:49 PM
****************************************/
void  paw_print_loggrid_information(FILE *fp)
{
    fprintf(fp,"\n");
    fprintf(fp," Logarithmic grid information ( r(i)=r0*pow(a,i) ):\n");
    fprintf(fp,"\n");

    fprintf(fp,"   a    = %le\n", paw_amesh_LogGrid());
    fprintf(fp,"   N    = %d\n", paw_N_LogGrid());
    fprintf(fp,"   r0   = %le\n",paw_r0_LogGrid());
    fprintf(fp,"   rmax = %le\n",paw_r_LogGrid()[paw_N_LogGrid()-1]);
    fprintf(fp,"\n");

}


int paw_get_grid_index(float r)
{
    int i;

    if (r > Lmax)
    {

        printf("grid point is out of range\n");
        exit(1);

    }
    else
    {
        i = (int) floor( log(r/rgrid[0])/log_amesh );
    }

    return i;
}

float* paw_scratch_LogGrid()
{

    return scratch;

}



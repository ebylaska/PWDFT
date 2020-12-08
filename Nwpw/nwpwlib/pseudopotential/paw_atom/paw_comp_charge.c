/*
   $Id$
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "paw_loggrid.h"
#include "paw_my_constants.h"
#include "paw_bisect.h"
#include "paw_comp_charge.h"
#include "paw_error_function.h"

static float comp_charge_tolerance = 1.0e-10;
static float r_comp;
static float sigma_comp;
static float* comp_charge_density;
static float* V_comp;


/****************************************
 Function name	  : paw_init_comp_charge
 Description	    :
 Return type		  : void
 Argument         :  float a_Z
 Argument         : float a_r_sphere
 Author     		  : Marat Valiev
 Date & Time		  : 4/9/99 3:26:49 PM
****************************************/
void paw_init_comp_charge( float comp_radius)
{

    r_comp = comp_radius;

    sigma_comp = paw_find_sigma_comp();

    comp_charge_density = paw_alloc_LogGrid();

    V_comp  = paw_alloc_LogGrid();

}


/****************************************
 Function name	  : paw_boundary_function
 Description	    :
 Return type		  : float
 Argument         : float sigma
 Author     		  : Marat Valiev
 Date & Time		  : 4/9/99 3:28:36 PM
****************************************/
float paw_boundary_function(float sigma)
{

    float tmp;
    float Z;

    Z = 1.0;

    tmp = Z*4/(sqrt(PI)*pow(sigma,3.0))*exp(-pow(r_comp/sigma,2.0));

    tmp = tmp - comp_charge_tolerance;

    return tmp;

}


/****************************************
 Function name	  : paw_find_sigma_comp
 Description	    :
 Return type		  : float
 Author     		  : Marat Valiev
 Date & Time		  : 4/9/99 3:28:28 PM
****************************************/
float paw_find_sigma_comp()
{
    float tmp;
    float rc1;
    float rc2;

    rc1 = r_comp/10;
    rc2 = 4*r_comp;

    tmp = paw_bisection(paw_boundary_function, rc1, rc2,0.0001*comp_charge_tolerance);

    return tmp;

}



/****************************************
 Function name	  : paw_find_comp_charge_potential
 Description	    :
 Return type		  : void
 Argument         : float charge
 Argument         : float ps_charge
 Argument         : float *V_comp
 Author     		  : Marat Valiev
 Date & Time		  : 4/9/99 3:28:30 PM
****************************************/
float* paw_find_comp_charge_potential(float Zion, float charge, float ps_charge)
{
    int k;
    int Ngrid;
    float comp_charge;
    float *rgrid;

    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    comp_charge = - Zion + charge - ps_charge;

    for (k = 0; k <= Ngrid-1; k++)
        comp_charge_density[k] = comp_charge*4/(sqrt(PI)*pow(sigma_comp,3.0))*
                                 exp(-pow(rgrid[k]/sigma_comp,2.0));

    for (k = 0; k <= Ngrid-1; k++)
        V_comp[k]    = comp_charge*paw_my_erf(rgrid[k]/sigma_comp)/rgrid[k];

    return V_comp;

}

float paw_get_sigma_comp()
{

    return sigma_comp;

}

float paw_get_comp_charge_radius()
{

    return r_comp;

}

void  paw_print_comp_charge_information(FILE *fp)
{
    fprintf(fp,"\n");

    fprintf(fp," Compensation charge information :\n");
    fprintf(fp,"\n");

    fprintf(fp,"   compensation charge radius    = %le\n",
            paw_get_comp_charge_radius() );

    fprintf(fp,"   gaussian width (sigma)        = %le\n",
            paw_get_sigma_comp());

    fprintf(fp,"   compensation charge tolerance = %le\n",
            comp_charge_tolerance);

    fprintf(fp,"\n\n");

}

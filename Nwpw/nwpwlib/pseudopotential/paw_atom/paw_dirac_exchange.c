/*
   $Id$
*/

#include	<stdio.h>
#include  <math.h>
#include  "paw_loggrid.h"
#include	"paw_dirac_exchange.h"
#include  "paw_my_constants.h"

static	float onethird;

static	float	alpha=0.6666666666666666667;
static	float *ex_functional;
static	float *Vx;

/****************************************
 Function name	  : paw_init_dirac_exchange()
 Description	    :
****************************************/
void paw_init_dirac_exchange()
{

    Vx		        = paw_alloc_LogGrid();
    ex_functional = paw_alloc_LogGrid();

    /* define constants */
    onethird  = 1.0/3.0;

}

/****************************************
 Function name	  : paw_generate_exchange_pot_LDA(float *rho)
 Description	    :
****************************************/
void paw_generate_exchange_pot_LDA(float *rho)

{
    int	i;
    float n;
    float n_onethird;
    float ux_p;

    /* loggrid variables */
    int	   Ngrid;


    /* access the loggrid variables */
    Ngrid     = paw_N_LogGrid();


    for (i=0; i<=Ngrid-1; i++)
    {

        n = rho[i]/(4.0*PI);
        n_onethird = pow((3.0*n/PI),onethird);

        ux_p = -(3.0/2.0)*alpha*n_onethird;

        Vx[i] = ux_p ;

    } /*for i*/


}
float* paw_get_exchange_potential()
{

    return Vx;

}

/****************************************
 Function name	  : paw_get_exchange_energy_LDA(float *rho)
 Description	    :
****************************************/
float paw_get_exchange_energy_LDA(float *rho)

{
    int	i;
    float n;
    float n_onethird;
    float ex_p;
    float Ex;
    float *tmp;

    /* loggrid variables */
    int	   Ngrid;

    /* access the loggrid variables */
    Ngrid     = paw_N_LogGrid();
    tmp		    = paw_scratch_LogGrid();


    for (i=0; i<Ngrid; ++i)
    {

        n = rho[i]/(4.0*PI);
        n_onethird = pow((3.0*n/PI),onethird);


        ex_p = -(9.0/8.0)*alpha*n_onethird;

        ex_functional[i] = ex_p;

    } /*for i*/



    for (i=0; i<Ngrid; ++i)
        tmp[i] = rho[i]*ex_functional[i];
    Ex = paw_Integrate_LogGrid(tmp);

    return Ex;

}



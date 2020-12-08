/*
   $Id$
*/

#include  <stdio.h>
#include  <math.h>
#include  "paw_loggrid.h"
#include  "paw_pred_cor.h"
#include  "paw_hartree.h"


static float Zion;
static float Eion;
static float *Vion;

/****************************************
 Function name	  : paw_init_ion(float Z)
 Description	    :
****************************************/
void paw_init_ion(float Z)
{
    int i;
    int Ngrid;
    float *rgrid;


    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    Zion = Z;

    Vion  = paw_alloc_LogGrid();
    for (i=0; i<Ngrid; ++i)
        Vion[i] = -Z/rgrid[i];

}

/****************************************
 Function name	  : paw_get_ion_energy
 Description	    :
 Return type		  : float
 Argument         : float *dn
 Author     		  : Marat Valiev
 Date & Time		  : 3/31/99 2:12:25 PM
****************************************/
float paw_get_ion_energy(float *dn)
{
    int    k;
    int Ngrid;
    float *tmp;


    Ngrid = paw_N_LogGrid();
    tmp   = paw_scratch_LogGrid();

    for (k=0; k<Ngrid; ++k)
    {
        tmp[k] = Vion[k]*dn[k];
    }

    Eion = paw_Integrate_LogGrid(tmp);
    return Eion;
}

float* paw_get_ion_pot()
{
    return Vion;
}

float paw_get_ion_charge()
{

    return Zion;


}


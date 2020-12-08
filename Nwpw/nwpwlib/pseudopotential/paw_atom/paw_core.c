/*
   $Id$
*/

#include        <stdlib.h>
#include        <stdio.h>
#include        <string.h>
#include        <math.h>

#include        "paw_my_constants.h"
#include        "paw_loggrid.h"
#include        "paw_core.h"
#include        "paw_orbitals.h"
#include        "paw_basis.h"
#include        "paw_atom.h"
#include        "paw_sdir.h"



/*internal core data*/
static float core_charge;
static float ps_core_charge;
static float core_kin_energy;
static float *core_density;
static float *ps_core_density;


/****************************************
 Function name	  : init_core
 Description	    :
****************************************/
void  paw_init_core()

{

    core_density    = paw_alloc_LogGrid();
    ps_core_density = paw_alloc_LogGrid();

}

void  paw_set_core()
{

    int k;
    int Ngrid;
    float *full_density;
    float *paw_basis_density;
    FILE *fp;
    float *rgrid;
    char output[200];

    Ngrid = paw_N_LogGrid();

    core_kin_energy = paw_get_atom_kinetic_energy() -
                      paw_get_paw_kinetic_energy();


    full_density      = paw_get_pointer_density();
    paw_basis_density = paw_get_pointer_paw_density();

    for (k=0;k<=Ngrid-1;k++)
        core_density[k] = full_density[k] - paw_basis_density[k];

    core_charge    = paw_Integrate_LogGrid(core_density);

    paw_Zero_LogGrid(ps_core_density);
    ps_core_charge = paw_Integrate_LogGrid(ps_core_density);

    if (paw_debug()) printf("core_charge = %f\n",core_charge);
    if (paw_debug()) printf("ps_core_charge = %f\n",ps_core_charge);
    if (paw_debug()) printf("core kinetic energy = %f\n",core_kin_energy);

    if (paw_debug())
    {
        sprintf(output, "%s%s", paw_sdir(),"core");
        fp = fopen(output,"w");

        rgrid = paw_r_LogGrid();
        for (k=0;k<Ngrid;++k)
            fprintf(fp,"%f  %f\n",rgrid[k],core_density[k]);
        fclose(fp);
    }

}


float* paw_get_pointer_core_density()
{

    return core_density;

}

float* paw_get_pointer_ps_core_density()
{

    return ps_core_density;

}

float paw_get_core_charge()
{

    return core_charge;

}

float paw_get_ps_core_charge()
{

    return ps_core_charge;

}

float paw_get_core_kinetic_energy()
{

    return core_kin_energy;

}


void  paw_print_core_information(FILE *fp)
{
    fprintf(fp,"\n");

    fprintf(fp," Core information :\n");
    fprintf(fp,"\n");

    fprintf(fp,"   core charge         = %le\n",
            core_charge );

    fprintf(fp,"   pseudo core charge  = %le\n",
            ps_core_charge );

    fprintf(fp,"   core kinetic energy = %le\n",
            core_kin_energy);

    fprintf(fp,"\n\n");

}

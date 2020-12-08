/*
   $Id$
*/

/************************************
  REVISION LOG ENTRY
  Revision By: ...
  Revised on 3/30/99 3:52:12 PM
  Comments: created file
************************************/

#include  <stdlib.h>
#include  <stdio.h>
#include  <string.h>
#include  <math.h>

#include "paw_loggrid.h"
#include "paw_potential.h"
#include "paw_vosko.h"
#include "paw_dirac_exchange.h"
#include "paw_hartree.h"
#include "paw_ion.h"
#include "paw_utilities.h"
#include "paw_sdir.h"

float    *Vi;

extern char     *atom_name;

/* global loggrid data */
extern int     Ngrid;
extern float  *rgrid;

extern   int   Hartree_Type ;
extern   int   Exchange_Type;
extern   int   Correlation_Type;


/****************************************
 Function name	  : paw_init_potential
 Description	    :
 Return type		  : void
 Argument         : float Z
 Author     		  : Marat Valiev
 Date & Time		  : 3/30/99 5:05:57 PM
****************************************/
void paw_init_potential()
{

    Vi   = paw_alloc_LogGrid();


}



/****************************************
 Function name	  : paw_get_kohn_sham_potential
 Description	    :
 Return type		  : void
 Argument         : float *dn
 Argument         : float **rho
 Argument         : float **Vo
 Author     		  : Marat Valiev
 Date & Time		  : 3/30/99 5:32:13 PM
****************************************/
void paw_find_kohn_sham_potential(float *rho,float *V_ks)
{
    int k;
    float    *Vion;
    float    *Vh;
    float    *Vx;
    float    *Vc;

    Vion = paw_get_ion_pot();

    paw_generate_hartree_pot(rho);
    Vh = paw_get_hartree_pot();

    paw_generate_exchange_pot_LDA(rho);
    Vx = paw_get_exchange_potential();

    paw_generate_corr_pot_LDA(rho);
    Vc = paw_get_corr_pot_LDA();

    for (k=0; k <= Ngrid-1; k++)
        V_ks[k] = Vion[k] + Vh[k] + Vx[k] + Vc[k];
}



/****************************************
 Function name	  : paw_set_Kohn_Sham_potential
 Description	    :
 Return type		  : void
 Argument         : float *rho
 Author     		  : Marat Valiev
 Date & Time		  : 4/10/99 6:59:53 PM
****************************************/
void paw_set_kohn_sham_potential(float *rho)
{
    paw_find_kohn_sham_potential(rho,Vi);
}



/****************************************
 Function name	  : *paw_get_Kohn_Sham_potential
 Description	    :
 Return type		  : float
 Author     		  : Marat Valiev
 Date & Time		  : 4/10/99 7:00:50 PM
****************************************/
float *paw_get_kohn_sham_potential()
{

    return Vi;

}


/****************************************
 Function name	  :   paw_Thomas_Fermi
 Description	    :
 Return type		  : void
 Argument         : float Z
 Argument         : float *V
 Author     		  : Marat Valiev
 Date & Time		  : 4/10/99 3:41:34 PM
****************************************/
void    paw_Thomas_Fermi(float Z, float *V)
{
    int    i;
    float  x,t;

    for (i=0; i<Ngrid; ++i)
    {
        x = rgrid[i]* pow( (Z/0.69395656), (1.0/3.0));
        t = Z/(
                1.0+sqrt(x)*(0.02747 - x*(0.1486-0.007298*x))
                +           x*(1.243 + x*(0.2302+0.006944*x))
            );

        if (t < 1.0)
            t=1.0;

        V[i] = -t/rgrid[i];

    }

}


/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/
/*****************************************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "paw_basis.h"
#include "paw_orbitals.h"
#include "paw_loggrid.h"
#include "paw_utilities.h"
#include "paw_my_constants.h"
#include "paw_bisect.h"
#include "paw_potential.h"
#include "paw_potential.h"
#include "paw_hartree.h"
#include "paw_comp_charge.h"
#include "paw_dirac_exchange.h"
#include "paw_vosko.h"
#include "paw_ion.h"
#include "paw_core.h"


static int nbasis;
static int* i_r_potential;
static float* rc_pot;
static float rc_ref;
static float lambda=6;
static float r_function_for_rc_pot;

static float *c;
static float pot_tolerance = 1.0e-10;
static float* r_potential;
static float r_ref;
static float **fcut;
static float **V_paw;
static float *V_ref;
static float *V_pseudo;

/****************************************
 Function name    : paw_function_for_rc_pot
 Description        :
 Return type              : float
 Argument         : float r
 Author                   : Marat Valiev
 Date & Time              : 4/9/99 12:50:27 PM
****************************************/
float paw_function_for_rc_pot( float r)
{

    float tmp;

    tmp = exp(-pow((r_function_for_rc_pot/r), lambda)) - pot_tolerance;

    return tmp;

}

float paw_find_rc_pot(float rcut_in)
{

    float tmp;
    float rc1;
    float rc2;

    r_function_for_rc_pot = rcut_in;

    rc1 = rcut_in/100;
    rc2 = 2*rcut_in;
    tmp = paw_bisection(paw_function_for_rc_pot, rc1, rc2, 0.00001*pot_tolerance);

    return tmp;

}


/****************************************
 Function name    : paw_init_paw_potential
 Description        :
 Return type              : void
 Argument         : int a_nbasis
 Argument         : float a_r_sphere
 Author                   : Marat Valiev
 Date & Time              : 4/11/99 4:11:34 PM
****************************************/
void paw_init_paw_potential(int a_nbasis,
                            float c0,
                            float a_r_ref,
                            float* a_r_potential,
                            float* V_ks)
{

    int   i;
    int   k;
    int    Ngrid;
    int ic;
    float rc;
    float a,d,b;
    float V_prime;
    float *rgrid;

    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    nbasis      = a_nbasis;

    r_ref  = a_r_ref;
    rc_ref = paw_find_rc_pot(r_ref);

    i_r_potential = (int *) malloc(nbasis * sizeof(int));
    r_potential   = (float *) malloc(nbasis * sizeof(float));
    rc_pot        = (float *) malloc(nbasis * sizeof(float));
    c             = (float *) malloc(nbasis * sizeof(float));
    fcut          = (float **) malloc(nbasis * sizeof(float *));
    V_pseudo      = paw_alloc_LogGrid();
    V_ref         = paw_alloc_LogGrid();
    V_paw         = (float **) malloc(nbasis * sizeof(float *));

    for (i = 0; i <= nbasis-1; ++i)
    {
        V_paw[i]  = paw_alloc_LogGrid();
        fcut[i]  = paw_alloc_LogGrid();
    }

    for (i=0; i <= nbasis-1;i++)
    {
        i_r_potential[i] = paw_get_grid_index(a_r_potential[i]);
        /*r_potential[i]   = rgrid[i_r_potential[i]];*/
        r_potential[i] = a_r_potential[i];
        rc_pot[i] = paw_find_rc_pot(r_potential[i]);
    }


    /* set ref potential
    for (k = 0; k <= Ngrid-1; ++k)
      V_ref[k] = V_ks[k]*(1 - exp(-pow((rgrid[k]/rc_ref),lambda)))+
                 c0*exp(-pow((rgrid[k]/rc_ref),lambda));

    */

    ic = paw_get_grid_index(a_r_ref);
    rc = rgrid[ic];

    V_prime = 0.5*(V_ks[ic+1] - V_ks[ic-1])/(rc*paw_log_amesh_LogGrid());


    a = c0;

    d = (0.5*V_prime*rc + a - V_ks[ic])/(rc*rc*rc*rc);

    b = (0.5*V_prime*rc - 2.0*d*rc*rc*rc*rc)/(rc*rc);

    for (k=0;k<ic;++k)
        V_ref[k] = a + b*rgrid[k]*rgrid[k] + d*rgrid[k]*rgrid[k]*rgrid[k]*rgrid[k];

    for (k=ic;k<Ngrid;++k)
        V_ref[k] = V_ks[k];

    /* Form cutoff function  for PS potential*/
    for (i = 0; i <= nbasis-1; i++)
        for (k = 0; k <= Ngrid-1; k++)
            fcut[i][k] = exp(-pow((rgrid[k]/rc_pot[i]),lambda));

    /*Initialize coefficients for the construction of PS potential*/
    for (i=0; i<nbasis; ++i)
        c[i] = V_ks[i_r_potential[i]];
    /*
      for (i = 0; i < nbasis; ++i)
        c[i] = V_ks[paw_get_grid_index(r_potential[i])];
    */


    /* set initial paw potential*/
    for (i = 0; i < nbasis; ++i)
        for (k = 0; k <= Ngrid-1; k++)
            V_paw[i][k] = V_ref[k] + c[i] * fcut[i][k];

}


/****************************************
 Function name    : paw_get_paw_potential
 Description        :
 Return type              : float*
 Argument         : int i
 Author                   : Marat Valiev
 Date & Time              : 1/11/99 11:11:16 AM
****************************************/
float* paw_get_paw_potential(int i)
{

    return V_paw[i];

}


float* paw_get_pointer_pseudopotential()
{

    return V_pseudo;

}

/****************************************
 Function name    : paw_update_paw_potential
 Description        :
 Return type              : void
 Argument         : int conv_status
 Argument         : int i
 Argument         : float *w
 Author                   : Marat Valiev
 Date & Time              : 1/10/99 6:27:44 PM
****************************************/
void paw_update_paw_potential(int *conv_status, int i, float eig, float eig_ps,float *w)
{

    int     k;
    int Ngrid;
    float  sv;
    float  dcl;
    float  eps;
    float  *tmp;

    eps = 1.0e-7;

    Ngrid = paw_N_LogGrid();
    tmp   = paw_scratch_LogGrid();

    for (k = 0; k <= Ngrid-1; k++)
    {
        tmp[k]= fcut[i][k] * w[k] * w[k];
    }

    sv = paw_Def_Integr(0.0,tmp,0.0,Ngrid-1);

    dcl = (eig - eig_ps) / sv;

    /*In case of the overflow rescale dcl*/

    if (fabs(dcl)>100)
        dcl = 10*fabs(dcl)/dcl;


    c[i] = c[i] + dcl;
    *conv_status = (fabs(dcl) <= eps);


    /*Form new paw potential if necessary*/
    if (!(*conv_status))
        for (k = 0; k <= Ngrid-1; k++)
            V_paw[i][k] = V_ref[k] + c[i] * fcut[i][k];

}

/****************************************
 Function name    : paw_get_ref_pot
 Description        :
 Return type              : float*
 Author                   : Marat Valiev
 Date & Time              : 1/25/99 11:29:24 AM
 Modifications    : 1/26/99
****************************************/
float* paw_get_ref_pot()
{
    return V_ref;
}


/****************************************
 Function name    : paw_generate_pseudopot
 Description        :
 Return type              : void
 Author                   : Marat Valiev
 Date & Time              : 4/10/99 7:40:00 PM
****************************************/
void paw_generate_pseudopot()
{

    int   k;
    int   Ngrid;
    float charge;
    float ps_charge;
    float Z;
    float *Vh;
    float *Vx;
    float *Vc;
    float *rho;
    float *rho_ps;
    float *rho_core;
    float *rho_core_ps;
    float *full_density;
    float *full_ps_density;
    float* V_comp;
    float *rgrid;
    FILE *fp;
    char data_filename[300];

    if ( !(paw_projectors_are_done()) )
    {
        printf("error, pseudopotential cannot be generated ");
        printf(" because projectors have not been yet \n");
        exit(1);
    }

    Ngrid = paw_N_LogGrid();
    rgrid = paw_r_LogGrid();

    Vh      = paw_alloc_LogGrid();

    Vx      = paw_alloc_LogGrid();
    Vc      = paw_alloc_LogGrid();

    full_density = paw_alloc_LogGrid();
    full_ps_density = paw_alloc_LogGrid();

    Z = paw_get_ion_charge();

    paw_set_core();

    /*get densities*/
    rho         = paw_get_pointer_paw_density();
    rho_ps      = paw_get_pointer_paw_ps_density();
    rho_core    = paw_get_pointer_core_density();
    rho_core_ps = paw_get_pointer_ps_core_density();

    paw_Zero_LogGrid(full_density);
    paw_Zero_LogGrid(full_ps_density);

    for (k=0;k<=Ngrid-1;k++)
    {
        full_density[k]    = rho[k] + rho_core[k];
        full_ps_density[k] = rho_ps[k] + rho_core_ps[k];

    }

    charge   = paw_Integrate_LogGrid(full_density);
    ps_charge      = paw_Integrate_LogGrid(full_ps_density);


    V_comp = paw_find_comp_charge_potential(Z,charge,ps_charge);

    paw_generate_hartree_pot(full_ps_density);
    Vh = paw_get_hartree_pot();

    paw_generate_exchange_pot_LDA(full_ps_density);
    Vx = paw_get_exchange_potential();

    paw_generate_corr_pot_LDA(full_ps_density);
    Vc = paw_get_corr_pot_LDA();


    /*form pseudopotential*/
    for (k=0;k<=Ngrid-1;k++)
    {

        V_pseudo[k] = V_ref[k] - Vh[k]- V_comp[k]- Vc[k] - Vx[k];



    }

    if (paw_debug())
    {
        sprintf(data_filename,"%sdensity",paw_sdir());
        fp = fopen(data_filename,"w");

        for (k=0;k<Ngrid;++k)
            fprintf(fp,"%f  %f   %f\n",rgrid[k],rho[k] - rho_ps[k],V_pseudo[k]);
        fclose(fp);
    }

    paw_dealloc_LogGrid(full_density);
    paw_dealloc_LogGrid(full_ps_density);

}

float paw_get_potential_matching_radius()
{
    return r_ref;
}

void  paw_print_paw_potential_information(FILE *fp)
{
    int i;
    int *n;
    int *l;

    fprintf(fp,"\n");

    fprintf(fp," Paw potential information :\n");
    fprintf(fp,"\n");

    fprintf(fp,"   reference potential matching radius    = %le\n",
            paw_get_potential_matching_radius());

    fprintf(fp,"   reference potential rcut parameter     = %le\n",
            rc_ref);

    fprintf(fp,"   lambda parameter                       = %le\n",
            lambda);

    fprintf(fp,"   potential tolerance                    = %le\n",
            pot_tolerance);

    fprintf(fp,"\n");

    n = paw_get_pointer_paw_n_array();
    //l = paw_get_pointer_l_array();
    l = paw_get_pointer_paw_l_array();

    fprintf(fp,"   nl     c[i]        rcut[i] \n");


    for (i=0;i<=nbasis-1;i++)
    {

        fprintf(fp,"   %d%s    %f    %f \n",n[i],paw_spd_Name(l[i]),c[i],
                rc_pot[i]);

    }

    fprintf(fp,"\n\n");

}

void paw_print_paw_potential_to_file(char* atom_name)
{
    int i;
    int j;
    int k;
    int Ngrid;
    int *prin_n;
    int * orb_l;
    float *rgrid;
    char data_filename[300];
    char script_filename[300];
    char nl_name[20];
    FILE *fp;

    if (paw_debug())
    {
        Ngrid = paw_N_LogGrid();
        rgrid = paw_r_LogGrid();

        prin_n = paw_get_pointer_paw_n_array();
        orb_l  = paw_get_pointer_paw_l_array();

        sprintf(data_filename,"%s%s_pot.dat",paw_sdir(),atom_name);
        fp = fopen(data_filename,"w+");

        for (k=0; k<=Ngrid-1; k++)
        {
            fprintf(fp,"%le\t%le\t%le", rgrid[k], V_ref[k], V_pseudo[k]);

            for (i=0; i<=nbasis-1; i++)
                fprintf(fp,"\t%le ",V_paw[i][k]);



            fprintf(fp,"\n");

        }

        fclose(fp);

        sprintf(script_filename,"%s%s_ref_pot.plt",paw_sdir(),atom_name);

        fp = fopen(script_filename,"w+");

        fprintf(fp,"set style data lines \n");
        fprintf(fp,"set nolabel \n");
        fprintf(fp,"set autoscale \n");
        fprintf(fp,"set xr[0:%f] \n",2*r_ref);
        fprintf(fp,"set grid \n");
        fprintf(fp,"set nokey \n");
        fprintf(fp,"set nolabel \n");

        fprintf(fp,"set xlabel \"r (a0)\" \n");
        fprintf(fp,"set title \" %s reference potential\n",atom_name);

        fprintf(fp,"plot \"%s\" using 1:2 \n",data_filename);

        fprintf(fp,"\n");
        fprintf(fp,"pause -1\n");
        fclose(fp);


        sprintf(script_filename,"%s%s_loc_pot.plt",paw_sdir(),atom_name);

        fp = fopen(script_filename,"w+");

        fprintf(fp,"set style data lines \n");
        fprintf(fp,"set nolabel \n");
        fprintf(fp,"set autoscale \n");
        fprintf(fp,"set xr[0:%f] \n",2*r_ref);
        fprintf(fp,"set grid \n");
        fprintf(fp,"set nokey \n");
        fprintf(fp,"set nolabel \n");

        fprintf(fp,"set xlabel \"r (a0)\" \n");
        fprintf(fp,"set title \" %s local potential\n",atom_name);

        fprintf(fp,"plot \"%s\" using 1:3 \n",data_filename);

        fprintf(fp,"\n");
        fprintf(fp,"pause -1\n");
        fclose(fp);


        /*gnu script file */
        for (i=0,j=4; i<=nbasis-1; i++,j=j+1)
        {

            sprintf(nl_name,"%d%s",prin_n[i],paw_spd_Name(orb_l[i]));

            sprintf(script_filename,"%s%s_%s_pot.plt",paw_sdir(),atom_name,nl_name);

            fp = fopen(script_filename,"w+");

            fprintf(fp,"set style data lines \n");
            fprintf(fp,"set nolabel \n");
            fprintf(fp,"set autoscale \n");
            fprintf(fp,"set xr[0:%f] \n",1.5*r_potential[i]);
            fprintf(fp,"set grid \n");
            fprintf(fp,"set nokey \n");
            fprintf(fp,"set nolabel \n");

            fprintf(fp,"set xlabel \"r (a0)\" \n");
            fprintf(fp,"set title \" %s paw potential for %s orbital\" \n",
                    atom_name,nl_name);

            fprintf(fp,"plot \"%s\" using 1:%d \n",data_filename,j);

            fprintf(fp,"\n");
            fprintf(fp,"pause -1\n");
            fclose(fp);

        }
    }
}










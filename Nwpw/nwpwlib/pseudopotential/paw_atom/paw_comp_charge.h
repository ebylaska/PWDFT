#ifndef _PAW_COMP_CHARGE_H_
#define _PAW_COMP_CHARGE_H_
/*
   $Id$
*/


extern void paw_init_comp_charge( float );
extern float paw_boundary_function(float );
extern float paw_find_sigma_comp();
extern float* paw_find_comp_charge_potential(float,float , float);
extern float paw_get_sigma_comp();
extern float paw_get_comp_charge_radius();
extern void  paw_print_comp_charge_information(FILE *fp);

#endif


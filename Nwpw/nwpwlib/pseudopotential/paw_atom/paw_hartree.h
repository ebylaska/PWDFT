#ifndef _PAW_HARTREE_H_
#define _PAW_HARTREE_H_
/*
   $Id$
*/


extern void paw_init_hartree();
extern void paw_generate_hartree_pot(float *n);
extern float paw_get_hartree_energy(float *n);
extern float* paw_get_hartree_pot();
#endif



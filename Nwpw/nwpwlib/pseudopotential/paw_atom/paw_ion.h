#ifndef _PAW_ION_H_
#define _PAW_ION_H_
/*
   $Id$
*/


extern void   paw_init_ion(float Z);
extern float paw_get_ion_energy(float *dn);
extern float* paw_get_ion_pot();
extern float paw_get_ion_charge();

#endif



#ifndef _PAW_ION_H_
#define _PAW_ION_H_
/*
   $Id$
*/

extern void paw_init_ion(double);
extern double paw_get_ion_energy(double *);
extern double *paw_get_ion_pot();
extern double paw_get_ion_charge();

#endif

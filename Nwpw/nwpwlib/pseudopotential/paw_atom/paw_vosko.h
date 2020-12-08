#ifndef _PAW_VOSKO_H_
#define _PAW_VOSKO_H_
/*
   $Id$
*/


extern void paw_init_vosko();
extern void paw_generate_corr_pot(float **rho, float **Vc);
extern float paw_get_correlation_energy(float **rho);
extern void paw_generate_corr_pot_LDA(float *rho);
extern float* paw_get_corr_pot_LDA();
extern float paw_get_correlation_energy_LDA(float *rho);

#endif



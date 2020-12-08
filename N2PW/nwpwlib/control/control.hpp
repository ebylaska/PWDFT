#ifndef _CONTROL_H_
#define _CONTROL_H_
/* control.h
   Author - Eric Bylaska

*/


#include	"Parallel.hpp"
#include	"rtdb.hpp"

extern void control_read(RTDB&);
extern int control_mapping();
extern int control_mapping1d();
extern int control_balance();
extern int control_np_orbital();
extern int control_ngrid(const int);
extern int control_ewald_ngrid(const int);
extern int control_ewald_ncut();
extern int control_loop(const int);
extern int control_pfft3_qsize();

extern float control_ewald_rcut();
extern float control_unita(const int, const int);
extern float control_ecut();
extern float control_wcut();
extern float control_time_step();
extern float control_fake_mass();
extern float control_tolerances(const int);

extern char *control_input_movecs_filename();

#endif

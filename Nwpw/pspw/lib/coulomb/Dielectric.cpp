/* Dielectric.cpp
   Author - Eric Bylaska
*/

/*
#include	<string>
#include        <iostream>
#include        <cstdio>
#include        <stdio.h>
#include        <cstdlib>
*/
#include <cmath>
#include "util.hpp"
#include	"Dielectric.hpp"


namespace pwdft {

/********************************************
 *                                          *
 * Dielectric_Operator::Dielectric_Operator *
 *                                          *
 ********************************************/
Dielectric_Operator::Dielectric_Operator(Pneb *mygrid, Control2& control)
{
   mypneb = mygrid;

   n2ft3d = mypneb->n2ft3d;

   r_grid     = mypneb->r_nalloc(3);
   mypneb->generate_r_sym_grid(r_grid);

   epsilon    = mypneb->r_alloc();
   sw         = mypneb->r_alloc();
   depsilondx = mypneb->r_alloc();
   depsilondy = mypneb->r_alloc();
   depsilondz = mypneb->r_alloc();

   dielec = control.gpoisson_dielec();
   rho0   = control.gpoisson_rho0();
   beta   = control.gpoisson_beta();
}

/********************************************
 *                                          *
 *     Dielectric_Operator::rho_generate    *
 *                                          *
 ********************************************/
void Dielectric_Operator::rho_generate(double *rho)
{
   double scal1 = 1.0/((double) ((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));
   double *Gx = mypneb->Gpackxyz(0,0);
   double *Gy = mypneb->Gpackxyz(0,1);
   double *Gz = mypneb->Gpackxyz(0,2);

   util_fattebert_dielec(n2ft3d,dielec,beta,rho0,rho,epsilon);

   /* calculate rho sw=rho(g) */
   mypneb->rr_SMul(scal1,epsilon,sw);
   mypneb->rc_fft3d(sw);
   mypneb->c_pack(0,sw);

   /* calculate (depsilondx,depsilondy,depsilondz) = grad epsilon */
   mypneb->tcc_pack_iMul(0,Gx,sw,depsilondx);
   mypneb->tcc_pack_iMul(0,Gy,sw,depsilondy);
   mypneb->tcc_pack_iMul(0,Gz,sw,depsilondz);
   mypneb->c_unpack(0,depsilondx);
   mypneb->c_unpack(0,depsilondy);
   mypneb->c_unpack(0,depsilondz);
   mypneb->cr_fft3d(depsilondx);
   mypneb->cr_fft3d(depsilondy);
   mypneb->cr_fft3d(depsilondz);
}


}

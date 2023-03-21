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

   mypneb->initialize_r_grid();
   r_grid = mypneb->r_grid;

   epsilon = mypneb->r_alloc();
   sw      = mypneb->r_alloc();
   p       = mypneb->r_alloc();

   dielec = control.gpoisson_dielec();
   rho0   = control.gpoisson_rho0();
   beta   = control.gpoisson_beta();
}

/********************************************
 *                                          *
 *  Dielectric_Operator::generate_dielec    *
 *                                          *
 ********************************************/
void Dielectric_Operator::generate_dielec(const double *rho)
{
   double scal1 = 1.0/((double) ((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));
   double *Gx = mypneb->Gpackxyz(0,0);
   double *Gy = mypneb->Gpackxyz(0,1);
   double *Gz = mypneb->Gpackxyz(0,2);

   util_fattebert_dielec(n2ft3d,dielec,beta,rho0,rho,epsilon);

   ///* calculate rho sw=rho(g) */
   //mypneb->rr_SMul(scal1,epsilon,sw);
   //mypneb->rc_fft3d(sw);
   //mypneb->c_pack(0,sw);

   /* calculate p = laplacian(sqrt(epsilon))/sqrt(epsilon) */
   for (auto i=0; i<n2ft3d; ++i)
       sw[i] = std::sqrt(epsilon[i]);
   mypneb->rr_SMul(scal1,epsilon,sw);
   mypneb->rc_fft3d(sw);
   mypneb->c_pack(0,sw);
   int kk = 0;
   for (auto k=0; k<(mypneb->npack(0)); ++k)
   {
      p[kk]   = -sw[kk]  *(Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k]);
      p[kk+1] = -sw[kk+1]*(Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k]);
      kk+=2;
   }
   mypneb->c_unpack(0,p);
   mypneb->cr_fft3d(p);
   for (auto i=0; i<n2ft3d; ++i)
       p[i] /= std::sqrt(epsilon[i]);
}

/**********************************************
 *                                            *
 *    Dielectric_Operator::generate_scaled    *
 *                                            *
 **********************************************/
void Dielectric_Operator::generate_scaled(double *scaled_rho)
{
   for (auto i=0; i<n2ft3d; ++i)
       scaled_rho[i] /= std::sqrt(epsilon[i]);
}



}

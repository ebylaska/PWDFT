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
 *            generate_rho_ion              *
 *                                          *
 ********************************************/
/*static void generate_rho_ion(Pneb *mypneb, Strfac *mystrfac, double rc,  double *rho_ion)
{
    int nion       = mystrfac->myion->nion;
    double *gauss  = mypneb->t_pack_allocate(0);
    double *katm   = mystrfac->myion->katm;
    double *zv_psp = mystrfac->myion->zv_psp;

    mypneb->c_pack_zero(0,rho_ion);
    mypneb->t_pack_gaussian(0,rc,gauss);
    for (auto ii=0; ii<nion; ++ii)
    {
       mystrfac->strfac_pack(0,ii,exi)
       mypneb->tcc_pack_aMulAdd(0,zv_psp[katm[ii]],gauss,exi,rho_ion);
    }
}
*/

/********************************************
 *                                          *
 * Dielectric_Operator::Dielectric_Operator *
 *                                          *
 ********************************************/
Dielectric_Operator::Dielectric_Operator(Pneb *mygrid, Strfac *mystrfacin, Control2& control)
{
   mypneb   = mygrid;
   mystrfac = mystrfacin;

   n2ft3d = mypneb->n2ft3d;

   mypneb->initialize_r_grid();
   r_grid = mypneb->r_grid;

   epsilon  = mypneb->r_alloc();
   depsilon = mypneb->r_alloc();
   sw      = mypneb->r_alloc();
   p       = mypneb->r_alloc();

   epsilon_x = mypneb->r_alloc();
   epsilon_y = mypneb->r_alloc();
   epsilon_z = mypneb->r_alloc();

   w_x = mypneb->r_alloc();
   w_y = mypneb->r_alloc();
   w_z = mypneb->r_alloc();

   rho_induced = mypneb->r_alloc();
   rho_ion     = mypneb->r_alloc();
   //generate_rho_ion(mypneb,mystrfac,rho_ion);

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

   std::cout << "dielec=" << dielec << std::endl;
   std::cout << "rho0=  " << rho0 << std::endl;
   std::cout << "beta=  " << beta << std::endl;
   util_fattebert_dielec(n2ft3d,dielec,beta,rho0,rho,epsilon);
   //util_dfattebert_dielec(n2ft3d,dielec,beta,rho0,rho,depsilon);

   /* calculate fft of epsilon */
   mypneb->rr_SMul(scal1,epsilon,sw);
   mypneb->rc_fft3d(sw);
   mypneb->c_pack(0,sw);

   /* calculate <epsilon_x,epsilon_y,epsilon_z> = grad epsilon */
   mypneb->tcc_pack_iMul(0,Gx,sw,epsilon_x);
   mypneb->tcc_pack_iMul(0,Gy,sw,epsilon_y);
   mypneb->tcc_pack_iMul(0,Gz,sw,epsilon_z);
   mypneb->c_unpack(0,epsilon_x);
   mypneb->c_unpack(0,epsilon_y);
   mypneb->c_unpack(0,epsilon_z);
   mypneb->cr_fft3d(epsilon_x);
   mypneb->cr_fft3d(epsilon_y);
   mypneb->cr_fft3d(epsilon_z);

   /* calculate p = laplacian(sqrt(epsilon))/sqrt(epsilon) */
   for (auto i=0; i<n2ft3d; ++i)
       sw[i] = std::sqrt(epsilon[i]);
   mypneb->r_SMul(scal1,sw);
   mypneb->rc_fft3d(sw);
   mypneb->c_pack(0,sw);
   int kk = 0;
   for (auto k=0; k<(mypneb->npack(0)); ++k)
   {
      double gg = (Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k]);

      p[kk]   = -sw[kk]  *gg;
      p[kk+1] = -sw[kk+1]*gg;
      kk+=2;
   }
   std::cout << "p[0]=" << Efmt(20,15) <<  p[0] << " " << p[1] << std::endl;


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

/************************************************
 *                                              *
 *  Dielectric_Operator::generate_over_epsilon  *
 *                                              *
 ************************************************/
void Dielectric_Operator::generate_over_epsilon(double *scaled_rho)
{
   for (auto i=0; i<n2ft3d; ++i)
       scaled_rho[i] /= epsilon[i];
}

/**********************************************
 *                                            *
 *  Dielectric_Operator::generate_dpotential  *
 *                                            *
 **********************************************/
void Dielectric_Operator::generate_dpotential(const double *w, double *dwdeps)
{
   double omega  = mypneb->lattice->omega();
   double scal1  = 1.0/((double) ((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));
   double scal2  = 1.0/omega;
   double fourpi = 16.0*std::atan(1.0);
   double dv     = omega*scal1;
   double *Gx = mypneb->Gpackxyz(0,0);
   double *Gy = mypneb->Gpackxyz(0,1);
   double *Gz = mypneb->Gpackxyz(0,2);

   /* calculate <w_x,w_y,w_z> = grad w */
   mypneb->tcc_pack_iMul(0,Gx,w,w_x);
   mypneb->tcc_pack_iMul(0,Gy,w,w_y);
   mypneb->tcc_pack_iMul(0,Gz,w,w_z);
   mypneb->c_unpack(0,w_x);
   mypneb->c_unpack(0,w_y);
   mypneb->c_unpack(0,w_z);
   mypneb->cr_fft3d(w_x);
   mypneb->cr_fft3d(w_y);
   mypneb->cr_fft3d(w_z);
   for (auto i=0; i<n2ft3d; ++i)
      dwdeps[i] = ( w_x[i]*epsilon_x[i]
                  + w_y[i]*epsilon_y[i]
                  + w_z[i]*epsilon_z[i] )/(fourpi*epsilon[i]);
   mypneb->r_SMul(scal2,dwdeps);

   double sumdw = mypneb->r_dsum(dwdeps)*dv; 
   std::cout << "sumdw =" << sumdw << std::endl;

   mypneb->r_abs(w_x);
   mypneb->r_abs(w_y);
   mypneb->r_abs(w_z);
   double sumx = mypneb->r_dsum(w_x)*dv; 
   double sumy = mypneb->r_dsum(w_y)*dv; 
   double sumz = mypneb->r_dsum(w_y)*dv; 
   std::cout << "sum_xyz =" << sumx << " " << sumy << " " << sumz << std::endl;

   double esumx = mypneb->r_dsum(epsilon_x)*dv; 
   double esumy = mypneb->r_dsum(epsilon_y)*dv; 
   double esumz = mypneb->r_dsum(epsilon_y)*dv; 
   std::cout << "esum_xyz =" << esumx << " " << esumy << " " << esumz << std::endl;

   mypneb->rr_copy(epsilon_x,w_x);
   mypneb->rr_copy(epsilon_y,w_y);
   mypneb->rr_copy(epsilon_z,w_z);
   mypneb->r_abs(w_x);
   mypneb->r_abs(w_y);
   mypneb->r_abs(w_z);
   double asumx = mypneb->r_dsum(w_x)*dv; 
   double asumy = mypneb->r_dsum(w_y)*dv; 
   double asumz = mypneb->r_dsum(w_y)*dv; 
   std::cout << "aesum_xyz =" << asumx << " " << asumy << " " << asumz << std::endl;

}


}

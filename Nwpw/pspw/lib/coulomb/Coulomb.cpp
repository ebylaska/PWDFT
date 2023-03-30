/* Kinetic.C - 
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
#include	"Coulomb.hpp"


namespace pwdft {

/*******************************************
 *                                         *
 *     Coulomb_Operator::Coulomb_Operator  *
 *                                         *
 *******************************************/
Coulomb_Operator::Coulomb_Operator(Pneb *mygrid, Control2& control)
{
   int k,pzero,zero,taskid;
   double gg;
   double *Gx  = mygrid->Gxyz(0);
   double *Gy  = mygrid->Gxyz(1);
   double *Gz  = mygrid->Gxyz(2);
   vg          = new double [mygrid->npack(0)];
   double *tmp = new double [mygrid->nfft3d];
   double fourpi = 16.0*atan(1.0);
   //double kk0 = 0.0; // where k0 = the uniform background potential 

   mypneb = mygrid;

   if (control.gpoisson_on())
   {
      has_dielec = true;
      mydielec   = new (std::nothrow) Dielectric_Operator(mypneb,control);
    //  kk0 = 1.0;
   }


   taskid = mypneb->d3db::parall->taskid_i();
   pzero  = mypneb->ijktop(0,0,0);
   zero   = mypneb->ijktoindex(0,0,0);

   for (k=0; k<(mypneb->nfft3d); ++k)
   {
      gg     = Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k];
      if ((pzero==taskid)&&(k==zero))
         tmp[k] = 0.0;
      else
         tmp[k] = fourpi/gg;
   }
   mypneb->t_pack(0,tmp);
   mypneb->tt_pack_copy(0,tmp,vg);

   delete [] tmp;


}



/*******************************************
 *                                         *
 *        Coulomb_Operator::vcoulomb       *
 *                                         *
 *******************************************/
void Coulomb_Operator::vcoulomb(const double *dng, double *vcout)
{
   int k,k1,ksize;

   ksize = (mypneb->npack(0));
   k1 = 0;
   for (k=0; k<ksize; ++k)
   {
      vcout[k1]   = vg[k]*dng[k1]; 
      vcout[k1+1] = vg[k]*dng[k1+1];
      k1 += 2;
   }
}

/*******************************************
 *                                         *
 *        Coulomb_Operator::ecoulomb       *
 *                                         *
 *******************************************/
double Coulomb_Operator::ecoulomb(const double *dng)
{
   int k,k1,k2,n,nsize,ksize1,ksize2;
   double ave;

   ksize1 = (mypneb->nzero(0));
   ksize2 = (mypneb->npack(0));
   ave = 0.0;
   k1  = 0;
   for (k=0; k<ksize1; ++k)
   {
      ave += vg[k]*(dng[k1]*dng[k1] + dng[k1+1]*dng[k1+1]); 
      k1  += 2;
   }
   for (k=ksize1; k<ksize2; ++k)
   {
      ave += 2.0*vg[k]*(dng[k1]*dng[k1] + dng[k1+1]*dng[k1+1]); 
      k1  += 2;
   }
   ave = mypneb->d3db::parall->SumAll(1,ave);
   //ave *= 0.5*lattice_omega();
   ave *= 0.5*(mypneb->lattice->omega());

   return ave;
}

/*******************************************
 *                                         *
 *    Coulomb_Operator::vcoulomb_dielec    *
 *                                         *
 *******************************************/
void Coulomb_Operator::vcoulomb_dielec(const double *rho, double *vcout)
{
    /* genereate dielectric and p function*/
    mydielec->generate_dielec(rho);

    double omega  = mypneb->lattice->omega();
    double scal1  = 1.0/((double) ((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));
    double scal2  = 1.0/omega;
    double fourpi = 16.0*atan(1.0);
    double *q0 = mypneb->r_alloc();
    double *q  = mypneb->r_alloc();
    double *w  = mypneb->r_alloc();
    double *p0 = mypneb->r_alloc();
    double dv = omega*scal1;
    double fac = 1.0/omega;
    double *tmp = mypneb->r_alloc();

    mypneb->rr_SMul((1.0/fourpi),mydielec->p,p0);

    double sum = mypneb->r_dsum(rho)*dv;
    std::cout << "sum(rho)*dv=" << Ffmt(20,15) << sum << std::endl;

    mypneb->rr_copy(rho,q0);

    // convert rho to q, q = rho/sqrt(epsilon);
    mydielec->generate_scaled(q0);
    mypneb->rr_copy(q0,q);

    double sum1=0.0;
    double sum2=0.0;
    double sum9=0.0;
    double sum1old=0.0;
    double sum2old=0.0;
    double sum9old=0.0;
    for (auto it=0; it<20; ++it) {
       mypneb->r_SMul(scal1,q);
       mypneb->r_zero_ends(q);
       mypneb->rc_fft3d(q);
       mypneb->c_pack(0,q);

       // set initial w in G-space - w = (4*pi/G2)*q
       this->vcoulomb(q,w);

       //Solve inhomogeneous Helmholtz equation in real-space
       //  -laplacian[ w(k) ] + FFT[p(r)*w(r)] = q(k)
       //mypneb->c_pack_SMul(0,fourpi,q);
       //mypneb->rcc_solve_Helmholtz(0,p,q,w);

       mypneb->c_unpack(0,w);
       mypneb->cr_fft3d(w);
       mypneb->r_zero_ends(w);

       //q = q0 + p*w
       mypneb->rr_copy(q0,q);
       //mypneb->rrr_Mul2Add(p0,w,q);
       mypneb->rrr_Mul(p0,w,tmp);
       mypneb->rr_daxpy(fac,tmp,q);
       std::cout <<  "fac=" << fac << std::endl;
       

       sum1old = sum1;
       sum2old = sum2;
       sum9old = sum9;
       sum1 = mypneb->r_dsum(q)*dv;
       sum2 = mypneb->rr_dot(rho,w)*dv;
       sum9 = mypneb->rr_dot(p0,w)*dv;
       std::cout << "sum(q)*dv    =" << Ffmt(20,15) << sum1 << " " << sum1-sum1old << std::endl;
       std::cout << "dot(rho,w)*dv=" << Ffmt(20,15) << sum2 << " " << sum2-sum2old << std::endl;
       std::cout << "dot(po,w)*dv =" << Ffmt(20,15) << sum9 << " " << sum9-sum9old << std::endl;
    }

    // convert w to vcout, vcout = w/sqrt(epsilon);
    mydielec->generate_scaled(w);

    double sum3 = mypneb->rr_dot(rho,w)*dv;
    std::cout << std::endl;
    std::cout << "Final dot(rho,w)*dv=" << Ffmt(20,15) << sum3 << std::endl;


    // Fourier transform from real-space to reciprocal space
    mypneb->r_SMul(scal1,w);
    mypneb->rc_fft3d(w);
    mypneb->c_pack(0,w);

    mypneb->cc_pack_copy(0,w,vcout);

    mypneb->r_dealloc(tmp);
    mypneb->r_dealloc(p0);
    mypneb->r_dealloc(w);
    mypneb->r_dealloc(q);
    mypneb->r_dealloc(q0);
}

/*******************************************
 *                                         *
 *    Coulomb_Operator::vcoulomb_dielec2   *
 *                                         *
 *******************************************/
void Coulomb_Operator::vcoulomb_dielec2(const double *rho, double *dng, double *vcout)
{
    /* genereate dielectric and p function*/
    mydielec->generate_dielec(rho);

    double omega  = mypneb->lattice->omega();
    double scal1  = 1.0/((double) ((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));
    double scal2  = 1.0/omega;
    double fourpi = 16.0*atan(1.0);
    double dv     = omega*scal1;
    double *q0 = mypneb->r_alloc();
    double *q  = mypneb->r_alloc();
    double *w  = mypneb->r_alloc();
    double *dwdeps_r = mypneb->r_alloc();

    double sum1=0.0;
    double sum2=0.0;
    double sum3=0.0;
    double sum1old=0.0;
    double sum2old=0.0;
    double sum3old=0.0;


    // convert rho to q, q = rho/epsilon;
    mypneb->rr_copy(rho,q0);
    mydielec->generate_over_epsilon(q0);
    mypneb->rr_copy(q0,q);

    for (auto it=0; it<3; ++it) {
       mypneb->r_SMul(scal1,q);
       mypneb->r_zero_ends(q);
       mypneb->rc_fft3d(q);
       mypneb->c_pack(0,q);

       // set initial w in G-space - w = (4*pi/G2)*q
       this->vcoulomb(q,w);

       //q = q0 + deps*dw
       mydielec->generate_dpotential(w,dwdeps_r);
       mypneb->rrr_Sum(q0,dwdeps_r,q);

       sum1old = sum1;
       sum2old = sum2;
       sum3old = sum2;
       sum1 = mypneb->r_dsum(q)*dv;
       sum2 = mypneb->cc_pack_dot(0,dng,w)*omega;
       sum3 = mypneb->r_dsum(dwdeps_r)*dv;
       std::cout << "sum(q)*dv        =" << Ffmt(20,15) << sum1 << " " << sum1-sum1old << std::endl;
       std::cout << "sum(dwdeps_r)*dv =" << Ffmt(20,15) << sum3 << " " << sum3-sum3old << std::endl;
       std::cout << "dot(dng,w)*omega =" << Ffmt(20,15) << sum2 << " " << sum2-sum2old << std::endl;
    }

    mypneb->cc_pack_copy(0,w,vcout);

    mypneb->r_dealloc(dwdeps_r);
    mypneb->r_dealloc(w);
    mypneb->r_dealloc(q);
    mypneb->r_dealloc(q0);
}


}

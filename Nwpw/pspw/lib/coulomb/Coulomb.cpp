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
void Coulomb_Operator::vcoulomb_dielec(const double *dng, double *vcout)
{
    double scal1  = 1.0/((double) ((mypneb->nx)*(mypneb->ny)*(mypneb->nz)));
    double fourpi = 16.0*atan(1.0);
    double *q = mypneb->r_alloc();
    double *w = mypneb->r_alloc();
    double *p = mydielec->p;


    mypneb->cc_pack_copy(0,dng,q);
    mypneb->c_unpack(0,q);
    mypneb->cr_fft3d(q);

    mydielec->generate_dielec(q);

    mydielec->generate_scaled(q);
    mypneb->r_SMul(-fourpi,q);

    // set initial w in real-space
    vcoulomb(dng,w);
    mypneb->c_unpack(0,w);
    mypneb->cr_fft3d(w);

 
    //Solve inhomogeneous Helmholtz equation in real-space
    //  laplacian[ w(r) ] + p(r)*w(r) = q(r)
    mypneb->rrr_solve_Helmholtz(0,p,q,w);

    // convert w to vcout, vcout = w/sqrt(epsilon);
    mydielec->generate_scaled(w);


    // Fourier transform from real-space to reciprocal space
    mypneb->r_SMul(scal1,w);
    mypneb->rc_fft3d(w);
    mypneb->c_pack(0,w);
    mypneb->cc_pack_copy(0,w,vcout);


    mypneb->r_dealloc(w);
    mypneb->r_dealloc(q);
}


}

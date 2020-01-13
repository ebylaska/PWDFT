/* Kinetic.C - 
   Author - Eric Bylaska
*/

/*
using namespace std;

#include	<string.h>
#include        <iostream>
#include        <cstdio>
#include        <stdio.h>
#include        <cstdlib>
*/
#include        <cmath>

#include	"PGrid.hpp"
#include	"Coulomb.hpp"


/*******************************************
 *                                         *
 *     Coulomb_Operator::Coulomb_Operator  *
 *                                         *
 *******************************************/
Coulomb_Operator::Coulomb_Operator(Pneb *mygrid)
{
   int k,pzero,zero,taskid;
   double gg;
   double *Gx  = mygrid->Gxyz(0);
   double *Gy  = mygrid->Gxyz(1);
   double *Gz  = mygrid->Gxyz(2);
   vg          = new double [mygrid->npack(0)];
   double *tmp = new double [mygrid->nfft3d];
   double fourpi = 16.0*atan(1.0);

   mypneb = mygrid;

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



void Coulomb_Operator::vcoulomb(double *dng, double *vcout)
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

double Coulomb_Operator::ecoulomb(double *dng)
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
   ave *= 0.5*lattice_omega();

   return ave;
}

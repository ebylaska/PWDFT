/* Lattice.cpp
   Author - Eric Bylaska

*/

#include	<cmath>
#include	"Control2.hpp"
#include	"Lattice.hpp"


static void get_cube(double *unita, double *unitg, double *omega)
{
   double twopi = 8.0*atan(1.0);
   unitg[0] = unita[4]*unita[8] - unita[5]*unita[7];
   unitg[1] = unita[5]*unita[6] - unita[3]*unita[8];
   unitg[2] = unita[3]*unita[7] - unita[4]*unita[6];
   unitg[3] = unita[7]*unita[2] - unita[8]*unita[1];
   unitg[4] = unita[8]*unita[0] - unita[6]*unita[2];
   unitg[5] = unita[6]*unita[1] - unita[7]*unita[0];
   unitg[6] = unita[1]*unita[5] - unita[2]*unita[4];
   unitg[7] = unita[2]*unita[3] - unita[0]*unita[5];
   unitg[8] = unita[0]*unita[4] - unita[1]*unita[3];
   double volume = unita[0]*unitg[0] + unita[1]*unitg[1] + unita[2]*unitg[2];
   for (int i=0; i<9; ++i) 
      unitg[i] *= twopi/volume;
   *omega = fabs(volume);
}


/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

Lattice::Lattice(Control2& control)
{
   int nx,ny,nz,nxh,nyh,nzh;
   double gx,gy,gz,gg,gg1,gg2,gg3,ecut0,wcut0;

   ecut0 = control.ecut();
   wcut0 = control.wcut();
   punita[0] = control.unita(0,0);
   punita[1] = control.unita(1,0);
   punita[2] = control.unita(2,0);

   punita[3] = control.unita(0,1);
   punita[4] = control.unita(1,1);
   punita[5] = control.unita(2,1);

   punita[6] = control.unita(0,2);
   punita[7] = control.unita(1,2);
   punita[8] = control.unita(2,2);
   get_cube(punita,punitg,&pomega);


   nx = control.ngrid(0);
   ny = control.ngrid(1);
   nz = control.ngrid(2);
   nxh = nx/2;
   nyh = ny/2;
   nzh = nz/2;

   gx = punitg[0]*((double) nxh);
   gy = punitg[1]*((double) nxh);
   gz = punitg[2]*((double) nxh);
   gg1 = gx*gx + gy*gy + gz*gz;

   gx = punitg[3]*((double) nyh);
   gy = punitg[4]*((double) nyh);
   gz = punitg[5]*((double) nyh);
   gg2 = gx*gx + gy*gy + gz*gz;

   gx = punitg[6]*((double) nzh);
   gy = punitg[7]*((double) nzh);
   gz = punitg[8]*((double) nzh);
   gg3 = gx*gx + gy*gy + gz*gz;

   gg = gg1;
   if (gg2 < gg) gg=gg2;
   if (gg3 < gg) gg=gg3;

   pecut = 0.50*gg;
   if (ecut0 < pecut) pecut = ecut0;
   pwcut = pecut;
   if (wcut0 < pwcut) pwcut = wcut0;

}


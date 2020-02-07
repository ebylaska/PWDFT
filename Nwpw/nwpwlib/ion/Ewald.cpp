/* Ewald.C -
   Author - Eric Bylaska
*/

using namespace std;

#include	<string.h>


#include        <iostream>
#include        <cstdio>
#include        <stdio.h>
#include        <cmath>
#include        <cstdlib>

#include	"control.hpp"
#include	"Ewald.hpp"
//#include	"Pseudopotential.hpp"

/*************************************
 *                                   *
 *          mandelung_get            *
 *                                   *
 *************************************/

double  mandelung_get()
{
   int n1,n2,n3;
   double ax,ay,az,gx,gy,gz,gg;
   double alpha,alpha1,alpha2,sum,ea;
   double rc,rs,epsilon,pi;
   int N=40;

   pi = 4.0*atan(1.0);
   rs = pow((3.0*lattice_omega()/(4.0*pi)),(1.0/3.0));
   rc = rs;
   epsilon = 1.0/rc;

   sum = 0.0;
   for (n1=(-N+1); n1<=(N-1); ++n1)
   for (n2=(-N+1); n2<=(N-1); ++n2)
   for (n3=(-N+1); n3<=(N-1); ++n3)
      if ((n1!=0)||(n2!=0)||(n3!=0))
      {
         ax = n1*lattice_unita(0,0) + n2*lattice_unita(0,1) + n3*lattice_unita(0,2);
         ay = n1*lattice_unita(1,0) + n2*lattice_unita(1,1) + n3*lattice_unita(1,2);
         az = n1*lattice_unita(2,0) + n2*lattice_unita(2,1) + n3*lattice_unita(2,2);
         ea = sqrt(ax*ax + ay*ay + az*az);
         sum += erfc(epsilon*ea)/ea;
      }
   alpha1 = sum;

   /* calculate alpha2 */
   sum = 0.0;
   for (n1=(-N+1); n1<=(N-1); ++n1)
   for (n2=(-N+1); n2<=(N-1); ++n2)
   for (n3=(-N+1); n3<=(N-1); ++n3)
      if ((n1!=0)||(n2!=0)||(n3!=0))
      {
         gx = n1*lattice_unitg(0,0) + n2*lattice_unitg(0,1) + n3*lattice_unitg(0,2);
         gy = n1*lattice_unitg(1,0) + n2*lattice_unitg(1,1) + n3*lattice_unitg(1,2);
         gz = n1*lattice_unitg(2,0) + n2*lattice_unitg(2,1) + n3*lattice_unitg(2,2);
         gg = gx*gx + gy*gy + gz*gz;
         sum += (4.0*pi/gg)*exp(-gg*rc*rc/4.0);
      }
   alpha2 = sum/lattice_omega();

   sum = alpha1 + alpha2 - pi*rc*rc/lattice_omega() - 2.0*epsilon/sqrt(pi);
   alpha = -sum*rs;
   return alpha;
}


/* Constructors */

/*********************************
 *                               *
 *          Ewald::Ewald         *
 *                               *
 *********************************/
//Ewald::Ewald(Parallel *inparall, Ion *inion, Pseudopotential *inpsp)
Ewald::Ewald(Parallel *inparall, Ion *inion, double *inzv)
{
   int i,j,k,l,k1,k2,k3;
   int enxh,enyh,enzh,enpack0;
   int tnp,tid,dutask;
   double g1,g2,g3,gg1,gg2,gg3,gg,ggcut;
   double pi,pi4,rs,w,term;
   double q,z,zz;
   double eps=1.0e-12;

   ewaldparall = inparall;
   ewaldion    = inion;
   tnp = ewaldparall->np();
   tid = ewaldparall->taskid();

   for (j=0; j<3; ++j)
   for (i=0; i<3; ++i)
   {
      unitg[i+j*3] = lattice_unitg(i,j);
      unita[i+j*3] = lattice_unita(i,j);
   }
      
   enx=control_ewald_ngrid(0);
   eny=control_ewald_ngrid(1);
   enz=control_ewald_ngrid(2);
   enxh=enx/2;
   enyh=eny/2;
   enzh=enz/2;


   /* determine ggcut */
   g1 = unitg[0]*(enxh);
   g2 = unitg[1]*(enxh);
   g3 = unitg[2]*(enxh);
   gg1 = g1*g1 + g2*g2 + g3*g3;

   g1 = unitg[3]*(enyh);
   g2 = unitg[4]*(enyh);
   g3 = unitg[5]*(enyh);
   gg2 = g1*g1 + g2*g2 + g3*g3;

   g1 = unitg[6]*(enzh);
   g2 = unitg[7]*(enzh);
   g3 = unitg[8]*(enzh);
   gg3 = g1*g1 + g2*g2 + g3*g3;

   ggcut = gg1;
   if (gg2<ggcut) ggcut = gg2;
   if (gg3<ggcut) ggcut = gg3;
   if ((2.0*control_ecut())<ggcut) 
      ggcut=2.0*control_ecut();
   eecut = 0.5*ggcut;


   /* determine enpack */
   dutask= 0;
   enpack= 0;
   enida = 0;
   k1 = 0;
   k2 = 0;
   k3 = 0;
   g1=k1*unitg[0]+k2*unitg[3]+k3*unitg[6];
   g2=k1*unitg[1]+k2*unitg[4]+k3*unitg[7];
   g3=k1*unitg[2]+k2*unitg[5]+k3*unitg[8];
   gg=g1*g1+g2*g2+g3*g3;
   if ((gg-ggcut)<(-eps))
   {
      if (dutask==tid)
      { 
         enpack++;
         enida++;
      }
      dutask = (dutask+1)%tnp;
   }
   for (k=1; k<enzh; ++k)
   {
      k1 = 0;
      k2 = 0;
      k3 = k;
      g1=k1*unitg[0]+k2*unitg[3]+k3*unitg[6];
      g2=k1*unitg[1]+k2*unitg[4]+k3*unitg[7];
      g3=k1*unitg[2]+k2*unitg[5]+k3*unitg[8];
      gg=g1*g1+g2*g2+g3*g3;
      if ((gg-ggcut)<(-eps))
      {
         if (dutask==tid) enpack++;
         dutask =(dutask+1)%tnp;
      }
   }

   for (k=(-enzh+1); k<enzh; ++k)
   for (j=1; j<enyh; ++j)
   {
      k1 = 0;
      k2 = j;
      k3 = k;
      g1=k1*unitg[0]+k2*unitg[3]+k3*unitg[6];
      g2=k1*unitg[1]+k2*unitg[4]+k3*unitg[7];
      g3=k1*unitg[2]+k2*unitg[5]+k3*unitg[8];
      gg=g1*g1+g2*g2+g3*g3;
      if ((gg-ggcut)<(-eps))
      {
         if (dutask==tid) enpack++;
         dutask = (dutask+1)%tnp;
      }
   }
   

   for (k=(-enzh+1); k<enzh; ++k)
   for (j=(-enyh+1); j<enyh; ++j)
   for (i=1; i<enxh; ++i)
   {
      k1=i;
      k2=j;
      k3=k;
      g1=k1*unitg[0]+k2*unitg[3]+k3*unitg[6];
      g2=k1*unitg[1]+k2*unitg[4]+k3*unitg[7];
      g3=k1*unitg[2]+k2*unitg[5]+k3*unitg[8];
      gg=g1*g1+g2*g2+g3*g3;

      if ((gg-ggcut)< (-eps))
      {
         if (dutask==tid) enpack++;
         dutask = (dutask+1)%tnp;
      }
   }
   enpack_all = ewaldparall->ISumAll(0,enpack);



   encut = control_ewald_ncut();
   enshl3d = (2*encut+1)*(2*encut+1)*(2*encut+1);
   ercut = control_ewald_rcut();
   pi  = 4.00*atan(1.0);
   pi4 = 4.00*pi;
   if (encut<=0) encut=1;
   if (ercut<=0.00)
   {
      rs = unita[0]*unita[0] + unita[1]*unita[1] + unita[2]*unita[2];
      rs = sqrt(rs);
      ercut=rs/pi;

      rs = unita[3]*unita[3] + unita[4]*unita[4] + unita[5]*unita[5];
      rs = sqrt(rs);
      w=rs/pi;
      if (w<ercut) ercut = w;

      rs = unita[6]*unita[6] + unita[7]*unita[7] + unita[8]*unita[8];
      rs = sqrt(rs);
      w=rs/pi;
      if (w<ercut) ercut = w;
   }
   w = 0.25*ercut*ercut;


   /* allocate memory */
   eG  = new double [3*enpack];
   vg  = new double [enpack];
   ss  = new double [2*enpack];
   exi  = new double [2*enpack];
   tmp3 = new double [enpack];
   ftmp = new double [3*(ewaldion->nion)];
   vcx  = new double [enpack];
   rcell = new double [3*enshl3d];
   ewx1 = new double [2*(ewaldion->nion)*enx];
   ewy1 = new double [2*(ewaldion->nion)*eny];
   ewz1 = new double [2*(ewaldion->nion)*enz];
   zv   = new double [ewaldion->nkatm];
   i_indx = new int[enpack];
   j_indx = new int[enpack];
   k_indx = new int[enpack];

   /* determine eG */
   for (i=0; i<(3*enpack); ++i) eG[i] = 0.0;
   dutask= 0;
   enpack0= 0;
   k1 = 0;
   k2 = 0;
   k3 = 0;
   g1=k1*unitg[0]+k2*unitg[3]+k3*unitg[6];
   g2=k1*unitg[1]+k2*unitg[4]+k3*unitg[7];
   g3=k1*unitg[2]+k2*unitg[5]+k3*unitg[8];
   gg=g1*g1+g2*g2+g3*g3;
   if ((gg-ggcut)<(-eps))
   {
      if (dutask==tid)
      { 
         eG[enpack0]          = g1;
         eG[enpack0+  enpack] = g2;
         eG[enpack0+2*enpack] = g3;
         i=k1; if (i<0) i += enx;
         j=k2; if (j<0) j += eny;
         k=k3; if (k<0) k += enz;
         i_indx[enpack0] = i;
         j_indx[enpack0] = j;
         k_indx[enpack0] = k;
         enpack0++;
      }
      dutask = (dutask+1)%tnp;
   }
   k1 = 0;
   k2 = 0;
   for (k3=1; k3<enzh; ++k3)
   {
      g1=k1*unitg[0]+k2*unitg[3]+k3*unitg[6];
      g2=k1*unitg[1]+k2*unitg[4]+k3*unitg[7];
      g3=k1*unitg[2]+k2*unitg[5]+k3*unitg[8];
      gg=g1*g1+g2*g2+g3*g3;
      if ((gg-ggcut)<(-eps))
      {
         if (dutask==tid)
         {
            eG[enpack0]          = g1;
            eG[enpack0+  enpack] = g2;
            eG[enpack0+2*enpack] = g3;
            i=k1; if (i<0) i += enx;
            j=k2; if (j<0) j += eny;
            k=k3; if (k<0) k += enz;
            i_indx[enpack0] = i;
            j_indx[enpack0] = j;
            k_indx[enpack0] = k;
            enpack0++;
         }
         dutask =(dutask+1)%tnp;
      }
   }

   k1 = 0;
   for (k3=(-enzh+1); k3<enzh; ++k3)
   for (k2=1; k2<enyh; ++k2)
   {
      g1=k1*unitg[0]+k2*unitg[3]+k3*unitg[6];
      g2=k1*unitg[1]+k2*unitg[4]+k3*unitg[7];
      g3=k1*unitg[2]+k2*unitg[5]+k3*unitg[8];
      gg=g1*g1+g2*g2+g3*g3;
      if ((gg-ggcut)<(-eps))
      {
         if (dutask==tid)
         {
            eG[enpack0]          = g1;
            eG[enpack0+  enpack] = g2;
            eG[enpack0+2*enpack] = g3;
            i=k1; if (i<0) i += enx;
            j=k2; if (j<0) j += eny;
            k=k3; if (k<0) k += enz;
            i_indx[enpack0] = i;
            j_indx[enpack0] = j;
            k_indx[enpack0] = k;
            enpack0++;
         }
         dutask =(dutask+1)%tnp;
      }
   }


   for (k3=(-enzh+1); k3<enzh; ++k3)
   for (k2=(-enyh+1); k2<enyh; ++k2)
   for (k1=1; k1<enxh; ++k1)
   {
      g1=k1*unitg[0]+k2*unitg[3]+k3*unitg[6];
      g2=k1*unitg[1]+k2*unitg[4]+k3*unitg[7];
      g3=k1*unitg[2]+k2*unitg[5]+k3*unitg[8];
      gg=g1*g1+g2*g2+g3*g3;

      if ((gg-ggcut)< (-eps))
      {
         if (dutask==tid)
         {
            eG[enpack0]          = g1;
            eG[enpack0+  enpack] = g2;
            eG[enpack0+2*enpack] = g3;
            i=k1; if (i<0) i += enx;
            j=k2; if (j<0) j += eny;
            k=k3; if (k<0) k += enz;
            i_indx[enpack0] = i;
            j_indx[enpack0] = j;
            k_indx[enpack0] = k;
            enpack0++;
         }
         dutask = (dutask+1)%tnp;
      }
   }

   /* find vg and vcx */
   for (i=0; i<enpack; ++i) vg[i]  = 0.0;
   for (i=0; i<enpack; ++i) vcx[i] = 0.0;

   for (k=enida; k<enpack; ++k)
   {
      g1 = eG[k];
      g2 = eG[k +   enpack];
      g3 = eG[k + 2*enpack];
      gg = g1*g1 + g2*g2 + g3*g3;
      term = pi4/gg;
      vcx[k] = term;
      vg[k] = term*exp(-w*gg);
   }

   /* set the mandelung constant */
   alpha = mandelung_get();

   /* set the ion charges */
   for (i=0; i<(ewaldion->nkatm); ++i)
      zv[i] = inzv[i];
      //zv[i] = inpsp->zv[i];

   /* ewald summation */
   rs = pow(3.0*lattice_omega()/pi4, 1.0/3.0);
   zz = 0.0;
   z  = 0.0;
   for (i=0; i<(ewaldion->nion); ++i)
   { 
      q = zv[ewaldion->katm[i]];
      zz += q*q;
      z  += q;
   }
   cewald = 0.0;
   for (i=0; i<enpack; ++i) cewald += vg[i];
   cewald *= 2.0;
   if (tnp>1) cewald = ewaldparall->SumAll(0,cewald);

   cewald = -0.50*zz*(alpha/rs + cewald/lattice_omega())
          -  0.50*(z*z-zz)*ercut*ercut*pi/lattice_omega();

   /* set rcell */
   l = 0;
   rcell[l]             = 0.0;
   rcell[l +   enshl3d] = 0.0;
   rcell[l + 2*enshl3d] = 0.0;
   for (k=-encut; k<=encut; ++k)
   for (j=-encut; j<=encut; ++j)
   for (i=-encut; i<=encut; ++i)
      if ( !( (i==0)&&(j==0)&&(k==0) ) )
      {
         ++l;
         rcell[l]             = i*unita[0] + j*unita[3] + k*unita[6];
         rcell[l +   enshl3d] = i*unita[1] + j*unita[4] + k*unita[7];
         rcell[l + 2*enshl3d] = i*unita[2] + j*unita[5] + k*unita[8];
      }

}


/*********************************
 *                               *
 *          Ewald::phafac        *
 *                               *
 *********************************/
void Ewald::phafac()
{
   int i,k,enxh,enyh,enzh;
   double a,b,sw1,sw2,sw3,pi;
   double cw1x,cw2x,cw3x;
   double cw1y,cw2y,cw3y;

   pi  = 4.00*atan(1.0);
   enxh = enx/2;
   enyh = eny/2;
   enzh = enz/2;

   for (i=0; i<(ewaldion->nion); ++i)
   {
      sw1 = unitg[0]*ewaldion->rion1[0+3*i]
          + unitg[1]*ewaldion->rion1[1+3*i]
          + unitg[2]*ewaldion->rion1[2+3*i]+pi;
      sw2 = unitg[3]*ewaldion->rion1[0+3*i]
          + unitg[4]*ewaldion->rion1[1+3*i]
          + unitg[5]*ewaldion->rion1[2+3*i]+pi;
      sw3 = unitg[6]*ewaldion->rion1[0+3*i]
          + unitg[7]*ewaldion->rion1[1+3*i]
          + unitg[8]*ewaldion->rion1[2+3*i]+pi;

      cw1x=cos(sw1); cw1y=-sin(sw1);
      cw2x=cos(sw2); cw2y=-sin(sw2);
      cw3x=cos(sw3); cw3y=-sin(sw3);
      
      ewx1[2*i*enx] = 1.0; ewx1[2*i*enx+1] = 0.0;
      ewy1[2*i*eny] = 1.0; ewy1[2*i*eny+1] = 0.0;
      ewz1[2*i*enz] = 1.0; ewz1[2*i*enz+1] = 0.0;
      for (k=1; k<=enxh; ++k)
      {
         a = ewx1[2*(k-1 + i*enx)];
         b = ewx1[2*(k-1 + i*enx)+1];
         ewx1[2*(k + i*enx)]   = a*cw1x - b*cw1y;
         ewx1[2*(k + i*enx)+1] = a*cw1y + b*cw1x;
         ewx1[2*(enx-k + i*enx)]   =  ewx1[2*(k + i*enx)];
         ewx1[2*(enx-k + i*enx)+1] = -ewx1[2*(k + i*enx)+1];
      }
      for (k=1; k<=enyh; ++k)
      {
         a = ewy1[2*(k-1 + i*eny)];
         b = ewy1[2*(k-1 + i*eny)+1];
         ewy1[2*(k + i*eny)]   = a*cw2x - b*cw2y;
         ewy1[2*(k + i*eny)+1] = a*cw2y + b*cw2x;
         ewy1[2*(eny-k + i*eny)]   =  ewy1[2*(k + i*eny)];
         ewy1[2*(eny-k + i*eny)+1] = -ewy1[2*(k + i*eny)+1];
      }
      for (k=1; k<=enzh; ++k)
      {
         a = ewz1[2*(k-1 + i*enz)];
         b = ewz1[2*(k-1 + i*enz)+1];
         ewz1[2*(k + i*enz)]   = a*cw3x - b*cw3y;
         ewz1[2*(k + i*enz)+1] = a*cw3y + b*cw3x;
         ewz1[2*(enz-k + i*enz)]   = ewz1[2*(k + i*enz)];
         ewz1[2*(enz-k + i*enz)+1] = -ewz1[2*(k + i*enz)+1];
      }

      ewx1[2*(enxh+i*enx)] = 0.0; ewx1[2*(enxh+i*enx)+1] = 0.0;
      ewy1[2*(enyh+i*eny)] = 0.0; ewy1[2*(enyh+i*eny)+1] = 0.0;
      ewz1[2*(enzh+i*enz)] = 0.0; ewz1[2*(enzh+i*enz)+1] = 0.0;
      

   }
}

void ewald_strfac_add_sub(const int npack,
                    const int indxi[],
                    const int indxj[],
                    const int indxk[],
                    const double exi[],
                    const double exj[],
                    const double exk[],
                    const double alpha,
                    double strx[])
{
   int i;
   double ai,aj,ak,c,d;
   double bi,bj,bk;
   for (i=0; i<npack; ++i)
   {
      ai = exi[2*indxi[i]]; bi = exi[2*indxi[i]+1];
      aj = exj[2*indxj[i]]; bj = exj[2*indxj[i]+1];
      ak = exk[2*indxk[i]]; bk = exk[2*indxk[i]+1];
      c  = aj*ak - bj*bk;
      d  = aj*bk + ak*bj;
      strx[2*i]   += alpha*(ai*c - bi*d);
      strx[2*i+1] += alpha*(ai*d + bi*c);
   }
}

/*********************************
 *                               *
 *          Ewald::energy        *
 *                               *
 *********************************/
double Ewald::energy()
{
   int i,j,k,l,nion,tnp,tid,dutask;
   double x,y,z,dx,dy,dz,zz,r,w;
   double etmp1,etmp2,eall;

   tnp = ewaldparall->np();
   tid = ewaldparall->taskid();
   nion = ewaldion->nion;

   for (k=0; k<(2*enpack); ++k) ss[k] = 0.0;

   for (i=0; i<nion; ++i)
   {
      ewald_strfac_add_sub(enpack,i_indx,j_indx,k_indx,
                     &ewx1[2*i*enx],
                     &ewy1[2*i*eny],
                     &ewz1[2*i*enz],
                     zv[ewaldion->katm[i]],
                     ss);
   }
   etmp1 = 0.0;
   for (k=0; k<(enpack); ++k)
   {
      x = ss[2*k]*ss[2*k]; 
      y = ss[2*k+1]*ss[2*k+1];
      etmp1 += (x+y)*vg[k];
   }
   if (tnp>1) etmp1 = ewaldparall->SumAll(0,etmp1);
   etmp1 = etmp1/lattice_omega() + cewald;

   dutask = 0;
   etmp2  = 0.0;
   for (i=0;   i<(nion-1); ++i)
   for (j=i+1; j<nion;     ++j)
   {
      if (dutask==tid)
      {
         dx = ewaldion->rion1[3*i]   - ewaldion->rion1[3*j];
         dy = ewaldion->rion1[3*i+1] - ewaldion->rion1[3*j+1];
         dz = ewaldion->rion1[3*i+2] - ewaldion->rion1[3*j+2];
         zz = zv[ewaldion->katm[i]]*zv[ewaldion->katm[j]];
         for (l=0; l<enshl3d; ++l)
         {
            x = rcell[l]             + dx;
            y = rcell[l +   enshl3d] + dy;
            z = rcell[l + 2*enshl3d] + dz;
            r = sqrt(x*x + y*y + z*z);
            w = r/ercut;
            etmp2 += zz*erfc(w)/r;
         }
      }
      dutask =(dutask+1)%tnp;
   }
   if (tnp>1) etmp2 = ewaldparall->SumAll(0,etmp2);

   eall = etmp1+etmp2;
   return eall;
}

void ewald_strfac_sub(const int npack,
                const int indxi[],
                const int indxj[],
                const int indxk[],
                const double exi[],
                const double exj[],
                const double exk[],
                double strx[])
{
   int i;
   double ai,aj,ak,c,d;
   double bi,bj,bk;
   for (i=0; i<npack; ++i)
   {
      ai = exi[2*indxi[i]]; bi = exi[2*indxi[i]+1];
      aj = exj[2*indxj[i]]; bj = exj[2*indxj[i]+1];
      ak = exk[2*indxk[i]]; bk = exk[2*indxk[i]+1];
      c  = aj*ak - bj*bk;
      d  = aj*bk + ak*bj;
      strx[2*i]   = (ai*c - bi*d);
      strx[2*i+1] = (ai*d + bi*c);
   }
}

void ewald_f_tmp3_sub(const int n, 
                      const double e[],
                      const double s[],
                      const double v[],
                      double t[])
{
   for (int i=0; i<n; ++i) 
   {
      t[i] = v[i]*(e[2*i]*s[2*i+1] - e[2*i+1]*s[2*i]);
   }
}
double ewald_f_ddot_sub(const int n, const double a[], const double b[])
{
   double sum = 0.0;
   for (int i=0; i<n; ++i) sum += a[i]*b[i];

   return sum;
}



/*********************************
 *                               *
 *          Ewald::force         *
 *                               *
 *********************************/
void Ewald::force(double *fion)
{
   int i,j,k,l,nion,tnp,tid,dutask;
   double x,y,z,dx,dy,dz,zz,r,w,zi,f;
   double scal2,sw1,sw2,sw3;
   double cerfc=1.128379167;

   scal2 = 1.0/lattice_omega();
   tnp = ewaldparall->np();
   tid = ewaldparall->taskid();
   nion = ewaldion->nion;


   for (k=0; k<(2*enpack); ++k) ss[k] = 0.0;
   for (i=0; i<nion; ++i)
   {
      ewald_strfac_add_sub(enpack,i_indx,j_indx,k_indx,
                     &ewx1[2*i*enx],
                     &ewy1[2*i*eny],
                     &ewz1[2*i*enz],
                     zv[ewaldion->katm[i]],
                     ss);
   }
   for (i=0; i<nion; ++i)
   {
      zi =  zv[ewaldion->katm[i]];
      ewald_strfac_sub(enpack,i_indx,j_indx,k_indx,
                     &ewx1[2*i*enx],
                     &ewy1[2*i*eny],
                     &ewz1[2*i*enz],
                     exi);
      ewald_f_tmp3_sub(enpack,exi,ss,vg,tmp3);

      ftmp[3*i]   = 2.0*scal2*zi*ewald_f_ddot_sub(enpack,eG,tmp3);
      ftmp[3*i+1] = 2.0*scal2*zi*ewald_f_ddot_sub(enpack,&eG[enpack],tmp3);
      ftmp[3*i+2] = 2.0*scal2*zi*ewald_f_ddot_sub(enpack,&eG[2*enpack],tmp3);
   }

   dutask = 0;
   for (i=0; i<(nion-1); ++i) 
   for (j=i+1; j<nion; ++j)
   {
      if (dutask==tid)
      {
         dx = ewaldion->rion1[3*i]   - ewaldion->rion1[3*j];
         dy = ewaldion->rion1[3*i+1] - ewaldion->rion1[3*j+1];
         dz = ewaldion->rion1[3*i+2] - ewaldion->rion1[3*j+2];
         zz = zv[ewaldion->katm[i]]*zv[ewaldion->katm[j]];
         sw1 = 0.0;
         sw2 = 0.0;
         sw3 = 0.0;
         for (l=0; l<enshl3d; ++l)
         {
            x = rcell[l]             + dx;
            y = rcell[l +   enshl3d] + dy;
            z = rcell[l + 2*enshl3d] + dz;
            r = sqrt(x*x + y*y + z*z);
            w = r/ercut;
            f = zz*(erfc(w)+cerfc*w*exp(-w*w))/(r*r*r);
            sw1 += x*f;
            sw2 += y*f;
            sw3 += z*f;
         }
         ftmp[3*i]   += sw1;
         ftmp[3*i+1] += sw2;
         ftmp[3*i+2] += sw3;
         ftmp[3*j]   -= sw1;
         ftmp[3*j+1] -= sw2;
         ftmp[3*j+2] -= sw3;
      }
      dutask =(dutask+1)%tnp;
   }
   if (tnp>1) ewaldparall->Vector_SumAll(0,3*nion,ftmp);
   for (i=0; i<3*nion; ++i) fion[i] += ftmp[i];
}

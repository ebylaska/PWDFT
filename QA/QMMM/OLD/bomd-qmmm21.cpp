#include        <iomanip>
#include        <cstring>
#include 	<cmath>
#include 	<iostream>
#include 	<fstream>
#include 	<cstdio>
#include 	<string>
#include	<unistd.h>
#include        <sys/stat.h>
#include	<string>
#include	<vector>
#include	<set>
#include	<map>
#include 	<algorithm>

#include "mpi.h"

#include "parsestring.hpp"

using namespace std;

extern int  c_lammps_pspw_qmmm_minimizer_filename(MPI_Comm,double*,double*,double*,double*,double*,bool,bool,const char*);
extern void c_lammps_pspw_input_filename(MPI_Comm,const char*,const char*);

#define	MASTER 		0
#define ANGTOBOHR	1.88972687777


#define Efmt(w,p) std::right << std::setw(w) << std::setprecision(p)  << std::scientific
#define Ffmt(w,p) std::right << std::setw(w) << std::setprecision(p)  << std::fixed
#define Ifmt(w)   std::right << std::setw(w) 



#define E124    std::right << std::setw(12) << std::setprecision(4)  << std::scientific
#define E156    std::right << std::setw(15) << std::setprecision(6)  << std::scientific
#define E1910   std::right << std::setw(19) << std::setprecision(10) << std::scientific
#define F206    std::right << std::setw(20) << std::setprecision(6)  << std::fixed
#define F105    std::right << std::setw(10) << std::setprecision(5)  << std::fixed

#define hxyzstream(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z) E124 << (a1x) << E124 << (a1y) << E124 << (a1z) << E124 << (a2x) <<  E124 << (a2y) << E124 << (a2z) << E124 << (a3x) << E124 << (a3y) << E124 << (a3z) 

#define xyzstream(S,X,Y,Z)     std::left << std::setw(3) << (S) << E124 << (X) << E124 << (Y) << E124 << (Z) 

void printxyz(std::ofstream *xyz, const int nion, const std::string symbol[], const double *unita, const double *rion) 
{
   double AACONV = 0.529177;

   double a1x = unita[0]*AACONV;
   double a1y = unita[1]*AACONV;
   double a1z = unita[2]*AACONV;
   double a2x = unita[3]*AACONV;
   double a2y = unita[4]*AACONV;
   double a2z = unita[5]*AACONV;
   double a3x = unita[6]*AACONV;
   double a3y = unita[7]*AACONV;
   double a3z = unita[8]*AACONV;

   *xyz << nion << std::endl << hxyzstream(a1x,a1y,a1z,a2x,a2y,a2z,a3x,a3y,a3z) << std::endl;

   for (auto ii=0; ii<nion; ++ii)
      *xyz << xyzstream(symbol[ii],rion[3*ii]*AACONV,rion[3*ii+1]*AACONV,rion[3*ii+2]*AACONV) << std::endl;
}

void printemotion(std::ofstream *emotion, const double current_time,
                  const double Etotal, const double E, 
                  const double KE, const double Eqm, const double Ecoul, const double ELJ)
{
   *emotion << E1910 << current_time
         << E1910 << Etotal
         << E1910 << E
         << E1910 << KE
         << E1910 << Eqm 
         << E1910 << Ecoul 
         << E1910 << ELJ
         << std::endl;
}

double symboltomass(std::string symbol)
{
   double mass = 1.008;
   if (mystring_contains(symbol,"H")) mass = 1.008;
   if (mystring_contains(symbol,"B")) mass = 11.00931;
   if (mystring_contains(symbol,"C")) mass = 12.0;
   if (mystring_contains(symbol,"N")) mass = 14.00307;
   if (mystring_contains(symbol,"O")) mass = 15.99491;
   if (mystring_contains(symbol,"F")) mass = 18.9984;
   if (mystring_contains(symbol,"P")) mass = 30.97376;
   if (mystring_contains(symbol,"S")) mass = 31.97207;
   if (mystring_contains(symbol,"K")) mass = 38.96371;

   if (mystring_contains(symbol,"He")) mass = 4.0026;
   if (mystring_contains(symbol,"Li")) mass = 7.016;
   if (mystring_contains(symbol,"Be")) mass = 9.01218;
   if (mystring_contains(symbol,"Ne")) mass = 19.99244;

   if (mystring_contains(symbol,"Na")) mass = 22.9898;
   if (mystring_contains(symbol,"Mg")) mass = 23.98504;
   if (mystring_contains(symbol,"Al")) mass = 26.98154;
   if (mystring_contains(symbol,"Si")) mass = 27.97693;
   if (mystring_contains(symbol,"Cl")) mass = 34.96885;
   return mass;
}


double symboltoepsilon(std::string symbol)
{
   double epsilon = 0.044;
   if (mystring_contains(symbol,"H")) epsilon = 0.044;
   if (mystring_contains(symbol,"B")) epsilon = 0.180;
   if (mystring_contains(symbol,"C")) epsilon = 0.105;
   if (mystring_contains(symbol,"N")) epsilon = 0.069;
   if (mystring_contains(symbol,"O")) epsilon = 0.060;
   if (mystring_contains(symbol,"F")) epsilon = 0.050;
   if (mystring_contains(symbol,"P")) epsilon = 0.305;
   if (mystring_contains(symbol,"S")) epsilon = 0.274;
   if (mystring_contains(symbol,"K")) epsilon = 0.035;

   if (mystring_contains(symbol,"He")) epsilon = 0.056;
   if (mystring_contains(symbol,"Li")) epsilon = 0.025;
   if (mystring_contains(symbol,"Be")) epsilon = 0.085;
   if (mystring_contains(symbol,"Ne")) epsilon = 0.042;

   if (mystring_contains(symbol,"Na")) epsilon = 0.030;
   if (mystring_contains(symbol,"Mg")) epsilon = 0.111;
   if (mystring_contains(symbol,"Al")) epsilon = 0.505;
   if (mystring_contains(symbol,"Si")) epsilon = 0.402;
   if (mystring_contains(symbol,"Cl")) epsilon = 0.227;
   return epsilon;
}

double symboltosigma(std::string symbol)
{
   double sigma = 2.886;
   if (mystring_contains(symbol,"H")) sigma = 2.886;
   if (mystring_contains(symbol,"B")) sigma = 4.083;
   if (mystring_contains(symbol,"C")) sigma = 3.851;
   if (mystring_contains(symbol,"N")) sigma = 3.66;
   if (mystring_contains(symbol,"O")) sigma = 3.50;
   if (mystring_contains(symbol,"F")) sigma = 3.364;
   if (mystring_contains(symbol,"P")) sigma = 4.147;
   if (mystring_contains(symbol,"S")) sigma = 4.035;
   if (mystring_contains(symbol,"K")) sigma = 3.812;

   if (mystring_contains(symbol,"He")) sigma = 2.362;
   if (mystring_contains(symbol,"Li")) sigma = 2.451;
   if (mystring_contains(symbol,"Be")) sigma = 2.745;
   if (mystring_contains(symbol,"Ne")) sigma = 3.243;

   if (mystring_contains(symbol,"Na")) sigma = 2.983;
   if (mystring_contains(symbol,"Mg")) sigma = 3.021;
   if (mystring_contains(symbol,"Al")) sigma = 4.499;
   if (mystring_contains(symbol,"Si")) sigma = 4.295;
   if (mystring_contains(symbol,"Cl")) sigma = 3.947;
   return sigma;
}

/**************************************************
 *                                                *
 *                spring_bond_frag                *
 *                                                *
 **************************************************/
double spring_bond_frag(const int nbond, const int indx[], const double Kr0[], const double rion[])
{
   double E = 0.0;
   for (auto b=0; b<nbond; ++b)
   {
      int ii = indx[0+2*b];
      int jj = indx[1+2*b];

      double x = rion[3*ii]   - rion[3*jj];
      double y = rion[3*ii+1] - rion[3*jj+1];
      double z = rion[3*ii+2] - rion[3*jj+2];
      double r = std::sqrt(x*x + y*y + z*z);
      double dr = r - Kr0[1+2*b];
      E += Kr0[2*b]*dr*dr;

   }
   return E;
}


/**************************************************
 *                                                *
 *                spring_bond_force_frag          *
 *                                                *
 **************************************************/
void spring_bond_force_frag(const int nbond, const int indx[], const double Kr0[], const double rion[], double fion[])
{
   for (auto b=0; b<nbond; ++b)
   {

      int ii = indx[0+2*b];
      int jj = indx[1+2*b];

      double x = rion[3*ii]   - rion[3*jj];
      double y = rion[3*ii+1] - rion[3*jj+1];
      double z = rion[3*ii+2] - rion[3*jj+2];
      double r = std::sqrt(x*x + y*y + z*z);
      double dr = r - Kr0[1+2*b];
      double dE = 2.0*Kr0[2*b]*dr/r;

      fion[3*ii]   -= x*dE;
      fion[3*ii+1] -= y*dE;
      fion[3*ii+2] -= z*dE;

      fion[3*jj]   += x*dE;
      fion[3*jj+1] += y*dE;
      fion[3*jj+2] += z*dE;
   }
}


/**************************************************
 *                                                *
 *               spring_angle_frag                *
 *                                                *
 **************************************************/
double spring_angle_frag(const int nangle, const int indx[], const double Kr0[], const double rion[])
{
   double E = 0.0;
   for (auto a=0; a<nangle; ++a)
   {
      int ii = indx[0+3*a];
      int jj = indx[1+3*a];
      int kk = indx[2+3*a];

      double x1 = rion[3*ii]   - rion[3*jj];
      double y1 = rion[3*ii+1] - rion[3*jj+1];
      double z1 = rion[3*ii+2] - rion[3*jj+2];
      double r1 = std::sqrt(x1*x1 + y1*y1 + z1*z1);

      double x2 = rion[3*kk]   - rion[3*jj];
      double y2 = rion[3*kk+1] - rion[3*jj+1];
      double z2 = rion[3*kk+2] - rion[3*jj+2];
      double r2 = std::sqrt(x2*x2 + y2*y2 + z2*z2);

      double denom = r1*r2;
      if (denom>1.0e-11)
      {
         double ctheta = (x1*x2+y1*y2+z1*z2)/(denom);
         if (ctheta > 1.0)  ctheta = 1.0;
         if (ctheta < -1.0) ctheta = -1.0;
         double q = std::acos(ctheta)-Kr0[1+2*a];
         E += Kr0[2*a]*q*q;
      }
   }
   return E;
}


/**************************************************
 *                                                *
 *              spring_angle_force_frag           *
 *                                                *
 **************************************************/
void spring_angle_force_frag(const int nangle, const int indx[], const double Kr0[], const double rion[], double fion[])
{
   for (auto a=0; a<nangle; ++a)
   {
      int ii = indx[0+3*a];
      int jj = indx[1+3*a];
      int kk = indx[2+3*a];

      double x1 = rion[3*ii]   - rion[3*jj];
      double y1 = rion[3*ii+1] - rion[3*jj+1];
      double z1 = rion[3*ii+2] - rion[3*jj+2];
      double r1sq = x1*x1 + y1*y1 + z1*z1;
      double r1 = std::sqrt(r1sq);

      double x2 = rion[3*kk]   - rion[3*jj];
      double y2 = rion[3*kk+1] - rion[3*jj+1];
      double z2 = rion[3*kk+2] - rion[3*jj+2];
      double r2sq = x2*x2 + y2*y2 + z2*z2;
      double r2 = std::sqrt(r2sq);

      double denom = r1*r2;
      if (denom>1.0e-11)
      {
         double ctheta = (x1*x2+y1*y2+z1*z2)/(denom);
         if (ctheta > 1.0)  ctheta = 1.0;
         if (ctheta < -1.0) ctheta = -1.0;
         double stheta = std::sqrt(1.0 - ctheta*ctheta);
         if (stheta < 0.001) stheta = 0.001;
         stheta = 1.0/stheta;

         double q  = std::acos(ctheta)-Kr0[1+2*a];
         double tk = Kr0[2*a]*q;

         double  aa     = 2.0*tk*stheta;
         double  a11    =  aa*ctheta/r1sq;
         double  a12    = -aa/(denom);
         double  a22    =  aa*ctheta/r2sq;

         double  vx1 = a11*x1 + a12*x2;
         double  vx2 = a22*x2 + a12*x1;

         double  vy1 = a11*y1 + a12*y2;
         double  vy2 = a22*y2 + a12*y1;

         double  vz1 = a11*z1 + a12*z2;
         double  vz2 = a22*z2 + a12*z1;

         fion[3*ii]   -= vx1;
         fion[3*ii+1] -= vy1;
         fion[3*ii+2] -= vz1;

         fion[3*jj]   += vx1 + vx2;
         fion[3*jj+1] += vy1 + vy2;
         fion[3*jj+2] += vz1 + vz2;

         fion[3*kk]   -= vx2;
         fion[3*kk+1] -= vy2;
         fion[3*kk+2] -= vz2;
      }
   }
}


/**************************************************
 *                                                *
 *                spring_Energy                   *
 *                                                *
 **************************************************/
double spring_Energy(const int nfrag, const int indxfrag_start[], const int kfrag[], 
                     const int nbond_frag[],  const int bond_start_frag[],  const int  bond_indx[], const double Krb[],
                     const int nangle_frag[], const int angle_start_frag[], const int angle_indx[], const double Kra[],
                     const double rion[])
{
   double   E = 0.0;
   for (auto w1=0; w1<nfrag; ++w1)
   {
      int ks1 = indxfrag_start[w1];
      int ia  = kfrag[w1];
      int nbs = nbond_frag[ia];
      int nas = nangle_frag[ia];
      if (nbs>0)
      {
         int bs1 = bond_start_frag[ia];
         E += spring_bond_frag(nbs,&bond_indx[2*bs1],&Krb[2*bs1],&rion[3*ks1]);
      }
      if (nas>0)
      {
         int ba1 = angle_start_frag[ia];
         E += spring_angle_frag(nas,&angle_indx[3*ba1],&Kra[2*ba1],&rion[3*ks1]);
      }
   }
   return E;
}

/**************************************************
 *                                                *
 *                spring_Force                    *
 *                                                *
 **************************************************/
void spring_Force(const int nfrag, const int indxfrag_start[], const int kfrag[], 
                  const int nbond_frag[],  const int bond_start_frag[],  const int  bond_indx[], const double Krb[],
                  const int nangle_frag[], const int angle_start_frag[], const int angle_indx[], const double Kra[],
                  const double rion[], double fion[])
{
   for (auto w1=0; w1<nfrag; ++w1)
   {
      int ks1 = indxfrag_start[w1];
      int ia  = kfrag[w1];
      int nbs = nbond_frag[ia];
      int nas = nangle_frag[ia];
      if (nbs>0)
      {
         int bs1 = bond_start_frag[ia];
         spring_bond_force_frag(nbs,&bond_indx[2*bs1],&Krb[2*bs1],&rion[3*ks1],&fion[3*ks1]);
      }
      if (nas>0)
      {
         int ba1 = angle_start_frag[ia];
         spring_angle_force_frag(nas,&angle_indx[3*ba1],&Kra[2*ba1],&rion[3*ks1],&fion[3*ks1]);
      }
   }
}







double LJ_energy(const double epsilon12, const double sigma12,
                 const double r1[], const double r2[])
{

   double x = r2[0] - r1[0];
   double y = r2[1] - r1[1];
   double z = r2[2] - r1[2];
   double r = std::sqrt(x*x + y*y + z*z);
   double u = (sigma12/r);
   double u6  = u*u*u*u*u*u;
   double u12 = u6*u6;

   return (4.0*epsilon12*(u12-u6));
}

void LJ_force(const double epsilon12, const double sigma12,
                const double r1[], double f1[],
                const double r2[], double f2[])
{
   double x = r2[0] - r1[0];
   double y = r2[1] - r1[1];
   double z = r2[2] - r1[2];
   double r = std::sqrt(x*x + y*y + z*z);
   double u = (sigma12/r);
   double u6  = u*u*u*u*u*u;
   double u12 = u6*u6;

   double dVLJ = -(4.00*epsilon12/r)*(12.0*u12-6.0*u6);

   f1[0] += (x/r)*dVLJ;
   f1[1] += (y/r)*dVLJ;
   f1[2] += (z/r)*dVLJ;

   f2[0] -= (x/r)*dVLJ;
   f2[1] -= (y/r)*dVLJ;
   f2[2] -= (z/r)*dVLJ;
}

double QMMM_LJ_energy(const int nion_qm, const int nion,
                      const double epsilon[], const double sigma[], const double rion[])
{
   double E = 0.0;
   // QMQM LJ == 0
   // QMMM LJ
   // qm and mm interactions
   for (auto qm=0; qm<nion_qm; ++qm)
      for (auto mm=nion_qm; mm<nion; ++mm)
      {
         double epsilon12 = std::sqrt(epsilon[qm]*epsilon[mm]);
         double sigma12 = 0.5*(sigma[qm] + sigma[mm]);
         double elj = LJ_energy(epsilon12,sigma12,&rion[3*qm],&rion[3*mm]);
         E += elj;
      }

   // mm interactions
   /*for (auto mm1=nion_qm; mm1<(nion-1); ++mm1)
      for (auto mm2=mm1+1; mm2<nion; ++mm2)
      {
         double epsilon12 = std::sqrt(epsilon[mm1]*epsilon[mm2]);
         double sigma12 = 0.5*(sigma[mm1] + sigma[mm2]);
         double elj = LJ_energy(epsilon12,sigma12,&rion[3*mm1],&rion[3*mm2]);
         E += elj;
      }
   */

   return E;
}

void QMMM_LJ_force(const int nion_qm, const int nion,
                   const double epsilon[], const double sigma[], const double rion[], double fion[])
{
   // qm and mm interactions
   for (auto qm=0; qm<nion_qm; ++qm)
      for (auto mm=nion_qm; mm<nion; ++mm) {
         double epsilon12 = std::sqrt(epsilon[qm]*epsilon[mm]);
         double sigma12 = 0.5*(sigma[qm] + sigma[mm]);
         LJ_force(epsilon12,sigma12,&rion[3*qm],&fion[3*qm],&rion[3*mm],&fion[3*mm]);
      }

   // mm interactions
   /*
   for (auto mm1=nion_qm; mm1<(nion-1); ++mm1)
      for (auto mm2=mm1+1; mm2<nion; ++mm2)
      {
         double epsilon12 = std::sqrt(epsilon[mm1]*epsilon[mm2]);
         double sigma12 = 0.5*(sigma[mm1] + sigma[mm2]);
         LJ_force(epsilon12,sigma12,&rion[3*mm1],&fion[3*mm1],&rion[3*mm2],&fion[3*mm2]);
      }
   */
}

double MMMM_LJ_energy(const int nion_qm, const int nion,
                      const double epsilon[], const double sigma[], const double rion[])
{
   double E = 0.0;

   // mm interactions
   for (auto mm1=nion_qm; mm1<(nion-1); ++mm1)
      for (auto mm2=mm1+1; mm2<nion; ++mm2)
      {
         double epsilon12 = std::sqrt(epsilon[mm1]*epsilon[mm2]);
         double sigma12 = 0.5*(sigma[mm1] + sigma[mm2]);
         double elj = LJ_energy(epsilon12,sigma12,&rion[3*mm1],&rion[3*mm2]);
         E += elj;
      }
   return E;
}

void MMMM_LJ_force(const int nion_qm, const int nion,
                   const double epsilon[], const double sigma[], const double rion[], double fion[])
{
   // mm interactions
   for (auto mm1=nion_qm; mm1<(nion-1); ++mm1)
      for (auto mm2=mm1+1; mm2<nion; ++mm2)
      {
         double epsilon12 = std::sqrt(epsilon[mm1]*epsilon[mm2]);
         double sigma12 = 0.5*(sigma[mm1] + sigma[mm2]);
         LJ_force(epsilon12,sigma12,&rion[3*mm1],&fion[3*mm1],&rion[3*mm2],&fion[3*mm2]);
      }
}





/**************************************************
 *                                                *
 *               Q_Switching                      *
 *                                                *
 **************************************************/
double Q_Switching(const double Rin, const double Rout, const double r)
{
   double s;

   if (r<=Rin)  return 0.0;
   if (r>=Rout) return 1.0;

   double c1 = (Rout-Rin);
   double c3 = c1*c1*c1;
   double c2 = 3.0*Rout-Rin;
   return(((r-Rin)*(r-Rin)) * (c2-2.0*r)/c3);
}

/**************************************************
 *                                                *
 *               Q_dSwitching                     *
 *                                                *
 **************************************************/
void Q_dSwitching(const double Rin, const double Rout, const double r,
                  double *S, double *dS)
{
   if (r<=Rin)  { *S = 0.0; *dS = 0.0; return; }
   if (r>=Rout) { *S = 1.0; *dS = 0.0; return; }
   
   double c1 = (Rout-Rin);
   double c3 = c1*c1*c1;
   double c2 = 3.0*Rout-Rin;
  
   *S  = ((r-Rin)*(r-Rin)) * (c2-2.0*r)/c3;
   *dS = 2.0*(r-Rin)*((c2-2.0*r)/c3) - 2.0*((r-Rin)*(r-Rin)) /c3;
   return;
}




double Q_Electrostatic_potential(const double r1[], const double q1,
                                 const double r2[])
{
   double x = r2[0] - r1[0];
   double y = r2[1] - r1[1];
   double z = r2[2] - r1[2];
   double r = std::sqrt(x*x + y*y + z*z);

   return (q1/r);
}

double Q_Electrostatic_self(const double r1[], const double q1,
                            const double r2[], const double q2)
{
   double x = r2[0] - r1[0];
   double y = r2[1] - r1[1];
   double z = r2[2] - r1[2];
   double r = std::sqrt(x*x + y*y + z*z);

   return (q1*q2/r);
}

void Q_Electrostatic_Force_self(const double r1[], const double q1, double f1[],
                                const double r2[], const double q2, double f2[])
{
   double x = r2[0] - r1[0];
   double y = r2[1] - r1[1];
   double z = r2[2] - r1[2];
   double rr = (x*x + y*y + z*z);
   double r  = std::sqrt(rr);

   double der = -q1*q2/rr;
   double fx = -(x/r)*der;
   double fy = -(y/r)*der;
   double fz = -(z/r)*der;

   f2[0] += fx;
   f2[1] += fy;
   f2[2] += fz;

   f1[0] -= fx;
   f1[1] -= fy;
   f1[2] -= fz;
}

/****************************************************
 *                                                  *
 *        QMMM_electrostatic_potential              *
 *                                                  *
 ****************************************************/
void QMMM_electrostatic_potential(const int nion_qm, const int nion,
                                  const double qion[], const double rion[], double uion[])
{
   for (auto qm=0; qm<nion_qm; ++qm)
   {
      uion[qm] = 0.0;
      for (auto mm=nion_qm; mm<nion; ++mm)
         uion[qm] += Q_Electrostatic_potential(&rion[3*mm],qion[mm], &rion[3*qm]);
   }
}


/****************************************************
 *                                                  *
 *        QMQM_electrostatic_energy                 *
 *                                                  *
 ****************************************************/
double QMQM_electrostatic_energy(const int nion_qm, const int nion,
                                 const double qion[], const double rion[])
{
   double E=0.0;

   // total qmqm interactions
   for (auto ii=0; ii<(nion_qm-1); ++ii)
      for (auto jj=ii+1; jj<nion_qm; ++jj)
         E += Q_Electrostatic_self(&rion[3*ii],qion[ii], 
                                   &rion[3*jj],qion[jj]);

   return E;
}


/****************************************************
 *                                                  *
 *        QMQM_electrostatic_force                  *
 *                                                  *
 ****************************************************/
void QMQM_electrostatic_force(const int nion_qm, const int nion,
                              const double qion[], const double rion[], double fion[])
{
   // total qmqm interactions
   for (auto ii=0; ii<(nion_qm-1); ++ii)
      for (auto jj=ii+1; jj<nion_qm; ++jj)
         Q_Electrostatic_Force_self(&rion[3*ii],qion[ii],&fion[3*ii],
                                    &rion[3*jj],qion[jj],&fion[3*jj]);

}


/****************************************************
 *                                                  *
 *        QMMM_electrostatic_energy                 *
 *                                                  *
 ****************************************************/
double QMMM_electrostatic_energy(const int nion_qm, const int nion,
                                 const double qion[], const double rion[])
{
   double E=0.0;

   // total qmmm interactions
   for (auto ii=0; ii<(nion_qm); ++ii)
      for (auto jj=nion_qm; jj<nion; ++jj)
         E += Q_Electrostatic_self(&rion[3*ii],qion[ii], 
                                   &rion[3*jj],qion[jj]);

   return E;
}

/****************************************************
 *                                                  *
 *        QMMM_electrostatic_force                  *
 *                                                  *
 ****************************************************/
void QMMM_electrostatic_force(const int nion_qm, const int nion,
                              const double qion[], const double rion[], double fion[])
{
   // total qmmm interactions
   for (auto ii=0; ii<(nion_qm); ++ii)
      for (auto jj=nion_qm; jj<nion; ++jj)
         Q_Electrostatic_Force_self(&rion[3*ii],qion[ii],&fion[3*ii],
                                    &rion[3*jj],qion[jj],&fion[3*jj]);

}


/****************************************************
 *                                                  *
 *                    Q_cm                          *
 *                                                  *
 ****************************************************/
void Q_cm(const int n, const int ks, const double amass[], const double rion[], 
         double rcm[])
{
   rcm[0] = 0.0; rcm[1] = 0.0; rcm[2] = 0.0;
   double m = 0.0;
   int kk = ks;
   for (auto k=0; k<n; ++k)
   {
      rcm[0] += amass[kk]*rion[3*kk];
      rcm[1] += amass[kk]*rion[3*kk+1];
      rcm[2] += amass[kk]*rion[3*kk+1];
      m += amass[kk];
      ++kk;
   }
   rcm[0] /= m;
   rcm[1] /= m;
   rcm[2] /= m;
}



/****************************************************
 *                                                  *
 *        MMMM_electrostatic_energy                 *
 *                                                  *
 ****************************************************/
double MMMM_electrostatic_energy(const int nfrag,       const int indx_frag_start[], 
                                 const int size_frag[], const int kfrag[],
                                 const double switch_Rin[], const double switch_Rout[],
                                 const double amass[],
                                 const double qion[], const double rion[])
{
   double E=0.0;
   double rw1_cm[3],rw2_cm[3];

   // total mmmm interactions
   for (auto w1=0; w1<(nfrag-1); ++ w1)
   {
      int ks1 = indx_frag_start[w1];
      int n1  = size_frag[w1];
      double Rin1  = switch_Rin[kfrag[w1]];
      double Rout1 = switch_Rin[kfrag[w1]];
      Q_cm(n1,ks1,amass,rion,rw1_cm);

      for (auto w2=w1+1; w2<nfrag; ++w2)
      {
         int ks2 = indx_frag_start[w2];
         int n2  = size_frag[w2];
         double Rin2  = switch_Rin[kfrag[w2]];
         double Rout2 = switch_Rin[kfrag[w2]];
         Q_cm(n2,ks2,amass,rion,rw2_cm);

         double Rin  = 0.5*(Rin1 +Rin2);
         double Rout = 0.5*(Rout1+Rout2);
      }
   }

/*   for (auto f1=0; f1<nfrag-1; ++f1)
      for (auto f2=f1+1; f2<nfrag; ++f2)
         E += Q_Electrostatic_self(&rion[3*ii],qion[ii], 
                                   &rion[3*jj],qion[jj]);

*/
   return E;
}





int main(int argc, char* argv[])
{
   std::ofstream *xyzfile;
   std::ofstream *emotionfile;

   //std::string nwfilename = "w2.nw";
   std::string nwfilename = "water2.nw";
   if (argc>1) nwfilename = argv[1];

   std::string debugfilename = "";

   int ierr,np,taskid;

   // Initialize MPI
   ierr = MPI_Init(&argc,&argv);
   ierr += MPI_Comm_size(MPI_COMM_WORLD,&np);
   ierr += MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

   if (taskid==MASTER) {
      debugfilename = "qmmm20.nwout";
      std::cout << "Hello world" << std::endl;
      std::cout << "np=" << np << " taskid=" << taskid << std::endl;
      std::cout << "argc=" << argc << std::endl;
      std::cout << "nwfilename=" << nwfilename << std::endl;
      std::cout << "debugfilename=" << debugfilename << std::endl;
   }

   // copy to c strings
   char cnwfilename[nwfilename.size()+1],cdebugfilename[debugfilename.size()+1];
   std::strcpy(cnwfilename,  nwfilename.c_str());
   std::strcpy(cdebugfilename,debugfilename.c_str());


   // read input deck
   int nwinput_size;
   std::string nwinput;
   //if (taskid==MASTER)
   {
      std::string line;

      // prepend print option to nwinput, options=none,off,low,medium,high,debug
      //nwinput = "\nprint off\n\n"; 
      nwinput = "";
      std::ifstream nwfile(nwfilename);
      if (nwfile.good())
      {
         while (getline(nwfile,line))
            nwinput += line + "\n";
      }
      nwfile.close();

      nwinput_size = nwinput.size();
   }

   // set up geometry sizes and indexes: nion, nkatm, katm[nion] 
   std::string geomblock;
   if (mystring_contains(mystring_lowercase(nwinput),"geometry"))
   {
      geomblock  = mystring_rtrim(mystring_ltrim(mystring_split(mystring_split(nwinput,"nocenter")[1],"end")[0]));
   }
   std::vector<std::string> geomlines = mystring_split(geomblock,"\n");
   int nion  = geomlines.size();
   int nkatm = 0;
   int katm[nion];
   {  std::set<std::string> anameset;
      for (auto & line: geomlines)
      {
        std::vector<std::string> ss = mystring_split0(line);
        anameset.emplace(ss[0]);
      }
      nkatm = anameset.size();
   }
   std::string aname[nkatm];

   {  int ia=0;
      int ii=0;
      int n = sizeof(aname)/sizeof(aname[0]);
      std::set<std::string> anameset;
      for (auto & line: geomlines)
      {
         std::vector<std::string> ss = mystring_split0(line);
         if (anameset.count(ss[0])==0)
         {
            anameset.emplace(ss[0]);
            aname[ia] = ss[0];
            ++ia;
         }
         auto itr = find(aname, aname + n, ss[0]);
         katm[ii] = std::distance(aname,itr);
         ++ii;
      }
   }
   if (taskid==MASTER) {
      std::cout << "katm= ";
      for (auto ii=0; ii<nion; ++ii)
         std::cout << " " << katm[ii];
      std::cout << std::endl;
  

      std::cout << std::endl;
      std::cout << "aname size=" << nkatm << std::endl;
      for (auto ia=0; ia<nkatm; ++ia)
          std::cout << "aname =" << ia << " " << aname[ia] << std::endl;
      std::cout << std::endl;
   }



   //set up qmmm fragments and sizes
   std::vector<std::string> fragments;
   double qmmm_lmbda = 1.0;
   int    nkfrag = 0;
   int    nfrag  = 0;
   int    nbond  = 0;
   int    nangle = 0;

   if  (mystring_contains(mystring_lowercase(nwinput),"fragment"))
   {
      fragments = mystring_split(nwinput,"fragment");
      fragments.erase(fragments.begin());
      nkfrag = fragments.size();
      for (auto & frag: fragments)
      {
         int nfrag0=0;
         frag = mystring_rtrim(mystring_ltrim(mystring_split(frag,"end")[0]));

         if (mystring_contains(mystring_lowercase(frag),"index_start"))
            nfrag0 = mystring_split0(mystring_trim(mystring_split(mystring_split(frag,"index_start")[1],"\n")[0])).size();
         
         if (mystring_contains(mystring_lowercase(frag),"bond_spring"))
            nbond += (mystring_split(frag,"bond_spring").size() - 1);

         if (mystring_contains(mystring_lowercase(frag),"angle_spring"))
            nangle += (mystring_split(frag,"angle_spring").size() - 1);

         nfrag += nfrag0;
      }
   }
   if (taskid==MASTER) {
      std::cout << "nkfrag=" << nkfrag << std::endl;
      std::cout << "nfrag =" << nfrag  << std::endl;
      std::cout << "nbond =" << nbond  << std::endl;
      std::cout << "nangle =" << nangle  << std::endl;
      std::cout << std::endl;
   }

   //set up qmmm lj stuff
   std::map<std::string,std::vector<double>> lj_qm_data,lj_mm_data;

   if  (mystring_contains(mystring_lowercase(nwinput),"lj_qm_parameters"))
   {
      std::vector<std::string> ss;
      std::vector<std::string> lj_fragments = mystring_split(nwinput,"lj_qm_parameters");
      lj_fragments.erase(lj_fragments.begin());
      for (auto & lj_frag: lj_fragments)
      {
         ss = mystring_split0(lj_frag);
         lj_qm_data[ss[0]].push_back(std::stod(ss[1]));
         lj_qm_data[ss[0]].push_back(std::stod(ss[2]));
      }
   }
   if  (mystring_contains(mystring_lowercase(nwinput),"lj_mm_parameters"))
   {
      std::vector<std::string> ss;
      std::vector<std::string> lj_fragments = mystring_split(nwinput,"lj_mm_parameters");
      lj_fragments.erase(lj_fragments.begin());
      for (auto & lj_frag: lj_fragments)
      {
         ss = mystring_split0(lj_frag);
         lj_mm_data[ss[0]].push_back(std::stod(ss[1]));
         lj_mm_data[ss[0]].push_back(std::stod(ss[2]));
      }
   }


   // Initialize lammps_pspw interface 
   c_lammps_pspw_input_filename(MPI_COMM_WORLD,cnwfilename,NULL);

   if (taskid==MASTER) std::cout << "out lammps inpute " << nangle  << std::endl;


   double unita[9] = {26.0,  0.0,  0.0,
                       0.0, 26.0,  0.0,
                       0.0,  0.0, 26.0};
   double dt = 5.0;
   double h  = 1.0/(2.0*dt);
   double rion0[3*nion],rion1[3*nion],rion2[3*nion],fion[3*nion]; 
   double mass[nion],dti[nion],KE,uion[nion],qion[nion];
   double Espring,ELJ,Eqm,Eqq;

   double epsilon[nion],sigma[nion];

   // switching parameters
   double switch_Rin[nkfrag]; //  = { (2.0160/0.529177) };
   double switch_Rout[nkfrag]; // = { (3.1287/0.529177) };

   int size_frag[nkfrag];
   int nfrag_frag[nkfrag];
   int nbond_frag[nkfrag];
   int nangle_frag[nkfrag];
   int bond_start_frag[nkfrag];
   int angle_start_frag[nkfrag];

   int indxfrag_start[nfrag];
   int kfrag[nfrag];

   // bond and angle spring parameters
   int    bond_indx[2*nbond],angle_indx[3*nangle];
   double Krb[2*nbond];
   double Kra[2*nangle];
   double pi = 4.0*std::atan(1.0);

   int bond_start_sum  = 0;
   int angle_start_sum = 0;
   int ntf = 0;
   int ia  = 0;
   for (auto & frag: fragments)
   {
      std::cout << "start" << std::endl;
      frag = mystring_rtrim(mystring_ltrim(mystring_split(frag,"end")[0]));

      int tsize = std::stoi(mystring_split0(mystring_split(frag,"size")[1])[0]);
      std::cout << "tsize=" << tsize << std::endl;

      size_frag[ia] = tsize;

      if (mystring_contains(mystring_lowercase(frag),"index_start"))
      {
         std::vector<std::string> ss = mystring_split0(mystring_trim(mystring_split(mystring_split(frag,"index_start")[1],"\n")[0]));
         nfrag_frag[ia] = ss.size();
         for (auto i=0; i<ss.size(); ++i)
         {
            indxfrag_start[ntf] = std::stoi(ss[i])-1;
            kfrag[ntf] = ia;
            ++ntf;
         }
      }

      if (mystring_contains(mystring_lowercase(frag),"bond_spring"))
      {
         std::string line;
         std::vector<std::string> ss;
         nbond_frag[ia]      = (mystring_split(frag,"bond_spring").size() - 1);
         bond_start_frag[ia] = bond_start_sum;
         bond_start_sum     += nbond_frag[ia];
         for (auto b=0; b<nbond_frag[ia]; ++b)
         {
            int b1 = b + bond_start_frag[ia];
            line = mystring_split(frag,"bond_spring")[b+1];
            ss   = mystring_split0(line);
            int i = std::stoi(ss[0]);
            int j = std::stoi(ss[1]);
            double k = std::stod(ss[2]);
            double r = std::stod(ss[3]);

            bond_indx[2*b1  ] = i-1;
            bond_indx[2*b1+1] = j-1;
            Krb[2*b1  ] = k/27.2116/23.06/ANGTOBOHR/ANGTOBOHR;
            Krb[2*b1+1] = r*ANGTOBOHR;
         }  
      }

      if (mystring_contains(mystring_lowercase(frag),"angle_spring"))
      {
         std::string line;
         std::vector<std::string> ss;
         nangle_frag[ia]      = (mystring_split(frag,"angle_spring").size() - 1);
         angle_start_frag[ia] = angle_start_sum;
         angle_start_sum     += nangle_frag[ia];
         for (auto a=0; a<nangle_frag[ia]; ++a)
         {
            int a1 = a + angle_start_frag[ia];
            line = mystring_split(frag,"angle_spring")[a+1];
            ss   = mystring_split0(line);
            int i = std::stoi(ss[0]);
            int j = std::stoi(ss[1]);
            int k = std::stoi(ss[2]);
            double ks    = std::stod(ss[3]);
            double theta = std::stod(ss[4]);
            std::cout << "i=" << i << " j=" << j << " k=" << k 
                      << " ks =" << ks << " theta=" << theta << std::endl;
            angle_indx[3*a1  ] = i-1;
            angle_indx[3*a1+1] = j-1;
            angle_indx[3*a1+2] = k-1;
            Kra[2*a1  ] = ks/27.2116/23.06;
            Kra[2*a1+1] = theta*pi/180.0;
         }
      }
      ++ia;
   }



   std::string symbol[nion];
   memset(rion0,0,3*nion*sizeof(double));
   memset(rion1,0,3*nion*sizeof(double));
   memset(rion2,0,3*nion*sizeof(double));
   memset(fion,0, 3*nion*sizeof(double));
   memset(uion,0,   nion*sizeof(double));
   memset(qion,0,   nion*sizeof(double));


   // get the qm atoms
   int nion_qm = 0;
   for (auto ii=0; ii<nion; ++ii)
   {
      if (!mystring_contains(mystring_lowercase(geomlines[ii]),"#"))
      {
         std::vector<std::string> ss = mystring_split0(geomlines[ii]);
         symbol[nion_qm]  = ss[0];
         if (lj_qm_data.count(symbol[nion_qm]) > 0)
         {
            epsilon[nion_qm] = lj_qm_data[symbol[nion_qm]][0]/23.06/27.2116;
            sigma[nion_qm]   = lj_qm_data[symbol[nion_qm]][1]*ANGTOBOHR;
         }
         else
         {
            epsilon[nion_qm] = symboltoepsilon(symbol[nion_qm])/23.06/27.2116;
            sigma[nion_qm]   = symboltosigma(symbol[nion_qm])*ANGTOBOHR;
         }
         mass[nion_qm]    = symboltomass(symbol[nion_qm])*1822.89;
         rion1[3*nion_qm]   = std::stod(ss[1])*ANGTOBOHR;
         rion1[3*nion_qm+1] = std::stod(ss[2])*ANGTOBOHR;
         rion1[3*nion_qm+2] = std::stod(ss[3])*ANGTOBOHR;
         dti[nion_qm] = dt*dt/mass[nion_qm];
         ++nion_qm;
      }
   }

   // get the mm atoms
   int nion_mm = 0;
   for (auto ii=0; ii<nion; ++ii)
   {
      if (mystring_contains(mystring_lowercase(geomlines[ii]),"#"))
      {
         std::vector<std::string> ss = mystring_split0(geomlines[ii]);
         std::string str = ss[0]; 
         str.erase(std::remove(str.begin(),str.end(),'#'),str.end());
         symbol[(nion_qm+nion_mm)]  = str;
         if (lj_mm_data.count(str) > 0)
         {
            epsilon[(nion_qm+nion_mm)] = lj_mm_data[symbol[(nion_qm+nion_mm)]][0]/23.06/27.2116;
            sigma[(nion_qm+nion_mm)]   = lj_mm_data[symbol[(nion_qm+nion_mm)]][1]*ANGTOBOHR;
         }
         else 
         {
            epsilon[(nion_qm+nion_mm)] = 0.0;
            sigma[(nion_qm+nion_mm)]   = 1.0;
         }
         mass[(nion_qm+nion_mm)]      = symboltomass(symbol[(nion_qm+nion_mm)])*1822.89;
         rion1[3*(nion_qm+nion_mm)]   = std::stod(ss[1])*ANGTOBOHR;
         rion1[3*(nion_qm+nion_mm)+1] = std::stod(ss[2])*ANGTOBOHR;
         rion1[3*(nion_qm+nion_mm)+2] = std::stod(ss[3])*ANGTOBOHR;
         dti[(nion_qm+nion_mm)]       = dt*dt/mass[(nion_qm+nion_mm)];
         qion[(nion_qm+nion_mm)]      = std::stod(ss[5]);
         ++nion_mm;
      }
   }

   if (taskid==MASTER) std::cout << "nion =" << nion << std::endl;
   if (taskid==MASTER) std::cout << "nion_mm =" << nion_mm << std::endl;
   if (taskid==MASTER) std::cout << "nion_qm =" << nion_qm << std::endl << std::endl;

   // kinetic energy
   Espring = spring_Energy(nfrag,indxfrag_start,kfrag, 
                           nbond_frag,bond_start_frag,bond_indx,Krb,
                           nangle_frag,angle_start_frag,angle_indx,Kra,rion1);


   std::cout << "Espring=" << Espring << std::endl;

   // kinetic energy
   KE = 0.0;
   for (auto ii=0; ii<nion; ++ii) {
      double vx = rion0[3*ii];
      double vy = rion0[3*ii+1]; 
      double vz = rion0[3*ii+2];
      KE += 0.5*mass[ii]*(vx*vx + vy*vy + vz*vz);
   }

   QMMM_electrostatic_potential(nion_qm,nion,qion,rion1,uion);

   if (taskid==MASTER)
   {
      std::cout << std::endl << std::endl;
      for (auto ii=0; ii<nion; ++ii)
            std::cout << "@ ii=" << ii << " " << symbol[ii] << "\trion: " << Ffmt(12,6) << rion1[3*ii] << " " << Ffmt(12,6) << rion1[3*ii+1] << " " << Ffmt(12,6) << rion1[3*ii+2]
                      << " mass = "  << mass[ii] << " uion = " << uion[ii] << " qion = " << qion[ii] 
                      << " epsilon =" << epsilon[ii] << " sigma=" << sigma[ii] << std::endl;
      std::cout << std::endl;
      std::cout << "@ Initial Kinetic Energy = " << Efmt(20,15) << KE << std::endl;
   }

   MPI_Barrier(MPI_COMM_WORLD);

   if (taskid==MASTER) {
      xyzfile = new (std::nothrow) std::ofstream;
      xyzfile->open("testqmm.xyz", std::ios::app);

      emotionfile = new (std::nothrow) std::ofstream;
      emotionfile->open("testqmm.emotion", std::ios::app);
   }


   // QM energy and forces
   ierr += c_lammps_pspw_qmmm_minimizer_filename(MPI_COMM_WORLD,rion1,uion,fion,qion,&Eqm,true,true,NULL);


   // Eqq = Electrostatic energy and forces
   Eqq = QMQM_electrostatic_energy(nion_qm,nion,qion,rion1) 
       + QMMM_electrostatic_energy(nion_qm,nion,qion,rion1);
       //+ MMMM_electrostatic_energy(nion_qm,nion,qion,rion1);
   QMQM_electrostatic_force(nion_qm,nion,qion,rion1,fion);
   QMMM_electrostatic_force(nion_qm,nion,qion,rion1,fion);
  // MMMM_electrostatic_force(nion_qm,nion,qion,rion1,fion);
   Eqm += Eqq;

   // ELJ = Lennard-Jones energy and forces
   ELJ = QMMM_LJ_energy(nion_qm,nion,epsilon,sigma,rion1);
   QMMM_LJ_force(nion_qm,nion,epsilon,sigma,rion1,fion);
   Eqm += ELJ;


   // Espring = MM spring energy and forces
   Espring = spring_Energy(nfrag,indxfrag_start,kfrag, 
                           nbond_frag,bond_start_frag,bond_indx,Krb,
                           nangle_frag,angle_start_frag,angle_indx,Kra,rion1);
   spring_Force(nfrag,indxfrag_start,kfrag, 
                nbond_frag,bond_start_frag,bond_indx,Krb,
                nangle_frag,angle_start_frag,angle_indx,Kra,rion1,fion);
   Eqm += Espring;



   if (taskid==MASTER) {
      std::cout << "@ Initial Energy = " << Efmt(20,15) << Eqm << " Initial LJ Energy = " << Efmt(20,15) << ELJ 
                << " Initial Spring Energy = " << Efmt(20,15) << Espring << std::endl;
      std::cout << "@" << std::endl;
      std::cout << "@ Initial Forces" << std::endl;
      for (auto ii=0; ii<nion; ++ii)
            std::cout << "@ ii=" << ii << " " << symbol[ii] << "\tfion: " << Ffmt(12,6) << fion[3*ii] << " " << Ffmt(12,6) << fion[3*ii+1] << " " << Ffmt(12,6) << fion[3*ii+2]
                      << " mass = "  << mass[ii] << " uion = " << uion[ii] << std::endl;
      std::cout << "@" << std::endl;
   }

   // take a Newton step
   for (auto ii=0; ii<nion; ++ii)
   {
      rion2[3*ii]   = rion1[3*ii]   + dt*rion0[3*ii]   + 0.5*dti[ii]*fion[3*ii];
      rion2[3*ii+1] = rion1[3*ii+1] + dt*rion0[3*ii+1] + 0.5*dti[ii]*fion[3*ii+1];
      rion2[3*ii+2] = rion1[3*ii+2] + dt*rion0[3*ii+2] + 0.5*dti[ii]*fion[3*ii+2];
   }
   MPI_Bcast(rion2,3*nion,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD);

   for (auto it=0; it<200; ++it)
   {

      memcpy(rion0,rion1,3*nion*sizeof(double));
      memcpy(rion1,rion2,3*nion*sizeof(double));
      memset(fion,0, 3*nion*sizeof(double));

      QMMM_electrostatic_potential(nion_qm,nion,qion,rion1,uion);

      ierr += c_lammps_pspw_qmmm_minimizer_filename(MPI_COMM_WORLD,rion1,uion,fion,qion,&Eqm,true,true,NULL);

      // Eqq = Electrostatic energy and forces
      Eqq = QMQM_electrostatic_energy(nion_qm,nion,qion,rion1);
      QMQM_electrostatic_force(nion_qm,nion,qion,rion1,fion);
      Eqm += Eqq;

      // ELJ = Lenard-Johnes energy and forces
      ELJ = QMMM_LJ_energy(nion_qm,nion,epsilon,sigma,rion1);
      QMMM_LJ_force(nion_qm,nion,epsilon,sigma,rion1,fion);
      Eqm += ELJ;

      // Espring = MM spring energy and forces
      Espring = spring_Energy(nfrag,indxfrag_start,kfrag, 
                              nbond_frag,bond_start_frag,bond_indx,Krb,
                              nangle_frag,angle_start_frag,angle_indx,Kra,rion1);
      spring_Force(nfrag,indxfrag_start,kfrag, 
                   nbond_frag,bond_start_frag,bond_indx,Krb,
                   nangle_frag,angle_start_frag,angle_indx,Kra,rion1,fion);
      Eqm += Espring;


      // take a position Verlet step
      for (auto ii=0; ii<nion; ++ii)
      {
         rion2[3*ii]   = 2.0*rion1[3*ii]   - rion0[3*ii]   + dti[ii]*fion[3*ii];
         rion2[3*ii+1] = 2.0*rion1[3*ii+1] - rion0[3*ii+1] + dti[ii]*fion[3*ii+1];
         rion2[3*ii+2] = 2.0*rion1[3*ii+2] - rion0[3*ii+2] + dti[ii]*fion[3*ii+2];
      }
      MPI_Bcast(rion2,3*nion,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD);

      // kinetic energy
      for (auto i=0; i<3*nion; ++i) rion0[i] = h*(rion2[i] - rion0[i]);
      KE = 0.0;
      for (auto ii=0; ii<nion; ++ii)
      {
         double vx = rion0[3*ii];
         double vy = rion0[3*ii+1];
         double vz = rion0[3*ii+2];
         KE += mass[ii]*(vx*vx + vy*vy + vz*vz);
      }
      KE *= 0.5;
      if (taskid==MASTER) {
         std::cout << "@ it = " << Ifmt(5) << it+1 << " Energy+Kinetic = " << Efmt(20,15) << Eqm+KE 
                   << " Energy = " << Efmt(20,15) << Eqm << " Kinetic Energy = " << Efmt(20,15) << KE 
                   <<  " QQ Energy = " << Efmt(20,15) << Eqq 
                   <<  " LJ Energy = " << Efmt(20,15) << ELJ 
                   <<  " Spring Energy = " << Efmt(20,15) << Espring << std::endl;
         std::cout << "@@ "     << Ifmt(5) << it+1 << " " << Efmt(20,15) << dt*(it+1) << " " << Efmt(20,15) << Eqm+KE << " " 
                   << Efmt(20,15) << Eqm << " " << Efmt(20,15) << KE 
                   << " " << Efmt(20,15) << Eqq 
                   << " " << Efmt(20,15) << ELJ 
                   << " " << Efmt(20,15) << Espring << std::endl;
      }


   }



   //for (auto ii=0; ii<nion; ++ii) dti[ii] = (dt*dt)/mass[ii];

   if (taskid==MASTER)
   {
      xyzfile->close();
      emotionfile->close();
      delete xyzfile;
      delete emotionfile;
   }



}

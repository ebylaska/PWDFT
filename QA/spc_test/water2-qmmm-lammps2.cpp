
#include        <iomanip>
#include        <cstring>
#include 	<cmath>
#include 	<iostream>
#include 	<fstream>
#include 	<cstdio>
#include 	<string>
#include	<unistd.h>
#include        <sys/stat.h>

#include "mpi.h"

using namespace std;

extern int  lammps_pspw_qmmm_minimizer(MPI_Comm, double*, double*, double*, double*, double*);
extern void lammps_pspw_input(MPI_Comm, std::string&);

#define	MASTER 		0
#define ANGTOBOHR	1.88972687777


/**************************************************
 *                                                *
 *                spring_bond                     *
 *                                                *
 **************************************************/
double spring_bond(const int nbond, const int indx[], const double Kr0[], const double rion[])
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
 *                spring_bond_force               *
 *                                                *
 **************************************************/
void spring_bond_force(const int nbond, const int indx[], const double Kr0[], const double rion[], double fion[])
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
 *               spring_angle                     *
 *                                                *
 **************************************************/
double spring_angle(const int nangle, const int indx[], const double Kr0[], const double rion[])
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
 *              spring_angle_force                *
 *                                                *
 **************************************************/
void spring_angle_force(const int nangle, const int indx[], const double Kr0[], const double rion[], double fion[])
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

   double dVLJ = -(4.00*epsilon12/(r*r))*(12.0*u12-6.0*u6);

   f1[0] += (x)*dVLJ;
   f1[1] += (y)*dVLJ;
   f1[2] += (z)*dVLJ;

   f2[0] -= (x)*dVLJ;
   f2[1] -= (y)*dVLJ;
   f2[2] -= (z)*dVLJ;
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


/****************************************************
 *                                                  *
 *                 QMMM_LJ_energy                   *
 *                                                  *
 ****************************************************/
double QMMM_LJ_energy(const int nion_qm, const int nion,
                      const double epsilon[], const double sigma[], const double rion[])
{
   double E = 0.0;
   for (auto qm=0; qm<nion_qm; ++qm)
      for (auto mm=nion_qm; mm<nion; ++mm)
      {
         double epsilon12 = std::sqrt(epsilon[qm]*epsilon[mm]);
         double sigma12 = 0.5*(sigma[qm] + sigma[mm]);

         double x = rion[3*qm]   - rion[3*mm];
         double y = rion[3*qm+1] - rion[3*mm+1];
         double z = rion[3*qm+2] - rion[3*mm+2];
         double r = std::sqrt(x*x + y*y + z*z);

         double u = (sigma12/r);
         double u6  = u*u*u*u*u*u;
         double u12 = u6*u6;

         E += (4.0*epsilon12*(u12-u6));

      }
   return E;
}

/****************************************************
 *                                                  *
 *                 QMMM_LJ_force                    *
 *                                                  *
 ****************************************************/
void QMMM_LJ_force(const int nion_qm, const int nion,
                   const double epsilon[], const double sigma[], const double rion[], double fion[])
{
   for (auto qm=0; qm<nion_qm; ++qm)
      for (auto mm=nion_qm; mm<nion; ++mm) 
      {
         double epsilon12 = std::sqrt(epsilon[qm]*epsilon[mm]);
         double sigma12 = 0.5*(sigma[qm] + sigma[mm]);

         double x = rion[3*qm]   - rion[3*mm];
         double y = rion[3*qm+1] - rion[3*mm+1];
         double z = rion[3*qm+2] - rion[3*mm+2];
         double r = std::sqrt(x*x + y*y + z*z);

         double u = (sigma12/r);
         double u6  = u*u*u*u*u*u;
         double u12 = u6*u6;

         double dVLJ = -(4.00*epsilon12/(r*r))*(12.0*u12-6.0*u6);

         fion[3*mm]   += (x)*dVLJ;
         fion[3*mm+1] += (y)*dVLJ;
         fion[3*mm+2] += (z)*dVLJ;

         fion[3*qm]   -= (x)*dVLJ;
         fion[3*qm+1] -= (y)*dVLJ;
         fion[3*qm+2] -= (z)*dVLJ;
      }
}


/****************************************************
 *                                                  *
 *        QMMM_electrostatic_potential              *
 *                                                  *
 ****************************************************/
void QMMM_electrostatic_potential(const int nion_qm, const int nion, 
                                  const double qion[], const double rion[], double uion[]) 
{
   memset(uion,0,nion_qm*sizeof(double));
   for (auto qm=0; qm<nion_qm; ++qm)
      for (auto mm=nion_qm; mm<nion; ++mm)
      {
         double x = rion[3*qm]   - rion[3*mm];
         double y = rion[3*qm+1] - rion[3*mm+1];
         double z = rion[3*qm+2] - rion[3*mm+2];
         double r = std::sqrt(x*x + y*y + z*z);

         uion[qm] += qion[mm]/r;
      }
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
   for (auto qm=0; qm<nion_qm; ++qm)
      for (auto mm=nion_qm; mm<nion; ++mm)
      {
         double x = rion[3*qm]   - rion[3*mm];
         double y = rion[3*qm+1] - rion[3*mm+1];
         double z = rion[3*qm+2] - rion[3*mm+2];
         double rr = (x*x + y*y + z*z);
         double r  = std::sqrt(rr);

         E += qion[qm]*qion[mm]/r;
      }

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
   for (auto qm=0; qm<nion_qm; ++qm)
   {
      for (auto mm=nion_qm; mm<nion; ++mm)
      {
         double x = rion[3*qm]   - rion[3*mm];
         double y = rion[3*qm+1] - rion[3*mm+1];
         double z = rion[3*qm+2] - rion[3*mm+2];
         double rr = (x*x + y*y + z*z);
         double rrr  = std::sqrt(rr)*rr;
         double mder = -qion[qm]*qion[mm]/(rrr);

         fion[3*mm]   += x*mder;
         fion[3*mm+1] += y*mder;
         fion[3*mm+2] += z*mder;

         fion[3*qm]   -= x*mder;
         fion[3*qm+1] -= y*mder;
         fion[3*qm+2] -= z*mder;
      }
   }
}

/****************************************************
 *                                                  *
 *        QMQM_electrostatic_energy                 *
 *                                                  *
 ****************************************************/
double  QMQM_electrostatic_energy(const int nion,const double qion[], const double rion[])
{

   double e1 = 0.0;
   for (auto jj=0; jj<(nion-1); ++jj)
      for (auto ii=jj+1; ii<nion; ++ii)
      {
         double x = rion[3*jj]  -rion[3*ii];
         double y = rion[3*jj+1]-rion[3*ii+1];
         double z = rion[3*jj+2]-rion[3*ii+2];
         double r = std::sqrt(x*x + y*y + z*z);
         e1 += qion[ii]*qion[jj]/r;
      }
   return e1;
}
/****************************************************
 *                                                  *
 *        QMQM_electrostatic_force                  *
 *                                                  *
 ****************************************************/
void  QMQM_electrostatic_force(const int nion,const double qion[], const double rion[], double fion[])
{
   for (auto jj=0; jj<(nion-1); ++jj)
      for (auto ii=jj+1; ii<nion; ++ii)
      {
         double x = rion[3*jj]  -rion[3*ii];
         double y = rion[3*jj+1]-rion[3*ii+1];
         double z = rion[3*jj+2]-rion[3*ii+2];
         double rr = (x*x + y*y + z*z);
         double rrr = std::sqrt(rr)*rr;
         double v = -qion[ii]*qion[jj]/(rrr);

         fion[3*ii]   += x*v;
         fion[3*ii+1] += y*v;
         fion[3*ii+2] += z*v;

         fion[3*jj]   -= x*v;
         fion[3*jj+1] -= y*v;
         fion[3*jj+2] -= z*v;
      }
}


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
                  const double KE, const double Eqm, const double Ecoul, const double ELJ, const double Espring)
{
   *emotion << E1910 << current_time
         << E1910 << Etotal
         << E1910 << E
         << E1910 << KE
         << E1910 << Eqm 
         << E1910 << Ecoul 
         << E1910 << ELJ
         << E1910 << Espring
         << std::endl;
}





int main(int argc, char* argv[])
{
   std::ofstream *xyzfile;
   std::ofstream *emotionfile;

   //std::string nwfilename = "w2.nw";
   std::string nwfilename = "water1.nw";
   if (argc>1) nwfilename = argv[1];

   int ierr,np,taskid;


   // Initialize MPI
   ierr = MPI_Init(&argc,&argv);
   ierr += MPI_Comm_size(MPI_COMM_WORLD,&np);
   ierr += MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

   if (taskid==MASTER) {
      std::cout << "Hello world" << std::endl;
      std::cout << "np=" << np << " taskid=" << taskid << std::endl;
      std::cout << "argc=" << argc << std::endl;
      std::cout << "nwfilename=" << nwfilename << std::endl;

      xyzfile = new (std::nothrow) std::ofstream;
      xyzfile->open("testqmm.xyz", std::ios::app);

      emotionfile = new (std::nothrow) std::ofstream;
      emotionfile->open("testqmm.emotion", std::ios::app);
   }


   int nion_qm = 3;
   int nwater  = 1;
   int nion    = nion_qm + 3*nwater;
   double rion0[3*nion],rion1[3*nion],rion2[3*nion];
   double uion[nion],qion[nion],fion[3*nion];
   double epsilon[nion],sigma[nion],mass[nion],dti[nion];
   double E,KE,Eqm,Ecoul,ELJ,Espring;

   std::string symbol[nion] = {"O","H","H","O","H","H"};

   double pi = 4.0*std::atan(1.0);
   int nbond  = 2;
   int nangle = 1;
   int bond_indx[4] = {3,4,3,5};
   int angle_indx[3] = {4,3,5};
   double Krb[2] = {554.25/27.2116/23.06/ANGTOBOHR/ANGTOBOHR,1.0*ANGTOBOHR};
   double Kra[2] = {47.744/27.2116/23.06, 109.4*pi/180.0};

   lammps_pspw_input(MPI_COMM_WORLD, nwfilename);


   // first water molecule - QM water
   rion1[0] =  0.021259*ANGTOBOHR; rion1[1] =  0.506771*ANGTOBOHR; rion1[2] =  2.831278*ANGTOBOHR;
   rion1[3] = -0.721039*ANGTOBOHR; rion1[4] =  1.083100*ANGTOBOHR; rion1[5] =  2.758378*ANGTOBOHR;
   rion1[6] =  0.158220*ANGTOBOHR; rion1[7] =  0.181883*ANGTOBOHR; rion1[8] =  1.745696*ANGTOBOHR;
   qion[0] = -0.8476; qion[1] = 0.4238; qion[2] = 0.4238;
   uion[0] =  0.0000; uion[1] = 0.0000; uion[2] = 0.0000;
   sigma[0] = 3.165558*ANGTOBOHR; epsilon[0] = 0.155394/(23.06*27.2116);
   sigma[1] = 0.700000*ANGTOBOHR; epsilon[1] = 0.044000/(23.06*27.2116);
   sigma[2] = 0.700000*ANGTOBOHR; epsilon[2] = 0.044000/(23.06*27.2116);
   mass[0] = 16.0*1822.89; mass[1]= 2.0*1822.89; mass[2]=2.0*1822.89;

   // second water molecule - MM water
   rion1[9]  =  0.161560*ANGTOBOHR; rion1[10] = -0.052912*ANGTOBOHR; rion1[11] =  0.033173*ANGTOBOHR;
   rion1[12] =  0.803054*ANGTOBOHR; rion1[13] =  0.369132*ANGTOBOHR; rion1[14] = -0.511660*ANGTOBOHR;
   rion1[15] = -0.325571*ANGTOBOHR; rion1[16] = -0.669574*ANGTOBOHR; rion1[17] = -0.488560*ANGTOBOHR;
   qion[3]  = -0.8476; qion[4] = 0.4238; qion[5] = 0.4238;
   uion[3]  =  0.0000; uion[4] = 0.0000; uion[5] = 0.0000;
   sigma[3] = 3.165558*ANGTOBOHR; epsilon[3] = 0.155394/(23.06*27.2116);
   sigma[4] = 0.700000*ANGTOBOHR; epsilon[4] = 0.0;
   sigma[5] = 0.700000*ANGTOBOHR; epsilon[5] = 0.0;
   mass[3] = 16.0*1822.89; mass[4]= 2.0*1822.89; mass[5]=2.0*1822.89;



   //lj_ion_parameters O 3.165558 0.155394
   //lj_ion_parameters H 0.700000 0.044000
   //lj_ion_parameters O^ 3.165558 0.155394a

   //shake units angstroms 1 2 3 cyclic 1.0 1.632993125 1.0
   double unita[9] = {20.0,  0.0,  0.0,
                       0.0, 20.0,  0.0,
                       0.0,  0.0, 20.0};
   double dsq[3] = {1.0*ANGTOBOHR, 1.632993125*ANGTOBOHR, 1.0*ANGTOBOHR}; for (auto i=0; i<3; ++i) dsq[i] = dsq[i]*dsq[i];
   int mm_water_indx[3] = {3, 4, 5};
   double dt = 5.0;
   double h  = 1.0/(2.0*dt);
   for (auto ii=0; ii<nion; ++ii) dti[ii] = (dt*dt)/mass[ii];

   //memcpy(rion2,rion1,3*nion*sizeof(double));
   //memset(rion0,0,3*nion*sizeof(double));
   //shake_chain(3,mm_water_indx,3,1.0e-6,55,dsq,mass,unita,rion1,rion2); 

   // zero potentials and forces 
   memset(uion,0,nion*sizeof(double));
   memset(fion,0,3*nion*sizeof(double));

   // calculate electrostatic potential on QM atoms, u(ii)
   QMMM_electrostatic_potential(nion_qm,nion,qion,rion1,uion);

   // Eqm = QM energy and forces 
   ierr += lammps_pspw_qmmm_minimizer(MPI_COMM_WORLD,rion1,uion,fion,qion,&Eqm);

   // Ecoul = QMQM Electrostatic energy and forces
   Ecoul = QMQM_electrostatic_energy(nion_qm,qion,rion1);
   QMQM_electrostatic_force(nion_qm,qion,rion1,fion);
   //Ecoul = QMQM_electrostatic_energy(nion,qion,rion1);
   //QMQM_electrostatic_force(nion,qion,rion1,fion);

   // Ecoul += QMMM Electrostatic energy and forces
   double EAPC = 0.0;
   for (auto ii=0; ii<nion_qm; ++ii)
      EAPC += qion[ii]*uion[ii];
   double eqmma = QMMM_electrostatic_energy(nion_qm,nion,qion,rion1);
   QMMM_electrostatic_force(nion_qm,nion,qion,rion1,fion);
   Ecoul += EAPC;
   //Ecoul += eqmma;
   std::cout << "EAPC=" << EAPC << " " << eqmma << std::endl;

   // ELJ = QMMM Lenard-Jones energy and forces
   ELJ = QMMM_LJ_energy(nion_qm,nion,epsilon,sigma,rion1);
   QMMM_LJ_force(nion_qm,nion,epsilon,sigma,rion1,fion);

   // Espring = MM spring energy and forces
   Espring = spring_bond(nbond,bond_indx,Krb,rion1) + spring_angle(nangle,angle_indx,Kra,rion1);
   spring_bond_force(nbond,bond_indx,Krb,rion1,fion);
   spring_angle_force(nangle,angle_indx,Kra,rion1,fion);

   // take a Newton step
   for (auto ii=0; ii<nion; ++ii)
   {
      rion2[3*ii]   = rion1[3*ii]   + dt*rion0[3*ii]   + 0.5*dti[ii]*fion[3*ii];
      rion2[3*ii+1] = rion1[3*ii+1] + dt*rion0[3*ii+1] + 0.5*dti[ii]*fion[3*ii+1];
      rion2[3*ii+2] = rion1[3*ii+2] + dt*rion0[3*ii+2] + 0.5*dti[ii]*fion[3*ii+2];
   }

   // water
   MPI_Bcast(rion2,3*nion,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD);
   

   // kinetic energy
   KE = 0.0;
   for (auto ii=0; ii<nion; ++ii)
   {
      double vx = rion0[3*ii];
      double vy = rion0[3*ii+1];
      double vz = rion0[3*ii+2];
      KE += mass[ii]*(vx*vx + vy*vy + vz*vz);
   }
   KE *= 0.5;

   // Current step output
   if (taskid==MASTER)
   {
      std::cout << std::endl << std::endl;
      for (auto ii=0; ii<nion; ++ii)
         std::cout << "ii=" << ii << " rion: "  << rion1[3*ii]   << " " << rion1[3*ii+1]   << " " << rion1[3*ii+2]   << " uion=" << uion[ii]   << std::endl;

      std::cout << "KE+energy=" << KE+Eqm+Ecoul+ELJ+Espring << " energy=" << Eqm+Ecoul+ELJ+Espring << " KE=" << KE << " Eqm=" << Eqm << " Ecoul=" << Ecoul << " ELJ=" << ELJ << " Espring=" << Espring << std::endl;

      for (auto ii=0; ii<nion; ++ii)
         std::cout << "ii=" << ii << " fion: "  << fion[3*ii]   << " " << fion[3*ii+2]   << " " << fion[3*ii+2]   << " qion=" << qion[ii]   << std::endl;
      std::cout << std::endl << std::endl;

      //printxyz(xyzfile,nion,symbol,unita,rion1) ;
      //printemotion(emotionfile,0.0,KE+Eqm+Ecoul+ELJ+Espring,Eqm+Ecoul+ELJ,KE,Eqm,Ecoul,ELJ,Espring);
   }




   // Verlet Iterations 
   int nsteps = 100;
   for(auto it=0; it<nsteps; ++it)
   {
      //for (auto ii=0; ii<(3*nion); ++ii) rion0[ii] = rion1[ii];
      //for (auto ii=0; ii<(3*nion); ++ii) rion1[ii] = rion2[ii];
      memcpy(rion0,rion1,3*nion*sizeof(double));
      memcpy(rion1,rion2,3*nion*sizeof(double));

      // zero potentials and forces 
      memset(uion,0,nion*sizeof(double));
      memset(fion,0,3*nion*sizeof(double));

      // calculate electrostatic potential on QM atoms, u(ii)
      QMMM_electrostatic_potential(nion_qm,nion,qion,rion1,uion);

      // QM energy and forces
      ierr += lammps_pspw_qmmm_minimizer(MPI_COMM_WORLD,rion1,uion,fion,qion,&Eqm);

      // Ecoul = QMQM Electrostatic energy and forces
      //Ecoul = 0.0;
      Ecoul = QMQM_electrostatic_energy(nion_qm,qion,rion1);
      QMQM_electrostatic_force(nion_qm,qion,rion1,fion);
      //Ecoul = QMQM_electrostatic_energy(nion,qion,rion1);
      //QMQM_electrostatic_force(nion,qion,rion1,fion);
      std::cout << "EQMQM = " << Ecoul << std::endl;

      // QMMM Electrostatic energy and forces
      double eqmm = 0.0;
      for (auto ii=0; ii<nion_qm; ++ii)
         eqmm += qion[ii]*uion[ii];
      double eqmm0 = QMMM_electrostatic_energy(nion_qm,nion,qion,rion1);
      std::cout << "Ecoul_qmm = " << eqmm << " " << eqmm0 << std::endl;
     //Ecoul += eqmm0;
      Ecoul += eqmm;
      QMMM_electrostatic_force(nion_qm,nion,qion,rion1,fion);

      // QMMM Electrostatic energy and forces
      ELJ = QMMM_LJ_energy(nion_qm,nion,epsilon,sigma,rion1);
      QMMM_LJ_force(nion_qm,nion,epsilon,sigma,rion1,fion);

      // Espring = MM spring energy and forces
      Espring = spring_bond(nbond,bond_indx,Krb,rion1) + spring_angle(nangle,angle_indx,Kra,rion1);
      spring_bond_force(nbond,bond_indx,Krb,rion1,fion);
      spring_angle_force(nangle,angle_indx,Kra,rion1,fion);

      // take a position Verlet step
      for (auto ii=0; ii<nion; ++ii)
      {
         rion2[3*ii]   = 2.0*rion1[3*ii]   - rion0[3*ii]   + dti[ii]*fion[3*ii];
         rion2[3*ii+1] = 2.0*rion1[3*ii+1] - rion0[3*ii+1] + dti[ii]*fion[3*ii+1];
         rion2[3*ii+2] = 2.0*rion1[3*ii+2] - rion0[3*ii+2] + dti[ii]*fion[3*ii+2];
      }

      // shake MM water
      MPI_Bcast(rion2,3*nion,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD);


      // kinetic energy
      for (auto i=0; i<3*nion; ++i) rion0[i] = h*(rion2[i] - rion0[i]);
      
      KE = 0.0;
      for (auto ii=0; ii<nion; ++ii)
      {
         double vx = rion0[3*ii];
         double vy = rion0[3*ii+1];
         double vz = rion0[3*ii+1];
         KE += mass[ii]*(vx*vx + vy*vy + vz*vz);
      }
      KE *= 0.5;

      if (taskid==MASTER)
      {
         std::cout << std::endl << std::endl;
         for (auto ii=0; ii<nion; ++ii)
            std::cout << "ii=" << ii << " rion: "  << rion1[3*ii]   << " " << rion1[3*ii+1]   << " " << rion1[3*ii+2]   
                                     << " vion: "  << rion0[3*ii]   << " " << rion0[3*ii+1]   << " " << rion0[3*ii+2]   
                      << " uion=" << uion[ii]   << std::endl;

         std::cout << "@ KE+energy=" << it << " " << KE+Eqm+Ecoul+ELJ+Espring << " energy=" << Eqm+Ecoul+ELJ+Espring 
                   << " KE=" << KE << " Eqm=" << Eqm << " Ecoul=" << Ecoul << " ELJ=" << ELJ << " Espring=" << Espring <<  std::endl;

         for (auto ii=0; ii<nion; ++ii)
            std::cout << "ii=" << ii << " fion: "  << fion[3*ii]   << " " << fion[3*ii+2]   << " " << fion[3*ii+2]   << " qion=" << qion[ii]   << std::endl;
         std::cout << std::endl << std::endl;

         printxyz(xyzfile,nion,symbol,unita,rion1) ;
         printemotion(emotionfile,it*5.0,KE+Eqm+Ecoul+ELJ+Espring,Eqm+Ecoul+ELJ+Espring,KE,Eqm,Ecoul,ELJ,Espring);
      }
   }

   if (taskid==MASTER)
   {
      xyzfile->close();
      emotionfile->close();
      delete xyzfile;
      delete emotionfile;
   }



}

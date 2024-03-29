
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
 *        QMMM_QMQM_electrostatic_energy            *
 *                                                  *
 ****************************************************/
double QMMM_QMQM_electrostatic_energy(const int nion_qm, const int nion,
                                      const double qion[], const double rion[])
{
   double E=0.0;

   // total qmqm and qmmm and mm interactions
   for (auto ii=0; ii<(nion-1); ++ii)
      for (auto jj=ii+1; jj<nion; ++jj)
         E += Q_Electrostatic_self(&rion[3*ii],qion[ii], 
                                   &rion[3*jj],qion[jj]);

   return E;
}


/****************************************************
 *                                                  *
 *        QMMM_QMQM_electrostatic_force             *
 *                                                  *
 ****************************************************/
void QMMM_QMQM_electrostatic_force(const int nion_qm, const int nion,
                                  const double qion[], const double rion[], double fion[])
{
   // total qmqm and qmmm and mm interactions
   for (auto ii=0; ii<(nion-1); ++ii)
      for (auto jj=ii+1; jj<nion; ++jj)
         Q_Electrostatic_Force_self(&rion[3*ii],qion[ii],&fion[3*ii],
                                    &rion[3*jj],qion[jj],&fion[3*jj]);

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

   std::string geomblock;
   if (mystring_contains(mystring_lowercase(nwinput),"geometry"))
   {
      geomblock  = mystring_rtrim(mystring_ltrim(mystring_split(mystring_split(nwinput,"nocenter")[1],"end")[0]));
   }
   std::vector<std::string> geomlines = mystring_split(geomblock,"\n");
   int nion = geomlines.size();
   //for (auto ii=0; ii<nion; ++ii)
   //   std::cout << ii << " geomline:" << geomlines[ii] << std::endl;

   int nbond = 0;
   if (mystring_contains(mystring_lowercase(nwinput),"bond_spring"))
   {
      nbond = mystring_split(nwinput,"bond_spring").size() - 1;
   }
   int nangle = 0;
   if (mystring_contains(mystring_lowercase(nwinput),"bond_spring"))
   {
      nangle = mystring_split(nwinput,"angle_spring").size() - 1;
   }


   //int nbond  = 2;
   //int nangle = 1;
   //int bond_indx[4] = {3,4,3,5};
   //int angle_indx[3] = {4,3,5};
   //double Krb[2] = {554.25/27.2116/23.06/ANGTOBOHR/ANGTOBOHR,1.0*ANGTOBOHR};
   //double Kra[2] = {47.744/27.2116/23.06, 109.4*pi/180.0};

   

   if (taskid==MASTER) std::cout << "nwinput = " << nwinput << std::endl << std::endl;
   if (taskid==MASTER) std::cout << "geomblock = " << geomblock << std::endl;
   if (taskid==MASTER) std::cout << "nbond =" << nbond << std::endl;
   if (taskid==MASTER) std::cout << "nangle =" << nangle << std::endl;

   // Initialize lammps_pspw interface 
   c_lammps_pspw_input_filename(MPI_COMM_WORLD,cnwfilename,NULL);


   double unita[9] = {26.0,  0.0,  0.0,
                       0.0, 26.0,  0.0,
                       0.0,  0.0, 26.0};
   double dt = 5.0;
   double h  = 1.0/(2.0*dt);
   double rion0[3*nion],rion1[3*nion],rion2[3*nion],fion[3*nion]; 
   double mass[nion],dti[nion],KE,uion[nion],qion[nion];
   double Espring,ELJ,Eqm,Eqq;

   double epsilon[nion],sigma[nion];

   int    bond_indx[2*nbond],angle_indx[3*nangle];
   double Krb[2*nbond];
   double Kra[2*nangle];


   if (mystring_contains(mystring_lowercase(nwinput),"bond_spring"))
   {
      std::string line;
      std::vector<std::string> ss;
      for (auto b=0; b<nbond; ++b)
      {
         line = mystring_split(nwinput,"bond_spring")[b+1];
         ss   = mystring_split0(line);
         int i = std::stoi(ss[0]);
         int j = std::stoi(ss[1]);
         double k = std::stod(ss[2]);
         double r = std::stod(ss[3]);
         bond_indx[2*b]   = i-1;
         bond_indx[2*b+1] = j-1;
         Krb[2*b]   = k/27.2116/23.06/ANGTOBOHR/ANGTOBOHR;
         Krb[2*b+1] = r*ANGTOBOHR;
      }
   }

   double pi = 4.0*std::atan(1.0);
   if (mystring_contains(mystring_lowercase(nwinput),"angle_spring"))
   {
      std::string line;
      std::vector<std::string> ss;
      for (auto a=0; a<nangle; ++a)
      {
         line = mystring_split(nwinput,"angle_spring")[a+1];
         ss   = mystring_split0(line);
         int i = std::stoi(ss[0]);
         int j = std::stoi(ss[1]);
         int k = std::stoi(ss[2]);
         double ks    = std::stod(ss[3]);
         double theta = std::stod(ss[4]);
         std::cout << "i=" << i << " j=" << j << " k=" << k 
                   << " ks =" << ks << " theta=" << theta << std::endl;
         angle_indx[3*a]   = i-1;
         angle_indx[3*a+1] = j-1;
         angle_indx[3*a+2] = k-1;
         Kra[2*a]   = ks/27.2116/23.06;
         Kra[2*a+1] = theta*pi/180.0;
      }
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
         epsilon[nion_qm] = symboltoepsilon(symbol[nion_qm])/23.06/27.2116;
         sigma[nion_qm]   = symboltosigma(symbol[nion_qm])*ANGTOBOHR;
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
         epsilon[(nion_qm+nion_mm)] = symboltoepsilon(symbol[(nion_qm+nion_mm)])/23.06/27.2116;
         sigma[(nion_qm+nion_mm)]     = symboltosigma(symbol[(nion_qm+nion_mm)])*ANGTOBOHR;
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
   Espring = spring_bond(nbond,bond_indx,Krb,rion1) + spring_angle(nangle,angle_indx,Kra,rion1);
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
   Eqq = QMMM_QMQM_electrostatic_energy(nion_qm,nion,qion,rion1);
   QMMM_QMQM_electrostatic_force(nion_qm,nion,qion,rion1,fion);
   Eqm += Eqq;

   // ELJ = Lenard-Johnes energy and forces
   ELJ = QMMM_LJ_energy(nion_qm,nion,epsilon,sigma,rion1);
   QMMM_LJ_force(nion_qm,nion,epsilon,sigma,rion1,fion);
   Eqm += ELJ;


   // Espring = MM spring energy and forces
   Espring = spring_bond(nbond,bond_indx,Krb,rion1) + spring_angle(nangle,angle_indx,Kra,rion1);
   spring_bond_force(nbond,bond_indx,Krb,rion1,fion);
   spring_angle_force(nangle,angle_indx,Kra,rion1,fion);
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
      Eqq = QMMM_QMQM_electrostatic_energy(nion_qm,nion,qion,rion1);
      QMMM_QMQM_electrostatic_force(nion_qm,nion,qion,rion1,fion);
      Eqm += Eqq;

      // ELJ = Lenard-Johnes energy and forces
      ELJ = QMMM_LJ_energy(nion_qm,nion,epsilon,sigma,rion1);
      QMMM_LJ_force(nion_qm,nion,epsilon,sigma,rion1,fion);
      Eqm += ELJ;

      // Espring = MM spring energy and forces
      Espring = spring_bond(nbond,bond_indx,Krb,rion1) + spring_angle(nangle,angle_indx,Kra,rion1);
      spring_bond_force(nbond,bond_indx,Krb,rion1,fion);
      spring_angle_force(nangle,angle_indx,Kra,rion1,fion);
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

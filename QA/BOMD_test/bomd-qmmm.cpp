
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

#include "mpi.h"

#include "parsestring.hpp"

using namespace std;

extern int  lammps_pspw_qmmm_minimizer(MPI_Comm, double*, double*, double*, double*, double*);
extern void lammps_pspw_input(MPI_Comm, std::string&);

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








int main(int argc, char* argv[])
{
   std::ofstream *xyzfile;
   std::ofstream *emotionfile;

   //std::string nwfilename = "w2.nw";
   std::string nwfilename = "water0.nw";
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
   }

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
   

   if (taskid==MASTER) std::cout << "nwinput = " << nwinput << std::endl << std::endl;
   if (taskid==MASTER) std::cout << "geomblock = " << geomblock << std::endl;

   lammps_pspw_input(MPI_COMM_WORLD, nwfilename);



   double unita[9] = {26.0,  0.0,  0.0,
                       0.0, 26.0,  0.0,
                       0.0,  0.0, 26.0};
   double dt = 6.0;
   double h  = 1.0/(2.0*dt);
   double rion0[3*nion],rion1[3*nion],rion2[3*nion],fion[3*nion]; 
   double mass[nion],dti[nion],KE,uion[nion],qion[nion];
   double Eqm;
   std::string symbol[nion];
   memset(rion0,0,3*nion*sizeof(double));
   memset(rion1,0,3*nion*sizeof(double));
   memset(rion2,0,3*nion*sizeof(double));
   memset(fion,0, 3*nion*sizeof(double));
   memset(uion,0,   nion*sizeof(double));
   memset(qion,0,   nion*sizeof(double));

   for (auto ii=0; ii<nion; ++ii)
   {
      std::vector<std::string> ss = mystring_split0(geomlines[ii]);
      symbol[ii] = ss[0];
      mass[ii]   = symboltomass(symbol[ii])*1822.89;
      rion1[3*ii]   = std::stod(ss[1])*ANGTOBOHR;
      rion1[3*ii+1] = std::stod(ss[2])*ANGTOBOHR;
      rion1[3*ii+2] = std::stod(ss[3])*ANGTOBOHR;
      dti[ii] = dt*dt/mass[ii];
   }

   // kinetic energy
   KE = 0.0;
   for (auto ii=0; ii<nion; ++ii) {
      double vx = rion0[3*ii];
      double vy = rion0[3*ii+1]; 
      double vz = rion0[3*ii+1];
      KE += 0.5*mass[ii]*(vx*vx + vy*vy + vz*vz);
   }

   if (taskid==MASTER)
   {
      std::cout << std::endl << std::endl;
      for (auto ii=0; ii<nion; ++ii)
            std::cout << "@ ii=" << ii << " " << symbol[ii] << "\trion: " << Ffmt(12,6) << rion1[3*ii] << " " << Ffmt(12,6) << rion1[3*ii+1] << " " << Ffmt(12,6) << rion1[3*ii+2]
                      << " mass = "  << mass[ii] << std::endl;
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
   ierr += lammps_pspw_qmmm_minimizer(MPI_COMM_WORLD,rion1,uion,fion,qion,&Eqm);

   if (taskid==MASTER) {
      std::cout << "@ Initial Energy = " << Efmt(20,15) << Eqm << std::endl;
      std::cout << "@" << std::endl;
      std::cout << "@ Initial Forces" << std::endl;
      for (auto ii=0; ii<nion; ++ii)
            std::cout << "@ ii=" << ii << " " << symbol[ii] << "\tfion: " << Ffmt(12,6) << fion[3*ii] << " " << Ffmt(12,6) << fion[3*ii+1] << " " << Ffmt(12,6) << fion[3*ii+2]
                      << " mass = "  << mass[ii] << std::endl;
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
      memset(uion,0,   nion*sizeof(double));
      memset(qion,0,   nion*sizeof(double));
      ierr += lammps_pspw_qmmm_minimizer(MPI_COMM_WORLD,rion1,uion,fion,qion,&Eqm);


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
         double vz = rion0[3*ii+1];
         KE += mass[ii]*(vx*vx + vy*vy + vz*vz);
      }
      KE *= 0.5;
      if (taskid==MASTER) {
         std::cout << "@ it = " << it+1 << " Energy+Kinetic = " << Efmt(20,15) << Eqm+KE << " Energy = " << Efmt(20,15) << Eqm << " Kinetic Energy = " << Efmt(20,15) << KE << std::endl;
         std::cout << "@@ " << it+1 << " " << Efmt(20,15) << Eqm+KE << " " 
                   << Efmt(20,15) << Eqm << " " << Efmt(20,15) << KE << std::endl;
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

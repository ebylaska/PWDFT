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

#include "qmmm.hpp"

using namespace std;

extern int  c_lammps_pspw_qmmm_minimizer_filename(MPI_Comm,double*,double*,double*,double*,double*,bool,bool,const char*);
extern void c_lammps_pspw_input_filename(MPI_Comm,const char*,const char*);

extern int  c_lammps_pspw_qmmm_nominimizer_filename(MPI_Comm,double*,double*,double*,double*,double*,bool,bool,const char*);

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
   std::string nwfilename = "ftest.nw";
   if (argc>1) nwfilename = argv[1];

   std::string debugfilename = "";

   int ierr,np,taskid;

   // Initialize MPI
   ierr = MPI_Init(&argc,&argv);
   ierr += MPI_Comm_size(MPI_COMM_WORLD,&np);
   ierr += MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

   if (taskid==MASTER) {
      debugfilename = "ftest.nwout";
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


   if (taskid==MASTER) std::cout << "nwinput = " << nwinput << std::endl << std::endl;
   if (taskid==MASTER) std::cout << "geomblock = " << geomblock << std::endl;
   if (taskid==MASTER) std::cout << "nion      = " << nion << std::endl;



   // Initialize lammps_pspw interface 
   c_lammps_pspw_input_filename(MPI_COMM_WORLD,cnwfilename,NULL);


   double unita[9] = {26.0,  0.0,  0.0,
                       0.0, 26.0,  0.0,
                       0.0,  0.0, 26.0};
   double dt = 5.0;
   double h  = 1.0/(2.0*dt);
   double rion0[3*nion],rion1[3*nion],rion2[3*nion],fion[3*nion],fion0[3*nion],fionm[3*nion],fionp[3*nion],rionm[3*nion],rionp[3*nion],fiont[3*nion]; 
   double dti[nion],KE,uion[nion],qion[nion];
   double Espring,ELJ,Eqm,Eqq,Etotal,Eqm0,Eqm1,Eqm2;
   double uionp[nion],uionm[nion],qionp[nion],qionm[nion];
   double uiont[nion],qiont[nion];


   std::string symbol[nion];
   std::memset(rion0,0,3*nion*sizeof(double));
   std::memset(rion1,0,3*nion*sizeof(double));
   std::memset(rion2,0,3*nion*sizeof(double));
   std::memset(fion,0, 3*nion*sizeof(double));
   std::memset(uion,0,   nion*sizeof(double));
   std::memset(qion,0,   nion*sizeof(double));

   std::memset(fion0,0, 3*nion*sizeof(double));

   if (taskid==MASTER) std::cout << "nion =" << nion << std::endl;

   for (auto ii=0; ii<nion; ++ii)
   {
      std::vector<std::string> ss = mystring_split0(geomlines[ii]);
      symbol[ii] = ss[0];
      rion1[3*ii]   = std::stod(ss[1])*ANGTOBOHR;
      rion1[3*ii+1] = std::stod(ss[2])*ANGTOBOHR;
      rion1[3*ii+2] = std::stod(ss[3])*ANGTOBOHR;
   }



   // QM energy and forces
   uion[0] = 0.01;
   uion[3] = -0.01;
   ierr += c_lammps_pspw_qmmm_minimizer_filename(MPI_COMM_WORLD,rion1,uion,fion,qion,&Eqm,true,false,cdebugfilename);
   Etotal = Eqm;
   std::cout << "EQM1 = " << Eqm << std::endl;

   std::memset(uionp,0,nion*sizeof(double));
   std::memset(qionp,0,nion*sizeof(double));
   double Ep,Em;
   double deltaxyz = 0.0001;
   std::memset(fiont,0, 3*nion*sizeof(double));
   std::memset(uiont,0,   nion*sizeof(double));
   std::memset(qiont,0,   nion*sizeof(double));
   //ierr += c_lammps_pspw_qmmm_nominimizer_filename(MPI_COMM_WORLD,rion1,uiont,fiont,qiont,&Espring,false,false,cdebugfilename);
   //std::cout << "EQM2 = " << Espring << std::endl;
   for (auto ii=0; ii<nion; ++ii)
      std::cout << " ii=" << ii << "\tqion=" << Ffmt(12,6) << qion[ii] << " qiont=" << Ffmt(12,6) << qiont[ii] << std::endl;
   for (auto ii=0; ii<nion; ++ii)
      std::cout << " ii=" << ii << "\tuion=" << Ffmt(12,6) << uion[ii] << " uiont=" << Ffmt(12,6) << uiont[ii] << std::endl;

   for (auto ii=0; ii<nion; ++ii)
   {
      for (auto ijk=0; ijk<3; ++ijk)
      {
         std::memcpy(rionp,rion1,3*nion*sizeof(double));
         rionp[3*ii+ijk] += deltaxyz;
         std::memcpy(uionp,uion,nion*sizeof(double));
         std::memcpy(qionp,qion,nion*sizeof(double));
         ierr += c_lammps_pspw_qmmm_nominimizer_filename(MPI_COMM_WORLD,rionp,uionp,fionp,qionp,&Ep,true,false,cdebugfilename);

         std::memcpy(rionm,rion1,3*nion*sizeof(double));
         rionm[3*ii+ijk] -= deltaxyz;
         std::memcpy(uionm,uion,nion*sizeof(double));
         std::memcpy(qionm,qion,nion*sizeof(double));
         ierr += c_lammps_pspw_qmmm_nominimizer_filename(MPI_COMM_WORLD,rionm,uionm,fionm,qionm,&Em,true,false,cdebugfilename);

         fion0[3*ii+ijk] = -(Ep-Em)/(2.0*deltaxyz);
         double ee0 = fion0[3*ii+ijk] - fion[3*ii+ijk];
         if (taskid==MASTER) {
            std::cout << "fioncheck: "     << Ifmt(5) << 3*ii+ijk << " " 
                      << " fd fion=" <<  Efmt(23,15) << fion0[3*ii+ijk] 
                      << " fion="    << Efmt(23,15)  << fion[3*ii+ijk] 
                      << " error="   << Efmt(23,15)  << ee0
                      << "   Ep="      << Efmt(22,15)  << Ep 
                      <<   " Em="      << Efmt(22,15)  << Em 
                      << "   dE="      << Efmt(13,6)  << Ep-Em << std::endl;
         }
      }
   }

}

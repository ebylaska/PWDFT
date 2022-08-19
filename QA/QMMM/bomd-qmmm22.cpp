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


   //set up qmmm
   QMMM_Operator qmmm(nwinput);

   int nion  = qmmm.nion;
   int nion_qm  = qmmm.nion_qm;
   int nion_mm  = qmmm.nion_mm;
   int nkatm = qmmm.nkatm;

   std::cout << "nion=" << nion << std::endl;
   std::cout << "nion_qm=" << nion_qm << std::endl;
   std::cout << "nion_mm=" << nion_mm << std::endl;
   std::cout << "nkatm="   << nkatm << std::endl;

   //if (taskid==MASTER) {
   //   std::cout << "katm= ";
   //   for (auto ii=0; ii<nion; ++ii)
   //      std::cout << " " << qmmm.katm[ii];
   //   std::cout << std::endl;
  

    //  std::cout << std::endl;
    //  std::cout << "aname size=" << nkatm << std::endl;
    //  for (auto ia=0; ia<nkatm; ++ia)
    //      std::cout << "aname =" << ia << " " << qmmm.aname[ia] << std::endl;
    //  std::cout << std::endl;
   //}




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


   std::string symbol[nion];
   std::memset(rion0,0,3*nion*sizeof(double));
   std::memcpy(rion1,qmmm.rion,3*nion*sizeof(double));
   std::memset(rion2,0,3*nion*sizeof(double));
   std::memset(fion,0, 3*nion*sizeof(double));
   std::memset(uion,0,   nion*sizeof(double));
   std::memset(qion,0,   nion*sizeof(double));


   if (taskid==MASTER) std::cout << "nion =" << nion << std::endl;
   if (taskid==MASTER) std::cout << "nion_mm =" << nion_mm << std::endl;
   if (taskid==MASTER) std::cout << "nion_qm =" << nion_qm << std::endl << std::endl;

   // kinetic energy
   Espring = qmmm.spring_Energy(rion1);
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
            std::cout << "@ ii=" << ii << " " << qmmm.symbol[ii] << "\trion: " << Ffmt(12,6) << rion1[3*ii] << " " << Ffmt(12,6) << rion1[3*ii+1] << " " << Ffmt(12,6) << rion1[3*ii+2]
                      << " mass = "  << qmmm.mass[ii] << " uion = " << uion[ii] << " qion = " << qion[ii] 
                      << " epsilon =" << qmmm.epsilon[ii] << " sigma=" << qmmm.sigma[ii] << std::endl;
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
   ELJ = qmmm.QMMM_LJ_energy(rion1);
   qmmm.QMMM_LJ_force(rion1,fion);
   Eqm += ELJ;


   // Espring = MM spring energy and forces
   Espring = qmmm.spring_Energy(rion1);
   qmmm.spring_Force(rion1,fion);
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
      ELJ = qmmm.QMMM_LJ_energy(rion1);
      qmmm.QMMM_LJ_force(rion1,fion);
      Eqm += ELJ;

      // Espring = MM spring energy and forces
      Espring = qmmm.spring_Energy(rion1);
      qmmm.spring_Force(rion1,fion);
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

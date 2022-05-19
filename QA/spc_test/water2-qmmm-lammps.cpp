
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
 *               shake_chain                      *
 *                                                *
 **************************************************/
/*
   Simple implementation of standard shake algorithm.
   Entry - n: number of atoms in molecule
           indx: index to atoms in r1 and r2 array
           nb: number of bonds (constraints) in molecule
           tol, maxit: tolerance and maxiteration
           dsq: array of squared distances, i.e., shake/bond distance constaints
           mass: array of masses
           unita: lattice vectors
           r1,r2: atom array
   Exit - r1,r2: atom arrays adjusted to satisfy bond distance constraints
*/
static void shake_chain(const int n, const int indx[], const int nb,
                 const double tol, const int maxit,
                 const double dsq[], double mass[], const double unita[],
                 double r1[], double r2[]) 
{
   bool   moving[n],moved[n]; 
   double da[3],ua[9];
   double rxi[n],ryi[n],rzi[n];
   double pxi[n],pyi[n],pzi[n];
   double rptol = 1.0e-6;

   // Determine the unit lattice vectors and distances
   for (auto j=0; j<3; ++j)
   {
      da[j] = std::sqrt(unita[3*j]*unita[3*j] + unita[3*j+1]*unita[3*j+1] + unita[3*j+2]*unita[3*j+2]);
      ua[3*j]   = unita[3*j]  /da[j];
      ua[3*j+1] = unita[3*j+1]/da[j];
      ua[3*j+2] = unita[3*j+2]/da[j];
   }

   double tol2 = 2.0*tol;
   for (auto a=0; a<n; ++a)
   {
      rxi[a] = r1[3*indx[a]];
      ryi[a] = r1[3*indx[a]+1];
      rzi[a] = r1[3*indx[a]+2];

      pxi[a] = r2[3*indx[a]];
      pyi[a] = r2[3*indx[a]+1];
      pzi[a] = r2[3*indx[a]+2];
      moving[a] = false;
      moved[a]  = true;
   }

   // shake iteration
   bool done = false;
   int it = 0;
   while ((!done) && (it<=maxit)) {
      done = true;

      for (auto a=0; a<nb; ++a) {
         auto b = (a+1)%n;
         if (moved[a] || moved[b]) {
            double pxab = pxi[a] - pxi[b];
            double pyab = pyi[a] - pyi[b];
            double pzab = pzi[a] - pzi[b];
            double cp1 = pxab*ua[0] + pyab*ua[1] + pzab*ua[2];
            double cp2 = pxab*ua[3] + pyab*ua[4] + pzab*ua[5];
            double cp3 = pxab*ua[6] + pyab*ua[7] + pzab*ua[8];
            cp1 -= round(cp1/da[0])*da[0];
            cp2 -= round(cp2/da[1])*da[1];
            cp3 -= round(cp3/da[2])*da[2];
            pxab = ua[0]*cp1 + ua[3]*cp2 + ua[6]*cp3;
            pyab = ua[1]*cp1 + ua[4]*cp2 + ua[7]*cp3;
            pzab = ua[2]*cp1 + ua[5]*cp2 + ua[8]*cp3;

            double pabsq  = pxab*pxab + pyab*pyab + pzab*pzab;
            double rabsq  = dsq[a];
            double diffsq = rabsq - pabsq;
            std::cout << "a=" << a << " b=" << b << " diffsq=" << diffsq << " dsq=" << dsq[a] << " pabsq=" << pabsq << std::endl;


            if (fabs(diffsq) > (rabsq*tol2)) {
               double rxab = rxi[a] - rxi[b];
               double ryab = ryi[a] - ryi[b];
               double rzab = rzi[a] - rzi[b];
               double cr1 = rxab*ua[0] + ryab*ua[1] + rzab*ua[2];
               double cr2 = rxab*ua[3] + ryab*ua[4] + rzab*ua[5];
               double cr3 = rxab*ua[6] + ryab*ua[7] + rzab*ua[8];
               cr1 -= round(cr1/da[0])*da[0];
               cr2 -= round(cr2/da[1])*da[1];
               cr3 -= round(cr3/da[2])*da[2];
               rxab = ua[0]*cr1 + ua[3]*cr2 + ua[6]*cr3;
               ryab = ua[1]*cr1 + ua[4]*cr2 + ua[7]*cr3;
               rzab = ua[2]*cr1 + ua[5]*cr2 + ua[8]*cr3;
               

               double rpab = rxab*pxab + ryab*pyab + rzab*pzab;

               if (rpab < (rabsq*rptol)) std::cout << "SHAKE CONSTRAINT FAILURE" << std::endl;
               
               double rma = 1.0/mass[a];
               double rmb = 1.0/mass[b];
               double gab = diffsq/(2.0*(rma+rmb) * rpab);
               double dx  = rxab*gab;
               double dy  = ryab*gab;
               double dz  = rzab*gab;

               pxi[a] += rma*dx;
               pyi[a] += rma*dy;
               pzi[a] += rma*dz;
               pxi[b] -= rmb*dx;
               pyi[b] -= rmb*dy;
               pzi[b] -= rmb*dz;

               moving[a] = true;
               moving[b] = true;
               done = false;
            }
         }
      }

      for (auto a=0; a<n; ++a) {
         moved[a] = moving[a];
         moving[a] = false;
      }
      ++it;
   }

   for (auto a=0; a<n; ++a) {
      r2[3*indx[a]]   = pxi[a];
      r2[3*indx[a]+1] = pyi[a];
      r2[3*indx[a]+2] = pzi[a];

      r1[3*indx[a]]   = rxi[a];
      r1[3*indx[a]+1] = ryi[a];
      r1[3*indx[a]+2] = rzi[a];
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

   f1[0] -= (x/r)*dVLJ;
   f1[1] -= (y/r)*dVLJ;
   f1[2] -= (z/r)*dVLJ;
   f2[0] += (x/r)*dVLJ;
   f2[1] += (y/r)*dVLJ;
   f2[2] += (z/r)*dVLJ;
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

   f2[0] -= fx;
   f2[1] -= fy;
   f2[2] -= fz;

   f1[0] += fx;
   f1[1] += fy;
   f1[2] += fz;
}

double QMMM_LJ_energy(const int nion_qm, const int nion,
                      const double epsilon[], const double sigma[], const double rion[])
{
   double E = 0.0;
   // QMQM LJ == 0
   // QMMM LJ
   for (auto qm=0; qm<nion_qm; ++qm)
      for (auto mm=nion_qm; mm<nion; ++mm)
      {
         double epsilon12 = std::sqrt(epsilon[qm]*epsilon[mm]);
         double sigma12 = 0.5*(sigma[qm] + sigma[mm]);
         E += LJ_energy(epsilon12,sigma12,&rion[3*qm],&rion[3*mm]);
      }
   return E;
}

void QMMM_LJ_force(const int nion_qm, const int nion,
                   const double epsilon[], const double sigma[], const double rion[], double fion[])
{
   for (auto qm=0; qm<nion_qm; ++qm)
      for (auto mm=nion_qm; mm<nion; ++mm) {
         double epsilon12 = std::sqrt(epsilon[qm]*epsilon[mm]);
         double sigma12 = 0.5*(sigma[qm] + sigma[mm]);
         LJ_force(epsilon12,sigma12,&rion[3*qm],&fion[3*qm],&rion[3*mm],&fion[3*mm]);
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
   for (auto qm=0; qm<nion_qm; ++qm)
   {
      uion[qm] = 0.0;
      for (auto mm=nion_qm; mm<nion; ++mm)
     {
         uion[qm] += Q_Electrostatic_potential(&rion[3*mm], qion[mm], &rion[3*qm]);
         std::cout << qm << " " << mm << " uion=" << uion[qm] << std::endl;
     }
     std::cout << std::endl;
   }
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
         double r  = std::sqrt(rr);
         double mder = qion[qm]*qion[mm]/(rr);
         fion[3*qm]   += mder*(x/r);
         fion[3*qm+1] += mder*(y/r);
         fion[3*qm+2] += mder*(z/r);

         fion[3*mm]   -= mder*(x/r);
         fion[3*mm+1] -= mder*(y/r);
         fion[3*mm+2] -= mder*(z/r);
      }

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





int main(int argc, char* argv[])
{
   std::ofstream *xyzfile;

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
      xyzfile->open("water2qmm.xyz", std::ios::app);
   }


   int nion_qm = 3;
   int nwater  = 1;
   int nion    = nion_qm + 3*nwater;
   double rion0[3*nion],rion1[3*nion],rion2[3*nion];
   double uion[nion],qion[nion],fion[3*nion];
   double epsilon[nion],sigma[nion],mass[nion];
   double E,KE,Eqm,Ecoul,ELJ;

   std::string symbol[nion] = {"O","H","H","O","H","H"};

   lammps_pspw_input(MPI_COMM_WORLD, nwfilename);


   // first water molecule - QM water
   memset(rion0,0,3*nion*sizeof(double));
   memset(rion2,0,3*nion*sizeof(double));
   rion1[0] =  0.021259*ANGTOBOHR; rion1[1] =  0.506771*ANGTOBOHR; rion1[2] =  2.831278*ANGTOBOHR;
   rion1[3] = -0.721039*ANGTOBOHR; rion1[4] =  1.083100*ANGTOBOHR; rion1[5] =  2.758378*ANGTOBOHR;
   rion1[6] =  0.158220*ANGTOBOHR; rion1[7] =  0.181883*ANGTOBOHR; rion1[8] =  1.945696*ANGTOBOHR;
   qion[0] = -0.8476; qion[1] = 0.4238; qion[2] = 0.4238;
   uion[0] =  0.0000; uion[1] = 0.0000; uion[2] = 0.0000;
   sigma[0] = 3.165558*ANGTOBOHR; epsilon[0] = 0.155394/(23.06*27.2116);
   sigma[1] = 0.700000*ANGTOBOHR; epsilon[1] = 0.044000/(23.06*27.2116);
   sigma[2] = 0.700000*ANGTOBOHR; epsilon[2] = 0.044000/(23.06*27.2116);
   mass[0] = 16.0*1822.89; mass[1]= 1.0*1822.89; mass[2]=1.0*1822.90;

   // second water molecule - MM water
   rion1[9]  =  0.161560*ANGTOBOHR; rion1[10] = -0.052912*ANGTOBOHR; rion1[11] =  0.033173*ANGTOBOHR;
   rion1[12] =  0.803054*ANGTOBOHR; rion1[13] =  0.369132*ANGTOBOHR; rion1[14] = -0.511660*ANGTOBOHR;
   rion1[15] = -0.325571*ANGTOBOHR; rion1[16] = -0.669574*ANGTOBOHR; rion1[17] = -0.488560*ANGTOBOHR;
   qion[3]  = -0.8476; qion[4] = 0.4238; qion[5] = 0.4238;
   uion[3]  =  0.0000; uion[4] = 0.0000; uion[5] = 0.0000;
   sigma[3] = 3.165558*ANGTOBOHR; epsilon[3] = 0.155394/(23.06*27.2116);
   sigma[4] = 0.700000*ANGTOBOHR; epsilon[4] = 0.0;
   sigma[5] = 0.700000*ANGTOBOHR; epsilon[5] = 0.0;
   mass[3] = 16.0*1822.89; mass[4]= 1.0*1822.89; mass[5]=1.0*1822.90;


   //lj_ion_parameters O 3.165558 0.155394
   //lj_ion_parameters H 0.700000 0.044000
   //lj_ion_parameters O^ 3.165558 0.155394a

   //shake units angstroms 1 2 3 cyclic 1.0 1.632993125 1.0
   double unita[9] = {20.0,  0.0,  0.0,
                       0.0, 20.0,  0.0,
                       0.0,  0.0, 20.0};
   double dsq[3] = {1.0*ANGTOBOHR, 1.632993125*ANGTOBOHR, 1.0*ANGTOBOHR};
   int mm_water_indx[3] = {3, 4, 5};
   double dt = 5.0;

   // zero potentials and forces 
   memset(uion,0,nion*sizeof(double));
   memset(fion,0,3*nion*sizeof(double));

   // calculate electrostatic potential on QM atoms, u(ii)
   QMMM_electrostatic_potential(nion_qm,nion,qion,rion1,uion);

   // QM energy and forces
   ierr += lammps_pspw_qmmm_minimizer(MPI_COMM_WORLD,rion1,uion,fion,qion,&Eqm);

   // QMMM Electrostatic energy and forces
   Ecoul = 0.0;
   for (auto ii=0; ii<nion_qm; ++ii)
      Ecoul += qion[ii]*uion[ii];
   QMMM_electrostatic_force(nion_qm,nion,qion,rion1,fion);

   // QMMM Electrostatic energy and forces
   ELJ = QMMM_LJ_energy(nion_qm,nion,epsilon,sigma,rion1);
   QMMM_LJ_force(nion_qm,nion,epsilon,sigma,rion1,fion);

   // take a Newton step
   for (auto ii=0; ii<nion; ++ii)
   {
      rion2[3*ii]   = rion1[3*ii]   + dt*rion0[3*ii]   + 0.5*(dt*dt/mass[ii])*fion[3*ii];
      rion2[3*ii+1] = rion1[3*ii+1] + dt*rion0[3*ii+1] + 0.5*(dt*dt/mass[ii])*fion[3*ii+1];
      rion2[3*ii+2] = rion1[3*ii+2] + dt*rion0[3*ii+2] + 0.5*(dt*dt/mass[ii])*fion[3*ii+2];
   }

   // shake MM water
   shake_chain(3,mm_water_indx,3,1.0e-6,25,dsq,mass,unita,rion1,rion2); 

   // kinetic energy
   KE = 0.0;
   for (auto ii=0; ii<nion; ++ii)
      KE += 0.5*mass[ii]*(rion0[3*ii]*rion0[3*ii] + rion0[3*ii+1]*rion0[3*ii+1] + rion0[3*ii+2]*rion0[3*ii+2]);

   // Current step output
   if (taskid==MASTER)
   {
      std::cout << std::endl << std::endl;
      for (auto ii=0; ii<nion; ++ii)
         std::cout << "ii=" << ii << " rion: "  << rion1[3*ii]   << " " << rion1[3*ii+1]   << " " << rion1[3*ii+2]   << " uion=" << uion[ii]   << std::endl;

      std::cout << "KE+energy=" << KE+Eqm+Ecoul+ELJ << " energy=" << Eqm+Ecoul+ELJ << " KE=" << KE << " Eqm=" << Eqm << " Ecoul=" << Ecoul << " ELJ=" << ELJ << std::endl;

      for (auto ii=0; ii<nion; ++ii)
         std::cout << "ii=" << ii << " fion: "  << fion[3*ii]   << " " << fion[3*ii+2]   << " " << fion[3*ii+2]   << " qion=" << qion[ii]   << std::endl;
      std::cout << std::endl << std::endl;

      printxyz(xyzfile,nion,symbol,unita,rion1) ;
   }




   // Verlet Iterations 
   for(auto it=0; it<1; ++it)
   {
      for (auto ii=0; ii<(3*nion); ++ii) rion0[ii] = rion1[ii];
      for (auto ii=0; ii<(3*nion); ++ii) rion1[ii] = rion2[ii];

      // zero potentials and forces 
      memset(uion,0,nion*sizeof(double));
      memset(fion,0,3*nion*sizeof(double));

      // calculate electrostatic potential on QM atoms, u(ii)
      QMMM_electrostatic_potential(nion_qm,nion,qion,rion1,uion);

      // QM energy and forces
      ierr += lammps_pspw_qmmm_minimizer(MPI_COMM_WORLD,rion1,uion,fion,qion,&Eqm);

      // QMMM Electrostatic energy and forces
      Ecoul = 0.0;
      for (auto ii=0; ii<nion_qm; ++ii)
         Ecoul += qion[ii]*uion[ii];
      QMMM_electrostatic_force(nion_qm,nion,qion,rion1,fion);

      // QMMM Electrostatic energy and forces
      ELJ = QMMM_LJ_energy(nion_qm,nion,epsilon,sigma,rion1);
      QMMM_LJ_force(nion_qm,nion,epsilon,sigma,rion1,fion);


      // take a position Verlet step
      for (auto ii=0; ii<nion; ++ii)
      {
         rion2[3*ii]   = 2.0*rion1[3*ii]   - rion0[3*ii]   + (dt*dt/mass[ii])*fion[3*ii];
         rion2[3*ii+1] = 2.0*rion1[3*ii+1] - rion0[3*ii+1] + (dt*dt/mass[ii])*fion[3*ii+1];
         rion2[3*ii+2] = 2.0*rion1[3*ii+2] - rion0[3*ii+2] + (dt*dt/mass[ii])*fion[3*ii+2];
      }

      // shake MM water
      shake_chain(3,mm_water_indx,3,1.0e-6,25,dsq,mass,unita,rion1,rion2); 


      // kinetic energy
      for (auto ii=0; ii<(3*nion); ++ii) 
         rion0[ii] = (rion2[ii]-rion0[ii])/(2.0*dt);
      KE = 0.0;
      for (auto ii=0; ii<nion; ++ii)
      {
         KE += 0.5*mass[ii]*(rion0[3*ii]*rion0[3*ii] + rion0[3*ii+1]*rion0[3*ii+1] + rion0[3*ii+2]*rion0[3*ii+2]);
         std::cout << ii << " mass, ke=" << mass[ii] << " " << KE << std::endl;
      }

      if (taskid==MASTER)
      {
         std::cout << std::endl << std::endl;
         for (auto ii=0; ii<nion; ++ii)
            std::cout << "ii=" << ii << " rion: "  << rion1[3*ii]   << " " << rion1[3*ii+1]   << " " << rion1[3*ii+2]   
                                     << " vion: "  << rion0[3*ii]   << " " << rion0[3*ii+1]   << " " << rion0[3*ii+2]   
                      << " uion=" << uion[ii]   << std::endl;

         std::cout << "KE+energy=" << KE+Eqm+Ecoul+ELJ << " energy=" << Eqm+Ecoul+ELJ << " KE=" << KE << " Eqm=" << Eqm << " Ecoul=" << Ecoul << " ELJ=" << ELJ << std::endl;

         for (auto ii=0; ii<nion; ++ii)
            std::cout << "ii=" << ii << " fion: "  << fion[3*ii]   << " " << fion[3*ii+2]   << " " << fion[3*ii+2]   << " qion=" << qion[ii]   << std::endl;
         std::cout << std::endl << std::endl;

         printxyz(xyzfile,nion,symbol,unita,rion1) ;
      }
   }

   if (taskid==MASTER)
   {
      xyzfile->close();
      delete xyzfile;
   }



}

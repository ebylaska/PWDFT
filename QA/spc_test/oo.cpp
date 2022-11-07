
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



double Q_r_energy(const double q12, const double r)
{
   return q12/r;
}

void Q_fd_force(const double q12, const double r12[], double f12[])
{
   double e0,epx,emx,epy,emy,epz,emz,r;

   double dx = r12[0];
   double dy = r12[1];
   double dz = r12[2];

   r = std::sqrt(dx*dx + dy*dy + dz*dz); 
   e0 = Q_r_energy(q12,r);

   dx = r12[0] + 0.001; r = std::sqrt(dx*dx + dy*dy + dz*dz); epx = Q_r_energy(q12,r);
   dx = r12[0] - 0.001; r = std::sqrt(dx*dx + dy*dy + dz*dz); emx = Q_r_energy(q12,r);
   dx = r12[0];
   f12[0] -= (epx-emx)/0.002;

   dy = r12[1] + 0.001; r = std::sqrt(dx*dx + dy*dy + dz*dz); epy = Q_r_energy(q12,r);
   dy = r12[1] - 0.001; r = std::sqrt(dx*dx + dy*dy + dz*dz); emy = Q_r_energy(q12,r);
   dy = r12[1];
   f12[1] -= (epy-emy)/0.002;

   dz = r12[2] + 0.001; r = std::sqrt(dx*dx + dy*dy + dz*dz); epz = Q_r_energy(q12,r);
   dz = r12[2] - 0.001; r = std::sqrt(dx*dx + dy*dy + dz*dz); emz = Q_r_energy(q12,r);
   dz = r12[2];
   f12[2] -= (epz-emz)/0.002;
}

double lj_r_energy(const double epsilon12, const double sigma12, const double r)
{
   double u = (sigma12/r);
   double u6  = u*u*u*u*u*u;
   double u12 = u6*u6;

   return (4.0*epsilon12*(u12-u6));
}

void lj_fd_force(const double epsilon12, const double sigma12, const double r12[], double f12[])
{
   double e0,epx,emx,epy,emy,epz,emz,r;

   double dx = r12[0];
   double dy = r12[1];
   double dz = r12[2];

   r = std::sqrt(dx*dx + dy*dy + dz*dz); 
   e0 = lj_r_energy(epsilon12,sigma12,r);

   dx = r12[0] + 0.001; r = std::sqrt(dx*dx + dy*dy + dz*dz); epx = lj_r_energy(epsilon12,sigma12,r);
   dx = r12[0] - 0.001; r = std::sqrt(dx*dx + dy*dy + dz*dz); emx = lj_r_energy(epsilon12,sigma12,r);
   dx = r12[0];
   f12[0] -= (epx-emx)/0.002;

   dy = r12[1] + 0.001; r = std::sqrt(dx*dx + dy*dy + dz*dz); epy = lj_r_energy(epsilon12,sigma12,r);
   dy = r12[1] - 0.001; r = std::sqrt(dx*dx + dy*dy + dz*dz); emy = lj_r_energy(epsilon12,sigma12,r);
   dy = r12[1];
   f12[1] -= (epy-emy)/0.002;

   dz = r12[2] + 0.001; r = std::sqrt(dx*dx + dy*dy + dz*dz); epz = lj_r_energy(epsilon12,sigma12,r);
   dz = r12[2] - 0.001; r = std::sqrt(dx*dx + dy*dy + dz*dz); emz = lj_r_energy(epsilon12,sigma12,r);
   dz = r12[2];
   f12[2] -= (epz-emz)/0.002;

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
                   const double epsilon[], const double sigma[], const double rion[], double fion[], double fion0[])
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

         double r12[3] = {x,y,z};
         double r21[3] = {-x,-y,-z};
         lj_fd_force(epsilon12,sigma12,r12,&fion0[3*qm]);
         lj_fd_force(epsilon12,sigma12,r21,&fion0[3*mm]);
      }
}

/****************************************************
 *                                                  *
 *        QMMM_electrostatic_Cmat                   *
 *                                                  *
 ****************************************************/
void QMMM_electrostatic_Cmat(const int nion_qm, const int nion, const double rion[], double Cmat[])
{
   memset(Cmat,0,nion_qm*nion_qm*sizeof(double));
   double *Cx = Cmat;
   double *Cy = &Cmat[nion_qm*nion_qm];
   double *Cz = &Cmat[2*nion_qm*nion_qm];

   for (auto k=0; k<nion_qm; ++k)
   for (auto j=0; j<nion_qm; ++j)
   {
      double tmpx = 0.0;
      double tmpy = 0.0;
      double tmpz = 0.0;
      for (auto m=nion_qm; m<nion; ++m)
      {
         double xkm = rion[3*k]   - rion[3*m];
         double ykm = rion[3*k+1] - rion[3*m+1];
         double zkm = rion[3*k+2] - rion[3*m+2];
         double rkm = std::sqrt(xkm*xkm + ykm*ykm + zkm*zkm);

         double xjm = rion[3*j]   - rion[3*m];
         double yjm = rion[3*j+1] - rion[3*m+1];
         double zjm = rion[3*j+2] - rion[3*m+2];
         double rjm = std::sqrt(xjm*xjm + yjm*yjm + zjm*zjm);
         tmpx += xkm/xjm* std::pow(rjm/rkm,3);
         tmpy += ykm/yjm* std::pow(rjm/rkm,3);
         tmpz += ykm/yjm* std::pow(rjm/rkm,3);
      }
      Cx[k+j*nion_qm] = tmpx;
      Cy[k+j*nion_qm] = tmpy;
      Cz[k+j*nion_qm] = tmpz;
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
                             const double qion[], const double rion[], double fion[], double fion0[]) 
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

         double r12[3] = {x,y,z};
         double r21[3] = {-x,-y,-z};
         Q_fd_force(qion[qm]*qion[mm],r12,&fion0[3*qm]);
         Q_fd_force(qion[qm]*qion[mm],r21,&fion0[3*mm]);
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





double get_QM_energy(MPI_Comm comm_world, const int nion_qm, const int nion, double rion[], double qion[])
{
   double fion[3*nion],uion[nion],Eqm;

   // zero potentials and forces 
   memset(uion,0,nion*sizeof(double));
   memset(fion,0,3*nion*sizeof(double));

   // calculate electrostatic potential on QM atoms, u(ii)
   QMMM_electrostatic_potential(nion_qm,nion,qion,rion,uion);

   // Eqm = QM energy and forces 
   int ierr = lammps_pspw_qmmm_minimizer(comm_world,rion,uion,fion,qion,&Eqm);

   return Eqm;
}

void get_QM_force(MPI_Comm comm_world, const int nion_qm, const int nion, double rion[], double qion[], double fion[])
{
   double ep,em;
   double dx = 0.001;

   for (auto ii=0; ii<nion; ++ii)
   {
      rion[3*ii] += 1*dx;   ep = get_QM_energy(comm_world, nion_qm, nion, rion, qion);
      rion[3*ii] -= 2*dx;   em = get_QM_energy(comm_world, nion_qm, nion, rion, qion);
      rion[3*ii] += 1*dx;
      fion[3*ii]   = -(ep-em)/(2.0*dx);

      rion[3*ii+1] += 1*dx; ep = get_QM_energy(comm_world, nion_qm, nion, rion, qion);
      rion[3*ii+1] -= 2*dx; em = get_QM_energy(comm_world, nion_qm, nion, rion, qion);
      rion[3*ii+1] += 1*dx;
      fion[3*ii+1] = -(ep-em)/(2.0*dx);

      rion[3*ii+2] += 1*dx; ep = get_QM_energy(comm_world, nion_qm, nion, rion, qion);
      rion[3*ii+2] -= 2*dx; em = get_QM_energy(comm_world, nion_qm, nion, rion, qion);
      rion[3*ii+2] += 1*dx;
      fion[3*ii+2] = -(ep-em)/(2.0*dx);
   }
}




int main(int argc, char* argv[])
{
   std::ofstream *xyzfile;
   std::ofstream *emotionfile;

   //std::string nwfilename = "w2.nw";
   //std::string nwfilename = "water0.nw";
   std::string nwfilename = "o1.nw";
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
      xyzfile->open("testooqmm.xyz", std::ios::app);

      emotionfile = new (std::nothrow) std::ofstream;
      emotionfile->open("testooqmm.emotion", std::ios::app);
   }


   int nion_qm = 1;
   int nion_mm = 1;
   int nion    = nion_qm + nion_mm;
   double rion0[3*nion],rion1[3*nion],rion2[3*nion];
   double uion[nion],qion[nion],fion[3*nion],fion1[3*nion],fion0[3*nion];
   double fionq1[3*nion],fionq0[3*nion];
   double epsilon[nion],sigma[nion],mass[nion],dti[nion];
   double E,KE,Eqm,Ecoul,ELJ;

   std::string symbol[nion] = {"O","O"};

   double pi = 4.0*std::atan(1.0);

   lammps_pspw_input(MPI_COMM_WORLD, nwfilename);


   // first water molecule - QM water
   rion1[0] =  0.161560*ANGTOBOHR; rion1[1] = -0.052912*ANGTOBOHR; rion1[2] =  0.033173*ANGTOBOHR;
   qion[0]  = 0.0;
   uion[0]  =  0.0000; 
   sigma[0] = 3.165558*ANGTOBOHR; epsilon[0] = 0.155394/(23.06*27.2116);
   mass[0] = 16.0*1822.89; 

   // O ion - MM 
   rion1[3] =  0.021259*ANGTOBOHR; rion1[4] =  0.506771*ANGTOBOHR; rion1[5] =  2.831278*ANGTOBOHR;
   qion[1]  = 0.0;
   uion[1]  = 0.0;
   sigma[1] = 3.165558*ANGTOBOHR; epsilon[1] = 0.155394/(23.06*27.2116);
   mass[1]  = 16.0*1822.89; 



   //lj_ion_parameters O 3.165558 0.155394
   //lj_ion_parameters H 0.700000 0.044000
   //lj_ion_parameters O^ 3.165558 0.155394a

   //shake units angstroms 1 2 3 cyclic 1.0 1.632993125 1.0
   double unita[9] = {20.0,  0.0,  0.0,
                       0.0, 20.0,  0.0,
                       0.0,  0.0, 20.0};
   double dt = 5.0;
   double h  = 1.0/(2.0*dt);
   for (auto ii=0; ii<nion; ++ii) dti[ii] = (dt*dt)/mass[ii];

   double cmat[3*nion_qm*nion_qm];
   QMMM_electrostatic_Cmat(nion_qm,nion,rion1,cmat);
   std::cout << "Cmatx = " << std::endl;
   for (auto i=0; i<nion_qm; ++i)
   {
      for (auto j=0; j<nion_qm; ++j)
         std::cout << cmat[i+j*nion_qm] << " ";
     std::cout << std::endl;
   }
   std::cout << std::endl << std::endl;
   std::cout << "Cmaty = " << std::endl;
   for (auto i=0; i<nion_qm; ++i)
   {
      for (auto j=0; j<nion_qm; ++j)
         std::cout << cmat[i+j*nion_qm+nion_qm*nion_qm] << " ";
     std::cout << std::endl;
   }
   std::cout << std::endl << std::endl;
   std::cout << "Cmatz = " << std::endl;
   for (auto i=0; i<nion_qm; ++i)
   {
      for (auto j=0; j<nion_qm; ++j)
         std::cout << cmat[i+j*nion_qm+2*nion_qm*nion_qm] << " ";
     std::cout << std::endl;
   }
   std::cout << std::endl << std::endl;

   //memcpy(rion2,rion1,3*nion*sizeof(double));
   //memset(rion0,0,3*nion*sizeof(double));
   //shake_chain(3,mm_water_indx,3,1.0e-6,55,dsq,mass,unita,rion1,rion2); 

   //double EE = get_QM_energy(MPI_COMM_WORLD, nion_qm, nion, rion1, qion,fion0);

   // zero potentials and forces 
   memset(uion,0,nion*sizeof(double));
   memset(fion,0,3*nion*sizeof(double));

   // calculate electrostatic potential on QM atoms, u(ii)
   QMMM_electrostatic_potential(nion_qm,nion,qion,rion1,uion);

   // Eqm = QM energy and forces 
   ierr += lammps_pspw_qmmm_minimizer(MPI_COMM_WORLD,rion1,uion,fion,qion,&Eqm);

   //std::cout << "EE=" << EE << " " << Eqm << std::endl;      

   // Ecoul = QMQM Electrostatic energy and forces
   //Ecoul = QMQM_electrostatic_energy(nion_qm,qion,rion1);
   //QMQM_electrostatic_force(nion_qm,qion,rion1,fion);
   Ecoul = 0.0;

   // Ecoul += QMMM Electrostatic energy and forces
   double EAPC = 0.0;
   for (auto ii=0; ii<nion_qm; ++ii)
      EAPC += qion[ii]*uion[ii];
   double eqmma = QMMM_electrostatic_energy(nion_qm,nion,qion,rion1);
   Ecoul += EAPC;

   memset(fionq0,0,3*nion*sizeof(double));
   memset(fionq1,0,3*nion*sizeof(double));
   QMMM_electrostatic_force(nion_qm,nion,qion,rion1,fionq1,fionq0);
   for (auto i=0; i<3*nion; ++i) fion[i] += fionq1[i];
   if (taskid==MASTER) std::cout << "EAPC=" << EAPC << " " << eqmma << std::endl;

   // ELJ = QMMM Lenard-Jones energy and forces
   ELJ = 0.0;
   //ELJ = QMMM_LJ_energy(nion_qm,nion,epsilon,sigma,rion1);
   memset(fion0,0,3*nion*sizeof(double));
   memset(fion1,0,3*nion*sizeof(double));
   //QMMM_LJ_force(nion_qm,nion,epsilon,sigma,rion1,fion1,fion0);
   //for (auto i=0; i<3*nion; ++i) fion[i] += fion1[i];
   

   // Espring = MM spring energy and forces

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

      std::cout << "KE+energy=" << KE+Eqm+Ecoul+ELJ << " energy=" << Eqm+Ecoul+ELJ << " KE=" << KE << " Eqm=" << Eqm << " Ecoul=" << Ecoul << " ELJ=" << ELJ << std::endl;

      for (auto ii=0; ii<nion; ++ii)
         std::cout << "ii=" << ii << " fion: "  << fion[3*ii]   << " " << fion[3*ii+1]   << " " << fion[3*ii+2]   << " qion=" << qion[ii]   << std::endl;
      std::cout << std::endl;

      for (auto ii=0; ii<nion; ++ii)
         std::cout << "ii=" << ii << " fion1:"  << fion1[3*ii]   << " " << fion1[3*ii+1]   << " " << fion1[3*ii+2]   << " qion=" << qion[ii]   << std::endl;
      std::cout << std::endl;
      for (auto ii=0; ii<nion; ++ii)
         std::cout << "ii=" << ii << " fion0:"  << fion0[3*ii]   << " " << fion0[3*ii+1]   << " " << fion0[3*ii+2]   << " qion=" << qion[ii]   << std::endl;
      std::cout << std::endl;
      for (auto ii=0; ii<nion; ++ii)
         std::cout << "ii=" << ii << " fionq1:"  << fionq1[3*ii]   << " " << fionq1[3*ii+1]   << " " << fionq1[3*ii+2]   << " qion=" << qion[ii]   << std::endl;
      std::cout << std::endl;
      for (auto ii=0; ii<nion; ++ii)
         std::cout << "ii=" << ii << " fionq0:"  << fionq0[3*ii]   << " " << fionq0[3*ii+1]   << " " << fionq0[3*ii+2]   << " qion=" << qion[ii]   << std::endl;
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
      //Ecoul = QMQM_electrostatic_energy(nion_qm,qion,rion1);
      //QMQM_electrostatic_force(nion_qm,qion,rion1,fion);
      //std::cout << "EQMQM = " << Ecoul << std::endl;
      Ecoul = 0.0;

      // QMMM Electrostatic energy and forces
      double eqmm = 0.0;
      for (auto ii=0; ii<nion_qm; ++ii)
         eqmm += qion[ii]*uion[ii];
      double eqmm0 = QMMM_electrostatic_energy(nion_qm,nion,qion,rion1);
      if (taskid==MASTER) std::cout << "Ecoul_qmm = " << eqmm << " " << eqmm0 << std::endl;
     //Ecoul += eqmm0;
      Ecoul += eqmm;
      memset(fionq0,0,3*nion*sizeof(double));
      memset(fionq1,0,3*nion*sizeof(double));
      QMMM_electrostatic_force(nion_qm,nion,qion,rion1,fionq1,fionq0);
      for (auto i=0; i<3*nion; ++i) fion[i] += fionq1[i];

      // QMMM Electrostatic energy and forces
      ELJ = 0.0;
      //ELJ = QMMM_LJ_energy(nion_qm,nion,epsilon,sigma,rion1);
      memset(fion0,0,3*nion*sizeof(double));
      memset(fion1,0,3*nion*sizeof(double));
      //QMMM_LJ_force(nion_qm,nion,epsilon,sigma,rion1,fion1,fion0);
      //for (auto i=0; i<3*nion; ++i) fion[i] += fion1[i];


      // Espring = MM spring energy and forces

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

         std::cout << "@ KE+energy=" << it << " " << KE+Eqm+Ecoul+ELJ << " energy=" << Eqm+Ecoul+ELJ 
                   << " KE=" << KE << " Eqm=" << Eqm << " Ecoul=" << Ecoul << " ELJ=" << ELJ <<   std::endl;

         for (auto ii=0; ii<nion; ++ii)
            std::cout << "ii=" << ii << " fion: "  << fion[3*ii]   << " " << fion[3*ii+1]   << " " << fion[3*ii+2]   << " qion=" << qion[ii]   << std::endl;
         std::cout << std::endl;

         for (auto ii=0; ii<nion; ++ii)
            std::cout << "ii=" << ii << " fion1:"  << fion1[3*ii]   << " " << fion1[3*ii+1]   << " " << fion1[3*ii+2]   << " qion=" << qion[ii]   << std::endl;
         std::cout << std::endl;
         for (auto ii=0; ii<nion; ++ii)
            std::cout << "ii=" << ii << " fion0:"  << fion0[3*ii]   << " " << fion0[3*ii+1]   << " " << fion0[3*ii+2]   << " qion=" << qion[ii]   << std::endl;
         std::cout << std::endl;
         for (auto ii=0; ii<nion; ++ii)
            std::cout << "ii=" << ii << " fionq1:"  << fionq1[3*ii]   << " " << fionq1[3*ii+1]   << " " << fionq1[3*ii+2]   << " qion=" << qion[ii]   << std::endl;
         std::cout << std::endl;
         for (auto ii=0; ii<nion; ++ii)
            std::cout << "ii=" << ii << " fionq0:"  << fionq0[3*ii]   << " " << fionq0[3*ii+1]   << " " << fionq0[3*ii+2]   << " qion=" << qion[ii]   << std::endl;

         std::cout << std::endl << std::endl;
         for (auto ii=0; ii<nion_qm; ++ii)
         {
            double dx = rion1[3*ii]   - rion1[9];
            double dy = rion1[3*ii+1] - rion1[10];
            double dz = rion1[3*ii+2] - rion1[11];
            double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
            double ee = qion[ii]*qion[3]/dist;
            std::cout << "ii="<< ii << " dist=" << dist << " ee=" << ee << std::endl;
         }
         printxyz(xyzfile,nion,symbol,unita,rion1) ;
         printemotion(emotionfile,it*5.0,KE+Eqm+Ecoul+ELJ,Eqm+Ecoul+ELJ,KE,Eqm,Ecoul,ELJ);
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

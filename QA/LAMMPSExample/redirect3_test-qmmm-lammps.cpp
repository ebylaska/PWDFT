#include 	<cmath>
#include 	<iostream>
#include 	<fstream>
#include 	<cstdio>
#include 	<string>
#include	<unistd.h>

#include "mpi.h"

#define MASTER 0



extern int  lammps_pspw_qmmm_minimizer(MPI_Comm, double*, double*, double*, double*, double*);
extern void lammps_pspw_input(MPI_Comm, std::string&);


int main(int argc, char* argv[])
{
   //std::string nwfilename = "w2.nw";
   std::string nwfilename = "h2-redirect.nw";
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

   double E;
   double uion[2], qion[2], rion[3*2],fion[3*2];

 
   lammps_pspw_input(MPI_COMM_WORLD, nwfilename);

   uion[0] = -0.01;
   uion[1] = 0.04;
   //rion[0] = 0.0; rion[1] = 0.0; rion[2] = -2.17;
   //rion[3] = 0.0; rion[4] = 0.0; rion[5] =  2.17;
   rion[0] = 0.0; rion[1] = 0.0; rion[2] = -0.7;
   rion[3] = 0.0; rion[4] = 0.0; rion[5] =  0.7;

   ierr += lammps_pspw_qmmm_minimizer(MPI_COMM_WORLD,rion,uion,fion,qion,&E);

   if (taskid==MASTER)
   {
      std::cout << std::endl << std::endl;
      std::cout << "rion: " << rion[0] << " " << rion[1] << " " << rion[2] << std::endl 
                << "      " << rion[3] << " " << rion[4] << " " << rion[5] << std::endl;
      std::cout << "uion: " << uion[0] << " " << uion[1] << std::endl;
      std::cout << std::endl;

      std::cout << "energy=" << E << std::endl;
      std::cout << "fion:  " << fion[0] << " " << fion[1] << " " << fion[2] << std::endl 
                << "       " << fion[3] << " " << fion[4] << " " << fion[5] << std::endl;
      std::cout << "qion:  " << qion[0] << " " << qion[1] << std::endl;
   }

}

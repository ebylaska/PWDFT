#include 	<cmath>
#include 	<iostream>
#include 	<fstream>
#include 	<cstdio>
#include 	<string>
#include	<unistd.h>

#include "mpi.h"



extern int  lammps_pspw_qmmm_minimizer(MPI_Comm, double*, double*, double*, double*, double*);
extern void lammps_pspw_input(MPI_Comm, std::string&);


#define	MASTER 0

// MACROS for redirecting stdout - note REDIRECT_INIT has to be executed 
//   in the same scope of both REDIRECT_ON and REDIRECT_OFF.
#define REDIRECT_INIT	 int stdout_fd;fpos_t stdout_pos;

#define REDIRECT_ON(OUTPUTFILE)	{ fgetpos(stdout,&stdout_pos); stdout_fd = dup(fileno(stdout)); freopen(OUTPUTFILE,"w",stdout); }

#define REDIRECT_OFF()	{ fflush(stdout); dup2(stdout_fd,fileno(stdout)); close(stdout_fd); clearerr(stdout); fsetpos(stdout,&stdout_pos); }


int main(int argc, char* argv[])
{
   //std::string nwfilename = "w2.nw";
   std::string nwfilename = "h2-noprint.nw";
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
   REDIRECT_INIT;

   // turn on redirect output 
   if (taskid==MASTER) {
      std::cout << "turning on redirected output" << std::endl;
      REDIRECT_ON("lammp_pspw_input_debug.output");
      std::cout << "redirected output turned on" << std::endl;
   }
 
   lammps_pspw_input(MPI_COMM_WORLD, nwfilename);

   // turn off redirect output
   if (taskid==MASTER) { 
      std::cout << "turning off redirected output" << std::endl;
      REDIRECT_OFF();
      std::cout << "redirected output turned off" << std::endl;
   }


   uion[0] = -0.01;
   uion[1] = 0.04;
   //rion[0] = 0.0; rion[1] = 0.0; rion[2] = -2.17;
   //rion[3] = 0.0; rion[4] = 0.0; rion[5] =  2.17;
   rion[0] = 0.0; rion[1] = 0.0; rion[2] = -0.7;
   rion[3] = 0.0; rion[4] = 0.0; rion[5] =  0.7;

   // turn on redirect output 
   if (taskid==MASTER) {
      std::cout << "turning on redirected output" << std::endl;
      REDIRECT_ON("lammp_pspw_qmmm_minimizer_debug.output");
      std::cout << "redirected output turned on" << std::endl;
   }

   ierr += lammps_pspw_qmmm_minimizer(MPI_COMM_WORLD,rion,uion,fion,qion,&E);

   // turn off redirect output
   if (taskid==MASTER) { 
      std::cout << "turning off redirected output" << std::endl;
      REDIRECT_OFF();
      std::cout << "redirected output turned off" << std::endl;
   }

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

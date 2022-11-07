#include 	<cmath>
#include 	<iostream>
#include 	<fstream>
#include 	<cstdio>
#include 	<string>

#include "mpi.h"



extern int  lammps_pspw_qmmm_minimizer(MPI_Comm, double*, double*, double*, double*, double*);
extern void lammps_pspw_input(MPI_Comm, std::string&);


#define	MASTER 0


int main(int argc, char* argv[])
{
   //std::string nwfilename = "w2.nw";
   std::string nwfilename = "h2-noprint.nw";
   if (argc>1) nwfilename = argv[1];
   

   // Initialize MPI
   int ierr,np,taskid;

   ierr = MPI_Init(&argc,&argv);
   ierr += MPI_Comm_size(MPI_COMM_WORLD,&np);
   ierr += MPI_Comm_rank(MPI_COMM_WORLD,&taskid);


   // Intializing files and stream buffers
   std::ofstream input_fs,running_fs;
   if (taskid==MASTER) {
      input_fs.open("lammps_pspw_init_debug.out");
      running_fs.open("lammps_pspw_qmmm_minimizer_debug.out");
   }
   std::streambuf* stream_buffer_cout = std::cout.rdbuf();


   if (taskid==MASTER) {
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11111!!" << std::endl;
      std::cout << "WARNING WARNING WARNING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      std::cout << "THIS DOES NOT WORK BECAUSE printf's ARE BEEING USED!!!" << std::endl;
      std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11111!!" << std::endl << std::endl;
      std::cout << "Hello world" << std::endl;
      std::cout << "np=" << np << " taskid=" << taskid << std::endl;
      std::cout << "argc=" << argc << std::endl;
      std::cout << "nwfilename=" << nwfilename << std::endl;
   }

   double E;
   double uion[2], qion[2], rion[3*2],fion[3*2];

   // turn on redirect output 
   if (taskid==MASTER) {
      std::cout << "turning on redirected output" << std::endl;
      std::cout.rdbuf(input_fs.rdbuf());
      std::cout << "redirected output turned on" << std::endl;
   }
 
   lammps_pspw_input(MPI_COMM_WORLD, nwfilename);

   // turn off redirect output
   if (taskid==MASTER) { 
      std::cout << "turning off redirected output" << std::endl;
      std::cout.rdbuf(stream_buffer_cout);
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
      std::cout.rdbuf(running_fs.rdbuf());
      std::cout << "redirected output turned on" << std::endl;
   }

   ierr += lammps_pspw_qmmm_minimizer(MPI_COMM_WORLD,rion,uion,fion,qion,&E);

   // turn off redirect output
   if (taskid==MASTER) { 
      std::cout << "turning off redirected output" << std::endl;
      std::cout.rdbuf(stream_buffer_cout);
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

   // close files
   if (taskid==MASTER) { input_fs.close(); running_fs.close(); }

}

#include 	<cmath>
#include 	<iostream>
#include 	<fstream>
#include 	<cstdio>
#include 	<string>

#include "mpi.h"

using namespace std;

//extern int dummy_pwdft();


extern int  lammps_pspw_minimizer(MPI_Comm, double*, double*, double*, double*, double*);
extern void lammps_pspw_input(MPI_Comm, std::string&);



int main(int argc, char* argv[])
{


   int ierr,np,taskid;

// Initialize MPI
   ierr = MPI_Init(&argc,&argv);
   ierr += MPI_Comm_size(MPI_COMM_WORLD,&np);
   ierr += MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

   if (taskid==0) {
      std::cout << "Hello world" << std::endl;
      std::cout << "np=" << np << " taskid=" << taskid << std::endl;
   }

   std::string nwfilename = "w2.nw";
 
   double E;
   double uion[2], qion[2], rion[3*2],fion[3*2];
   lammps_pspw_input(MPI_COMM_WORLD, nwfilename);

   uion[0] = 0.0;
   uion[1] = 0.0;
   rion[0] = 0.0; rion[1] = 0.0; rion[2] = -2.17;
   rion[3] = 0.0; rion[4] = 0.0; rion[5] =  2.17;
   ierr += lammps_pspw_minimizer(MPI_COMM_WORLD,rion,uion,rion,qion,&E);

   ierr += MPI_Finalize();
}

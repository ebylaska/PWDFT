// A simple program that computes the square root of a number
#include <cmath>
//#include <cstdlib>
#include <iostream>
#include        <cstdio>
#include        <stdio.h>
#include <string>
#include "NwpwConfig.h"
#include "mpi.h"


//extern int cpsd(int argc, char *argv[]);
extern int cpsd(MPI_Comm comm_world0);

int main(int argc, char* argv[])
{
  int ierr=0;
  int taskid,np;
  int MASTER=0;

  // Initialize MPI
  ierr = MPI_Init(&argc,&argv);
  ierr += MPI_Comm_size(MPI_COMM_WORLD,&np);
  ierr += MPI_Comm_rank(MPI_COMM_WORLD,&taskid);


  if (taskid==0)
  {
     std::cout << argv[0] << " (NWChemEx) - Version " << Nwpw_VERSION_MAJOR << "."
               << Nwpw_VERSION_MINOR << std::endl;
  }


  //int ijk = cpsd(argc,argv);
  ierr += cpsd(MPI_COMM_WORLD);


  // Finalize MPI
  ierr += MPI_Finalize();


  return ierr;
}

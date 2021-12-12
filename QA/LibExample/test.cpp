#include 	<cmath>
#include 	<iostream>
#include 	<fstream>
#include 	<cstdio>
#include 	<string>

#include "mpi.h"

using namespace std;

//extern int dummy_pwdft();

namespace pwdft {
using namespace pwdft;

extern char *util_date();
extern void seconds(double *);
extern int pspw_geovib(MPI_Comm, std::string&);
extern int pspw_minimizer(MPI_Comm, std::string&);

}


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
      std::cout <<  "dummy=" << pwdft::util_date() << std::endl;
   }

   string line,nwinput;
   ifstream nwfile("w2-pspw.json");
   if (nwfile.good())
   {
      while (getline(nwfile,line))
         nwinput += line + "\n";
   }
   nwfile.close();


   if (taskid==0) {
      std::cout << "rtdbstr=" << nwinput << std::endl;
   }

   ierr += pwdft::pspw_minimizer(MPI_COMM_WORLD,nwinput);

   if (taskid==0) {
      std::cout << "output rtdbstr=" << nwinput << std::endl;
   }


   ierr += MPI_Finalize();
}

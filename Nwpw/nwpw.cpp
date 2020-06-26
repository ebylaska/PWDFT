
////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                //
//                              NWChemEx PWDFT NWPW  (MPI code)                                   //
//                                                                                                //
//    This is a developing PW DFT parallel code for NWChemEX. This code uses pseudopotentials     //
// and plane-wave basis sets to perform Density Functional Theory calculations.  A significant    //
// advantage of this code is its ability to simulate dynamics on a ground state potential surface //
// directly at run-time using the Car-Parrinello algorithm. This method's efficiency and accuracy //
// make it a desirable first principles method of simulation in the study of complex molecular,   //
// liquid, and solid state systems. Applications for this first principles method include the     //
// calculation of free energies, search for global minima, explicit simulation of solvated        //
// molecules, and simulations of complex vibrational modes that cannot be described within the    //
// harmonic approximation.                                                                        //
//                                                                                                //
// The code is a collection of two methods.                                                       //
//                                                                                                //
// PSPW - (PSeudopotential Plane-Wave) A gamma point code for calculating molecules, liquids,     //
//        crystals, and surfaces.                                                                 //
// Band - A band structure code for calculating crystals and surfaces with small band gaps (e.g.  //
//        semi-conductors and metals).                                                            //
//                                                                                                //
// The PSPW and Band  methdods can be used to compute the energy and optimize the geometry. Both  //
// the methods can also be used to find saddle points, and compute numerical second derivatives.  //
//                                                                                                //
////////////////////////////////////////////////////////////////////////////////////////////////////

#include <cmath>
//#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <stdio.h>
#include <string>

#include "NwpwConfig.h"
#include "mpi.h"

#include "util_date.hpp"
#include "parse_pwdft.hpp"

#include "nwpw.hpp"

using namespace std;


int main(int argc, char* argv[])
{
  int ierr=0;
  int taskid,np;
  int MASTER=0;
  string line,nwinput,nwfilename;

  // Initialize MPI
  ierr = MPI_Init(&argc,&argv);
  ierr += MPI_Comm_size(MPI_COMM_WORLD,&np);
  ierr += MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

  bool oprint = (taskid==MASTER);

  /* set the pseudopotential library directory */
  const char *nwpw_libraryps = Nwpw_LIBRARYPS_Default;
  if (const char* libraryps0 = std::getenv("NWPW_LIBRARY"))
     nwpw_libraryps = libraryps0;
  else if (const char*  libraryps0 = std::getenv("NWCHEM_NWPW_LIBRARY"))
     nwpw_libraryps = libraryps0;


  if (oprint) 
  {
     std::cout << argv[0] << " (NWChemEx) - Version " << Nwpw_VERSION_MAJOR << "."
               << Nwpw_VERSION_MINOR << std::endl;
     if (argc>1) 
     {
        nwfilename = argv[1];
        nwinput = "";
        ifstream nwfile(argv[1]);
        if (nwfile.good())
        {
           while (getline(nwfile,line))
              nwinput += line + "\n";
        }
        nwfile.close();
     }
     else
     {
        nwfilename = "stdin";
        while (getline(std::cin,line))
              nwinput += line + "\n";
     }
     std::cout << std::endl;
     std::cout << "============================== echo of input deck ==============================\n";
     std::cout << nwinput;
     std::cout << "================================================================================\n\n";
     std::cout << "              NorthwestEx Computational Chemistry Package 1.0.0\n";
     std::cout << "           --------------------------------------------------------\n\n";
     std::cout << "                  Pacific Northwest National Laboratory\n";
     std::cout << "                           Richland, WA 99354\n\n";
     std::cout << "                         Copyright (c) 2020\n";
     std::cout << "                  Pacific Northwest National Laboratory\n";
     std::cout << "                       Battelle Memorial Institute\n\n";

     std::cout << "        NWChemEx is an open-source computational chemistry package\n";
     std::cout << "                   distributed under the terms of the\n";
     std::cout << "                 Educational Community License (ECL) 2.0\n";
     std::cout << "        A copy of the license is included with this distribution\n";
     std::cout << "                         in the LICENSE.TXT file\n\n";
     std::cout << "                             ACKNOWLEDGMENT\n";
     std::cout << "                             --------------\n\n";
     std::cout << "       This software and its documentation were developed at the\n";
     std::cout << "       Pacific Northwest National Laboratory, a multiprogram\n";
     std::cout << "       national laboratory, operated for the U.S. Department of Energy\n";
     std::cout << "       by Battelle under Contract Number DE-AC05-76RL01830. Support\n";
     std::cout << "       for this work was provided by the Department of Energy \n";
     std::cout << "       Office of Advanced Scientific Computing and the Office of Basic\n";
     std::cout << "       Energy Sciences.\n\n";
     std::cout << "       Job information\n";
     std::cout << "       ---------------\n";
     std::cout << "       program          = pwdft (NWChemEx)\n";
     std::cout << "       build configured = " << Nwpw_COMPILE_TIMESTAMP << std::endl;
     std::cout << "       source           = " << Nwpw_TOP << std::endl;
     std::cout << "       version          = " << Nwpw_VERSION_MAJOR << "." << Nwpw_VERSION_MINOR << std::endl;
     std::cout << "       psp libraries    = " << nwpw_libraryps << std::endl << std::endl;
     std::cout << "       date             = " << util_date() << std::endl;
     std::cout << "       nproc            = " << np << std::endl;
     std::cout << "       input            = " << nwfilename << std::endl << std::endl;
     std::cout << std::endl << std::endl;

  }
 
  // Broadcast nwinput across MPI tasks 
  if (np>1)
  {
      int nwinput_size = nwinput.size();
      MPI_Bcast(&nwinput_size,1,MPI_INT,MASTER,MPI_COMM_WORLD);
      if (taskid != MASTER)
         nwinput.resize(nwinput_size);
      MPI_Bcast(const_cast<char*>(nwinput.data()),nwinput_size,MPI_CHAR,MASTER,MPI_COMM_WORLD);
  }

  string rtdbstr  = parse_nwinput(nwinput);
  int task = parse_task(rtdbstr);
  if (oprint) std::cout << "rtdbstr=" << rtdbstr << std::endl;
  if (oprint) std::cout << "task0=" << task << std::endl;
  while (task>0)
  {
     if (task==5) ierr += cpsd(MPI_COMM_WORLD,rtdbstr);
     if (task==6) ierr += cpmd(MPI_COMM_WORLD,rtdbstr);
    
     rtdbstr = parse_rtdbstring(rtdbstr);
     task    = parse_task(rtdbstr);
     if (oprint) std::cout << "rtdbstr=" << rtdbstr << std::endl;
     if (oprint) std::cout << "task =" << task << std::endl;
  }

  //int ijk = cpsd(argc,argv);
  //ierr += cpsd(MPI_COMM_WORLD,rtdbstr);


  if (taskid==MASTER) parse_write(rtdbstr);

  // Finalize MPI
  ierr += MPI_Finalize();


  return ierr;
}

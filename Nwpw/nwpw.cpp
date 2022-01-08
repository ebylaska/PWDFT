
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
#include <string>

#include "NwpwConfig.h"
#include "mpi.h"

#include "util_date.hpp"
#include "parse_pwdft.hpp"
#include "psp_library.hpp"

#include "nwpw.hpp"

#include "json.hpp"
using json = nlohmann::json;

using namespace std;
using namespace pwdft;

static std::string fortran_rtdbstring;

extern "C" void pspw_fortran_minimizer_(int *comm_world, double *rion, double *uion, double *E, double *fion, double *qion)
{
   std::string line,nwinput;

   std::cout << "HELLO from pspw_minimizer_fortran_   comm= " << *comm_world <<  std::endl;
   std::cout << " " << MPI_COMM_WORLD << std::endl;

   auto fortran_rtdbjson =  json::parse(fortran_rtdbstring);

   int nion = fortran_rtdbjson["geometries"]["geometry"]["nion"];

   fortran_rtdbjson["geometries"]["geometry"]["coords"] = std::vector<double>(rion,&rion[3*nion]);
   fortran_rtdbjson["nwpw"]["apc"]["u"] = std::vector<double>(uion,&uion[nion]);
   fortran_rtdbstring    = fortran_rtdbjson.dump();

   
   std::cout << "input fortran_rtdbstring= " << fortran_rtdbstring << std::endl;;
   int  ierr = pwdft::pspw_minimizer(MPI_COMM_WORLD, fortran_rtdbstring);
   std::cout << "output fortran_rtdbstring= " << fortran_rtdbstring << std::endl;;


   fortran_rtdbjson =  json::parse(fortran_rtdbstring);
   *E = fortran_rtdbjson["pspw"]["energy"];

   std::vector<double> v = fortran_rtdbjson["pspw"]["fion"];
   std::copy(v.begin(),v.end(), fion);

   std::vector<double> vv = fortran_rtdbjson["nwpw"]["apc"]["q"];
   std::copy(vv.begin(),vv.end(), qion);
}


extern "C" void pspw_fortran_input_(char *filename, int *flen)
{
   int taskid,np,ierr,nwinput_size;;
   int MASTER=0;
   ierr = MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
   ierr = MPI_Comm_size(MPI_COMM_WORLD,&np);

   std::cout << "filename=" << filename << " len=" << *flen << std::endl;

   std::string nwinput;

   MPI_Barrier(MPI_COMM_WORLD);
   if (taskid==MASTER)
   {
      std::string line;
      std::string nwfilename(filename);
      //std::string nwfilename("w2.nw");

      std::cout << "nwfilename =" << nwfilename << std::endl;
      nwinput = "";
      std::ifstream nwfile(nwfilename);
      if (nwfile.good())
      {
         while (getline(nwfile,line))
            nwinput += line + "\n";
      }
      nwfile.close();
      std::cout << "nwinput =" << nwinput << std::endl;;

      nwinput_size = nwinput.size();
   }
   std::cout << "taskid=" << taskid << std::endl;
   MPI_Barrier(MPI_COMM_WORLD);
   std::cout << "new taskid=" << taskid << " np=" << np << std::endl;

   // Broadcast nwinput across MPI tasks 
   if (np>1)
   {
      MPI_Bcast(&nwinput_size,1,MPI_INT,MASTER,MPI_COMM_WORLD);
      if (taskid != MASTER)
         nwinput.resize(nwinput_size);
      MPI_Bcast(const_cast<char*>(nwinput.data()),nwinput_size,MPI_CHAR,MASTER,MPI_COMM_WORLD);
   }

   MPI_Barrier(MPI_COMM_WORLD);
   fortran_rtdbstring = pwdft::parse_nwinput(nwinput);

}



static std::string lammps_rtdbstring;

extern int lammps_pspw_minimizer(MPI_Comm comm_world, double *rion, double *uion, double *fion, double *qion, double *E)
{ 
   std::string line,nwinput;

   auto lammps_rtdbjson = json::parse(lammps_rtdbstring);
   int nion             = lammps_rtdbjson["geometries"]["geometry"]["nion"];

   lammps_rtdbjson["geometries"]["geometry"]["coords"] = std::vector<double>(rion,&rion[3*nion]);
   lammps_rtdbjson["nwpw"]["apc"]["u"] = std::vector<double>(uion,&uion[nion]);
   lammps_rtdbstring    = lammps_rtdbjson.dump();

   int  ierr = pwdft::pspw_minimizer(comm_world, lammps_rtdbstring);

   lammps_rtdbjson =  json::parse(lammps_rtdbstring);
   *E = lammps_rtdbjson["pspw"]["energy"];

   std::vector<double> v = lammps_rtdbjson["pspw"]["fion"];
   std::copy(v.begin(),v.end(), fion);

   std::vector<double> vv = lammps_rtdbjson["nwpw"]["apc"]["q"];
   std::copy(vv.begin(),vv.end(), qion);

   return ierr;
}


extern void lammps_pspw_input(MPI_Comm comm_world, std::string& nwfilename)
{
   int taskid,np,ierr,nwinput_size;
   int MASTER=0;
   ierr = MPI_Comm_rank(comm_world,&taskid);
   ierr = MPI_Comm_size(comm_world,&np);

   std::string nwinput;

   MPI_Barrier(comm_world);
   if (taskid==MASTER)
   {  
      std::string line;
      
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
   MPI_Barrier(comm_world);

   // Broadcast nwinput across MPI tasks 
   if (np>1)
   {
      MPI_Bcast(&nwinput_size,1,MPI_INT,MASTER,comm_world);
      if (taskid != MASTER)
         nwinput.resize(nwinput_size);
      MPI_Bcast(const_cast<char*>(nwinput.data()),nwinput_size,MPI_CHAR,MASTER,comm_world);
   }

   MPI_Barrier(comm_world);
   lammps_rtdbstring = pwdft::parse_nwinput(nwinput);
}




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

  /* Fetch  the pseudopotential library directory */
  const char *nwpw_libraryps = Nwpw_LIBRARYPS_Default;
  if (const char *libraryps0 = std::getenv("NWPW_LIBRARY"))
     nwpw_libraryps = libraryps0;
  else if (const char *libraryps0 = std::getenv("NWCHEM_NWPW_LIBRARY"))
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
     std::cout << "       program               = pwdft (NWChemEx)\n";
     std::cout << "       build configured      = " << Nwpw_COMPILE_TIMESTAMP << std::endl;
     std::cout << "       source                = " << Nwpw_TOP << std::endl;
     std::cout << "       version               = " << Nwpw_VERSION_MAJOR << "." << Nwpw_VERSION_MINOR << std::endl;
     //std::cout << "       psp libraries    = " << nwpw_libraryps << std::endl << std::endl;
     std::cout << "       default psp libraries = " << psp_library("").nwpw_libraryps_dir << std::endl << std::endl;
     std::cout << "       date                  = " << util_date() << std::endl;
     std::cout << "       nproc                 = " << np << std::endl;
     std::cout << "       input                 = " << nwfilename << std::endl << std::endl;
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

  MPI_Barrier(MPI_COMM_WORLD);
  string rtdbstr  = parse_nwinput(nwinput);
  int task = parse_task(rtdbstr);
  MPI_Barrier(MPI_COMM_WORLD);

  if (oprint) std::cout << "First rtdbstr=" << rtdbstr << std::endl;
  if (oprint) std::cout << "First task=" << task << std::endl << std::endl;


  // Initialize wavefunction
  {
     std::string input_wavefunction_filename = parse_input_wavefunction_filename(rtdbstr);
     int wfound=0; if (taskid==MASTER) { ifstream wfile(input_wavefunction_filename); if (wfile.good()) wfound=1; wfile.close(); }
     MPI_Bcast(&wfound,1,MPI_INT,MASTER,MPI_COMM_WORLD);
     if (!wfound) 
     {
        auto lowlevel_rtdbstrs =  parse_gen_lowlevel_rtdbstrs(rtdbstr);
        for (const auto & elem: lowlevel_rtdbstrs) {
           if (oprint) std::cout << std::endl << "Running staged energy optimization - lowlevel_rtdbstr = " << elem << std::endl << std::endl;
           std::string dum_rtdbstr  = elem;
           ierr += pspw_minimizer(MPI_COMM_WORLD,dum_rtdbstr);
        }
     }
  }

  // Tasks
  while (task>0)
  {
     /* Energy or Gradient task */
     if ((task==1) || (task==2))
     {
        MPI_Barrier(MPI_COMM_WORLD);
        ierr += pwdft::pspw_minimizer(MPI_COMM_WORLD,rtdbstr);
     }

     /* Optimize task */
     if (task==3) 
     {
        if (oprint) std::cout << std::endl << "Running geometry optimization calculation - rtdbstr = " << rtdbstr << std::endl << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        ierr += pwdft::pspw_geovib(MPI_COMM_WORLD,rtdbstr);
     }

     /* Frequency task */
     if (task==4) 
     {
        if (oprint) std::cout << std::endl << "Running frequency calculation - rtdbstr = " << rtdbstr << std::endl << std::endl;
     }

     /* Steepest descent task */
     if (task==5) 
     { 
        ierr += pwdft::cpsd(MPI_COMM_WORLD,rtdbstr); /* Steepest_Descent task */
     }

     /* Car-Parrinello task */
     if (task==6) 
     {
        ierr += pwdft::cpmd(MPI_COMM_WORLD,rtdbstr); /* Car-Parrinello task */
     }
     MPI_Barrier(MPI_COMM_WORLD);


     // parse json string 
     rtdbstr = parse_rtdbstring(rtdbstr);
     MPI_Barrier(MPI_COMM_WORLD);

     // Find next task 
     task    = parse_task(rtdbstr);
     if (oprint) std::cout << std::endl << "Next rtdbstr=" << rtdbstr << std::endl;
     if (oprint) std::cout << "Next task =" << task << std::endl << std::endl;
     MPI_Barrier(MPI_COMM_WORLD);
  }

  //int ijk = cpsd(argc,argv);

  // DEBUG
  //MPI_Barrier(MPI_COMM_WORLD);
  //ierr += cpsd_debug(MPI_COMM_WORLD,rtdbstr);
  //MPI_Barrier(MPI_COMM_WORLD);

  if (taskid==MASTER) parse_write(rtdbstr);
  MPI_Barrier(MPI_COMM_WORLD);


  // Finalize MPI
  ierr = MPI_Finalize();



  return ierr;
}


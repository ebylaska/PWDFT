
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
#include <unistd.h>
#include <iomanip>

#include "NwpwConfig.h"
#include "mpi.h"

#include "util_date.hpp"
#include "parse_pwdft.hpp"
#include "psp_library.hpp"
#include "ion_ion.hpp"

#include "nwpw.hpp"

#include "json.hpp"
using json = nlohmann::json;

using namespace pwdft;

#define MASTER 0


// *************************************************************************
//
// Fortran Interface Routines
//
//
// *************************************************************************

// MACROS for redirecting stdout - note REDIRECT_INIT has to be executed
//   in the same scope of both REDIRECT_ON and REDIRECT_OFF.
#define REDIRECT_INIT    int stdout_fd;fpos_t stdout_pos;

#define REDIRECT_ON(OUTPUTFILE) { fgetpos(stdout,&stdout_pos); stdout_fd = dup(fileno(stdout)); freopen(OUTPUTFILE,"w",stdout); }

#define REDIRECT_OFF()  { fflush(stdout); dup2(stdout_fd,fileno(stdout)); close(stdout_fd); clearerr(stdout); fsetpos(stdout,&stdout_pos); }





static std::string fortran_rtdbstring;

extern "C" void pspw_fortran_minimizer_(MPI_Fint *fcomm_world, double *rion, double *uion, double *E, double *fion, double *qion)
{
   MPI_Comm comm_world = MPI_Comm_f2c(*fcomm_world);

   std::string line,nwinput;

   //std::cout << "HELLO from pspw_minimizer_fortran_   comm= " <<  comm_world <<  std::endl;
   //std::cout << " MPI_COMM_WORLD= " << MPI_COMM_WORLD << std::endl;

   auto fortran_rtdbjson =  json::parse(fortran_rtdbstring);

   int nion = fortran_rtdbjson["geometries"]["geometry"]["nion"];

   fortran_rtdbjson["geometries"]["geometry"]["coords"] = std::vector<double>(rion,&rion[3*nion]);
   fortran_rtdbjson["nwpw"]["apc"]["u"] = std::vector<double>(uion,&uion[nion]);
   fortran_rtdbjson["current_task"]  = "gradient";

   fortran_rtdbstring    = fortran_rtdbjson.dump();

   
   //std::cout << "input fortran_rtdbstring= " << fortran_rtdbstring << std::endl;;

   int  ierr = pwdft::pspw_minimizer(comm_world, fortran_rtdbstring, std::cout);

   //std::cout << "output fortran_rtdbstring= " << fortran_rtdbstring << std::endl;;


   fortran_rtdbjson =  json::parse(fortran_rtdbstring);
   *E = fortran_rtdbjson["pspw"]["energy"];

   std::vector<double> v = fortran_rtdbjson["pspw"]["fion"];
   std::copy(v.begin(),v.end(), fion);

   std::vector<double> vv = fortran_rtdbjson["nwpw"]["apc"]["q"];
   std::copy(vv.begin(),vv.end(), qion);
}


extern "C" void pspw_fortran_input_(MPI_Fint *fcomm_world, char *filename, int *flen)
{
   MPI_Comm comm_world = MPI_Comm_f2c(*fcomm_world);
   int taskid,np,ierr,nwinput_size;;
   ierr = MPI_Comm_rank(comm_world,&taskid);
   ierr = MPI_Comm_size(comm_world,&np);

   //std::cout << "taskid=" << taskid << " np=" << np << std::endl;
   //std::cout << "filename=" << filename << " len=" << *flen << std::endl;

   std::string nwinput;

   MPI_Barrier(comm_world);
   if (taskid==MASTER)
   {
      std::string line;
      char tfilename[*flen]; for (auto i=0; i<*flen; ++i) tfilename[i] = filename[i];
      std::string nwfilename(tfilename);
      //std::string nwfilename("w2.nw");

      // preappend print option to nwinput, options=none,off,low,medium,high,debug
      std::cout << "nwfilename =" << nwfilename << std::endl;
      //nwinput = "\nprint off\n\n";
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
   //std::cout << "taskid=" << taskid << std::endl;
   MPI_Barrier(comm_world);
   //std::cout << "new taskid=" << taskid << " np=" << np << std::endl;

   // Broadcast nwinput across MPI tasks 
   if (np>1)
   {
      MPI_Bcast(&nwinput_size,1,MPI_INT,MASTER,comm_world);
      if (taskid != MASTER)
         nwinput.resize(nwinput_size);
      MPI_Bcast(const_cast<char*>(nwinput.data()),nwinput_size,MPI_CHAR,MASTER,comm_world);
   }

   MPI_Barrier(comm_world);
   fortran_rtdbstring = pwdft::parse_nwinput(nwinput);

   // Initialize wavefunction if it doesn't exist
   if (parse_initialize_wvfnc(fortran_rtdbstring,true))
   {
     bool wvfnc_initialize = parse_initialize_wvfnc(fortran_rtdbstring,false);
     if (wvfnc_initialize) fortran_rtdbstring = parse_initialize_wvfnc_set(fortran_rtdbstring,false);

     int taskid;
     int ierr = MPI_Comm_rank(comm_world,&taskid);
     bool oprint = (taskid==MASTER) && false;
     std::string input_wavefunction_filename = pwdft::parse_input_wavefunction_filename(fortran_rtdbstring);
     int wfound=0; if (taskid==MASTER) { std::ifstream wfile(input_wavefunction_filename); if (wfile.good()) wfound=1; wfile.close(); }
     MPI_Bcast(&wfound,1,MPI_INT,MASTER,comm_world);

     if ((!wfound) || (wvfnc_initialize))
     {
        auto lowlevel_rtdbstrs = pwdft::parse_gen_lowlevel_rtdbstrs(fortran_rtdbstring);
        for (const auto & elem: lowlevel_rtdbstrs) {
           std::string dum_rtdbstr  = elem;
           if (wvfnc_initialize) {
              dum_rtdbstr = parse_initialize_wvfnc_set(dum_rtdbstr,true);
              wvfnc_initialize = false;
           }
           if (oprint) std::cout << std::endl << "Running staged energy optimization - lowlevel_rtdbstr = " << dum_rtdbstr << std::endl << std::endl;
           ierr += pwdft::pspw_minimizer(comm_world,dum_rtdbstr,std::cout);
        }
     }
   }


}



// *************************************************************************
//
// LAMMPS Interface Routines - output sent to io stream
//
// *************************************************************************

static std::string lammps_output_filename;
static std::string lammps_rtdbstring;
static bool printqmmm = false;


extern int lammps_pspw_aimd_minimizer(MPI_Comm comm_world, double *rion, double *fion, double *E, std::ostream& coutput)
{
   int ierr;
   auto lammps_rtdbjson = json::parse(lammps_rtdbstring);
   int nion             = lammps_rtdbjson["geometries"]["geometry"]["nion"];
   lammps_rtdbjson["geometries"]["geometry"]["coords"] = std::vector<double>(rion,&rion[3*nion]);

   lammps_rtdbstring = lammps_rtdbjson.dump();

   ierr = pwdft::pspw_minimizer(comm_world,lammps_rtdbstring,coutput);

   lammps_rtdbjson   =  json::parse(lammps_rtdbstring);

   // fetch output - energy and forces
   double ee       = lammps_rtdbjson["pspw"]["energy"];
   *E = ee;

   std::vector<double> v = lammps_rtdbjson["pspw"]["fion"];
   std::copy(v.begin(),v.end(), fion);

   return ierr;
}


extern int lammps_pspw_qmmm_minimizer(MPI_Comm comm_world, double *rion, double *uion, double *fion, double *qion, double *E, 
                                      bool removeqmmmcoulomb, bool removeqmqmcoulomb, std::ostream& coutput)
{ 
   int ierr;
   //int taskid,np,ierr;
   //ierr = MPI_Comm_rank(comm_world,&taskid);
   //ierr = MPI_Comm_size(comm_world,&np);

   auto lammps_rtdbjson = json::parse(lammps_rtdbstring);

   int nion             = lammps_rtdbjson["geometries"]["geometry"]["nion"];
   lammps_rtdbjson["geometries"]["geometry"]["coords"] = std::vector<double>(rion,&rion[3*nion]);

   if (lammps_rtdbjson["nwpw"].is_null()) { json nwpw; lammps_rtdbjson["nwpw"] = nwpw; }
   if (lammps_rtdbjson["nwpw"]["apc"].is_null()) { json apc; lammps_rtdbjson["nwpw"]["apc"] = apc; }
   lammps_rtdbjson["nwpw"]["apc"]["on"] = true;
   lammps_rtdbjson["nwpw"]["apc"]["u"]  = std::vector<double>(uion,&uion[nion]);
   lammps_rtdbjson["current_task"]  = "gradient";

   lammps_rtdbstring = lammps_rtdbjson.dump();

   // run minimizer
   //if (printqmmm) std::cout << "lammps_rtdbstring = " << lammps_rtdbstring << std::endl;

   ierr = pwdft::pspw_minimizer(comm_world,lammps_rtdbstring,coutput);


   lammps_rtdbjson   =  json::parse(lammps_rtdbstring);

   // fetch output - energy, forces, and apc charges
   double ee   = lammps_rtdbjson["pspw"]["energy"];
   double eapc = lammps_rtdbjson["pspw"]["energies"][51];

   if (removeqmmmcoulomb)
      *E = (ee-eapc);
   else
      *E = (ee);

   std::vector<double> v = lammps_rtdbjson["pspw"]["fion"];
   std::copy(v.begin(),v.end(), fion);

   std::vector<double> vv = lammps_rtdbjson["nwpw"]["apc"]["q"];
   std::copy(vv.begin(),vv.end(), qion);

   //coutput << " pwdft EAPC=" << std::fixed << std::setw(15) << std::setprecision(11) << eapc << std::endl;


   // pre-remove qm/qm electrostatic interactions
   if (removeqmqmcoulomb)
   {
      double ecoul = pwdft::ion_ion_e(nion,qion,rion);
      //std::cout << " pwdft QMQM Ecoul=" << std::fixed << std::setw(15) << std::setprecision(11) << ecoul 
      //          << " qion=" << lammps_rtdbjson["nwpw"]["apc"]["q"]  
      //          << std::endl;
      *E -= ecoul;
      pwdft::ion_ion_m_f(nion,qion,rion,fion);
   }

   return ierr;
}


extern int lammps_pspw_qmmm_nominimizer(MPI_Comm comm_world, double *rion, double *uion, double *fion, double *qion, double *E, 
                                        bool removeqmmmcoulomb, bool removeqmqmcoulomb, std::ostream& coutput)
{
   int ierr;
   //int taskid,np,ierr;
   //ierr = MPI_Comm_rank(comm_world,&taskid);
   //ierr = MPI_Comm_size(comm_world,&np);

   auto lammps_rtdbjson = json::parse(lammps_rtdbstring);

   int nion             = lammps_rtdbjson["geometries"]["geometry"]["nion"];
   lammps_rtdbjson["geometries"]["geometry"]["coords"] = std::vector<double>(rion,&rion[3*nion]);

   if (lammps_rtdbjson["nwpw"].is_null()) { json nwpw; lammps_rtdbjson["nwpw"] = nwpw; }
   if (lammps_rtdbjson["nwpw"]["apc"].is_null()) { json apc; lammps_rtdbjson["nwpw"]["apc"] = apc; }
   lammps_rtdbjson["nwpw"]["apc"]["on"] = true;
   lammps_rtdbjson["nwpw"]["apc"]["u"]  = std::vector<double>(uion,&uion[nion]);
   lammps_rtdbjson["current_task"]  = "noit_gradient";

   lammps_rtdbstring = lammps_rtdbjson.dump();

   // run minimizer
   //if (printqmmm) std::cout << "qmmm_nominimizer:lammps_rtdbstring = " << lammps_rtdbstring << std::endl;

   ierr = pwdft::pspw_minimizer(comm_world,lammps_rtdbstring,coutput);


   lammps_rtdbjson   =  json::parse(lammps_rtdbstring);

   // fetch output - energy, forces, and apc charges
   double ee   = lammps_rtdbjson["pspw"]["energy"];
   double eapc = lammps_rtdbjson["pspw"]["energies"][51];

   if (removeqmmmcoulomb)
      *E = (ee-eapc);
   else
      *E = (ee);

   std::vector<double> v = lammps_rtdbjson["pspw"]["fion"];
   std::copy(v.begin(),v.end(), fion);

   //std::vector<double> vv = lammps_rtdbjson["nwpw"]["apc"]["q"];
   //std::copy(vv.begin(),vv.end(), qion);

   //coutput << " pwdft EAPC=" << std::fixed << std::setw(15) << std::setprecision(11) << eapc << std::endl;


   // pre-remove qm/qm electrostatic interactions
   if (removeqmqmcoulomb)
   {
      double ecoul = pwdft::ion_ion_e(nion,qion,rion);
      //std::cout << " pwdft QMQM Ecoul=" << std::fixed << std::setw(15) << std::setprecision(11) << ecoul << " qion=" << lammps_rtdbjson["nwpw"]["apc"]["q"] <<  std::endl;
      *E -= ecoul;
      pwdft::ion_ion_m_f(nion,qion,rion,fion);
   }

   return ierr;
}





extern void lammps_pspw_input(MPI_Comm comm_world, std::string& nwfilename, std::ostream& coutput)
{
   int taskid,np,ierr,nwinput_size;
   ierr = MPI_Comm_rank(comm_world,&taskid);
   ierr = MPI_Comm_size(comm_world,&np);


   std::string nwinput;

   MPI_Barrier(comm_world);
   if (taskid==MASTER)
   {  
      printqmmm = true;
      std::string line;
      
      // prepend print option to nwinput, options=none,off,low,medium,high,debug
      //nwinput = "\nprint off\n\n"; 
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


   auto lammps_rtdbjson = json::parse(lammps_rtdbstring);


   // Initialize wavefunction if it doesn't exist or if initialize_wavefunction set
   if (parse_initialize_wvfnc(lammps_rtdbstring,true))
   {
      bool wvfnc_initialize = parse_initialize_wvfnc(lammps_rtdbstring,false);
      if (wvfnc_initialize) lammps_rtdbstring = parse_initialize_wvfnc_set(lammps_rtdbstring,false);

      int taskid;
      int ierr = MPI_Comm_rank(comm_world,&taskid);
      bool oprint = (taskid==MASTER) && false;
      std::string input_wavefunction_filename = pwdft::parse_input_wavefunction_filename(lammps_rtdbstring);
      int wfound=0; if (taskid==MASTER) { std::ifstream wfile(input_wavefunction_filename); if (wfile.good()) wfound=1; wfile.close(); }
      MPI_Bcast(&wfound,1,MPI_INT,MASTER,comm_world);
     
      if ((!wfound) || (wvfnc_initialize))
      {
        auto lowlevel_rtdbstrs = pwdft::parse_gen_lowlevel_rtdbstrs(lammps_rtdbstring);
        for (const auto & elem: lowlevel_rtdbstrs) {
           std::string dum_rtdbstr  = elem;
           if (wvfnc_initialize) {
              dum_rtdbstr = parse_initialize_wvfnc_set(dum_rtdbstr,true);
              wvfnc_initialize = false;
           }
           if (oprint) coutput << std::endl << "Running staged energy optimization - lowlevel_rtdbstr = " << dum_rtdbstr << std::endl << std::endl;
           ierr += pwdft::pspw_minimizer(comm_world,dum_rtdbstr,coutput);
        }
      }

   }
}

extern int lammps_pspw_cpmd_start(MPI_Comm comm_world,std::ostream& coutput) {
   auto ierr = pwdft::ctask_cpmd_start(comm_world,lammps_rtdbstring,coutput);
   return ierr;
}

extern int lammps_pspw_cpmd_run(MPI_Comm comm_world, double *rion, double *uion, 
                                double *fion, double *qion, double *E,
                                bool removeqmmmcoulomb, bool removeqmqmcoulomb, std::ostream& coutput) {
   double Etot,Eapc;
   auto ierr = ctask_cpmd_run(comm_world,rion,uion,qion,fion,&Etot,&Eapc,coutput);
   if (removeqmmmcoulomb)
      *E = Etot - Eapc;
   else
      *E = Etot;
   
   if (removeqmqmcoulomb)
   {
   }


   return ierr;
}

extern int lammps_pspw_cpmd_stop(MPI_Comm comm_world,std::ostream& coutput) {
   auto ierr = pwdft::ctask_cpmd_stop(comm_world,coutput);
   return ierr;
}



// *************************************************************************
//
// LAMMPS Interface filename Routines - output is appended to filename
//
// *************************************************************************
class NullBuffer : public std::streambuf
{
public:
  int overflow(int c) { return c; }
};

extern int lammps_pspw_aimd_minimizer_filename(MPI_Comm comm_world, double *rion, double *fion, double *E, std::string& filename)
{
   int ierr;
   if (filename.empty())
   {
      NullBuffer null_buffer;
      std::ostream null_stream(&null_buffer);
      ierr = lammps_pspw_aimd_minimizer(comm_world,rion,fion,E,null_stream);
   }
   else
   {
      std::ofstream nwout(filename,std::ios_base::app);
      ierr = lammps_pspw_aimd_minimizer(comm_world,rion,fion,E,nwout);
   }
   return ierr;
}
extern int lammps_pspw_qmmm_minimizer_filename(MPI_Comm comm_world, double *rion, double *uion, double *fion, double *qion, double *E, 
                                               bool removeqmmmcoulomb, bool removeqmqmcoulomb, std::string& filename)
{
   int ierr;
   if (filename.empty())
   {
      NullBuffer null_buffer;
      std::ostream null_stream(&null_buffer);
      ierr = lammps_pspw_qmmm_minimizer(comm_world,rion,uion,fion,qion,E,removeqmmmcoulomb,removeqmqmcoulomb,null_stream);
   }
   else
   {
      std::ofstream nwout(filename,std::ios_base::app);
      ierr = lammps_pspw_qmmm_minimizer(comm_world,rion,uion,fion,qion,E,removeqmmmcoulomb,removeqmqmcoulomb,nwout);
   }
   return ierr;
}

extern int lammps_pspw_qmmm_nominimizer_filename(MPI_Comm comm_world, double *rion, double *uion, double *fion, double *qion, double *E,
                                                 bool removeqmmmcoulomb, bool removeqmqmcoulomb, std::string& filename)
{
   int ierr;
   if (filename.empty())
   {
      NullBuffer null_buffer;
      std::ostream null_stream(&null_buffer);
      ierr = lammps_pspw_qmmm_nominimizer(comm_world,rion,uion,fion,qion,E,removeqmmmcoulomb,removeqmqmcoulomb,null_stream);
   }
   else
   {
      std::ofstream nwout(filename,std::ios_base::app);
      ierr = lammps_pspw_qmmm_nominimizer(comm_world,rion,uion,fion,qion,E,removeqmmmcoulomb,removeqmqmcoulomb,nwout);
   }
   return ierr;
}




extern void lammps_pspw_input_filename(MPI_Comm comm_world, std::string& nwfilename, std::string& filename)
{
   if (filename.empty())
   {
      NullBuffer null_buffer;
      std::ostream null_stream(&null_buffer);
      lammps_pspw_input(comm_world,nwfilename,null_stream);
   }
   else
   {
      std::ofstream nwout(filename,std::ios_base::app);
      lammps_pspw_input(comm_world,nwfilename,nwout);
   }
}



// *************************************************************************
//
// LAMMPS c - Interface filename Routines - output is appended to filename
//
// *************************************************************************
static std::string convertcstring(const char *cstring)
{
   if (cstring)
   {
      std::string rstring(cstring);
      return rstring;
   }
   else
   {
      std::string rstring("");
      return rstring;
   }
}
extern int c_lammps_pspw_aimd_minimizer_filename(MPI_Comm comm_world, double *rion, double *fion, double *E, const char *cfilename)
{
   std::string filename = convertcstring(cfilename);
   return lammps_pspw_aimd_minimizer_filename(comm_world,rion,fion,E,filename);
}
extern int c_lammps_pspw_qmmm_minimizer_filename(MPI_Comm comm_world, double *rion, double *uion, double *fion, double *qion, double *E, 
                                                 bool removeqmmmcoulomb, bool removeqmqmcoulomb, const char *cfilename)
{
   std::string filename = convertcstring(cfilename);
   return lammps_pspw_qmmm_minimizer_filename(comm_world,rion,uion,fion,qion,E,removeqmmmcoulomb,removeqmqmcoulomb,filename);
} 
extern int c_lammps_pspw_qmmm_nominimizer_filename(MPI_Comm comm_world, double *rion, double *uion, double *fion, double *qion, double *E, 
                                                   bool removeqmmmcoulomb, bool removeqmqmcoulomb, const char *cfilename)
{
   std::string filename = convertcstring(cfilename);
   return lammps_pspw_qmmm_nominimizer_filename(comm_world,rion,uion,fion,qion,E,removeqmmmcoulomb,removeqmqmcoulomb,filename);
} 

extern void c_lammps_pspw_input_filename(MPI_Comm comm_world, const char *cnwfilename, const char *cfilename)
{
   std::string nwfilename = convertcstring(cnwfilename);
   std::string filename   = convertcstring(cfilename);
   lammps_pspw_input_filename(comm_world,nwfilename,filename);
}




int main(int argc, char* argv[])
{
  int ierr=0;
  int taskid,np;
  std::string line,nwinput,nwfilename;

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
        std::ifstream nwfile(argv[1]);
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
  std::string rtdbstr  = parse_nwinput(nwinput);
  int task = parse_task(rtdbstr);
  MPI_Barrier(MPI_COMM_WORLD);

  if (oprint) std::cout << "First rtdbstr=" << rtdbstr << std::endl;
  if (oprint) std::cout << "First task=" << task << std::endl << std::endl;


  // Initialize wavefunction
  if (parse_initialize_wvfnc(rtdbstr,true))
  {
     bool wvfnc_initialize = parse_initialize_wvfnc(rtdbstr,false);
     if (wvfnc_initialize) rtdbstr = parse_initialize_wvfnc_set(rtdbstr,false);

     std::string input_wavefunction_filename = parse_input_wavefunction_filename(rtdbstr);
     int wfound=0; if (taskid==MASTER) { std::ifstream wfile(input_wavefunction_filename); if (wfile.good()) wfound=1; wfile.close(); }
     MPI_Bcast(&wfound,1,MPI_INT,MASTER,MPI_COMM_WORLD);

     if ((!wfound) || (wvfnc_initialize))
     {
        auto lowlevel_rtdbstrs =  parse_gen_lowlevel_rtdbstrs(rtdbstr);
        for (const auto & elem: lowlevel_rtdbstrs) {
           // add initiale to el
           std::string dum_rtdbstr  = elem;
           if (wvfnc_initialize) {
              dum_rtdbstr = parse_initialize_wvfnc_set(dum_rtdbstr,true);
              wvfnc_initialize = false;
           }
           if (oprint) std::cout << std::endl << "Running staged energy optimization - lowlevel_rtdbstr = " << dum_rtdbstr << std::endl << std::endl;
           ierr += pspw_minimizer(MPI_COMM_WORLD,dum_rtdbstr,std::cout);
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
        ierr += pwdft::pspw_minimizer(MPI_COMM_WORLD,rtdbstr,std::cout);
     }

     /* Optimize task */
     if (task==3) 
     {
        if (oprint) std::cout << std::endl << "Running geometry optimization calculation - rtdbstr = " << rtdbstr << std::endl << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        ierr += pwdft::pspw_geovib(MPI_COMM_WORLD,rtdbstr,std::cout);
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

     /* Born-Oppenheimer task */
     if (task==7) 
     {
        if (oprint) std::cout << std::endl << "Running Born-Oppenheimer molecular dynamics - rtdbstr = " << rtdbstr << std::endl << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        ierr += pwdft::pspw_bomd(MPI_COMM_WORLD,rtdbstr,std::cout);
     }

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


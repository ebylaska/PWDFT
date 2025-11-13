


#include <cmath>

#include "vdw.hpp"
#include "Control2.hpp"
#include "compressed_io.hpp"

#include "NwpwLibraryVdwConfig.hpp"

#include "vdw_DF.hpp"

#include <iostream>
#include "iofmt.hpp"
#include <cstring>



namespace pwdft {


/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
/**
 * @brief Constructor for the `d3db` class.
 *
 * This constructor initializes an instance of the `d3db` class with the given parameters.
 *
 * @param inparall A pointer to a `Parallel` object.
 * @param inmaptype An integer specifying the mapping type.
 * @param nx The number of grid points in the x-direction.
 * @param ny The number of grid points in the y-direction.
 * @param nz The number of grid points in the z-direction.
 */
vdw_DF::vdw_DF(PGrid *inmygrid, Control2 &control)
{
   mygrid   = inmygrid;
   myparall = mygrid->d3db::parall;
   has_vdw = true;

   nfft3d = mygrid->nfft3d;
   npack0 = mygrid->npack(0);
   n2ft3d = mygrid->n2ft3d;

   bool oprint = (myparall->is_master() && control.print_level("medium"));


   const std::string nwpw_vdw_qmesh = std::string(Nwpw_LIBRARYVDW_Default) + "/vdw_qmesh.dat";
   //const char *nwpw_libraryps = Nwpw_LIBRARYPS_Default + "/VDW/vdw_qmesh.dat";

   char datafile[256];
   strcpy(datafile, "vdw_kernels.dat");
   control.add_permanent_dir(datafile);
   //std::cout << "nwpw_vdw = " << nwpw_vdw_qmesh << " datafile=" << datafile << std::endl;

   int ifound = cfileexists(datafile);
   if (ifound == 0)
   {
      if (oprint) std::cout << "Generating VDW kernel filename:" << datafile << std::endl;
      //vdw_DF_kernel_gen_data(datafile)
      vdw_DF_kernel_gen_data(myparall,datafile,nwpw_vdw_qmesh.c_str());
   }

   if (myparall->is_master())
   {
      openfile(5,datafile,"r");
      iread(5,&Nqs,1);
      iread(5,&nk,1);
      dread(5,&kmax,1);
   }
   myparall->Brdcst_iValue(0,MASTER,&Nqs);
   myparall->Brdcst_iValue(0,MASTER,&nk);
   myparall->Brdcst_Values(0,MASTER,1,&kmax);
   nk1 = nk + 1;


   double *qmesh  = new (std::nothrow) double[Nqs]();
   double *ya     = new (std::nothrow) double[Nqs*Nqs]();
   double *ya2    = new (std::nothrow) double[Nqs*Nqs]();
   double *gphi   = new (std::nothrow) double[nk1]();
   double *phi    = new (std::nothrow) double[nk1*Nqs*(Nqs+1)]();
   double *theta  = new (std::nothrow) double[Nqs*n2ft3d]();
   double *ufunc  = new (std::nothrow) double[Nqs*n2ft3d]();
   double *xcp    = new (std::nothrow) double[2*n2ft3d]();
   double *xce    = new (std::nothrow) double[2*n2ft3d]();
   double *xxp    = new (std::nothrow) double[2*n2ft3d]();
   double *xxe    = new (std::nothrow) double[2*n2ft3d]();
   double *rho    = new (std::nothrow) double[2*n2ft3d]();
   double *Gpack  = new (std::nothrow) double[npack0]();
   int    *nxpack = new (std::nothrow) int[npack0]();


   if (myparall->is_master())
   {
      dread(5,qmesh,Nqs);
      dread(5,phi,nk1*Nqs*(Nqs+1));
   }

      


};


}





#include	<iostream>
#include	<cstdio>
#include	<cmath>
#include	<cstdlib>
#include	<string>
#include	<vector>
//
#include "iofmt.hpp"
#include "util_linesearch.hpp"
#include	"Parallel.hpp"
//#include	"control.hpp"
#include	"Control2.hpp"
#include	"Lattice.hpp"
#include	"util_date.hpp"
#include	"PGrid.hpp"
#include	"Pneb.hpp"
#include	"Ion.hpp"
#include	"Ewald.hpp"
#include	"Kinetic.hpp"
#include	"exchange_correlation.hpp"
#include	"Pseudopotential.hpp"
#include	"Electron.hpp"
#include	"psi.hpp"
//#include	"rtdb.hpp"
#include	"mpi.h"


#include	"psp_library.hpp"
#include	"psp_file_check.hpp"
#include	"nwpw_timing.hpp"
#include "gdevice.hpp"

#include "nwpw_dplot.hpp"


#include "json.hpp"
using json = nlohmann::json;




namespace pwdft {

/******************************************
 *                                        *
 *              pspw_dplot                *
 *                                        *
 ******************************************/
int pspw_dplot(MPI_Comm comm_world0,std::string& rtdbstring,std::ostream& coutput)
{
   Parallel myparallel(comm_world0);

   int ne[2];
   char date[26];
   double *psi1,*psi_r,*dn,*rho;


   Control2 control(myparallel.np(),rtdbstring);
   int flag =  control.task();

   bool hprint = (myparallel.is_master() && control.print_level("high"));
   bool oprint = (myparallel.is_master() && control.print_level("medium"));
   bool lprint = (myparallel.is_master() && control.print_level("low"));

   /* reset Parallel base_stdio_print = lprint */
   myparallel.base_stdio_print = lprint;


   if (oprint)
   {
      std::ios_base::sync_with_stdio();
      coutput << "          *****************************************************\n";
      coutput << "          *                                                   *\n";
      coutput << "          *                 PWDFT PSPW dplot                  *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *  [              C++ implementation             ]  *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *              version #1.00   03/15/23             *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *    This code was developed by Eric J. Bylaska,    *\n";
      coutput << "          *    Abhishek Bagusetty, David H. Bross, ...        *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *****************************************************\n";
      coutput << "          >>> job started at       " << util_date() << " <<<\n";
   }

   /* initialize processor grid structure */
   myparallel.init2d(control.np_orbital(),control.pfft3_qsize());

   /* initialize lattice */
   Lattice mylattice(control);

   /* read in ion structure */
   //Ion myion(myrtdb);
   Ion myion(rtdbstring,control);

   /* Check for and generate psp files                       */
   /* - this routine also sets the valence charges in myion, */
   /*   and total_ion_charge and ne in control               */
   psp_file_check(&myparallel,&myion,control,coutput);
   MPI_Barrier(comm_world0);


   /* initialize parallel grid structure */
   Pneb mygrid(&myparallel,&mylattice,control,control.ispin(),control.ne_ptr());


   /* initialize psi1 */
   int ispin = control.ispin();
   int neall = mygrid.neq[0] + mygrid.neq[1];
   int shift1 = 2*(mygrid.npack(1));
   int shift2 = (mygrid.n2ft3d);
   int n2ft3d = (mygrid.n2ft3d);

   double omega = mylattice.omega();
   double scal1 = 1.0/((double) ((mygrid.nx)*(mygrid.ny)*(mygrid.nz)));
   double scal2 = 1.0/omega;
   double dv    = omega*scal1;
   ne[0] = control.ne(0);
   ne[1] = control.ne(1);
   psi1  = mygrid.g_allocate(1);
   psi_r = mygrid.h_allocate();
   dn    = mygrid.r_nalloc(ispin);
   rho   = mygrid.r_alloc();
   gdevice_psi_alloc(mygrid.npack(1),mygrid.neq[0]+mygrid.neq[1],control.tile_factor());


   bool newpsi = psi_read(&mygrid,control.input_movecs_filename(),
                                  control.input_movecs_initialize(),psi1,coutput);
   MPI_Barrier(comm_world0);

   /* initialize nwpw_dplot */
   nwpw_dplot mydplot(&myion,&mygrid,control);

   /* convert psi(G) to psi(r) - Expensive */
   int indx1 = 0;
   int indx2 = 0;
   for (auto i=0; i<neall; ++i)
   {
      mygrid.cc_pack_copy(1,psi1+indx1,psi_r+indx2);
      mygrid.c_unpack(1,psi_r+indx2);
      mygrid.cr_fft3d(psi_r+indx2);

      indx1 += shift1;
      indx2 += shift2;
   }

   /* generate dn */
   mygrid.hr_aSumSqr(scal2,psi_r,dn);

   //std::cout << "number of cubefiles = " << control.number_cubefiles() << std::endl;
   if (oprint) coutput << std::endl;
   for (auto i=0; i<control.number_cubefiles(); ++i)
   {
      int cubetype = control.cubetype_cubefiles(i);
      if (cubetype!=0)
      {
         std::string cubename = control.cubename_cubefiles(i);
         std::string cubekey  = control.cubekey_cubefiles(i);
         std::string cube_comment;
         // orbital
         if (cubetype>0)
         {
            if (oprint) coutput << " generating cubefile - orbital     " << Ifmt(5) << cubetype 
                                << " - cubefilename = " << cubename 
                                << " (json key=" << cubekey << ")" << std::endl;
             int ishift = (cubetype-1)*n2ft3d;;
             std::memcpy(rho,psi_r+ishift,n2ft3d*sizeof(double));
             cube_comment = "SCF Molecular Orbitals";
         }
         // total density
         else if (cubetype==-1)
         {
            if (oprint) coutput << " generating cubefile - total density    "
                                << " - cubefilename = " << cubename 
                                << " (json key=" << cubekey << ")" << std::endl;
            mygrid.rrr_Sum(dn,dn+(ispin-1)*n2ft3d,rho);
            cube_comment = "SCF Total Density";
         }
         // difference density
         else if (cubetype==-2)
         {
            if (oprint) coutput << " generating cubefile - diff density     "
                                << " - cubefilename = " << cubename 
                                << " (json key=" << cubekey << ")" << std::endl;
            std::memcpy(rho,dn,n2ft3d*sizeof(double));
            mygrid.rrr_Minus(dn,dn+(ispin-1)*n2ft3d,rho);
            cube_comment = "SCF Spin Density";

         }
         // alpha density
         else if (cubetype==-3)
         {
            if (oprint) coutput << " generating cubefile - alpha density    "
                                << " - cubefilename = " << cubename 
                                << " (json key=" << cubekey << ")" << std::endl;
            std::memcpy(rho,dn,n2ft3d*sizeof(double));
            cube_comment = "SCF Alpha Density";
         }
         // beta density
         else if (cubetype==-4)
         {
            if (oprint) coutput << " generating cubefile - beta density     "
                                << " - cubefilename = " << cubename 
                                << " (json key=" << cubekey << ")" << std::endl;
            std::memcpy(rho,dn+(ispin-1)*n2ft3d,n2ft3d*sizeof(double));
            cube_comment = "SCF Beta Density";
         }
         // laplacian density
         else if (cubetype==-5)
         {
            if (oprint) coutput << " generating cubefile - laplacian density"
                                << " - cubefilename = " << cubename 
                                << " (json key=" << cubekey << ")" << std::endl;
            mygrid.rrr_Sum(dn,dn+(ispin-1)*n2ft3d,rho);
            mygrid.r_SMul(scal1,rho);
            mygrid.rc_fft3d(rho);
            mygrid.r_zero_ends(rho);
            double *Gx  = mygrid.Gxyz(0);
            double *Gy  = mygrid.Gxyz(1);
            double *Gz  = mygrid.Gxyz(2);
            for (auto k=0; k<(mygrid.nfft3d); ++k)
               rho[k] *= -0.5*(Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k]);
            mygrid.cr_fft3d(rho);
            cube_comment = "SCF Laplacian Density";
         }
         // potential density
         else if (cubetype==-6)
         {
            if (oprint) coutput << " generating cubefile - potential density"
                                << " - cubefilename = " << cubename 
                                << " (json key=" << cubekey << ")" << std::endl;
            cube_comment = "SCF Potential Density";
         }

         mydplot.gcube_write(cubename,cubetype,cube_comment,rho);
      }


   }


//                 |**************************|
// *****************   report consumed time   **********************
//                 |**************************|

   /* deallocate memory */
   mygrid.g_deallocate(psi1);
   mygrid.h_deallocate(psi_r);
   mygrid.r_dealloc(dn);
   mygrid.r_dealloc(rho);
   gdevice_psi_dealloc();

   MPI_Barrier(comm_world0);

   return 0;
}

}



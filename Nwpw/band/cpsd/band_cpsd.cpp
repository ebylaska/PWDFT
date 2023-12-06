#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "Parallel.hpp"
#include "iofmt.hpp"
//#include	"control.hpp"
#include "Control2.hpp"

#include "Ewald.hpp"
#include "Ion.hpp"
#include "Lattice.hpp"
#include "Brillouin.hpp"
#include "cKinetic.hpp"
#include "cCoulomb.hpp"
#include "cExchange_Correlation.hpp"
#include "Cneb.hpp"
#include "CPseudopotential.hpp"
#include "cpsi.hpp"

#include "CStrfac.hpp"

#include "util_date.hpp"
#include "mpi.h"

//#include "gdevice.hpp"
#include "nwpw_timing.hpp"
#include "psp_file_check.hpp"
#include "psp_library.hpp"

#include "json.hpp"
using json = nlohmann::json;

namespace pwdft {

/******************************************
 *                                        *
 *                band_cpsd               *
 *                                        *
 ******************************************/
int band_cpsd(MPI_Comm comm_world0, std::string &rtdbstring)
{
   Parallel myparallel(comm_world0);

   int version, nfft[3], ne[2], ispin, nbrillq, nbrillouin;
   int i, ii, ia, nn, ngrid[3], matype, nelem, icount, done;
   char date[26];
   double sum1, sum2, ev, zv;
   double cpu1, cpu2, cpu3, cpu4;
   double E[80], deltae, deltac, deltar, viral, unita[9], en[2];
   double *psi1, *psi2, *Hpsi, *psi_r;
   double *dn;
   double *hml, *lmbda, *eig;
 
   Control2 control(myparallel.np(), rtdbstring);

   bool hprint = (myparallel.is_master() && control.print_level("high"));
   bool oprint = (myparallel.is_master() && control.print_level("medium"));
   bool lprint = (myparallel.is_master() && control.print_level("low"));
 
   /* reset Parallel base_stdio_print = lprint */
   myparallel.base_stdio_print = lprint;
 
   for (ii = 0; ii < 70; ++ii)
      E[ii] = 0.0;
 
   if (myparallel.is_master())
      seconds(&cpu1);
   if (oprint) {
      std::ios_base::sync_with_stdio();
      std::cout << "          *****************************************************\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *     Car-Parrinello solid-state calculation        *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *     [     steepest descent minimization   ]       *\n";
      std::cout << "          *     [          C++ implementation         ]       *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *            version #7.00   11/02/23               *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *    This code was developed by Eric J. Bylaska     *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *****************************************************\n";
      std::cout << "          >>> job started at       " << util_date() << " <<<\n";
   }
 
    /* initialize processor grid structure */
   myparallel.init3d(control.np_dimensions(1),control.np_dimensions(2),control.pfft3_qsize());
   MPI_Barrier(comm_world0);

   /* initialize lattice */
   Lattice mylattice(control);

   /* read in ion structure */
   // Ion myion(myrtdb);
   Ion myion(rtdbstring, control);
   MPI_Barrier(comm_world0);


   /* Check for and generate psp files                       */
   /* - this routine also sets the valence charges in myion, */
   /*   and total_ion_charge and ne in control               */
   psp_file_check(&myparallel, &myion, control, std::cout);
   MPI_Barrier(comm_world0);


   /* read in Brillouin zone */
   Brillouin mybrillouin(rtdbstring,&mylattice,control);
   //control.set_total_ion_charge(8);

   //std::cout << "ispin=" << control.ispin() << " ne=" << control.ne_ptr()[0] << " " << control.ne_ptr()[1] 
   //          << " nbrillioun=" << mybrillouin.nbrillouin << std::endl;
   /* initialize parallel grid structure */
   Cneb mygrid(&myparallel, &mylattice, control, control.ispin(),control.ne_ptr(),&mybrillouin);

   /* initialize psi1 and psi2 */
   //ispin = control.ispin(); 
   nbrillouin = mygrid.nbrillouin; ispin = mygrid.ispin; ne[0] = mygrid.ne[0]; ne[1] = mygrid.ne[1]; nbrillq = mygrid.nbrillq;
   psi1 = mygrid.g_allocate_nbrillq_all();
   psi2 = mygrid.g_allocate_nbrillq_all();
   Hpsi = mygrid.g_allocate_nbrillq_all();
   psi_r = mygrid.h_allocate_nbrillq_all();
   dn    = mygrid.r_nalloc(ispin);
   hml   = mygrid.w_allocate_nbrillq_all();
   lmbda = mygrid.w_allocate_nbrillq_all();
   eig   = new double[nbrillq*(ne[0]+ne[1])];


   // psi_read(&mygrid,&version,nfft,unita,&ispin,ne,nbrill,psi2,control.input_movecs_filename());
   bool newpsi = cpsi_read(&mygrid,control.input_movecs_filename(),control.input_movecs_initialize(),psi2,std::cout);
   mygrid.gg_copy(psi2,psi1);
   MPI_Barrier(comm_world0);

   /* setup structure factor */
   CStrfac mystrfac(&myion,&mygrid);
   mystrfac.phafac();

   /* initialize operators */
   cKinetic_Operator mykin(&mygrid);
   cCoulomb_Operator mycoulomb(&mygrid, control);
   cXC_Operator myxc(&mygrid, control);

   /* initialize psps */
   CPseudopotential mypsp(&myion,&mygrid,&mystrfac,control,std::cout);

   /* setup ewald */
   Ewald myewald(&myparallel,&myion,&mylattice,control,mypsp.zv);
   myewald.phafac();



   if (oprint)
   {
      std::cout << std::endl;
      std::cout << "     ===================  summary of input  =======================" << std::endl;
      std::cout << "\n input psi filename: " << control.input_movecs_filename() << std::endl;
      std::cout << std::endl;
      std::cout << " number of processors used: " << myparallel.np() << std::endl;
      std::cout << " processor grid           : " << myparallel.np_i() << " x " << myparallel.np_j() <<  " x " << myparallel.np_k() << std::endl;
      if (mygrid.maptype == 1) std::cout << " parallel mapping         : 1d-slab" << std::endl;
      if (mygrid.maptype == 2) std::cout << " parallel mapping         : 2d-hilbert" << std::endl;
      if (mygrid.maptype == 3) std::cout << " parallel mapping         : 2d-hcurve" << std::endl;
      if (mygrid.isbalanced())
         std::cout << " parallel mapping         : balanced" << std::endl;
      else
         std::cout << " parallel mapping         : not balanced" << std::endl;
      if (mygrid.staged_gpu_fft_pipeline) std::cout << " parallel mapping         : staged gpu fft" << std::endl;
      if (control.tile_factor() > 1)
         std::cout << " GPU tile factor          : " << control.tile_factor() << std::endl;

      std::cout << std::endl;
      std::cout << " options:" << std::endl;
      std::cout << "   ion motion           = ";
      if (control.geometry_optimize())
         std::cout << "yes" << std::endl;
      else
         std::cout << "no" << std::endl;
      std::cout << "   boundary conditions  = periodic" << std::endl;
      std::cout << "   electron spin        = ";
      if (ispin == 1)
         std::cout << "restricted" << std::endl;
      else
         std::cout << "unrestricted" << std::endl;
      std::cout << myxc;


      std::cout << mypsp.print_pspall();

      std::cout << std::endl;
      std::cout << " atom composition:" << std::endl;
      for (ia = 0; ia < myion.nkatm; ++ia)
         std::cout << "   " << myion.atom(ia) << " : " << myion.natm[ia];
      std::cout << std::endl << std::endl;
      std::cout << " initial ion positions (au):" << std::endl;
      for (ii = 0; ii < myion.nion; ++ii)
         std::cout << Ifmt(4) << ii + 1 << " " << myion.symbol(ii) << "\t( "
                   << Ffmt(10,5) << myion.rion1[3*ii] << " "
                   << Ffmt(10,5) << myion.rion1[3*ii+1] << " "
                   << Ffmt(10,5) << myion.rion1[3*ii+2] << " ) - atomic mass = "
                   << Ffmt(6,3) << myion.amu(ii) << std::endl;
      std::cout << "   G.C.\t( "
                << Ffmt(10,5) << myion.gc(0) << " "
                << Ffmt(10,5) << myion.gc(1) << " "
                << Ffmt(10,5) << myion.gc(2) << " )" << std::endl;
      std::cout << " C.O.M.\t( "
                << Ffmt(10,5) << myion.com(0) << " "
                << Ffmt(10,5) << myion.com(1) << " "
                << Ffmt(10,5) << myion.com(2) << " )" << std::endl;

      std::cout << std::endl;
      std::cout << myion.print_symmetry_group();

      if (control.geometry_optimize())
         std::cout << std::endl << myion.print_constraints(0);

      std::cout << std::endl;
      std::cout << " number of electrons: spin up ="
                << Ifmt(6) << mygrid.ne[0] << " ("
                << Ifmt(4) << mygrid.neq[0] << " per task) down ="
                << Ifmt(6) << mygrid.ne[ispin-1] << " ("
                << Ifmt(4) << mygrid.neq[ispin-1] << " per task)" << std::endl;

      std::cout << std::endl;
      std::cout << " supercell:" << std::endl;
      std::cout << "      volume = " << Ffmt(10,2) << mylattice.omega()
                << std::endl;
      std::cout << "      lattice:    a1 = < "
                << Ffmt(8,3) << mylattice.unita(0,0) << " "
                << Ffmt(8,3) << mylattice.unita(1,0) << " "
                << Ffmt(8,3) << mylattice.unita(2,0) << " >\n";
      std::cout << "                  a2 = < "
                << Ffmt(8,3) << mylattice.unita(0,1) << " "
                << Ffmt(8,3) << mylattice.unita(1,1) << " "
                << Ffmt(8,3) << mylattice.unita(2,1) << " >\n";
      std::cout << "                  a3 = < "
                << Ffmt(8,3) << mylattice.unita(0,2) << " "
                << Ffmt(8,3) << mylattice.unita(1, 2) << " "
                << Ffmt(8,3) << mylattice.unita(2, 2) << " >\n";
      std::cout << "      reciprocal: b1 = < "
                << Ffmt(8,3) << mylattice.unitg(0, 0) << " "
                << Ffmt(8,3) << mylattice.unitg(1, 0) << " "
                << Ffmt(8,3) << mylattice.unitg(2,0) << " >\n";
      std::cout << "                  b2 = < "
                << Ffmt(8,3) << mylattice.unitg(0,1) << " "
                << Ffmt(8,3) << mylattice.unitg(1,1) << " "
                << Ffmt(8,3) << mylattice.unitg(2,1) << " >\n";
      std::cout << "                  b3 = < "
                << Ffmt(8,3) << mylattice.unitg(0,2) << " "
                << Ffmt(8,3) << mylattice.unitg(1,2) << " "
                << Ffmt(8,3) << mylattice.unitg(2,2) << " >\n";

      {
         double aa1, bb1, cc1, alpha1, beta1, gamma1;
         mylattice.abc_abg(&aa1, &bb1, &cc1, &alpha1, &beta1, &gamma1);
         std::cout << "      lattice:    a =    "
                   << Ffmt(8,3) << aa1 << " b =   "
                   << Ffmt(8,3) << bb1 << " c =    "
                   << Ffmt(8,3) << cc1 << std::endl;
         std::cout << "                  alpha ="
                   << Ffmt(8,3) << alpha1 << " beta ="
                   << Ffmt(8,3) << beta1 << " gamma ="
                   << Ffmt(8,3) << gamma1 << std::endl;
      }

      
      std::cout << "      density cutoff ="
                << Ffmt(7,3) << mylattice.ecut()
                << " fft =" << Ifmt(4) << mygrid.nx << " x "
                            << Ifmt(4) << mygrid.ny << " x "
                            << Ifmt(4) << mygrid.nz
                << "  (" << Ifmt(8) << mygrid.npack_all(0) << " waves "
                         << Ifmt(8) << mygrid.npack(0) << " per task)" << std::endl;
      std::cout << "      wavefnc cutoff ="
                << Ffmt(7,3) << mylattice.wcut()
                << " fft =" << Ifmt(4) << mygrid.nx << " x "
                            << Ifmt(4) << mygrid.ny << " x "
                            << Ifmt(4) << mygrid.nz
                << "  (" << Ifmt(8) << mygrid.npack_all(1) << " waves "
                         << Ifmt(8) << mygrid.npack(1) << " per task)" << std::endl;
    
      std::cout << "\n";
      std::cout << " Ewald parameters:\n";
      std::cout << "      energy cutoff = "
                << Ffmt(7,3) << myewald.ecut()
                << " fft= " << Ifmt(4) << myewald.nx() << " x "
                            << Ifmt(4) << myewald.ny() << " x "
                            << Ifmt(4) << myewald.nz()
                << "  (" << Ifmt(8) << myewald.npack_all() << " waves "
                         << Ifmt(8) << myewald.npack() << " per task)" << std::endl;
      std::cout << "      Ewald summation: cut radius = "
                << Ffmt(7,3) << myewald.rcut() << " and " << Ifmt(3) << myewald.ncut() << std::endl;
      std::cout << "                       Mandelung Wigner-Seitz ="
                << Ffmt(12,8) << myewald.mandelung()
                << " (alpha =" << Ffmt(12,8) << myewald.rsalpha()
                << " rs =" << Ffmt(12,8) << myewald.rs() << ")" << std::endl;

       /* print nbrillouin */
      std::cout << std::endl;
      std::cout << " brillouin zone:" << std::endl;
      std::cout << mybrillouin.print_zone();
      std::cout << std::endl;
   }

   MPI_Barrier(comm_world0);

   // write wavefunctions
   version = 5;
   nfft[0] = mygrid.nx;
   nfft[1] = mygrid.ny;
   nfft[2] = mygrid.nz;
   std::cout << "version=" << version<< std::endl;
   std::cout << "nfft=" << nfft[0] << " " << nfft[1] << " " << nfft[2] << std::endl;
   std::cout << "ispin=" << ispin << std::endl;
   std::cout << "ne=" << ne[0] << " " << ne[1] << std::endl;
   cpsi_write(&mygrid,&version,nfft,mylattice.unita_ptr(),&mygrid.ispin,mygrid.ne,&mygrid.nbrillouin,psi1,control.output_movecs_filename(),std::cout);


   /* deallocate memory */
   std::cout << "deallocate memory, taskid=" << myparallel.taskid() 
             << " taskid_i=" << myparallel.taskid_i() 
             << " taskid_j=" << myparallel.taskid_j() 
             << " taskid_k=" << myparallel.taskid_k() << std::endl;
   mygrid.g_deallocate(psi1);
   mygrid.g_deallocate(psi2);
   mygrid.g_deallocate(Hpsi);
   mygrid.h_deallocate(psi_r);
   mygrid.r_dealloc(dn);
   mygrid.w_deallocate(hml);
   mygrid.w_deallocate(lmbda);
   delete [] eig;
   //mygrid.d3db::mygdevice.psi_dealloc();


 
   
   MPI_Barrier(comm_world0);
   return 0;
}

} // namespace pwdft

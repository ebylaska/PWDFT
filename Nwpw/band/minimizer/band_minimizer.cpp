
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
//
#include "Parallel.hpp"
#include "iofmt.hpp"
#include "util_linesearch.hpp"
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
#include "Solid.hpp"
#include "cElectron.hpp"

#include "util_date.hpp"
//#include	"rtdb.hpp"
#include "mpi.h"

//#include "gdevice.hpp"
#include "nwpw_timing.hpp"
#include "psp_file_check.hpp"
#include "psp_library.hpp"

#include "band_cgsd_energy.hpp"
#include "band_cgsd_excited.hpp"

#include "json.hpp"
using json = nlohmann::json;

namespace pwdft {

/******************************************
 *                                        *
 *            band_minimizer              *
 *                                        *
 ******************************************/
int band_minimizer(MPI_Comm comm_world0, std::string &rtdbstring, std::ostream &coutput)
{
   // Parallel myparallel(argc,argv);
   Parallel myparallel(comm_world0);
   // RTDB myrtdb(&myparallel, "eric.db", "old");
 
   int version, nfft[3], ne[2], ispin;
   int i, ii, ia, nn, ngrid[3], matype, nelem, icount, done;
   char date[26];
   double sum1, sum2, ev, zv;
   double cpu1, cpu2, cpu3, cpu4;
   double E[70], deltae, deltac, deltar, viral, unita[9];
 
   // double *psi1,*psi2,*Hpsi,*psi_r;
   // double *dn;
   // double *hml,*lmbda,*eig;
 
   for (ii=0; ii<70; ++ii)
      E[ii] = 0.0;
 
   Control2 control(myparallel.np(), rtdbstring);
   int flag = control.task();

   bool hprint = (myparallel.is_master() && control.print_level("high"));
   bool oprint = (myparallel.is_master() && control.print_level("medium"));
   bool lprint = (myparallel.is_master() && control.print_level("low"));
  
   /* reset Parallel base_stdio_print = lprint */
   myparallel.base_stdio_print = lprint;
  
   if (myparallel.is_master()) seconds(&cpu1);
   if (oprint) 
   {
      std::ios_base::sync_with_stdio();
      coutput << "          *****************************************************\n";
      coutput << "          *                                                   *\n";
      coutput << "          *               PWDFT BAND Calculation              *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *  [      (bundled Grassmann/Stiefel manifold)   ]  *\n";
      coutput << "          *  [              C++ implementation             ]  *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *              version #7.00   07/19/24             *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *    This code was developed by Eric J. Bylaska,    *\n";
      coutput << "          *                                        ...        *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *****************************************************\n";
      coutput << "          >>> job started at       " << util_date() << " <<<\n";
   }
  
   // control_read(myrtdb);
   // control_read(myparallel.np(),rtdbstring);
  
   // initialize processor grid structure
   myparallel.init3d(control.np_dimensions(1),control.np_dimensions(2),control.pfft3_qsize());
   MPI_Barrier(comm_world0);
   
  
 
   // initialize lattice
   Lattice mylattice(control);
  
   // read in ion structure
   // Ion myion(myrtdb);
   Ion myion(rtdbstring, control);
  
   // Check for and generate psp files
   // - this routine also sets the valence charges in myion,
   //   and total_ion_charge and ne in control
   psp_file_check(&myparallel,&myion,control,coutput);
   MPI_Barrier(comm_world0);


   // fetch ispin and ne psi information from control
   ispin = control.ispin();
   ne[0] = control.ne(0);
   ne[1] = control.ne(1);
   nfft[0] = control.ngrid(0);
   nfft[1] = control.ngrid(1);
   nfft[2] = control.ngrid(2);
   unita[0] = mylattice.unita1d(0);
   unita[1] = mylattice.unita1d(1);
   unita[2] = mylattice.unita1d(2);
   unita[3] = mylattice.unita1d(3);
   unita[4] = mylattice.unita1d(4);
   unita[5] = mylattice.unita1d(5);
   unita[6] = mylattice.unita1d(6);
   unita[7] = mylattice.unita1d(7);
   unita[8] = mylattice.unita1d(8);
   version = control.version;




   /* read in Brillouin zone */
   Brillouin mybrillouin(rtdbstring,&mylattice,control);

   // initialize parallel grid structure
   Cneb mygrid(&myparallel, &mylattice, control, control.ispin(),control.ne_ptr(),&mybrillouin);
  
   // initialize gdevice memory
   mygrid.c3db::mygdevice.psi_alloc(mygrid.npack1_max(),mygrid.neq[0]+mygrid.neq[1],control.tile_factor());


  
   // setup structure factor
   CStrfac mystrfac(&myion,&mygrid);
   mystrfac.phafac();
   mystrfac.phafac_k();
  
   // initialize operators
   cKinetic_Operator mykin(&mygrid);
   cCoulomb_Operator mycoulomb(&mygrid,control);
   cXC_Operator myxc(&mygrid, control);

   /* initialize psps */
   CPseudopotential mypsp(&myion,&mygrid,&mystrfac,control,coutput);

   // initialize electron operators
   cElectron_Operators myelectron(&mygrid,&mykin,&mycoulomb,&myxc,&mypsp);
  

   /* setup ewald */
   Ewald myewald(&myparallel,&myion,&mylattice,control,mypsp.zv);
   myewald.phafac();



  
   // initialize solid
   Solid mysolid(control.input_movecs_filename(),
                 control.input_movecs_initialize(),&mygrid,&myion,
                 &mystrfac,&myewald,&myelectron,&mypsp,control,coutput);
  
   /* intialize the linesearch */
   util_linesearch_init();
  
   //                 |**************************|
   // *****************   summary of input data  **********************
   //                 |**************************|
  
   MPI_Barrier(comm_world0);

   if (oprint)
   {
      coutput << std::endl;
      coutput << "     ===================  summary of input  =======================" << std::endl;
      coutput << "\n input psi filename: " << control.input_movecs_filename() << std::endl;
      coutput << std::endl;
      coutput << " number of processors used: " << myparallel.np() << std::endl;
      coutput << " processor grid           : " << myparallel.np_i() << " x " << myparallel.np_j() <<  " x " << myparallel.np_k() << std::endl;
      if (mygrid.maptype == 1) coutput << " parallel mapping         : 1d-slab" << std::endl;
      if (mygrid.maptype == 2) coutput << " parallel mapping         : 2d-hilbert" << std::endl;
      if (mygrid.maptype == 3) coutput << " parallel mapping         : 2d-hcurve" << std::endl;
      if (mygrid.isbalanced())
         coutput << " parallel mapping         : balanced" << std::endl;
      else
         coutput << " parallel mapping         : not balanced" << std::endl;

      if (mygrid.c3db::mygdevice.has_gpu())
      {
         //coutput << " parallel mapping         : has GPUs" << std::endl;
         if (mygrid.c3db::mygdevice.type_gpu()==1) coutput << " parallel mapping         : CUDA GPUs" << std::endl;
         if (mygrid.c3db::mygdevice.type_gpu()==2) coutput << " parallel mapping         : SYCL GPUs" << std::endl;
         if (mygrid.c3db::mygdevice.type_gpu()==3) coutput << " parallel mapping         : HIP SYCL GPUs" << std::endl;
         if (mygrid.c3db::mygdevice.type_gpu()==4) coutput << " parallel mapping         : OpenCL GPUs" << std::endl;
         if (mygrid.staged_gpu_fft_pipeline) coutput << " parallel mapping         : staged GPU FFT" << std::endl;
         if (control.tile_factor() > 0)      coutput << " GPU tile factor          : " << control.tile_factor() << std::endl;
      }

      coutput << std::endl;
      coutput << " options:" << std::endl;
      coutput << "   boundary conditions  = periodic" << std::endl;
      coutput << "   electron spin        = ";
      if (ispin == 1)
         coutput << "restricted" << std::endl;
      else
         coutput << "unrestricted" << std::endl;
      coutput << myxc;

      coutput << mypsp.print_pspall();

      coutput << std::endl;
      coutput << " atom composition:" << std::endl;
      for (ia = 0; ia < myion.nkatm; ++ia)
         coutput << "   " << myion.atom(ia) << " : " << myion.natm[ia];
      coutput << std::endl << std::endl;
      coutput << " initial ion positions (au):" << std::endl;
      for (ii = 0; ii < myion.nion; ++ii)
         coutput << Ifmt(4) << ii + 1 << " " << myion.symbol(ii) << "\t( "
                   << Ffmt(10,5) << myion.rion1[3*ii] << " "
                   << Ffmt(10,5) << myion.rion1[3*ii+1] << " "
                   << Ffmt(10,5) << myion.rion1[3*ii+2] << " ) - atomic mass = "
                   << Ffmt(6,3) << myion.amu(ii) << std::endl;
      coutput << "   G.C.\t( "
                << Ffmt(10,5) << myion.gc(0) << " "
                << Ffmt(10,5) << myion.gc(1) << " "
                << Ffmt(10,5) << myion.gc(2) << " )" << std::endl;
      coutput << " C.O.M.\t( "
                << Ffmt(10,5) << myion.com(0) << " "
                << Ffmt(10,5) << myion.com(1) << " "
                << Ffmt(10,5) << myion.com(2) << " )" << std::endl;


      coutput << std::endl;
      coutput << myion.print_symmetry_group();

      if (control.geometry_optimize())
         coutput << std::endl << myion.print_constraints(0);

      coutput << std::endl;
      coutput << " number of electrons: spin up ="
                << Ifmt(6) << mygrid.ne[0] << " ("
                << Ifmt(4) << mygrid.neq[0] << " per task) down ="
                << Ifmt(6) << mygrid.ne[ispin-1] << " ("
                << Ifmt(4) << mygrid.neq[ispin-1] << " per task)" << std::endl;

      coutput << std::endl;
      coutput << " supercell:" << std::endl;
      coutput << "      volume = " << Ffmt(10,2) << mylattice.omega()
                << std::endl;
      coutput << "      lattice:    a1 = < "
                << Ffmt(8,3) << mylattice.unita(0,0) << " "
                << Ffmt(8,3) << mylattice.unita(1,0) << " "
                << Ffmt(8,3) << mylattice.unita(2,0) << " >\n";
      coutput << "                  a2 = < "
                << Ffmt(8,3) << mylattice.unita(0,1) << " "
                << Ffmt(8,3) << mylattice.unita(1,1) << " "
                << Ffmt(8,3) << mylattice.unita(2,1) << " >\n";
      coutput << "                  a3 = < "
                << Ffmt(8,3) << mylattice.unita(0,2) << " "
                << Ffmt(8,3) << mylattice.unita(1, 2) << " "
                << Ffmt(8,3) << mylattice.unita(2, 2) << " >\n";
      coutput << "      reciprocal: b1 = < "
                << Ffmt(8,3) << mylattice.unitg(0, 0) << " "
                << Ffmt(8,3) << mylattice.unitg(1, 0) << " "
                << Ffmt(8,3) << mylattice.unitg(2,0) << " >\n";
      coutput << "                  b2 = < "
                << Ffmt(8,3) << mylattice.unitg(0,1) << " "
                << Ffmt(8,3) << mylattice.unitg(1,1) << " "
                << Ffmt(8,3) << mylattice.unitg(2,1) << " >\n";
      coutput << "                  b3 = < "
                << Ffmt(8,3) << mylattice.unitg(0,2) << " "
                << Ffmt(8,3) << mylattice.unitg(1,2) << " "
                << Ffmt(8,3) << mylattice.unitg(2,2) << " >\n";


      coutput << "\n";
      coutput << " Ewald parameters:\n";
      coutput << "      energy cutoff = "
                << Ffmt(7,3) << myewald.ecut()
                << " fft= " << Ifmt(4) << myewald.nx() << " x "
                            << Ifmt(4) << myewald.ny() << " x "
                            << Ifmt(4) << myewald.nz()
                << "  (" << Ifmt(8) << myewald.npack_all() << " waves "
                         << Ifmt(8) << myewald.npack() << " per task)" << std::endl;
      coutput << "      Ewald summation: cut radius = "
                << Ffmt(7,3) << myewald.rcut() << " and " << Ifmt(3) << myewald.ncut() << std::endl;
      coutput << "                       Mandelung Wigner-Seitz ="
                << Ffmt(12,8) << myewald.mandelung()
                << " (alpha =" << Ffmt(12,8) << myewald.rsalpha()
                << " rs =" << Ffmt(12,8) << myewald.rs() << ")" << std::endl;

       /* print nbrillouin */
      coutput << std::endl;
      coutput << " brillouin zone:" << std::endl;
      coutput << mybrillouin.print_zone();

       /* print nbrillouin */
      coutput << std::endl;
      coutput << " computational grids:" << std::endl;
      coutput << "      density     cutoff ="
                << Ffmt(7,3) << mylattice.ecut()
                << " fft =" << Ifmt(4) << mygrid.nx << " x "
                            << Ifmt(4) << mygrid.ny << " x "
                            << Ifmt(4) << mygrid.nz
                << "  (" << Ifmt(8) << mygrid.npack_all(0) << " waves "
                         << Ifmt(8) << mygrid.npack(0) << " per task)" << std::endl;
      for (auto nb=0; nb<mygrid.nbrillouin; ++nb)
         coutput << "      wavefnc"
                   << Ifmt(4) << nb+1
                   << " cutoff ="
                   << Ffmt(7,3) << mylattice.wcut()
                   << " fft =" << Ifmt(4) << mygrid.nx << " x "
                               << Ifmt(4) << mygrid.ny << " x "
                               << Ifmt(4) << mygrid.nz
                   << "  (" << Ifmt(8) << mygrid.npack_all_print(nb) << " waves "
                            << Ifmt(8) << mygrid.npack_print(nb) << " per task)" << std::endl;
      coutput << std::endl;



      if (flag > 0)
      {
         coutput << std::endl;
         coutput << " technical parameters:\n";
        if (control.io_buffer()) coutput << "      using io buffer " << std::endl;
         coutput << "      fixed step: time step =" << Ffmt(12, 2)
                 << control.time_step() << "  ficticious mass =" << Ffmt(12, 2)
                 << control.fake_mass() << std::endl;
         coutput << "      tolerance =" << Efmt(12, 3) << control.tolerances(0)
                 << " (energy) " << Efmt(12, 3) << control.tolerances(1)
                 << " (density) " << Efmt(12, 3) << control.tolerances(2)
                 << " (ion)\n";
         coutput << "      max iterations = " << Ifmt(10)
                 << control.loop(0) * control.loop(1) << " (" << Ifmt(5)
                 << control.loop(0) << " inner " << Ifmt(5) << control.loop(1)
                 << " outer)\n";
         if (control.minimizer()==1) coutput << "      minimizer = Grassmann conjugate gradient\n";
         if (control.minimizer()==2) coutput << "      minimizer = Grassmann lmbfgs\n";
         if (control.minimizer()==4) coutput << "      minimizer = Stiefel conjugate gradient\n";
         if (control.minimizer()==5) coutput << "      minimizer = scf (potential)\n";
         if (control.minimizer()==7) coutput << "      minimizer = Stiefel lmbfgs\n";
         if (control.minimizer()==8) coutput << "      minimizer = scf (density)\n";
         if ((control.minimizer()==5) || (control.minimizer()==8))
         {
            coutput << std::endl;
            coutput << " Kohn-Sham scf parameters:\n";
            coutput << "     Kohn-Sham algorithm  = conjugate gradient\n";
            coutput << "     SCF algorithm        = simple mixing\n";
            coutput << "     SCF mixing parameter =    x.xxxx\n";
            coutput << "     Kohn-Sham iterations = xxxx\n";
            if (control.minimizer()==5) coutput << "     SCF mixing type      = potential\n";
            if (control.minimizer()==8) coutput << "     SCF mixing type      = density\n";
            coutput << "     Kerker damping       =    x.xxxx\n";
         }
      }
      else
      {
         coutput << std::endl;
         coutput << " technical parameters:\n";
         coutput << "      optimization of psi and densities turned off" << std::endl;
      }
      coutput << std::endl << std::endl << std::endl;


   }

   MPI_Barrier(comm_world0);
   if (myparallel.is_master()) seconds(&cpu2);
  
   //*                |***************************|
   //******************     call CG minimizer     **********************
   //*                |***************************|
  
   // calculate energy
   double EV = 0.0;
  
   if (flag < 0) 
   {
      EV = band_cgsd_noit_energy(mysolid, true, coutput);
   }
   else 
   {
      EV = band_cgsd_energy(control, mysolid, true, coutput);
   }
   if (myparallel.is_master()) seconds(&cpu3);
  
   // write energy results to the json
   auto rtdbjson = json::parse(rtdbstring);
   rtdbjson["band"]["energy"] = EV;
   rtdbjson["band"]["energies"] = mysolid.E;
   rtdbjson["band"]["eigenvalues"] = mysolid.eig_vector();

   // calculate fion
   if (flag == 2) 
   {
      // double *fion = new double[3*myion.nion];
      double fion[3*myion.nion];
      band_cgsd_energy_gradient(mysolid, fion);
      if (lprint) 
      {
         coutput << std::endl << " Ion Forces (au):" << std::endl;
         for (ii = 0; ii < myion.nion; ++ii)
            coutput << Ifmt(4) << ii + 1 << " " << myion.symbol(ii) << "\t( "
                    << Ffmt(10,5) << fion[3 *ii] << " "
                    << Ffmt(10,5) << fion[3*ii+1] << " "
                    << Ffmt(10,5) << fion[3*ii+2] << " )\n";
         coutput << std::endl << std::endl;
      }
      rtdbjson["band"]["fion"] = std::vector<double>(fion, &fion[3 * myion.nion]);
      for (ii=0; ii<(3*myion.nion); ++ii) fion[ii] *= -1.0;
      rtdbjson["band"]["gradient"] = std::vector<double>(fion, &fion[3 * myion.nion]);
     
      // delete [] fion;
   }  


   // calculate excited state orbitals 
   band_cgsd_excited(control, mysolid, true, coutput);



   // write psi
   // write psi
   if (flag > 0)
      mysolid.writepsi(control.output_movecs_filename(), coutput);
   MPI_Barrier(comm_world0);
  
   rtdbstring = rtdbjson.dump();
   myion.writejsonstr(rtdbstring);
   //                 |**************************|
   // *****************   report consumed time   **********************
   //                 |**************************|
   if (myparallel.is_master()) seconds(&cpu4);
   if (oprint)
   {
      double t1 = cpu2 - cpu1;
      double t2 = cpu3 - cpu2;
      double t3 = cpu4 - cpu3;
      double t4 = cpu4 - cpu1;
      //double av = t2 / ((double)control.loop(0) * icount);
      double av = t2 / ((double)myelectron.counter);
      // coutput.setf(ios::scientific);
      coutput << std::scientific;
      coutput << std::endl;
      coutput << " -----------------" << std::endl;
      coutput << " cputime in seconds" << std::endl;
      coutput << " prologue    : " << t1 << std::endl;
      coutput << " main loop   : " << t2 << std::endl;
      coutput << " epilogue    : " << t3 << std::endl;
      coutput << " total       : " << t4 << std::endl;
      coutput << " cputime/step: " << av << std::endl;
      coutput << std::endl;

      nwpw_timing_print_final(myelectron.counter, coutput);

      coutput << std::endl;
      coutput << " >>> job completed at     " << util_date() << " <<<" << std::endl;
   }




   mygrid.c3db::mygdevice.psi_dealloc();
   MPI_Barrier(comm_world0);
   return 0;
}

} // namespace pwdft

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "Parallel.hpp"
#include "iofmt.hpp"
//#include	"control.hpp"
#include "Control2.hpp"
#include "Coulomb12.hpp"
#include "DFPT.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "Kinetic.hpp"
#include "Lattice.hpp"
#include "PGrid.hpp"
#include "Pneb.hpp"
#include "Pseudopotential.hpp"
#include "Strfac.hpp"
#include "exchange_correlation.hpp"
#include "HFX.hpp"
#include "inner_loop.hpp"
#include "psi.hpp"
#include "util_date.hpp"
#include "nwpw_aimd_running_data.hpp"
//#include	"rtdb.hpp"
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
 *                cpsd                    *
 *                                        *
 ******************************************/
int cpsd(MPI_Comm comm_world0, std::string &rtdbstring)
{
   // Parallel myparallel(argc,argv);
   Parallel myparallel(comm_world0);
   // RTDB myrtdb(&myparallel, "eric.db", "old");

   int version, nfft[3], ne[2], nextra[2],  ispin;
   int i, ii, ia, nn, ngrid[3], matype, nelem, icount, done;
   char date[26];
   double sum1, sum2, ev, zv;
   double cpu1, cpu2, cpu3, cpu4;
   double E[80], deltae, deltac, deltar, viral, unita[9], en[2];
   double *psi1, *psi2, *Hpsi, *psi_r;
   double *occ1, *occ2;
   double *dn;
   double *hml, *lmbda, *eig, *eig_prev;

   // psi smearing block
   bool fractional;
   int smearoccupation, smeartype;
   double smearfermi[2], smearcorrection, smearkT, smearcorrection_old;
 
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
      std::cout << "          *     Car-Parrinello calculation for molecules,     *\n";
      std::cout << "          *       microclusters, liquids, and materials       *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *     [     steepest descent minimization   ]       *\n";
      std::cout << "          *     [          C++ implementation         ]       *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *            version #7.00   09/20/18               *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *    This code was developed by Eric J. Bylaska     *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *****************************************************\n";
      std::cout << "          >>> job started at       " << util_date() << " <<<\n";
   }
 
   // control_read(myrtdb);
   // control_read(myparallel.np(),rtdbstring);
 
   /* initialize processor grid structure */
   myparallel.init2d(control.np_orbital(), control.pfft3_qsize());
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
 
   /* debug output - print charge, ispin, and ne */
   /*
      if (myparallel.is_master())
      {
          std::cout << endl;
          std::cout << "total_ion_charge = " << myion.total_zv() << endl;
          std::cout << "ispin = " << control.ispin() << endl;
          std::cout << "ne = " << control.ne(0) << " " << control.ne(1) << endl;
          std::cout << "ne = " << control.ne_ptr()[0] << " " <<
      control.ne_ptr()[1] << endl;
      }
   */
 
   /* fetch fractional control variables */
   fractional = control.fractional();
   if (fractional)
   {
      nextra[0] = control.fractional_orbitals(0);
      if (control.ispin()==2)
         nextra[1] = control.fractional_orbitals(1);
      else
         nextra[1] = 0;
      smearcorrection = 0.0;
      smeartype = control.fractional_smeartype();
      smearkT   = control.fractional_kT();
   }
   else
   {
      nextra[0] = 0;
      nextra[1] = 0;
   }

   /* fetch ispin and ne psi information from control */
   ispin = control.ispin();
   ne[0] = control.ne(0) + nextra[0];
   ne[1] = control.ne(1) + nextra[1];

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
   version = 3;
 
   /* initialize parallel grid structure */
   //Pneb mygrid(&myparallel, &mylattice, control, control.ispin(),control.ne_ptr());
   Pneb mygrid(&myparallel,&mylattice,control,ispin,ne);
 
   /* initialize psi1 and psi2 */
   psi1 = mygrid.g_allocate(1);
   psi2 = mygrid.g_allocate(1);
   Hpsi = mygrid.g_allocate(1);
   psi_r = mygrid.h_allocate();
   dn = mygrid.r_nalloc(ispin);
   hml = mygrid.m_allocate(-1,1);
   lmbda = mygrid.m_allocate(-1,1);
   eig = new double[ne[0] + ne[1]];
   eig_prev = new double[ne[0] + ne[1]];
   if (fractional)
   {
      occ1 = mygrid.initialize_occupations_with_allocation(nextra); 
      occ2 = mygrid.initialize_occupations_with_allocation(nextra); 
   }

   mygrid.d3db::mygdevice.psi_alloc(mygrid.npack(1),mygrid.neq[0]+mygrid.neq[1],control.tile_factor());
 
   // psi_read(&mygrid,&version,nfft,unita,&ispin,ne,psi2,control.input_movecs_filename());
   bool newpsi = psi_read(&mygrid,control.input_movecs_filename(),
                          control.input_movecs_initialize(),psi2,
                          &smearoccupation,occ2,std::cout);
   if (fractional)
   {
      smearoccupation = 1;
      std::vector<double> filling = control.fractional_filling();
      if (filling.size() > 0)
      {
         int sz = filling.size();
         if (sz > (ne[0]+ne[1])) sz = ne[0]+ne[1];
         std::memcpy(occ2,filling.data(),sz*sizeof(double));
      }

      std::memcpy(occ1,occ2,(ne[0]+ne[1])*sizeof(double));
   }
   
   MPI_Barrier(comm_world0);
 
   /* setup structure factor */
   Strfac mystrfac(&myion,&mygrid);
   mystrfac.phafac();
 
   /* initialize operators */
   Kinetic_Operator mykin(&mygrid);
   Coulomb12_Operator mycoulomb12(&mygrid, control);
   mycoulomb12.initialize_dielectric(&myion,&mystrfac);
   XC_Operator myxc(&mygrid, control);
   HFX_Operator myhfx(&mygrid, mycoulomb12.has_coulomb2, mycoulomb12.mycoulomb2, control);
 
   /* initialize psps */
   Pseudopotential mypsp(&myion,&mygrid,&mystrfac,control,std::cout);
 
   /* setup ewald */
   Ewald myewald(&myparallel,&myion,&mylattice,control,mypsp.zv);
   myewald.phafac();
 
   /* setup dfpt */
   // DFPT_Operators mydfpt(&mygrid,&mykin,&mycoulomb12,&myxc,&mypsp,1.0e-3);
 
   //                 |**************************|
   // *****************   summary of input data  **********************
   //                 |**************************|
 
   if (oprint) 
   {
      std::cout << std::endl;
      std::cout << "     ===================  summary of input  =======================" << std::endl;
      std::cout << "\n input psi filename: " << control.input_movecs_filename() << std::endl;
      std::cout << std::endl;
      std::cout << " number of processors used: " << myparallel.np() << std::endl;
      std::cout << " processor grid           : " << myparallel.np_i() << " x " << myparallel.np_j() << std::endl;
      if (mygrid.maptype == 1) std::cout << " parallel mapping         : 1d-slab" << std::endl;
      if (mygrid.maptype == 2) std::cout << " parallel mapping         : 2d-hilbert" << std::endl;
      if (mygrid.maptype == 3) std::cout << " parallel mapping         : 2d-hcurve" << std::endl;
      if (mygrid.isbalanced())
         std::cout << " parallel mapping         : balanced" << std::endl;
      else
         std::cout << " parallel mapping         : not balanced" << std::endl;

      if (mygrid.d3db::mygdevice.has_gpu())
      {
         //std::cout << " parallel mapping         : has GPU" << std::endl;
         if (mygrid.d3db::mygdevice.type_gpu()==1) std::cout << " parallel mapping         : CUDA GPUs" << std::endl;
         if (mygrid.d3db::mygdevice.type_gpu()==2) std::cout << " parallel mapping         : SYCL GPUs" << std::endl;
         if (mygrid.d3db::mygdevice.type_gpu()==3) std::cout << " parallel mapping         : HIP SYCL GPUs" << std::endl;
         if (mygrid.d3db::mygdevice.type_gpu()==4) std::cout << " parallel mapping         : OpenCL GPUs" << std::endl;
         if (mygrid.staged_gpu_fft_pipeline) std::cout << " parallel mapping         : staged GPU FFT" << std::endl;
         if (control.tile_factor() > 0)      std::cout << " GPU tile factor          : " << control.tile_factor() << std::endl;
      }
      std::cout << " fft container size       : " << control.fft_container_size() << std::endl;
      
      std::cout << "\n options:\n";
      std::cout << "   ion motion           = ";
      if (control.geometry_optimize())
         std::cout << "yes\n";
      else
         std::cout << "no\n";
      std::cout << "   boundary conditions  = ";
      if (control.version == 3) std::cout << "periodic\n";
      if (control.version == 4) std::cout << "aperiodic\n";
      
      std::cout << "   electron spin        = ";
      if (ispin == 1)
         std::cout << "restricted\n";
      else
         std::cout << "unrestricted\n";

      std::cout << myxc;

      if (fractional) std::cout << "   using fractional" << std::endl;

      std::cout << myhfx;
     
      std::cout << mypsp.print_pspall();
      
      std::cout << "\n total charge =" << Ffmt(8,3) << control.total_charge() << std::endl;
      
      std::cout << "\n atom composition:" << std::endl;
      for (ia = 0; ia < myion.nkatm; ++ia)
         std::cout << "   " << myion.atom(ia) << " : " << myion.natm[ia];
      std::cout << "\n\n initial ion positions (au):" << std::endl;
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
      
      std::cout << mypsp.myefield->shortprint_efield();
      std::cout << mycoulomb12.shortprint_dielectric();
      
      std::cout << std::endl;
      if (fractional)
      {
         double oen[ispin];
         int n=0;
         for (auto ms=0; ms<ispin; ++ms)
         {
            oen[ms] = 0;
            for (auto i=0; i<ne[ms]; ++i)
            {
               oen[ms] += occ1[n];
               ++n;
            }
         }
         std::cout << " number of electrons: spin up ="  
                   << Ffmt(6,2)  << oen[0] << "  " 
                   << Ffmt(27,2) << oen[ispin-1] << " (   fractional)" << std::endl;
      }
      else
         std::cout << " number of electrons: spin up ="  
                   << Ifmt(6) << mygrid.ne[0] << " (" 
                   << Ifmt(4) << mygrid.neq[0] << " per task) down =" 
                   << Ifmt(6) << mygrid.ne[ispin-1] << " ("
                   << Ifmt(4) << mygrid.neq[ispin-1] << " per task)" << std::endl;

      std::cout << " number of orbitals:  spin up ="  
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
      
      std::cout << std::endl;
      std::cout << " technical parameters:\n";
      if (control.io_buffer()) std::cout << "      using io buffer " << std::endl;
      std::cout << "      fixed step: time step =" << Ffmt(12,2) << control.time_step() 
                << "  ficticious mass =" << Ffmt(12,2) << control.fake_mass() << std::endl;
      std::cout << "      tolerance =" << Efmt(12,3) << control.tolerances(0)
                << " (energy) " << Efmt(12,3) << control.tolerances(1)
                << " (density) " << Efmt(12,3) << control.tolerances(2)
                << " (ion)\n";
      std::cout << "      max iterations = " << Ifmt(10) << control.loop(0)*control.loop(1) 
                << " (" << Ifmt(5) << control.loop(0) << " inner " 
                        << Ifmt(5) << control.loop(1) << " outer)" << std::endl;
      if (!control.deltae_check()) std::cout << "      allow DeltaE > 0" << std::endl;

      if (control.fractional())
      {
         std::cout << std::endl;
         std::cout << " fractional smearing parameters:" << std::endl;
         std::cout << "    smearing algorithm = " << smeartype << std::endl;
         std::cout << "    smearing parameter = ";
         if (smeartype==-1) std::cout << "fixed_occupation" << std::endl;
         if (smeartype==0) std::cout << "step function" << std::endl;
         if (smeartype==1) std::cout << "Fermi-Dirac" << std::endl;
         if (smeartype==2) std::cout << "Gaussian" << std::endl;
         if (smeartype==3) std::cout << "Hermite" << std::endl;
         if (smeartype==4) std::cout << "Marazari-Vanderbilt" << std::endl;
         if (smeartype==5) std::cout << "Methfessel-Paxton" << std::endl;
         if (smeartype==6) std::cout << "Cold smearing" << std::endl;
         if (smeartype==7) std::cout << "Lorentzian" << std::endl;
         if (smeartype==8) std::cout << "step" << std::endl;
         if (smeartype>=0) 
         {
            std::cout <<  "    smearing parameter = " << Ffmt(9,3) << smearkT 
                                              << " (" << Ffmt(7,1) << control.fractional_temperature() << " K)" <<  std::endl;
            std::cout <<  "    mixing parameter(alpha)   = " << Ffmt(7,2) << control.fractional_alpha() << std::endl;
            std::cout <<  "    mixing parameter(alpha_min)   = " << Ffmt(7,2) << control.fractional_alpha_min() << std::endl;
            std::cout <<  "    mixing parameter(alpha_max)   = " << Ffmt(7,2) << control.fractional_alpha_max() << std::endl;
            std::cout <<  "    mixing parameter(beta)   = " << Ffmt(7,2) << control.fractional_beta() << std::endl;
            std::cout <<  "    mixing parameter(gamma)   = " << Ffmt(7,2) << control.fractional_gamma() << std::endl;
            std::cout <<  "    rmsd occupation tolerance   = " << Efmt(12,3) << control.fractional_rmsd_tolerance() << std::endl;
            if (ispin==2) 
               std::cout <<  "    extra orbitals     : up=" << nextra[0] << " down= " << nextra[1] << std::endl;
            else
               std::cout <<  "    extra orbitals     = " << Ifmt(7) << nextra[0] << std::endl;
            if (control.fractional_frozen()) 
               std::cout <<  "    frozen orbitals" << std::endl;
            //if (control.fractional_frozen())
            {
               size_t total_occupations = ne[0] + ne[1];
               std::ostringstream oss;
               oss << "    initial occupations = [";
               for (auto i=0; i<total_occupations-1; ++i) 
                  oss << occ1[i] << " ";
               oss << occ1[total_occupations - 1] << "]";
               std::cout << oss.str() << std::endl;
            }
         }

      }
   }
 
   //                 |**************************|
   // *****************   start iterations       **********************
   //                 |**************************|
   double rmsd_occupation = 0.0;
   bool occupation_update = fractional && !control.fractional_frozen();

   if (myparallel.is_master()) seconds(&cpu2);
   if (oprint) 
   {
      std::cout << std::endl << std::endl << std::endl;
      std::cout << "     ========================== iteration ==========================" << std::endl;
      std::cout << "          >>> iteration started at " << util_date() << "  <<<" << std::endl;
      std::cout << "     iter.             Energy       DeltaE     DeltaPsi     DeltaIon" << std::endl;
      std::cout << "     ---------------------------------------------------------------" << std::endl;
   }

   if (control.loop(1) > 0) {

      // Initialize AIMD running data
      nwpw_aimd_running_data mymotion_data(control,&myparallel,&mylattice,&myion,E,hml,psi1,dn);

      done   = 0;
      icount = 0;
      while (!done) 
      {
         ++icount;
         inner_loop(control, &mygrid, &myion, &mykin, &mycoulomb12, &myxc, &mypsp,
                    &mystrfac, &myewald, psi1, psi2, Hpsi, psi_r, dn, hml, lmbda, E,
                    &deltae, &deltac, &deltar,
                    fractional, occ1, occ2);
         //if (fractional) E[0] += smearcorrection;
        
         // mydfpt.start(psi1,psi_r
     
         // Update AIMD Running data
         if (control.geometry_optimize()) mymotion_data.update_iteration(icount);
     
        
         if (oprint)
            std::cout << Ifmt(10) << icount*control.loop(0) 
                      << Efmt(19,10) << E[0]
                      << Efmt(13,5) << deltae 
                      << Efmt(13,5) << deltac 
                      << Efmt(13,5) << deltar;
        
         /* check for competion */
         if ((!fractional) && ((deltae > 0.0) && (icount > 1) && control.deltae_check()))
         {
            done = 1;
            if (oprint) std::cout << std::endl << "         *** Energy going up. iteration terminated";
         } 
         else if ((std::fabs(deltae) < control.tolerances(0)) &&
                   (deltac < control.tolerances(1)) &&
                   (deltar < control.tolerances(2)) && 
                   (rmsd_occupation < control.fractional_rmsd_tolerance()))
         {
            done = 1;
            if (oprint) std::cout << std::endl << "         *** tolerance ok.    iteration terminated";
         } 
         else if (icount >= control.loop(1)) 
         {
            done = 1;
            if (oprint) std::cout << std::endl << "          *** arrived at the Maximum iteration.    terminated ***";
         }

         // update occupations
         if ((fractional) && (!control.fractional_frozen()))
         {
            mygrid.m_diagonalize(hml, eig);
         
            if ((icount < control.loop(1)) && (smeartype>=0) && (occupation_update) && (not done))
            {
               double alpha = control.fractional_alpha(); // alpha must be between 0 and 1
               double gamma = control.fractional_gamma(); // gamma must be between 0 and 1
         
               // Adaptive alpha parameters
               double alpha_min = control.fractional_alpha_min();
               double alpha_max = control.fractional_alpha_max();
               double beta = control.fractional_beta();
               double rmsd_threshold = control.fractional_rmsd_threshold();  
               //if (alpha < alpha_min) alpha_min =  alpha;
               //if (alpha > alpha_max) alpha_max =  alpha;
         
               if (icount == 0)  // Initialize eig_prev for the first iteration
                  std::memcpy(eig_prev,eig,(ne[0]+ne[1])*sizeof(double));
               else             // Smooth eigenvalues in subsequent iterations
                  for (size_t i=0; i<(ne[0]+ne[1]); ++i)
                     eig_prev[i] = (1.0-gamma)*eig_prev[i] + gamma*eig[i];
         
               //ZZ = myion.total_zv() - control.total_charge()
         
               // Define occupations based on smoothed eigenvalues
               smearcorrection_old = smearcorrection;
               mygrid.m_0define_occupation(-1.0, false,
                             control.multiplicity(),
                             myion.total_zv(), control.total_charge(),
                             eig_prev, hml, occ2,
                             smeartype,smearkT,smearfermi,&smearcorrection);
         
               // RMSD occupation computation
               double rmsd_occupation = 0.0;
               for (size_t i=0; i<(ne[0] + ne[1]); ++i)
               {
                   double delta_occ = occ2[i] - occ1[i];
                   rmsd_occupation += delta_occ * delta_occ;
               }
               rmsd_occupation = std::sqrt(rmsd_occupation / (ne[0] + ne[1]));
         
               // Adaptive alpha adjustment
              //std::cout << "alpha=" << alpha << " rmsd_occupation=" << rmsd_occupation << " rmsd_threshold" << rmsd_threshold << std::endl;
               if (rmsd_occupation < rmsd_threshold)  // Converging well
               {
                  alpha = std::min(alpha_max, alpha * (1.0 + beta));
                  //std::cout << "f alpha = " << alpha << " alpha_max =" << alpha_max <<  std::endl;
               }
               else  // Oscillations or divergence
                  alpha = std::max(alpha_min, alpha * (1.0 - beta));
         
               // Update occupations
               if (icount>0)
               {
                  for (size_t i=0; i<(ne[0]+ne[1]); ++i)
                     occ1[i] = (1.0-alpha)*occ1[i] + alpha*occ2[i];
               }
               else
               {
                  std::memcpy(occ1,occ2,(ne[0]+ne[1])*sizeof(double));
               }
            
         
                E[28] = smearcorrection;
                //E[0] += E[28];
         
                // Debugging output (optional)
                if (oprint)
                   std::cout << " Iteration: " << icount
                             << ", RMSD: " << rmsd_occupation
                             << ", Alpha: " << alpha
                             << ", Smear Correction: " << smearcorrection 
                             << ", Delta Smear Correction: " << smearcorrection - smearcorrection_old;
         
                //occupation_update  =  (rmsd_occupation > control.fractional_rmsd_tolerance());
            }
            else 
            {
               smearfermi[0]   =  mygrid.define_smearfermi(ne[0],eig,occ1);
               smearcorrection =  mygrid.add_smearcorrection(smeartype,ne[0],eig,occ1,smearfermi[0],smearkT);
               if (ispin==1)
               {
                  smearcorrection *= 2.0;
               }
               else
               {
                  smearfermi[1]   =  mygrid.define_smearfermi(ne[1],eig+ne[0],occ1+ne[0]);
                  smearcorrection +=  mygrid.add_smearcorrection(smeartype,ne[1],eig+ne[0],occ1+ne[0],smearfermi[1],smearkT);
               }
               E[28] = smearcorrection;
               //E[0] += E[28];
            }
         }
         if (oprint) std::cout << std::endl;
      }
   }
   if (myparallel.is_master()) seconds(&cpu3);
   if (oprint) { std::cout << "          >>> iteration ended at   " << util_date() << "  <<<\n"; }

   /* diagonalize the hamiltonian */
   if ((!fractional) || (control.fractional_frozen()))
      mygrid.m_diagonalize(hml, eig);
 
   /* calculate real-space number of electrons, en */
   {
      double omega = mylattice.omega();
      double scal1 = 1.0/((double)((mygrid.nx) * (mygrid.ny) * (mygrid.nz)));
      double dv = omega * scal1;
     
      en[0] = dv*mygrid.r_dsum(dn);
      en[1] = en[0];
      if (ispin>1)
         en[1] = dv*mygrid.r_dsum(&dn[mygrid.n2ft3d]);
   }
 
   // calculate APC if v_apc_on==false
   if (mypsp.myapc->apc_on)
     if (!(mypsp.myapc->v_apc_on))
       mypsp.myapc->dngen_APC(dn, false);
 
   mypsp.mydipole->gen_dipole(dn);
 
   //                  |***************************|
   // ****************** report summary of results **********************
   //                  |***************************|
   if (oprint)
   {
      std::cout << std::endl << std::endl;
      std::cout << "     ===================  summary of results ======================" << std::endl;
      std::cout << "\n final ion positions (au):" << std::endl;
      for (ii = 0; ii < myion.nion; ++ii)
        std::cout << Ifmt(4) << ii + 1 << " " << myion.symbol(ii) << "\t( "
                  << Ffmt(10,5) << myion.rion1[3*ii] << " " 
                  << Ffmt(10,5) << myion.rion1[3*ii+1] << " " 
                  << Ffmt(10,5) << myion.rion1[3*ii+2] << " ) - atomic mass = " 
                  << Ffmt(6,3) << myion.amu(ii) << std::endl;
      std::cout << "   G.C.\t( " << Ffmt(10,5) << myion.gc(0) << " "
                                 << Ffmt(10,5) << myion.gc(1) << " " 
                                 << Ffmt(10,5) << myion.gc(2) << " )" << std::endl;
      std::cout << " C.O.M.\t( " << Ffmt(10,5) << myion.com(0) << " "
                                 << Ffmt(10,5) << myion.com(1) << " " 
                                 << Ffmt(10,5) << myion.com(2) << " )" << std::endl;

      if (control.geometry_optimize())
      {
         std::cout << "\n final ion forces (au):" << std::endl;
         for (ii = 0; ii < myion.nion; ++ii)
            std::cout << Ifmt(4) << ii+1 << " " << myion.symbol(ii) << "\t( "
                      << Ffmt(10,5) << myion.fion1[3*ii]   << " " 
                      << Ffmt(10,5) << myion.fion1[3*ii+1] << " " 
                      << Ffmt(10,5) << myion.fion1[3*ii+2] << " )" << std::endl;

         std::cout << "   |F|/nion  = " << Ffmt(12,6) << myion.rms_fion(myion.fion1) << std::endl
                   << "   max|Fatom|= " << Ffmt(12,6) << myion.max_fion(myion.fion1) << "  ("
                                        << Ffmt(8,3) << myion.max_fion(myion.fion1) * (27.2116/0.529177)
                                        << " eV/Angstrom)" << std::endl << std::endl;
      }
 
      // if (mypsp.myapc->v_apc_on)
      //    std::cout << mypsp.myapc->shortprint_APC();

      std::cout << std::endl;
      std::cout << std::fixed << " number of electrons: spin up= " 
                << Ffmt(11,5) << en[0] << "  down= " 
                << Ffmt(11,5) << en[ispin-1] << " (real space)";
      std::cout << std::endl << std::endl;
      double E0 = E[0];
      if   (fractional) E0 += E[28];
      std::cout << " total     energy    : " << Efmt(19,10) << E0 << " (" 
                << Efmt(15,5) << E0/myion.nion << " /ion)" << std::endl;

      if (mypsp.myefield->efield_on) {
        std::cout << std::endl;
        std::cout << " QM Energies" << std::endl;
        std::cout << " -----------" << std::endl;
        std::cout << " total  QM energy    : " << Efmt(19,10) << (E[0]-E[48]-E[49]) 
                  << " (" << Efmt(15,5) << (E[0]-E[48]-E[49])/(mygrid.ne[0]+mygrid.ne[1])
                  << " /electron)" << std::endl;
       }

      std::cout << " total orbital energy: " 
                << Efmt(19,10) << E[1] << " ("
                << Efmt(15,5) << E[1]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      std::cout << " hartree energy      : " 
                << Efmt(19,10) << E[2] << " ("
                << Efmt(15,5) << E[2]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      std::cout << " exc-corr energy     : " 
                << Efmt(19,10) << E[3] << " ("
                << Efmt(15,5) << E[3]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;

      if (mycoulomb12.dielectric_on())
         std::cout << " dielectric energy   : " 
                   << Efmt(19,10) << E[61] << " ("
                   << Efmt(15,5) << E[61]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;

      if (mypsp.myapc->v_apc_on)
         std::cout << " APC energy          : " 
                   << Efmt(19,10) << E[51] << " ("
                   << Efmt(15,5) << E[51]/myion.nion << " /ion)" << std::endl;

      std::cout << " ion-ion energy      : " 
                << Efmt(19,10) << E[4] << " ("
                << Efmt(15,5) << E[4]/myion.nion << " /ion)" << std::endl;


      if ((fractional) && (!control.fractional_frozen()))
         std::cout << " smearing energy     : " 
                << Efmt(19,10) << E[28] << " ("
                << Efmt(15,5) << E[28]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;

      std::cout << std::endl;
      std::cout << " K.S. kinetic energy : " 
                << Efmt(19,10) << E[5] << " ("
                << Efmt(15,5) << E[5]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      std::cout << " K.S. V_l energy     : " 
                << Efmt(19,10) << E[6] << " (" << Efmt(15,5) << E[6]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      std::cout << " K.S. V_nl_energy    : " 
                << Efmt(19,10) << E[7] << " ("
                << Efmt(15,5) << E[7]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      std::cout << " K.S. V_Hart energy  : " 
                << Efmt(19,10) << E[8] << " ("
                << Efmt(15,5) << E[8]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      std::cout << " K.S. V_xc energy    : " 
                << Efmt(19,10) << E[9] << " ("
                << Efmt(15,5) << E[9]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;

      if (mycoulomb12.dielectric_on())
         std::cout << " K.S. V_dielec energy: " 
                   << Efmt(19,10) << E[62] << " ("
                   << Efmt(15,5) << E[62]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;

      if (mypsp.myapc->v_apc_on)
         std::cout << " K.S. V_APC energy   : " 
                   << Efmt(19,10) << E[52] << " ("
                   << Efmt(15,5) << E[52]/myion.nion << " /ion)" << std::endl;

      viral = (E[9] + E[8] + E[7] + E[6]) / E[5];
      std::cout << " Viral Coefficient   : " 
                << Efmt(19,10) << viral << std::endl;

      if ((fractional) && (!control.fractional_frozen()))
      {
         std::cout << std::endl;
         ev = 27.2116;
         if (ispin==1)
            std::cout << " Fermi energy        : " 
                      << Efmt(19,10) << smearfermi[0] << " (" << Ffmt(8,3) << smearfermi[0]*ev << " eV)" << std::endl;
         else
            std::cout << "  Fermi energy = " 
                      << Efmt(19,10) << smearfermi[0] << " (" << Ffmt(8,3) << smearfermi[0]*ev << " eV)"
                      << Efmt(19,10) << smearfermi[1] << " (" << Ffmt(8,3) << smearfermi[1]*ev << " eV)" << std::endl;
      }
      if (myion.has_ion_constraints())
      {
         std::cout << std::endl;
         if (myion.has_ion_bond_constraints())     
            std::cout << " spring bond         : " << Efmt(19,10) << E[70] << " ("
                                                   << Efmt(15,5)  << E[70]/myion.nion << " /ion)" << std::endl;
         if (myion.has_ion_bondings_constraints()) 
            std::cout << " spring bondings     : " << Efmt(19,10) << E[71] << " ("
                                                   << Efmt(15,5)  << E[71]/myion.nion << " /ion)" << std::endl;
      }

      if (mypsp.myefield->efield_on) {
         std::cout << std::endl;
         std::cout << " Electric Field Energies" << std::endl;
         std::cout << " -----------------------" << std::endl;
         std::cout << " - Electric Field Energy   : " << Efmt(19,10) << E[48] + E[49] << std::endl;
         std::cout << " - Electric Field/Electron : " << Efmt(19,10) << E[48] << std::endl;
         std::cout << " - Electric Field/Ion      : " << Efmt(19,10) << E[49] << std::endl;
      }

      std::cout << "\n orbital energies:\n";
      nn = ne[0] - ne[1];
      ev = 27.2116;
      for (i=0; i<nn; ++i) 
      {
         if (fractional)
            if ((occ1[i] < 1.e-3) && (occ1[i]>1.0e-12))
               std::cout << Efmt(18,7) << eig[i] << " (" << Ffmt(8,3) << eig[i] * ev << "eV) occ=" 
                         << Efmt(9,3) << occ1[i] << std::endl;
            else
               std::cout << Efmt(18,7) << eig[i] << " (" << Ffmt(8,3) << eig[i] * ev << "eV) occ=" 
                         << Ffmt(5,3) << occ1[i] << std::endl;
         else
            std::cout << Efmt(18,7) << eig[i] << " (" << Ffmt(8,3) << eig[i] * ev << "eV)" << std::endl;
      }
      for (i=0; i<ne[1]; ++i)
      {
         if (fractional) 
         {
            if ((occ1[i+nn] < 1.e-3) && (occ1[i+nn]>1.0e-12)) 
               std::cout << Efmt(18,7) << eig[i+nn]    << " (" 
                         << Ffmt(8,3)  << eig[i+nn]*ev << "eV)  occ=" 
                         << Efmt(9,3)  << occ1[i+nn]   << " ";
                      //<< Efmt(18,7) << eig[i+(ispin-1)*ne[0]]    << " (" 
                      //<< Ffmt(8,3)  << eig[i+(ispin-1)*ne[0]]*ev << "eV) occ="  
                      //<< Efmt(9,3)  << occ1[i + (ispin-1)*ne[0]] << std::endl;
            else
               std::cout << Efmt(18,7) << eig[i+nn]    << " (" 
                         << Ffmt(8,3)  << eig[i+nn]*ev << "eV)  occ=" 
                         << Ffmt(5,3)  << occ1[i+nn]   << " ";

            if ((occ1[i+(ispin-1)*ne[0]] < 1.e-3) && (occ1[i+(ispin-1)*ne[0]]>1.0e-12))
               std::cout << Efmt(18,7) << eig[i+(ispin-1)*ne[0]]    << " (" 
                         << Ffmt(8,3)  << eig[i+(ispin-1)*ne[0]]*ev << "eV) occ="  
                         << Efmt(9,3)  << occ1[i + (ispin-1)*ne[0]] << std::endl;
            else
               std::cout << Efmt(18,7) << eig[i+(ispin-1)*ne[0]]    << " (" 
                         << Ffmt(8,3)  << eig[i+(ispin-1)*ne[0]]*ev << "eV) occ="  
                         << Ffmt(5,3)  << occ1[i + (ispin-1)*ne[0]] << std::endl;
         }
         else
            std::cout << Efmt(18,7) << eig[i+nn] << " (" 
                      << Ffmt(8,3)  << eig[i + nn] * ev << "eV) " 
                      << Efmt(18,7) << eig[i+(ispin-1)*ne[0]] << " (" 
                      << Ffmt(8,3)  << eig[i+(ispin-1)*ne[0]]*ev << "eV)" << std::endl;
         // printf("%18.7le",eig[i+nn]); printf(" (");
         // printf("%8.3lf",eig[i+nn]*ev); printf("eV) ");
         // printf("%18.7le",eig[i+(ispin-1)*ne[0]]); printf(" (");
         // printf("%8.3lf",eig[i+(ispin-1)*ne[0]]*ev); printf("eV)\n");
      }
      std::cout << std::endl;

      // write geometry and constraints analysis
      if (control.geometry_optimize())
      {
         std::cout << myion.print_bond_angle_torsions();
         std::cout << std::endl << myion.print_constraints(1);
      }

      // write APC analysis
      if (mypsp.myapc->apc_on)
         std::cout << mypsp.myapc->print_APC(mypsp.zv);

      // write dipoles
      std::cout << mypsp.mydipole->shortprint_dipole();

   }

   psi_write(&mygrid,&version,nfft,unita,&ispin,ne,psi1,&smearoccupation,occ1,control.output_movecs_filename(),std::cout);
   MPI_Barrier(comm_world0);

   /* deallocate memory */
   mygrid.g_deallocate(psi1);
   mygrid.g_deallocate(psi2);
   mygrid.g_deallocate(Hpsi);
   mygrid.h_deallocate(psi_r);
   mygrid.r_dealloc(dn);
   mygrid.m_deallocate(hml);
   mygrid.m_deallocate(lmbda);
   delete[] eig;
   delete[] eig_prev;
   if (fractional)
   {
      delete[] occ1;
      delete[] occ2;
   }
   mygrid.d3db::mygdevice.psi_dealloc();

   // write results to the json
   auto rtdbjson = json::parse(rtdbstring);
   rtdbjson["pspw"]["energy"] = E[0];
   rtdbjson["pspw"]["energies"] = E;
   if (mypsp.myapc->apc_on) 
   {
      double qion[myion.nion];
      for (auto ii = 0; ii < myion.nion; ++ii)
         qion[ii] = -mypsp.myapc->Qtot_APC(ii) + mypsp.zv[myion.katm[ii]];
      rtdbjson["nwpw"]["apc"]["q"] = std::vector<double>(qion, &qion[myion.nion]);
   }
   
   if (mypsp.mydipole->dipole_on) 
   {
      double *mdipole = mypsp.mydipole->mdipole;
      double mu = std::sqrt(mdipole[0]*mdipole[0] + mdipole[1]*mdipole[1] + mdipole[2]*mdipole[2]);
      rtdbjson["nwpw"]["dipole"] = std::vector<double>(mdipole, &mdipole[3]);
      rtdbjson["nwpw"]["dipole_magnitude"] = mu;
   }
   
   // set rtdbjson initialize_wavefunction option to false
   if (rtdbjson["nwpw"]["initialize_wavefunction"].is_boolean())
      rtdbjson["nwpw"]["initialize_wavefunction"] = false;
   
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
      double av = t2 / ((double)control.loop(0) * icount);
      // std::cout.setf(ios::scientific);
      std::cout << std::scientific;
      std::cout << std::endl;
      std::cout << " -----------------" << std::endl;
      std::cout << " cputime in seconds" << std::endl;
      std::cout << " prologue    : " << t1 << std::endl;
      std::cout << " main loop   : " << t2 << std::endl;
      std::cout << " epilogue    : " << t3 << std::endl;
      std::cout << " total       : " << t4 << std::endl;
      std::cout << " cputime/step: " << av << std::endl;
      std::cout << std::endl;
    
      nwpw_timing_print_final(control.loop(0) * icount, std::cout);
    
      std::cout << std::endl;
      std::cout << " >>> job completed at     " << util_date() << " <<<" << std::endl;
   }
   
   MPI_Barrier(comm_world0);
   return 0;
}

} // namespace pwdft

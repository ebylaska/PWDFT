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
#include "band_inner_loop.hpp"
#include "nwpw_aimd_running_data.hpp"
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

   int version, nfft[3], ne[2], nextra[2], ispin, nbrillq, nbrillouin;
   int i, ii, ia, nn, ngrid[3], matype, nelem, icount, done;
   char date[26];
   double sum1, sum2, ev, zv;
   double cpu1, cpu2, cpu3, cpu4;
   double E[80], deltae, deltac, deltar, viral, unita[9], en[2];
   double *psi1, *psi2, *Hpsi, *psi_r;
   double *occ1 = nullptr;
   double *occ2 = nullptr;
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
   ispin = control.ispin();
   ne[0] = control.ne(0) + nextra[0];
   ne[1] = control.ne(1) + nextra[1];



   /* read in Brillouin zone */
   Brillouin mybrillouin(rtdbstring,&mylattice,control);
   //control.set_total_ion_charge(8);

   //std::cout << "ispin=" << control.ispin() << " ne=" << control.ne_ptr()[0] << " " << control.ne_ptr()[1] 
   //          << " nbrillioun=" << mybrillouin.nbrillouin << std::endl;
   /* initialize parallel grid structure */
   //Cneb mygrid(&myparallel, &mylattice, control, control.ispin(),control.ne_ptr(),&mybrillouin);
   Cneb mygrid(&myparallel, &mylattice, control, ispin,ne,&mybrillouin);

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
   eig      = new double[nbrillq*(ne[0]+ne[1])];
   eig_prev = new double[nbrillq*(ne[0]+ne[1])];
   if (fractional)
   {
      occ1 = mygrid.initialize_occupations_with_allocation(nextra);
      occ2 = mygrid.initialize_occupations_with_allocation(nextra);
   }

   mygrid.c3db::mygdevice.psi_alloc(mygrid.npack1_max(),mygrid.neq[0]+mygrid.neq[1],control.tile_factor());


   // psi_read(&mygrid,&version,nfft,unita,&ispin,ne,nbrill,psi2,control.input_movecs_filename());
   bool newpsi = cpsi_read(&mygrid,control.input_movecs_filename(),control.input_movecs_initialize(),
                           psi2,&smearoccupation,occ2,std::cout);
   mygrid.gg_copy(psi2,psi1);
   if (fractional)
   {
      smearoccupation = 1;
      std::vector<double> filling = control.fractional_filling();
      if (filling.size() > 0)
      {
         int sz = filling.size();
         if (sz > (ne[0]+ne[1])) sz = (ne[0]+ne[1]);
         for (auto nbq=0; nbq<nbrillq; ++nbq)
         {
            int ishift = (ne[0]+ne[1]);
            std::memcpy(occ2+ishift,filling.data(),sz*sizeof(double));
         }
      }

      std::memcpy(occ1,occ2,nbrillq*(ne[0]+ne[1])*sizeof(double));
   }
   MPI_Barrier(comm_world0);

   /* setup structure factor */
   CStrfac mystrfac(&myion,&mygrid);
   mystrfac.phafac();
   mystrfac.phafac_k();

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

      if (mygrid.c3db::mygdevice.has_gpu())
      {
         //std::cout << " parallel mapping         : has GPUs" << std::endl;
         if (mygrid.c3db::mygdevice.type_gpu()==1) std::cout << " parallel mapping         : CUDA GPUs" << std::endl;
         if (mygrid.c3db::mygdevice.type_gpu()==2) std::cout << " parallel mapping         : SYCL GPUs" << std::endl;
         if (mygrid.c3db::mygdevice.type_gpu()==3) std::cout << " parallel mapping         : HIP SYCL GPUs" << std::endl;
         if (mygrid.c3db::mygdevice.type_gpu()==4) std::cout << " parallel mapping         : OpenCL GPUs" << std::endl;
         if (mygrid.staged_gpu_fft_pipeline) std::cout << " parallel mapping         : staged GPU FFT" << std::endl;
         if (control.tile_factor() > 0)      std::cout << " GPU tile factor          : " << control.tile_factor() << std::endl;
      }

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

      if (control.fractional()) std::cout << "   using fractional" << std::endl;


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
   }

   double oen[2] = {0.0,0.0};
   if (control.fractional())
   {
      int n=0;

      //The outer loop nbq is split over dimension 3
      for (auto nbq=0; nbq<mygrid.nbrillq; ++nbq)
      {
         double weight =  mygrid.pbrill_weight(nbq);
         for (auto ms=0; ms<ispin; ++ms)
         {
            //The orbital loop index i is split over dimension 2
            for (auto i=0; i<mygrid.neq[ms]; ++i)
            {
               oen[ms] += weight*occ1[n];
               ++n;
            }
         }
     }
     myparallel.Vector_SumAll(2,2,oen); // Reduce over orbitals (orbital distribution)
     myparallel.Vector_SumAll(3,2,oen); // Then reduce over BZ points (k-point distribution)
     // myparallel.Reduce_Valuese(2,0,2,oen_distribute,oen);
     // myparallel.Reduce_Valuese(2,0,2,oen_distribute,oen);
   }

   if (oprint)
   {
      if (control.fractional())
         std::cout << " number of electrons: spin up ="
                   << Ffmt(6,2)  << oen[0] << "  "
                   << Ffmt(27,2) << oen[ispin-1] << " (   fractional)" << std::endl;
      else
         std::cout << " number of electrons: spin up =" << Ifmt(6) << mygrid.ne[0]
                   << " (" << Ifmt(4) << mygrid.neq[0]
                   << " per task) down =" << Ifmt(6) << mygrid.ne[ispin-1] << " ("
                   << Ifmt(4) << mygrid.neq[ispin - 1] << " per task)" << std::endl;

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

       /* print nbrillouin */
      std::cout << std::endl;
      std::cout << " computational grids:" << std::endl;
      std::cout << "      density     cutoff ="
                << Ffmt(7,3) << mylattice.ecut()
                << " fft =" << Ifmt(4) << mygrid.nx << " x "
                            << Ifmt(4) << mygrid.ny << " x "
                            << Ifmt(4) << mygrid.nz
                << "  (" << Ifmt(8) << mygrid.npack_all(0) << " waves "
                         << Ifmt(8) << mygrid.npack(0) << " per task)" << std::endl;
      for (auto nb=0; nb<mygrid.nbrillouin; ++nb)
         std::cout << "      wavefnc"  
                   << Ifmt(4) << nb+1
                   << " cutoff ="
                   << Ffmt(7,3) << mylattice.wcut()
                   << " fft =" << Ifmt(4) << mygrid.nx << " x "
                               << Ifmt(4) << mygrid.ny << " x "
                               << Ifmt(4) << mygrid.nz
                   << "  (" << Ifmt(8) << mygrid.npack_all_print(nb) << " waves "
                            << Ifmt(8) << mygrid.npack_print(nb) << " per task)" << std::endl;

      

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
   }

   if (fractional)
   {
      const int total_occ = ne[0] + ne[1];
      const int nbrillq = mygrid.nbrillq;
      const int ksize = myparallel.np_k();  // number of k-points
      const int kmy = myparallel.taskid_k();  // my k-point
      const int nlocal = nbrillq * total_occ;
                   
      std::vector<double> local_occ(nlocal);
      std::memcpy(local_occ.data(), occ1, nlocal * sizeof(double));
 
      std::vector<double> gathered_occ;
      if (myparallel.is_master_d(3))
      {       
         gathered_occ.resize(ksize * nlocal);
      }        
                    
      myparallel.Vector_GatherAll(3, nlocal, local_occ.data(), gathered_occ.data(), 0);
                   
      if (oprint)
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
         if (smeartype==4) std::cout << "Marzari-Vanderbilt" << std::endl;
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
         }
         if (ispin==2)
            std::cout <<  "    extra orbitals     : up=" << nextra[0] << " down= " << nextra[1] << std::endl;
         else
            std::cout <<  "    extra orbitals     = " << Ifmt(7) << nextra[0] << std::endl;
         if (control.fractional_frozen())
            std::cout <<  "    frozen orbitals" << std::endl;

         // Print occupations by k-point
         std::cout << "    initial occupations per Brillouin zone:" << std::endl;

         auto format_occ = [](double val) -> std::string {
            if (std::fabs(val - std::round(val)) < 1e-8)
               return std::to_string(static_cast<int>(std::round(val)));
            else {
               std::ostringstream oss;
               oss << std::fixed << std::setprecision(3) << val;
               return oss.str();
            }
         };

         for (auto nb=0; nb<mygrid.nbrillouin; ++nb)
         {
            std::cout << "        k-point " << nb+1 << ": [";
            for (int i=0; i<total_occ; ++i)
            {
                int idx = nb*total_occ + i;
                std::cout << format_occ(gathered_occ[idx]);
                if (i < total_occ * nbrillq - 1)
                    std::cout << " ";
            }
            std::cout << "]" << std::endl;
         }

         
      }
   }

   MPI_Barrier(comm_world0);

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
         band_inner_loop(control, &mygrid, &myion, &mykin, &mycoulomb, &myxc, &mypsp,
                         &mystrfac, &myewald, psi1, psi2, occ1, Hpsi, psi_r, dn, hml, lmbda, E,
                         &deltae, &deltac, &deltar);

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
         if ((!fractional) && (deltae > 0.0) && (icount > 1) && control.deltae_check())
         {
            done = 1;
            if (oprint) std::cout << "         *** Energy going up. iteration terminated\n";
         }
         else if ((std::fabs(deltae) < control.tolerances(0)) &&
                  (deltac < control.tolerances(1)) &&
                  (deltar < control.tolerances(2)))
         {
            done = 1;
            if (oprint) std::cout << "         *** tolerance ok.    iteration terminated\n";
         }
         else if (icount >= control.loop(1))
         {
            done = 1;
            if (oprint) std::cout << "          *** arrived at the Maximum iteration.    terminated ***\n";
         }


         // update occupations
         if ((fractional) && (!control.fractional_frozen()))
         {
            mygrid.w_diagonalize(hml, eig);
 
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
                   std::memcpy(eig_prev,eig,nbrillq*(ne[0]+ne[1])*sizeof(double));
                else             // Smooth eigenvalues in subsequent iterations
                   for (size_t i=0; i < nbrillq*(ne[0]+ne[1]); ++i)
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
                for (size_t i=0; i< nbrillq*(ne[0] + ne[1]); ++i)
                {
                    double delta_occ = occ2[i] - occ1[i];
                    rmsd_occupation += delta_occ * delta_occ;
                }
                rmsd_occupation = std::sqrt(rmsd_occupation / (nbrillq*(ne[0] + ne[1])) );
 
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
                   for (size_t i=0; i< nbrillq*(ne[0]+ne[1]); ++i)
                      occ1[i] = (1.0-alpha)*occ1[i] + alpha*occ2[i];
                }
                else
                {
                   std::memcpy(occ1,occ2,nbrillq*(ne[0]+ne[1])*sizeof(double));
                }
 
                 E[28] = smearcorrection;
                 //E[0] += E[28];
 
                 // Debugging output (optional)
                 if (oprint)
                    std::cout << " Iter.: " << icount
                              << ", RMSD: " << rmsd_occupation
                              << ", Alpha: " << alpha
                              << ", Smear Corr.: " << smearcorrection
                              << ", Delta Smear: " << smearcorrection - smearcorrection_old;
    
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
      mygrid.w_diagonalize(hml, eig);


   /* calculate real-space number of electrons, en */
   {        
      double omega = mylattice.omega();
      double scal1 = 1.0/((double)((mygrid.nx) * (mygrid.ny) * (mygrid.nz)));
      double dv = omega * scal1;

      en[0] = dv*mygrid.r_dsum(dn);
      en[1] = en[0];
      if (ispin>1)                      
         en[1] = dv*mygrid.r_dsum(dn + mygrid.nfft3d);
   }  


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

      std::cout << std::endl;
      std::cout << std::fixed << " number of electrons: spin up= "
                << Ffmt(11,5) << en[0] << "  down= "
                << Ffmt(11,5) << en[ispin-1] << " (real space)";
      std::cout << std::endl << std::endl;
      std::cout << " total     energy    : " << Efmt(19,10) << E[0] << " ("
                << Efmt(15,5) << E[0]/myion.nion << " /ion)" << std::endl;

      std::cout << " total orbital energy: "
                << Efmt(19,10) << E[1] << " ("
                << Efmt(15,5) << E[1]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      std::cout << " hartree energy      : "
                << Efmt(19,10) << E[2] << " ("
                << Efmt(15,5) << E[2]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      std::cout << " exc-corr energy     : "
                << Efmt(19,10) << E[3] << " ("
                << Efmt(15,5) << E[3]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      std::cout << " ion-ion energy      : "
                << Efmt(19,10) << E[4] << " ("
                << Efmt(15,5) << E[4]/myion.nion << " /ion)" << std::endl;
      if (std::fabs(E[28])>1.0e-9)
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

      viral = (E[9] + E[8] + E[7] + E[6]) / E[5];
      std::cout << " Virial Coefficient  : "
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
   }

   for (auto nb=0; nb<mygrid.nbrillouin; ++nb)
   {

      int nbq = mygrid.ktoindex(nb);
      int pk = mygrid.ktop(nb);
      int n = ne[0] + ne[1];
      double tmpeig[n];
      double tmpocc[n];
      std::memset(tmpeig,0,n*sizeof(double));
      std::memset(tmpocc,0,n*sizeof(double));

      if (pk==myparallel.taskid_k()) 
      {
         std::memcpy(tmpeig,eig+nbq*n, n*sizeof(double));
         std::memcpy(tmpocc,occ1+nbq*n,n*sizeof(double));
      }
      myparallel.Vector_SumAll(3,n,tmpeig);
      myparallel.Vector_SumAll(3,n,tmpocc);

      if (oprint)
      {
         std::cout << std::endl;
         std::cout << mybrillouin.print_zone_point(nb);

         std::cout << "\n orbital energies:\n";
         nn = ne[0] - ne[1];
         ev = 27.2116;

         for (int i=0; i<nn; ++i)
         {           
            //os << eig1stream(mymolecule.eig[i], mymolecule.eig[i] * ev);
            if (fractional)
            {
               if ((tmpocc[i] < 1.e-3) && (tmpocc[i]>1.0e-12))
                  std::cout << Efmt(18,7) << tmpeig[i] << " (" << Ffmt(8,3) << tmpeig[i] * ev << "eV) occ="
                     << Efmt(9,3) << tmpocc[i] << std::endl;
               else
                  std::cout << Efmt(18,7) << tmpeig[i] << " (" << Ffmt(8,3) << tmpeig[i] * ev << "eV) occ="
                     << Ffmt(5,3) << tmpocc[i] << std::endl;
            }
            else
               std::cout << Efmt(18,7) << tmpeig[i] << " (" << Ffmt(8,3) << tmpeig[i] * ev << "eV)" << std::endl;
         }
         for (int i=0; i<ne[1]; ++i)
         {
             if (fractional)
             {
                if ((tmpocc[i+nn] < 1.e-3) && (tmpocc[i+nn]>1.0e-12))
                   std::cout << Efmt(18,7) << eig[i+nn]    << " ("
                      << Ffmt(8,3)  << tmpeig[i+nn]*ev << "eV)  occ="
                      << Efmt(9,3)  << tmpocc[i+nn]   << " ";
                else
                   std::cout << Efmt(18,7) << eig[i+nn]    << " ("
                      << Ffmt(8,3)  << tmpeig[i+nn]*ev << "eV)  occ="
                      << Ffmt(5,3)  << tmpocc[i+nn]   << " "
                      << Efmt(18,7) << tmpeig[i+(ispin-1)*ne[0]]    << " ("
                      << Ffmt(8,3)  << tmpeig[i+(ispin-1)*ne[0]]*ev << "eV) occ="
                      << Ffmt(5,3)  << tmpocc[i+(ispin-1)*ne[0]] << std::endl;

                if ((tmpocc[i+(ispin-1)*ne[0]] < 1.e-3) && (tmpocc[i+(ispin-1)*ne[0]]>1.0e-12))
                   std::cout << Efmt(18,7) << tmpeig[i+(ispin-1)*ne[0]]    << " ("
                      << Ffmt(8,3)  << tmpeig[i+(ispin-1)*ne[0]]*ev << "eV) occ="
                      << Efmt(9,3)  << tmpocc[i + (ispin-1)*ne[0]] << std::endl;
                else
                   std::cout << Efmt(18,7) << tmpeig[i+(ispin-1)*ne[0]]    << " ("
                      << Ffmt(8,3)  << tmpeig[i+(ispin-1)*ne[0]]*ev << "eV) occ="
                      << Ffmt(5,3)  << tmpocc[i + (ispin-1)*ne[0]] << std::endl;
             }
             else
                std::cout << Efmt(18,7) << tmpeig[i+nn] << " ("
                   << Ffmt(8,3)  << tmpeig[i + nn] * ev << "eV) "
                   << Efmt(18,7) << tmpeig[i+(ispin-1)*ne[0]] << " ("
                   << Ffmt(8,3)  << tmpeig[i+(ispin-1)*ne[0]]*ev << "eV)" << std::endl;
         }
         std::cout << std::endl;
      }
   }
   if (oprint)
   {
      // write geometry and constraints analysis
      if (control.geometry_optimize())
      {
         std::cout << myion.print_bond_angle_torsions();
         std::cout << std::endl << myion.print_constraints(1);
      }
   }

   // write wavefunctions
   version = 5;
   nfft[0] = mygrid.nx;
   nfft[1] = mygrid.ny;
   nfft[2] = mygrid.nz;
   if (occ1)
      cpsi_write(&mygrid,&version,nfft,mylattice.unita_ptr(),&mygrid.ispin,mygrid.ne,&mygrid.nbrillouin,psi1,&smearoccupation,occ1,control.output_movecs_filename(),std::cout);
   else
      cpsi_write(&mygrid,&version,nfft,mylattice.unita_ptr(),&mygrid.ispin,mygrid.ne,&mygrid.nbrillouin,psi1,control.output_movecs_filename(),std::cout);


   /* deallocate memory */
   mygrid.g_deallocate(psi1);
   mygrid.g_deallocate(psi2);
   mygrid.g_deallocate(Hpsi);
   mygrid.h_deallocate(psi_r);
   mygrid.r_dealloc(dn);
   mygrid.w_deallocate(hml);
   mygrid.w_deallocate(lmbda);
   delete [] eig;
   delete [] eig_prev;
   if (fractional)
   {
      delete[] occ1;
      delete[] occ2;
   }
   mygrid.c3db::mygdevice.psi_dealloc();


  // write results to the json
   auto rtdbjson = json::parse(rtdbstring);
   rtdbjson["pspw"]["energy"] = E[0];
   rtdbjson["pspw"]["energies"] = E;

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


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
//#include	"control.hpp"
#include "Control2.hpp"
#include "Coulomb12.hpp"
#include "Electron.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "Kinetic.hpp"
#include "Lattice.hpp"
#include "Molecule.hpp"
#include "PGrid.hpp"
#include "Pneb.hpp"
#include "Pseudopotential.hpp"
#include "Strfac.hpp"
#include "exchange_correlation.hpp"
#include "inner_loop.hpp"
#include "psi.hpp"
#include "util_date.hpp"
//#include	"rtdb.hpp"
#include "mpi.h"

//#include "gdevice.hpp"
#include "nwpw_timing.hpp"
#include "psp_file_check.hpp"
#include "psp_library.hpp"

#include "cgsd_energy.hpp"
#include "cgsd_excited.hpp"

#include "json.hpp"
using json = nlohmann::json;

namespace pwdft {

/******************************************
 *                                        *
 *            pspw_minimizer              *
 *                                        *
 ******************************************/
int pspw_minimizer(MPI_Comm comm_world0, std::string &rtdbstring, std::ostream &coutput)
{
   // Parallel myparallel(argc,argv);
   Parallel myparallel(comm_world0);
   // RTDB myrtdb(&myparallel, "eric.db", "old");
 
   int version, nfft[3], ne[2], ispin, nextra[2];
   int i, ii, ia, nn, ngrid[3], matype, nelem, icount, done;
   char date[26];
   double sum1, sum2, ev, zv;
   double cpu1, cpu2, cpu3, cpu4;
   double E[70], deltae, deltac, deltar, viral, unita[9];

   bool fractional;
 
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
      coutput << "          *               PWDFT PSPW Calculation              *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *  [ (Grassmann/Stiefel manifold implementation) ]  *\n";
      coutput << "          *  [              C++ implementation             ]  *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *              version #7.00   02/27/21             *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *    This code was developed by Eric J. Bylaska,    *\n";
      coutput << "          *    Abhishek Bagusetty, David H. Bross, ...        *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *****************************************************\n";
      coutput << "          >>> job started at       " << util_date() << " <<<\n";
   }
  
   // control_read(myrtdb);
   // control_read(myparallel.np(),rtdbstring);
  
   // initialize processor grid structure
   myparallel.init2d(control.np_orbital(), control.pfft3_qsize());
  
 
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
   fractional = control.fractional();
   if (fractional)
   {
      nextra[0] = control.fractional_orbitals(0);
      if (control.ispin()==2)
         nextra[1] = control.fractional_orbitals(1);
      else
         nextra[1] = 0;
   }
   else
   {
      nextra[0] = 0;
      nextra[1] = 0;
   }

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
   version = control.version;
  
   // initialize parallel grid structure
   Pneb mygrid(&myparallel,&mylattice,control,control.ispin(),ne);
  
   // initialize gdevice memory
   mygrid.d3db::mygdevice.psi_alloc(mygrid.npack(1),mygrid.neq[0]+mygrid.neq[1],control.tile_factor());
  
   // setup structure factor
   Strfac mystrfac(&myion,&mygrid);
   mystrfac.phafac();
  
   // initialize operators
   Kinetic_Operator mykin(&mygrid);
   Coulomb12_Operator mycoulomb12(&mygrid,control);
   mycoulomb12.initialize_dielectric(&myion,&mystrfac);
  
   // initialize xc
   XC_Operator myxc(&mygrid,control);
  
   // initialize psp
   Pseudopotential mypsp(&myion,&mygrid,&mystrfac,control,coutput);
  
   // append Born information to rtdb for restarts
   if (mypsp.myapc->born_on)
      mypsp.myapc->myborn->writejsonstr(rtdbstring);
  
   // initialize electron operators
   Electron_Operators myelectron(&mygrid,&mykin,&mycoulomb12,&myxc,&mypsp);
  
   // setup ewald
   Ewald myewald(&myparallel,&myion,&mylattice,control,mypsp.zv);
   myewald.phafac();
  
   // initialize Molecule
   Molecule mymolecule(control.input_movecs_filename(),
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
      coutput << "\n";
      coutput << "     ===================  summary of input  =======================" << std::endl;
      coutput << "\n input psi filename: " << control.input_movecs_filename() << std::endl;
      coutput << "\n";
      coutput << " number of processors used: " << myparallel.np() << std::endl;
      coutput << " processor grid           : " << myparallel.np_i() << " x " << myparallel.np_j() << std::endl;
      if (mygrid.maptype==1) coutput << " parallel mapping         : 1d-slab" << std::endl;
      if (mygrid.maptype==2) coutput << " parallel mapping         : 2d-hilbert" << std::endl;
      if (mygrid.maptype==3) coutput << " parallel mapping         : 2d-hcurve" << std::endl;
      if (mygrid.isbalanced())
         coutput << " parallel mapping         : balanced" << std::endl;
      else
         coutput << " parallel mapping         : not balanced" << std::endl;

      if (mygrid.d3db::mygdevice.has_gpu())
      {
         //coutput << " parallel mapping         : has GPU" << std::endl;
         if (mygrid.d3db::mygdevice.type_gpu()==1) coutput << " parallel mapping         : CUDA GPUs" << std::endl;
         if (mygrid.d3db::mygdevice.type_gpu()==2) coutput << " parallel mapping         : SYCL GPUs" << std::endl;
         if (mygrid.d3db::mygdevice.type_gpu()==3) coutput << " parallel mapping         : HIP SYCL GPUs" << std::endl;
         if (mygrid.d3db::mygdevice.type_gpu()==4) coutput << " parallel mapping         : OpenCL GPUs" << std::endl;
         if (mygrid.staged_gpu_fft_pipeline) coutput << " parallel mapping         : staged GPU FFT" << std::endl;
         if (control.tile_factor() > 0)      coutput << " GPU tile factor          : " << control.tile_factor() << std::endl;
      }
      coutput << " fft container size       : " << control.fft_container_size() << std::endl;
     
      coutput << "\n options:\n";
      // coutput << "   geometry optimize    = ";
      // if (control.geometry_optimize() || flag==3)
      //    coutput << "yes\n";
      // else
      //    coutput << "no\n";
      coutput << "   boundary conditions  = ";
      if (control.version==3) coutput << "periodic\n";
      if (control.version==4) coutput << "aperiodic\n";
     
      coutput << "   electron spin        = ";
      if (ispin == 1)
         coutput << "restricted\n";
      else
         coutput << "unrestricted\n";
      coutput << myxc;

      if (fractional) std::cout << "   using fractional" << std::endl;
      
      coutput << mypsp.print_pspall();
     
      coutput << "\n total charge =" << Ffmt(8, 3) << control.total_charge() << std::endl;
     
      coutput << "\n atom composition:" << "\n";
      for (ia = 0; ia < myion.nkatm; ++ia)
         coutput << "   " << myion.atom(ia) << " : " << myion.natm[ia];
      coutput << "\n\n initial ion positions (au):"
              << "\n";
      for (ii = 0; ii < myion.nion; ++ii)
         coutput << Ifmt(4) << ii+1 << " " << myion.symbol(ii) << "\t( "
                 << Ffmt(10,5) << myion.rion1[3*ii] << " " 
                 << Ffmt(10,5) << myion.rion1[3*ii+1] << " " 
                 << Ffmt(10,5) << myion.rion1[3*ii+2] << " ) - atomic mass = " 
                 << Ffmt(6,3)  << myion.amu(ii) << std::endl;
      coutput << "   G.C.\t( " << Ffmt(10,5) << myion.gc(0) << " " 
                               << Ffmt(10,5) << myion.gc(1) << " " 
                               << Ffmt(10,5) << myion.gc(2) << " )" << std::endl;
      coutput << " C.O.M.\t( " << Ffmt(10,5) << myion.com(0) << " "
                               << Ffmt(10,5) << myion.com(1) << " " 
                               << Ffmt(10,5) << myion.com(2) << " )" << std::endl;

      coutput << std::endl;
      coutput << myion.print_symmetry_group();
      coutput << mypsp.myefield->shortprint_efield();
      coutput << mycoulomb12.shortprint_dielectric();
     
      coutput << "\n";
      if (fractional)
      {
         double oen[ispin];
         int n=0;
         for (auto ms=0; ms<ispin; ++ms)
         {
            oen[ms] = 0;
            for (auto i=0; i<ne[ms]; ++i)
            {
               oen[ms] += mymolecule.occ1[n];
               ++n;
            }
         }
         coutput << " number of electrons: spin up ="
                 << Ffmt(6,2)  << oen[0] << "  "
                 << Ffmt(27,2) << oen[ispin-1] << " (   fractional)" << std::endl;
      }
      else
         coutput << " number of electrons: spin up =" << Ifmt(6) << mygrid.ne[0]
                 << " (" << Ifmt(4) << mygrid.neq[0]
                 << " per task) down =" << Ifmt(6) << mygrid.ne[ispin-1] << " ("
                 << Ifmt(4) << mygrid.neq[ispin - 1] << " per task)" << std::endl;

      coutput << " number of orbitals:  spin up ="
              << Ifmt(6) << mygrid.ne[0] << " ("
              << Ifmt(4) << mygrid.neq[0] << " per task) down ="
              << Ifmt(6) << mygrid.ne[ispin-1] << " ("
              << Ifmt(4) << mygrid.neq[ispin-1] << " per task)" << std::endl;

      coutput << std::endl;
      coutput << " supercell:\n";
      coutput << "      volume = " << Ffmt(10,2) << mylattice.omega()
              << std::endl;
      coutput << "      lattice:    a1 = < " 
              << Ffmt(8,3) << mylattice.unita(0, 0) << " " 
              << Ffmt(8,3) << mylattice.unita(1, 0) << " " 
              << Ffmt(8,3) << mylattice.unita(2, 0) << " >\n";
      coutput << "                  a2 = < " 
              << Ffmt(8,3) << mylattice.unita(0, 1) << " " 
              << Ffmt(8,3) << mylattice.unita(1, 1) << " " 
              << Ffmt(8,3) << mylattice.unita(2, 1) << " >\n";
      coutput << "                  a3 = < " 
              << Ffmt(8,3) << mylattice.unita(0, 2) << " " 
              << Ffmt(8,3) << mylattice.unita(1, 2) << " " 
              << Ffmt(8,3) << mylattice.unita(2, 2) << " >\n";
      coutput << "      reciprocal: b1 = < " 
              << Ffmt(8,3) << mylattice.unitg(0, 0) << " " 
              << Ffmt(8,3) << mylattice.unitg(1, 0) << " " 
              << Ffmt(8,3) << mylattice.unitg(2, 0) << " >\n";
      coutput << "                  b2 = < " 
              << Ffmt(8,3) << mylattice.unitg(0, 1) << " " 
              << Ffmt(8,3) << mylattice.unitg(1, 1) << " " 
              << Ffmt(8,3) << mylattice.unitg(2, 1) << " >\n";
      coutput << "                  b3 = < " 
              << Ffmt(8,3) << mylattice.unitg(0, 2) << " " 
              << Ffmt(8,3) << mylattice.unitg(1, 2) << " " 
              << Ffmt(8,3) << mylattice.unitg(2, 2) << " >\n";
     
      {
        double aa1, bb1, cc1, alpha1, beta1, gamma1;
        mylattice.abc_abg(&aa1, &bb1, &cc1, &alpha1, &beta1, &gamma1);
        coutput << "      lattice:    a =    " << Ffmt(8,3) << aa1
                << " b =   " << Ffmt(8,3) << bb1 << " c =    " << Ffmt(8,3) << cc1 << std::endl;
        coutput << "                  alpha =" << Ffmt(8,3) << alpha1
                << " beta =" << Ffmt(8,3) << beta1 << " gamma =" << Ffmt(8,3) << gamma1 << std::endl;
      }
      coutput << "      density cutoff =" << Ffmt(7,3) << mylattice.ecut()
              << " fft =" << Ifmt(4) << mygrid.nx << " x " << Ifmt(4) << mygrid.ny << " x " << Ifmt(4) << mygrid.nz 
              << "  (" << Ifmt(8) << mygrid.npack_all(0) << " waves " << Ifmt(8) << mygrid.npack(0) << " per task)" 
              << std::endl;
      coutput << "      wavefnc cutoff =" << Ffmt(7,3) << mylattice.wcut()
              << " fft =" << Ifmt(4) << mygrid.nx << " x " << Ifmt(4) << mygrid.ny << " x " << Ifmt(4) << mygrid.nz 
              << "  (" << Ifmt(8) << mygrid.npack_all(1) << " waves " << Ifmt(8) << mygrid.npack(1) << " per task)" 
              << std::endl;
      coutput << std::endl;
      if (control.version==3)
      {
         coutput << " Ewald parameters:\n";
         coutput << "      energy cutoff = " << Ffmt(7,3) << myewald.ecut()
                 << " fft =" << Ifmt(4) << myewald.nx() << " x " << Ifmt(4) << myewald.ny() << " x " << Ifmt(4) << myewald.nz() 
                 << "  (" << Ifmt(8) << myewald.npack_all() << " waves " << Ifmt(8) << myewald.npack() << " per task)" 
                 << std::endl;
         coutput << "      Ewald summation: cut radius = " << Ffmt(7,3) << myewald.rcut() << " and " 
                 << Ifmt(3) << myewald.ncut() << std::endl;
         coutput << "                       Mandelung Wigner-Seitz =" 
                 << Ffmt(12, 8) << myewald.mandelung() << " (alpha =" << Ffmt(12,8) << myewald.rsalpha() 
                 << " rs =" << Ffmt(12,8) << myewald.rs() << ")" << std::endl;
      }
     
      if (flag > 0) 
      {
         coutput << std::endl;
         coutput << " technical parameters:\n";
         if (control.nolagrange()) coutput << "      disabling Lagrange multiplier " << std::endl;
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
            coutput << "     Kohn-Sham iterations = " << control.ks_maxit_orb() 
                                                      << " ( " << control.ks_maxit_orbs() << " outer)\n";
            if (control.scf_algorithm()==0) coutput << "     SCF algorithm        = simple mixing\n";
            if (control.scf_algorithm()==1) coutput << "     SCF algorithm        = Broyden mixing\n";
            if (control.scf_algorithm()==2) coutput << "     SCF algorithm        = Johnson-Pulay mixing"
                                                    << " (" << Ifmt(3) <<  control.diis_histories() << " histories)\n";
            if (control.scf_algorithm()==3) coutput << "     SCF algorithm        = Anderson mixing\n";
            if (control.scf_algorithm()==4) coutput << "     SCF algorithm        = Thomas-Fermi mixing\n";
            if (control.minimizer()==5) coutput << "     SCF mixing type      = potential\n";
            if (control.minimizer()==8) coutput << "     SCF mixing type      = density\n";
            if (control.scf_extra_rotate()) coutput << "     SCF extra rotate\n";
            if (control.scf_algorithm()==4) 
               coutput << "     SCF mixing parameters: alpha=" << control.scf_alpha() << " beta=" << control.scf_beta() << std::endl;
            else
               coutput << "     SCF mixing parameter: alpha= " << control.scf_alpha() << std::endl;
            if (control.kerker_g0()>0.0) coutput << "     Kerker damping       = " << control.kerker_g0() << std::endl;

            coutput << std::endl;

            if (control.fractional())
            {
               coutput <<  " fractional smearing parameters:" << std::endl;
               coutput <<  "    smearing algorithm = " << mymolecule.smeartype << std::endl;
               coutput <<  "    smearing parameter = ";
               if (mymolecule.smeartype==-1) coutput << "fixed_occupation" << std::endl;
               if (mymolecule.smeartype==0) coutput << "step function" << std::endl;
               if (mymolecule.smeartype==1) coutput << "Fermi-Dirac" << std::endl;
               if (mymolecule.smeartype==2) coutput << "Gaussian" << std::endl;
               if (mymolecule.smeartype==4) coutput << "Marazar-Vanderbilt" << std::endl;
               if (mymolecule.smeartype>=0)
               {
                  coutput <<  "    smearing parameter = " << Ffmt(9,3) << mymolecule.smearkT
                                                          << " (" << Ffmt(7,1) << control.fractional_temperature() << " K)" <<  std::endl;
                  coutput <<  "    mixing parameter   = " << Ffmt(7,1) << control.fractional_alpha() << std::endl;
                  if (ispin==2)
                     coutput <<  "    extra orbitals     : up=" << nextra[0] << " down= " << nextra[1] << std::endl;
                  else
                     coutput <<  "    extra orbitals     = " << Ifmt(7) << nextra[0] << std::endl;
               }
            }
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
      EV = cgsd_noit_energy(mymolecule, true, coutput);
   }
   else 
   {
      EV = cgsd_energy(control, mymolecule, true, coutput);
   }
   if (myparallel.is_master()) seconds(&cpu3);

   // write energy results to the json
   auto rtdbjson = json::parse(rtdbstring);
   rtdbjson["pspw"]["energy"] = EV;
   rtdbjson["pspw"]["energies"] = mymolecule.E;
   rtdbjson["pspw"]["eigenvalues"] = mymolecule.eig_vector();
  
   // calculate fion
   if (flag == 2) 
   {
      // double *fion = new double[3*myion.nion];
      double fion[3*myion.nion];
      cgsd_energy_gradient(mymolecule, fion);
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
      rtdbjson["pspw"]["fion"] = std::vector<double>(fion, &fion[3 * myion.nion]);
      for (ii=0; ii<(3*myion.nion); ++ii) fion[ii] *= -1.0;
      rtdbjson["pspw"]["gradient"] = std::vector<double>(fion, &fion[3 * myion.nion]);
     
      // delete [] fion;
   }

   // calculate excited state orbitals 
   cgsd_excited(control, mymolecule, true, coutput);

   // calculate oep orbitals 

   // calculate the spin contamination
  
   // APC analysis
   if (mypsp.myapc->apc_on) 
   {
      if (!(mypsp.myapc->v_apc_on))
         mypsp.myapc->gen_APC(mymolecule.dng1, false);
     
      // set qion
      double qion[myion.nion];
      for (auto ii = 0; ii < myion.nion; ++ii)
         qion[ii] = -mypsp.myapc->Qtot_APC(ii) + mypsp.zv[myion.katm[ii]];
      rtdbjson["nwpw"]["apc"]["q"] = std::vector<double>(qion, &qion[myion.nion]);
     
      if (lprint) 
      {
         coutput << mypsp.myapc->print_APC(mypsp.zv);
      }
   }
  
   // dipole analysis
   if (mypsp.mydipole->dipole_on) 
   {
      double *mdipole = mypsp.mydipole->mdipole;
      double mu = std::sqrt(mdipole[0]*mdipole[0] + mdipole[1]*mdipole[1] + mdipole[2]*mdipole[2]);
      rtdbjson["nwpw"]["dipole"] = std::vector<double>(mdipole,&mdipole[3]);
      rtdbjson["nwpw"]["dipole_magnitude"] = mu;
   }
  
   // write psi
   if (flag > 0)
      mymolecule.writepsi(control.output_movecs_filename(), coutput);
   MPI_Barrier(comm_world0);
  
   // set rtdbjson initialize_wavefunction option to false
   if (rtdbjson["nwpw"]["initialize_wavefunction"].is_boolean())
      rtdbjson["nwpw"]["initialize_wavefunction"] = false;
  
   // write rtdbjson
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
      double av = t2 / ((double)myelectron.counter);
      // cout.setf(ios::scientific);
      coutput << std::scientific;
      coutput << std::endl;
      coutput << " ------------------" << std::endl;
      coutput << " cputime in seconds" << std::endl;
      coutput << " prologue    : " << Efmt(9,3) << t1 << std::endl;
      coutput << " main loop   : " << Efmt(9,3) << t2 << std::endl;
      coutput << " epilogue    : " << Efmt(9,3) << t3 << std::endl;
      coutput << " total       : " << Efmt(9,3) << t4 << std::endl;
      coutput << " cputime/step: " << Efmt(9,3) << av << " ( "
              << myelectron.counter << " evaluations, "
              << util_linesearch_counter() << " linesearches)" << std::endl;;
      coutput << std::endl;
     
      nwpw_timing_print_final(myelectron.counter, coutput);
     
      coutput << std::endl;
      coutput << " >>> job completed at     " << util_date() << " <<<" << std::endl;;
   }
  
   // deallocate memory
   mygrid.d3db::mygdevice.psi_dealloc();
  
   MPI_Barrier(comm_world0);

   return 0;
}

} // namespace pwdft

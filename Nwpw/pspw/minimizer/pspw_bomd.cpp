
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
#include "nwpw_Nose_Hoover.hpp"
#include "nwpw_aimd_running_data.hpp"
#include "psi.hpp"
#include "util_date.hpp"
//#include	"rtdb.hpp"
#include "mpi.h"

//#include "gdevice.hpp"
#include "nwpw_timing.hpp"
#include "psp_file_check.hpp"
#include "psp_library.hpp"

#include "nwpw_lmbfgs.hpp"

#include "cgsd_energy.hpp"

#include "json.hpp"
using json = nlohmann::json;

namespace pwdft {

/******************************************
 *                                        *
 *              pspw_bomd                 *
 *                                        *
 ******************************************/
int pspw_bomd(MPI_Comm comm_world0,std::string &rtdbstring,std::ostream &coutput) 
{
   // Parallel myparallel(argc,argv);
   Parallel myparallel(comm_world0);
   // RTDB myrtdb(&myparallel, "eric.db", "old");
 
   bool verlet, SA;
   int version, nfft[3], ne[2], ispin;
   int i, ii, ia, nn, ngrid[3], matype, nelem, icount;
   char date[26];
   double sum1, sum2, ev, zv;
   double cpu1, cpu2, cpu3, cpu4, cpustep;
   double E[70], deltae, deltac, deltar, viral, unita[9];
   double sa_alpha[2], sa_decay[2], Te_init, Tr_init, Te_new, Tr_new;
   double kb = 3.16679e-6;
 
   Control2 control(myparallel.np(), rtdbstring);
   int flag = control.task();
 
   bool hprint = (myparallel.is_master() && control.print_level("high"));
   bool oprint = (myparallel.is_master() && control.print_level("medium"));
   bool lprint = (myparallel.is_master() && control.print_level("low"));
 
   /* reset Parallel base_stdio_print = lprint */
   myparallel.base_stdio_print = lprint;
 
   // double *psi1,*psi2,*Hpsi,*psi_r;
   // double *dn;
   // double *hml,*lmbda,*eig;
 
   for (ii=0; ii<70; ++ii)
      E[ii] = 0.0;
 
   if (myparallel.is_master()) seconds(&cpu1);
   if (oprint) 
   {
      std::ios_base::sync_with_stdio();
      coutput << "          *****************************************************\n";
      coutput << "          *                                                   *\n";
      coutput << "          *   PWDFT PSPW Born-Oppenheimer molecular dynamics  *\n";
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
 
   /* initialize processor grid structure */
   myparallel.init2d(control.np_orbital(), control.pfft3_qsize());
 
   /* initialize lattice */
   Lattice mylattice(control);
 
   /* read in ion structure */
   // Ion myion(myrtdb);
   Ion myion(rtdbstring, control);
 
   /* Check for and generate psp files                       */
   /* - this routine also sets the valence charges in myion, */
   /*   and total_ion_charge and ne in control               */
   psp_file_check(&myparallel, &myion, control, coutput);
 
   /* fetch ispin and ne psi information from control */
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
   version = 3;
 
   /* initialize parallel grid structure */
   Pneb mygrid(&myparallel, &mylattice, control, control.ispin(),
               control.ne_ptr());
 
   /* initialize gdevice memory */
   mygrid.d3db::mygdevice.psi_alloc(mygrid.npack(1), mygrid.neq[0] + mygrid.neq[1],
                     control.tile_factor());
 
   /* setup structure factor */
   Strfac mystrfac(&myion, &mygrid);
   mystrfac.phafac();
 
   /* initialize operators */
   Kinetic_Operator mykin(&mygrid);
   Coulomb12_Operator mycoulomb12(&mygrid, control);
   mycoulomb12.initialize_dielectric(&myion,&mystrfac);
 
   /* initialize xc */
   XC_Operator myxc(&mygrid, control);
 
   /* initialize psp */
   Pseudopotential mypsp(&myion, &mygrid, &mystrfac, control, coutput);
 
   /* initialize electron operators */
   Electron_Operators myelectron(&mygrid, &mykin, &mycoulomb12, &myxc, &mypsp);
 
   // setup ewald
   Ewald myewald(&myparallel, &myion, &mylattice, control, mypsp.zv);
   myewald.phafac();
 
   /* initialize thermostats */
   double w = 0.01;
   nwpw_Nose_Hoover mynose(myion, (mygrid.ne[0] + mygrid.ne[1]), w, control);
 
   /* initialize simulated annealing */
   SA = false;
   Tr_init = 0.0;
   sa_alpha[1] = 1.0;
   if (control.SA()) {
     sa_decay[1] = control.SA_decay(1);
     if (mynose.on()) {
       SA = true;
       Tr_init = mynose.Tr;
     } else {
       double dt = control.bo_time_step();
       SA = false;
       sa_alpha[1] = exp(-(dt / sa_decay[1]));
     }
   }
 
   // initialize Molecule
   Molecule mymolecule(control.input_movecs_filename(),
                       control.input_movecs_initialize(), &mygrid, &myion,
                       &mystrfac, &myewald, &myelectron, &mypsp, control, coutput);
 
   MPI_Barrier(comm_world0);
 
   // driver parameters
   int maxit = control.driver_maxiter();
   double tol_Gmax = control.driver_gmax();
   double tol_Grms = control.driver_grms();
   double tol_Xrms = control.driver_xrms();
   double tol_Xmax = control.driver_xmax();
   double trust = control.driver_trust();
   int lmbfgs_size = control.driver_lmbfgs_size();
 
   //                 |**************************|
   // *****************   summary of input data  **********************
   //                 |**************************|
 
   if (oprint) {
     coutput << std::endl;
     coutput << "     ===================  summary of input  =======================" << std::endl;
     coutput << "\n input psi filename: " << control.input_movecs_filename() << std::endl;
     coutput << std::endl;
     coutput << " number of processors used: " << myparallel.np() << std::endl;
     coutput << " processor grid           : " << myparallel.np_i() << " x" << myparallel.np_j() << std::endl;
     if (mygrid.maptype == 1) coutput << " parallel mapping         : 1d-slab" << std::endl;
     if (mygrid.maptype == 2) coutput << " parallel mapping         : 2d-hilbert" << std::endl;
     if (mygrid.maptype == 3) coutput << " parallel mapping         : 2d-hcurve" << std::endl;
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
        if (control.tile_factor() > 1)      coutput << " GPU tile factor          : " << control.tile_factor() << std::endl;
     }
     coutput << " fft container size       : " << control.fft_container_size() << std::endl;

     coutput << "\n options:\n";
     coutput << "   boundary conditions  = ";
     if (control.version == 3)
       coutput << "periodic\n";
     if (control.version == 4)
       coutput << "aperiodic\n";
 
     coutput << "   electron spin        = ";
     if (ispin == 1)
       coutput << "restricted\n";
     else
       coutput << "unrestricted\n";
     coutput << myxc;
 
     // coutput << "\n elements involved in the cluster:\n";
     // for (ia=0; ia<myion.nkatm; ++ia)
     //{
     //    printf("    %2d : %4s   core charge: %4.1lf  lmax=%1d\n",
     //            ia+1,myion.atom(ia),mypsp.zv[ia],mypsp.lmax[ia]);
     //    printf("           comment : %s\n",mypsp.comment[ia]);
     //    printf("           pseudopotential type            :
     //    %3d\n",mypsp.psp_type[ia]); printf("           highest angular
     //    component       : %3d\n",mypsp.lmax[ia]); printf("           local
     //    potential used            : %3d\n",mypsp.locp[ia]); printf(" number of
     //    non-local projections : %3d\n",mypsp.nprj[ia]); if
     //    (mypsp.semicore[ia])
     //       printf("           semicore corrections included   : %6.3lf
     //       (radius) %6.3lf (charge)\n",mypsp.rcore[ia],mypsp.ncore(ia));
     //    printf("           cutoff = ");
     //    for (ii=0; ii<=mypsp.lmax[ia]; ++ii)
     //       printf("%8.3lf",mypsp.rc[ia][ii]);
     //    printf("\n");
     // }
     coutput << mypsp.print_pspall();
 
     coutput << "\n total charge =" << Ffmt(8, 3) << control.total_charge()
             << std::endl;
 
     coutput << "\n atom composition:"
             << "\n";
     for (ia = 0; ia < myion.nkatm; ++ia)
       coutput << "   " << myion.atom(ia) << " : " << myion.natm[ia];
     coutput << "\n\n initial ion positions (au):"
             << "\n";
     for (ii = 0; ii < myion.nion; ++ii)
       coutput << Ifmt(4) << ii + 1 << " " << myion.symbol(ii) << "\t( "
               << Ffmt(10, 5) << myion.rion1[3 * ii] << " " << Ffmt(10, 5)
               << myion.rion1[3 * ii + 1] << " " << Ffmt(10, 5)
               << myion.rion1[3 * ii + 2] << " ) - atomic mass = " << Ffmt(6, 3)
               << myion.amu(ii) << std::endl;
     coutput << "   G.C.\t( " << Ffmt(10, 5) << myion.gc(0) << " " << Ffmt(10, 5)
             << myion.gc(1) << " " << Ffmt(10, 5) << myion.gc(2) << " )"
             << std::endl;
     coutput << " C.O.M.\t( " << Ffmt(10, 5) << myion.com(0) << " "
             << Ffmt(10, 5) << myion.com(1) << " " << Ffmt(10, 5) << myion.com(2)
             << " )" << std::endl;

     coutput << std::endl;
     //coutput << myion.print_symmetry_group();
     coutput << myion.print_symmetry_group(rtdbstring);
     coutput << mypsp.myefield->shortprint_efield();
     coutput << mycoulomb12.shortprint_dielectric();
 
     coutput << "\n";
     coutput << " number of electrons: spin up =" << Ifmt(6) << mygrid.ne[0]
             << " (" << Ifmt(4) << mygrid.neq[0]
             << " per task) down =" << Ifmt(6) << mygrid.ne[ispin - 1] << " ("
             << Ifmt(4) << mygrid.neq[ispin - 1] << " per task)" << std::endl;
 
     coutput << "\n";
     coutput << " supercell:\n";
     coutput << "      volume = " << Ffmt(10, 2) << mylattice.omega()
             << std::endl;
     coutput << "      lattice:    a1 = < " << Ffmt(8, 3)
             << mylattice.unita(0, 0) << " " << Ffmt(8, 3)
             << mylattice.unita(1, 0) << " " << Ffmt(8, 3)
             << mylattice.unita(2, 0) << " >\n";
     coutput << "                  a2 = < " << Ffmt(8, 3)
             << mylattice.unita(0, 1) << " " << Ffmt(8, 3)
             << mylattice.unita(1, 1) << " " << Ffmt(8, 3)
             << mylattice.unita(2, 1) << " >\n";
     coutput << "                  a3 = < " << Ffmt(8, 3)
             << mylattice.unita(0, 2) << " " << Ffmt(8, 3)
             << mylattice.unita(1, 2) << " " << Ffmt(8, 3)
             << mylattice.unita(2, 2) << " >\n";
     coutput << "      reciprocal: b1 = < " << Ffmt(8, 3)
             << mylattice.unitg(0, 0) << " " << Ffmt(8, 3)
             << mylattice.unitg(1, 0) << " " << Ffmt(8, 3)
             << mylattice.unitg(2, 0) << " >\n";
     coutput << "                  b2 = < " << Ffmt(8, 3)
             << mylattice.unitg(0, 1) << " " << Ffmt(8, 3)
             << mylattice.unitg(1, 1) << " " << Ffmt(8, 3)
             << mylattice.unitg(2, 1) << " >\n";
     coutput << "                  b3 = < " << Ffmt(8, 3)
             << mylattice.unitg(0, 2) << " " << Ffmt(8, 3)
             << mylattice.unitg(1, 2) << " " << Ffmt(8, 3)
             << mylattice.unitg(2, 2) << " >\n";
 
     {
       double aa1, bb1, cc1, alpha1, beta1, gamma1;
       mylattice.abc_abg(&aa1, &bb1, &cc1, &alpha1, &beta1, &gamma1);
       coutput << "      lattice:    a =    " << Ffmt(8, 3) << aa1
               << " b =   " << Ffmt(8, 3) << bb1 << " c =    " << Ffmt(8, 3)
               << cc1 << std::endl;
       coutput << "                  alpha =" << Ffmt(8, 3) << alpha1
               << " beta =" << Ffmt(8, 3) << beta1 << " gamma =" << Ffmt(8, 3)
               << gamma1 << std::endl;
     }
     coutput << "      density cutoff =" << Ffmt(7, 3) << mylattice.ecut()
             << " fft =" << Ifmt(4) << mygrid.nx << " x " << Ifmt(4) << mygrid.ny
             << " x " << Ifmt(4) << mygrid.nz << "  (" << Ifmt(8)
             << mygrid.npack_all(0) << " waves " << Ifmt(8) << mygrid.npack(0)
             << " per task)" << std::endl;
     coutput << "      wavefnc cutoff =" << Ffmt(7, 3) << mylattice.wcut()
             << " fft =" << Ifmt(4) << mygrid.nx << " x " << Ifmt(4) << mygrid.ny
             << " x " << Ifmt(4) << mygrid.nz << "  (" << Ifmt(8)
             << mygrid.npack_all(1) << " waves " << Ifmt(8) << mygrid.npack(1)
             << " per task)" << std::endl;
     coutput << "\n";
     coutput << " Ewald parameters:\n";
     coutput << "      energy cutoff = " << Ffmt(7, 3) << myewald.ecut()
             << " fft =" << Ifmt(4) << myewald.nx() << " x " << Ifmt(4)
             << myewald.ny() << " x " << Ifmt(4) << myewald.nz() << "  ("
             << Ifmt(8) << myewald.npack_all() << " waves " << Ifmt(8)
             << myewald.npack() << " per task)" << std::endl;
     coutput << "      Ewald summation: cut radius = " << Ffmt(7, 3)
             << myewald.rcut() << " and " << Ifmt(3) << myewald.ncut()
             << std::endl;
     coutput << "                       Mandelung Wigner-Seitz =" << Ffmt(12, 8)
             << myewald.mandelung() << " (alpha =" << Ffmt(12, 8)
             << myewald.rsalpha() << " rs =" << Ffmt(12, 8) << myewald.rs()
             << ")" << std::endl;
 
     if (flag > 0) 
     {
        coutput << std::endl;
        coutput << " technical parameters:\n";
        if (control.nolagrange()) coutput << "      disabling Lagrange multiplier" << std::endl;
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
        if (control.minimizer() == 1)
           coutput << "      minimizer = Grassmann conjugate gradient\n";
        if (control.minimizer() == 2)
           coutput << "      minimizer = Grassmann lmbfgs\n";
        if (control.minimizer() == 4)
           coutput << "      minimizer = Stiefel conjugate gradient\n";
        if (control.minimizer() == 5)
           coutput << "      minimizer = scf (potential)\n";
        if (control.minimizer() == 7)
           coutput << "      minimizer = Stiefel lmbfgs\n";
        if (control.minimizer() == 8)
           coutput << "      minimizer = scf (density)\n";
        if ((control.minimizer() == 5) || (control.minimizer() == 8)) 
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
           if (control.kerker_g0() > 0.0) coutput << "     Kerker damping       = " << control.kerker_g0() << std::endl;
           coutput << std::endl;
           if (control.fractional())
           { 
              coutput << " fractional smearing parameter:" << std::endl;
              if (control.fractional_smeartype()==-1) coutput << "     smearing algorithm = fixed occupation\n";
              if (control.fractional_smeartype()==0)  coutput << "     smearing algorithm = step function\n";
              if (control.fractional_smeartype()==1)  coutput << "     smearing algorithm = Fermi-Dirac\n";
              if (control.fractional_smeartype()==2)  coutput << "     smearing algorithm = Gaussian\n";
              if (control.fractional_smeartype()==4)  coutput << "     smearing algorithm = Marzari-Vanderbilt\n";
              if (control.fractional_smeartype()>=0)  coutput << "     smearing parameter = "
                                                              << control.fractional_kT() 
                                                              << " (" <<  Ffmt(8,1) <<  control.fractional_temperature() << " K)\n"
                                                              << "     mixing parameter   =  "
                                                              << control.fractional_alpha() << std::endl;
              coutput << std::endl;
           }
        }
     } else {
       coutput << std::endl;
       coutput << " technical parameters:\n";
       coutput << "      optimization of psi and densities turned off"
               << std::endl;
     }
     coutput << std::endl << std::endl << std::endl;
 
     coutput << " -----------------------------------------------------------------------------------\n";
     coutput << " ----------------------- Born-Oppenheimer molecular dynamics -----------------------\n";
     coutput << " -----------------------------------------------------------------------------------\n\n";
     coutput << "\n";
     coutput << " molecular dynamics parameters:\n";
     if (myion.fix_translation)
       coutput << "      translation constrained\n";
     if (myion.fix_rotation)
       coutput << "      rotation constrained\n";
 
     coutput << "      time step =" << Ffmt(11, 2) << control.bo_time_step()
             << " iterations = " << Ifmt(10)
             << control.bo_steps(0) * control.bo_steps(1) << " (" << Ifmt(5)
             << control.bo_steps(0) << " inner " << Ifmt(5)
             << control.bo_steps(1) << " outer)\n";
     if (control.bo_algorithm() == 0)
       coutput << "      integration algorithm = position Verlet\n";
     if (control.bo_algorithm() == 1)
       coutput << "      integration algorithm = velocity Verlet\n";
     if (control.bo_algorithm() == 2)
       coutput << "      integration algorithm = leap frog\n";
 
     coutput << std::endl;
     coutput << " scaling parameters:   " << std::endl;
     coutput << "      cooling/heating rate   =" << Efmt(12, 5)
             << control.ion_scaling() << " (ion)" << std::endl;
     coutput << "      initial kinetic energy =" << Efmt(12, 5) << myion.eki0
             << " (ion)" << std::endl;
     coutput << "                           " << Efmt(15, 5) << myion.ekg
             << " (C.O.M.)" << std::endl;
     coutput << "      after scaling          =" << Efmt(12, 5) << myion.eki1
             << " (ion)" << std::endl;
     coutput << "      increased energy       =" << Efmt(12, 5)
             << myion.eki1 - myion.eki0 << " (ion)" << std::endl;
     coutput << std::endl;
 
     if (mynose.on())
       coutput << mynose.inputprint();
     else
       coutput << " constant energy simulation" << std::endl;
 
     if (SA)
       coutput << "      SA decay rate = " << Ffmt(10, 3) << sa_decay[1]
               << "(ion)\n";
   }
   if (myparallel.is_master())
     seconds(&cpu2);
 
   //*                |***************************|
   //******************     call GeoVibminimizer  **********************
   //*                |***************************|
 
   /* intialize the linesearch */
   util_linesearch_init();
 
   /*  calculate energy and gradient */
   double g, gg, Gmax, Grms, Xrms, Xmax;
   double Eold = 0.0;
   double EV = 0.0;
   double eki = 0.0;
   double *Emol = mymolecule.E;
 
   int nfsize = 3 * myion.nion;
   int one = 1;
   double mrone = -1.0;
 
   // allocate temporary memory from stack
   double fion[3 * myion.nion];
 
   bool done = false;
 
   /*  calculate energy */
   if (oprint) {
 
     coutput << "\n\n";
     coutput << " -----------------------------------------------------------------------------------\n";
     coutput << " -----------------------------    Initial Geometry     -----------------------------\n";
     coutput << " -----------------------------------------------------------------------------------\n\n";
     coutput << " ---------------------------------\n";
     coutput << "         Initial Geometry         \n";
     coutput << " ---------------------------------\n";
     coutput << mymolecule.myion->print_bond_angle_torsions();
     coutput << "\n\n\n";
     coutput << "         ================ Born-Oppenheimer AIMD iteration ================\n";
     coutput << "     >>> iteration started at " << util_date() << " <<<\n";
     coutput << "     iter.          KE+Energy             Energy        KE_Ion   Temperature\n";
     coutput << "     -----------------------------------------------------------------------\n";
   }
 
   nwpw_aimd_running_data mymotion_data(control, &myparallel, &mylattice, &myion, E,
                                        mymolecule.hml, mymolecule.psi1,
                                        mymolecule.rho1);
 
   if (control.bo_steps(1) > 0) {
     bool nose = mynose.on();
     bool verlet = (control.bo_algorithm() == 0);
     bool vverlet = (control.bo_algorithm() == 1);
     bool leapfrog = (control.bo_algorithm() == 2);
     double dt = control.bo_time_step();
     int it_out = control.bo_steps(1);
     int it_in = control.bo_steps(0);
     double r = 1.0;
     double ssr = 1.0;
     double eki1 = 0.0;
     icount = 0;
 
     double fion1[3 * myion.nion];
 
     EV = cgsd_energy(control, mymolecule, false, coutput);
     cgsd_energy_gradient(mymolecule, fion);
 
     if (nose)
       r = (1.0 - 0.5 * dt * mynose.dXr());
     myion.Newton_step(fion, sa_alpha[1] * r);
 
     eki = myion.eki1;
     if (nose)
       mynose.Newton_step(0.01, eki);
 
     E[0] = EV + eki;
     E[1] = EV;
     E[2] = 0.0;
     E[3] = eki;
     for (auto i = 0; i < 56; ++i)
       E[i + 4] = Emol[i];
 
     done = false;
     // outer loop iteration
     while (!done) {
       ++icount;
 
       // inner loop iteration
       for (auto it = 0; it < it_in; ++it) {
         if (vverlet) {
           memcpy(fion1, fion, 3 * myion.nion * sizeof(double));
           myion.shift21();
         } else {
           myion.shift();
           if (nose)
             mynose.shift();
         }
 
         //  calculate the energy and gradient
         EV = cgsd_energy(control,mymolecule,false,coutput);
         cgsd_energy_gradient(mymolecule, fion);
 
         if (nose) 
         {
            ssr = mynose.ssr();
            myion.Nose_step(ssr, fion);
            eki = myion.eki1;
            mynose.Verlet_step(0.01, eki);
         } 
         else if (vverlet) 
         {
            myion.vVerlet_step(fion, fion1);
            myion.vshift();
            myion.Newton_step(fion, sa_alpha[1] * r);
            eki = myion.eki1;
         } 
         else 
         {
            myion.Verlet_step(fion, sa_alpha[1]);
            eki = myion.eki1;
         }
 
         E[0] = EV + eki;
         E[1] = EV;
         E[2] = 0.0;
         E[3] = eki;
         for (auto i = 0; i < 56; ++i)
           E[i + 4] = Emol[i];
 
       } // end inner loop
 
       // Update AIMD Running data
       mymotion_data.update_iteration(icount);
       // eave = mymotion_data.eave; evar = mymotion_data.evar;
       // have = mymotion_data.have; hvar = mymotion_data.hvar;
       // qave = mymotion_data.qave; qvar = mymotion_data.qvar;
 
       // Write out loop energies
       if (oprint) 
       {
          if (SA) 
          {
             if (mynose.on())
                coutput << Ifmt(10) << (icount*it_in) 
                        << Efmt(19,10) << EV+eki+mynose.r_energy() 
                        << Efmt(19,10) << EV
                        << Efmt(14,5)  << eki 
                        << Ffmt(9,1)   << Tr_new 
                        << std::endl;
             else
                coutput << Ifmt(10) << (icount*it_in) 
                        << Efmt(19,10) << EV+eki
                        << Efmt(19,10) << EV 
                        << Efmt(14,5)  << eki 
                        << Ffmt(9,1)   << Tr_new 
                        << std::endl;
          } 
          else 
          {
             if (mynose.on())
                coutput << Ifmt(10) << (icount*it_in) 
                        << Efmt(19,10) << EV+eki+mynose.r_energy() 
                        << Efmt(19,10) << EV
                        << Efmt(14,5)  << eki 
                        << Ffmt(14,2)  << myion.Temperature()
                        << std::endl;
             else
                coutput << Ifmt(10) << (icount*it_in) 
                        << Efmt(19,10) << EV+eki
                        << Efmt(19,10) << EV 
                        << Efmt(14,5)  << eki 
                        << Ffmt(14,2)  << myion.Temperature() << std::endl;
          }
       }
 
       // Check for running out of time
       if (control.out_of_time()) {
         done = true;
         if (oprint)
           coutput << "         *** out of time. iteration terminated."
                   << std::endl;
       }
       // Check for Completion
       else if (icount >= it_out) {
         done = true;
         if (oprint)
           coutput
               << "         *** arrived at the Maximum iteration.   terminated."
               << std::endl;
       }
     } // end outer loop
 
   } // end bomd iterations
 
   if (oprint) {
     coutput << "\n\n";
     coutput << " ---------------------------------\n";
     coutput << "  Final Geometry \n";
     coutput << " ---------------------------------\n";
     coutput << mymolecule.myion->print_bond_angle_torsions();
     coutput << "\n\n";
   }
   if (myparallel.is_master())
     seconds(&cpu3);
 
   //*******************************************************************
 
   // write energy results to the json
   auto rtdbjson = json::parse(rtdbstring);
   rtdbjson["pspw"]["energy"] = EV;
   rtdbjson["pspw"]["energies"] = mymolecule.E;
   rtdbjson["pspw"]["eigenvalues"] = mymolecule.eig_vector();
 
   // APC analysis
   if (mypsp.myapc->apc_on) {
     if (!(mypsp.myapc->v_apc_on))
       mypsp.myapc->gen_APC(mymolecule.dng1, false);
 
     // set qion
     double qion[myion.nion];
     for (auto ii = 0; ii < myion.nion; ++ii)
       qion[ii] = -mypsp.myapc->Qtot_APC(ii) + mypsp.zv[myion.katm[ii]];
     rtdbjson["nwpw"]["apc"]["q"] = std::vector<double>(qion, &qion[myion.nion]);
 
     if (oprint) {
       coutput << mypsp.myapc->print_APC(mypsp.zv);
     }
   }
 
   // set rtdbjson initialize_wavefunction option to false
   if (rtdbjson["nwpw"]["initialize_wavefunction"].is_boolean())
     rtdbjson["nwpw"]["initialize_wavefunction"] = false;
 
   MPI_Barrier(comm_world0);
 
   /* write psi */
   if (flag > 0)
     mymolecule.writepsi(control.output_movecs_filename(), coutput);
 
   /* write rtdbjson */
   rtdbstring = rtdbjson.dump();
   myion.writejsonstr(rtdbstring);
 
   //                 |**************************|
   // *****************   report consumed time   **********************
   //                 |**************************|
   if (myparallel.is_master())
     seconds(&cpu4);
   if (oprint) {
     double t1 = cpu2 - cpu1;
     double t2 = cpu3 - cpu2;
     double t3 = cpu4 - cpu3;
     double t4 = cpu4 - cpu1;
     double av = t2 / ((double)myelectron.counter);
     // coutput.setf(ios::scientific);
     coutput << std::scientific;
     coutput << "\n";
     coutput << " -----------------"
             << "\n";
     coutput << " cputime in seconds"
             << "\n";
     coutput << " prologue    : " << t1 << "\n";
     coutput << " main loop   : " << t2 << "\n";
     coutput << " epilogue    : " << t3 << "\n";
     coutput << " total       : " << t4 << "\n";
     coutput << " cputime/step: " << av << " ( " << myelectron.counter
             << " evaluations, " << util_linesearch_counter()
             << " linesearches)\n";
     coutput << "\n";
 
     nwpw_timing_print_final(myelectron.counter, coutput);
 
     coutput << "\n";
     coutput << " >>> job completed at     " << util_date() << " <<<\n";
   }
 
   /* deallocate memory */
   mygrid.d3db::mygdevice.psi_dealloc();
   // delete [] fion;
   // delete [] sion;
 
   MPI_Barrier(comm_world0);
   std::cout << "END pspw_bomd" << std::endl;
 
   return 0;
}

} // namespace pwdft

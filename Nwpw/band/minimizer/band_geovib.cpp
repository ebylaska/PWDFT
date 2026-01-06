
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

#include "nwpw_lmbfgs.hpp"

#include "json.hpp"
using json = nlohmann::json;

namespace pwdft {

/******************************************
 *                                        *
 *            band_geovib                 *
 *                                        *
 ******************************************/
int band_geovib(MPI_Comm comm_world0, std::string &rtdbstring, std::ostream &coutput)
{
   // Parallel myparallel(argc,argv);
   Parallel myparallel(comm_world0);
   // RTDB myrtdb(&myparallel, "eric.db", "old");
 
   int version, nfft[3], ne[2], nextra[2], ispin;
   int i, ii, ia, nn, ngrid[3], matype, nelem, icount;
   char date[26];
   double sum1, sum2, ev, zv;
   double cpu1, cpu2, cpu3, cpu4, cpustep;
   double E[80], deltae, deltac, deltar, viral, unita[9];

   bool fractional;
 
   // double *psi1,*psi2,*Hpsi,*psi_r;
   // double *dn;
   // double *hml,*lmbda,*eig;
 
   for (ii=0; ii<80; ++ii)
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
      coutput << "          *            PWDFT BAND GeoVib Calculation          *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *  [      (bundled Grassmann/Stiefel manifold)   ]  *\n";
      coutput << "          *  [              C++ implementation             ]  *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *              version #7.00   07/19/24             *\n";
      coutput << "          *                                                   *\n";
      coutput << "          *    This code was developed by Eric J. Bylaska,    *\n";
      coutput << "          *    Raymundo Hernandez Esparza, Duo Song,          *\n";
      coutput << "          *    David H. Bross, ...                            *\n";
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

   // fetch ispin and ne psi information from control
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


   /* read in Brillouin zone */
   Brillouin mybrillouin(rtdbstring,&mylattice,control);

   // initialize parallel grid structure
   Cneb mygrid(&myparallel, &mylattice, control, control.ispin(),ne,&mybrillouin);

  
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

      if (control.fractional()) coutput << "   using fractional" << std::endl;

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
      //coutput << myion.print_symmetry_group();
      coutput << myion.print_symmetry_group(rtdbstring);

      if (control.geometry_optimize())
         coutput << std::endl << myion.print_constraints(0);

      coutput << std::endl;

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
               oen[ms] += weight*mysolid.occ1[n];
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
         coutput << " number of electrons: spin up ="
                 << Ffmt(6,2)  << oen[0] << "  "
                 << Ffmt(27,2) << oen[ispin-1] << " (   fractional)" << std::endl;
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
      coutput << " Brillouin zone:" << std::endl;
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
         if ((control.minimizer()==3) || 
             (control.minimizer()==5) || 
             (control.minimizer()==6) || 
             (control.minimizer()==8))
         {
            coutput << std::endl;
            coutput << " Kohn-Sham scf parameters:\n";
            if (control.ks_algorithm()==-1) coutput << "      Kohn-Sham algorithm  = block conjugate gradient\n";
            if (control.ks_algorithm()==0) coutput << "      Kohn-Sham algorithm  = conjugate gradient\n";
            if (control.ks_algorithm()==1) coutput << "      Kohn-Sham algorithm  = rmm-diis\n";
            if (control.ks_algorithm()==2) coutput << "      Kohn-Sham algorithm  = Grassmann conjugate gradient\n";
            if (control.ks_algorithm()==3) coutput << "      Kohn-Sham algorithm  = Grassmann Stiefel\n";

            if ((control.ks_algorithm()==2) || (control.ks_algorithm()==3))
               coutput << "      Kohn-Sham iterations = " << "( " << control.loop(0) << " inner)" << std::endl;
            else
               coutput << "      Kohn-Sham iterations = " << control.ks_maxit_orb()
                                                          << " ( " << control.ks_maxit_orbs() << " outer)\n";

            if (control.scf_algorithm()==0) coutput << "      scf algorithm        = simple mixing\n";
            if (control.scf_algorithm()==1) coutput << "      scf algorithm        = Broyden mixing\n";
            if (control.scf_algorithm()==2) coutput << "      scf algorithm        = Johnson-Pulay mixing"
                                                    << " (" << Ifmt(3) <<  control.diis_histories() << " histories)\n";
            if (control.scf_algorithm()==3) coutput << "      scf algorithm        = Anderson mixing\n";
            if (control.scf_algorithm()==4) coutput << "      scf algorithm        = Thomas-Fermi mixing\n";

            if (control.minimizer()==5) coutput << "      scf mixing type      = potential\n";
            if (control.minimizer()==8) coutput << "      scf mixing type      = density\n";
            if (control.scf_extra_rotate()) coutput << "     scf extra rotate\n";
            if (control.scf_algorithm()==4)
               coutput << "      scf mixing parameters: alpha=" << control.scf_alpha() << " beta=" << control.scf_beta() << std::endl;
            else
               coutput << "      scf mixing parameter: alpha= " << control.scf_alpha() << std::endl;
            if (control.kerker_g0()>0.0) coutput << "      Kerker damping       = " << control.kerker_g0() << std::endl;
         }
      }
      else
      {
         coutput << std::endl;
         coutput << " technical parameters:\n";
         coutput << "      optimization of psi and densities turned off" << std::endl;
      }
   }

   if (control.fractional())
   {

      const int total_occ = ne[0] + ne[1];
      const int nbrillq = mygrid.nbrillq;
      const int ksize = myparallel.np_k();  // number of k-points
      const int kmy = myparallel.taskid_k();  // my k-point
      const int nlocal = nbrillq * total_occ;

      std::vector<double> local_occ(nlocal);
      std::memcpy(local_occ.data(), mysolid.occ1, nlocal * sizeof(double));

      std::vector<double> gathered_occ;
      if (myparallel.is_master_d(3)) 
       {
         gathered_occ.resize(ksize * nlocal);
      }

      myparallel.Vector_GatherAll(3, nlocal, local_occ.data(), gathered_occ.data(), 0);

      if (oprint) 
      {
         coutput <<  std::endl;
         coutput <<  " fractional smearing parameters:" << std::endl;
         
         coutput <<  "      smearing algorithm = " << mysolid.smeartype << std::endl;
         coutput <<  "      smearing parameter = ";
         if (mysolid.smeartype==-1) coutput << "fixed occupation" << std::endl;
         if (mysolid.smeartype==0) coutput << "step function" << std::endl;
         if (mysolid.smeartype==1) coutput << "Fermi-Dirac" << std::endl;
         if (mysolid.smeartype==2) coutput << "Gaussian" << std::endl;
         if (mysolid.smeartype==3) coutput << "Hermite" << std::endl;
         if (mysolid.smeartype==4) coutput << "Marzari-Vanderbilt" << std::endl;
         if (mysolid.smeartype==5) coutput << "Methfessel-Paxton" << std::endl;
         if (mysolid.smeartype==6) coutput << "Cold smearing" << std::endl;
         if (mysolid.smeartype==7) coutput << "Lorentzian" << std::endl;
         if (mysolid.smeartype==8) coutput << "step" << std::endl;
        
         if (mysolid.smeartype>=0)
         {
            coutput <<  "      smearing parameter = " << Ffmt(9,3) << mysolid.smearkT
                                                    << " (" << Ffmt(7,1) << control.fractional_temperature() << " K)" <<  std::endl;
            coutput <<  "      mixing parameter   = " << Ffmt(7,1) << control.fractional_alpha() << std::endl;
            coutput <<  "      mixing parameter(alpha)   = " << Ffmt(7,2) << control.fractional_alpha() << std::endl;
            coutput <<  "      mixing parameter(alpha_min)   = " << Ffmt(7,2) << control.fractional_alpha_min() << std::endl;
            coutput <<  "      mixing parameter(alpha_max)   = " << Ffmt(7,2) << control.fractional_alpha_max() << std::endl;
            coutput <<  "      mixing parameter(beta)   = " << Ffmt(7,2) << control.fractional_beta() << std::endl;
            coutput <<  "      mixing parameter(gamma)   = " << Ffmt(7,2) << control.fractional_gamma() << std::endl;
            coutput <<  "      rmsd occupation tolerance   = " << Efmt(12,3) << control.fractional_rmsd_tolerance() << std::endl;
          }
          {
            if (ispin==2)
                coutput <<  "      extra orbitals     : up=" << nextra[0] << " down= " << nextra[1] << std::endl;
            else
                coutput <<  "      extra orbitals     = " << Ifmt(7) << nextra[0] << std::endl;
       
            if (control.fractional_frozen()) coutput <<  "      frozen oribtals" << std::endl;


            // Print occupations by k-point
            coutput << "      initial occupations per Brillouin zone:" << std::endl;

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
               coutput << "          k-point " << nb+1 << ": [";
               for (int i=0; i<total_occ; ++i)
               {
                   int idx = nb*total_occ + i;
                   coutput << format_occ(gathered_occ[idx]);
                   if (i < total_occ * nbrillq - 1)
                       coutput << " ";
               }
               coutput << "]" << std::endl;
            }
         }
      }
   }

   if (oprint)
   {
     coutput << std::endl << std::endl << std::endl;

      coutput << " --------------------------------------------------------------"
                 "---------------------\n";
      coutput << " ----------------------------- Geometry Optimization "
                 "-------------------------------\n";
      coutput << " --------------------------------------------------------------"
                 "---------------------\n\n";
      coutput << " Optimization parameters:  " << Ifmt(5) << maxit
              << " (maxiter) " << Efmt(12, 3) << tol_Gmax << " (gmax) "
              << Efmt(12, 3) << tol_Grms << " (grms) " << Efmt(12, 3) << tol_Xrms
              << " (xrms) " << Efmt(12, 3) << tol_Xmax << " (xmax)\n";

      coutput << "    maximum number of steps       (maxiter) = " << Ifmt(4)
              << maxit << std::endl;
      coutput << "    maximum gradient threshold       (Gmax) = " << Efmt(12, 6)
              << tol_Gmax << std::endl;
      coutput << "    rms gradient threshold           (Grms) = " << Efmt(12, 6)
              << tol_Grms << std::endl;
      coutput << "    rms cartesian step threshold     (Xrms) = " << Efmt(12, 6)
              << tol_Xrms << std::endl;
      coutput << "    maximum cartesian step threshold (Xmax) = " << Efmt(12, 6)
              << tol_Xmax << std::endl;
      coutput << std::endl;
      coutput << "    fixed trust radius              (Trust) = " << Efmt(12, 6)
              << trust << std::endl;
      coutput << "    number lmbfgs histories   (lmbfgs_size) = " << Ifmt(4)
              << lmbfgs_size << std::endl;
   }

   MPI_Barrier(comm_world0);
   if (myparallel.is_master()) seconds(&cpu2);
  
   //*                |***************************|
   //******************     call CG minimizer     **********************
   //*                |***************************|

   /*  calculate energy and gradient */
   double g, gg, Gmax, Grms, Xrms, Xmax;
   double Eold = 0.0;
   double EV = 0.0;

   int nfsize = 3 * myion.nion;
   int one = 1;
   double mrone = -1.0;

   // allocate temporary memory from stack
   double fion[3 * myion.nion];
   double sion[3 * myion.nion];


   bool done = false;
   int it = 0;

   /*  calculate energy */
   if (oprint) {
     coutput << "\n\n";
     coutput << " --------------------------------------------------------------"
                "---------------------\n";
     coutput << " -----------------------------    Initial Geometry     "
                "-----------------------------\n";
     coutput << " --------------------------------------------------------------"
                "---------------------\n\n";
     coutput << " ---------------------------------\n";
     coutput << "         Initial Geometry         \n";
     coutput << " ---------------------------------\n";
     coutput << mysolid.myion->print_bond_angle_torsions();
     coutput << "\n\n\n";
     coutput << " ---------------------------------\n";
     coutput << "     Calculate Initial Energy     \n";
     coutput << " ---------------------------------\n\n";
   }
   EV = band_cgsd_energy(control, mysolid, true, coutput);

   /*  calculate the gradient */
   if (oprint) {
     coutput << "\n";
     coutput << " ---------------------------------\n";
     coutput << "    Calculate Initial Gradient    \n";
     coutput << " ---------------------------------\n\n";
   }
   band_cgsd_energy_gradient(mysolid, fion);
   if (oprint)
   {
      coutput << " ion forces (au):" << std::endl;
      for (auto ii = 0; ii < mysolid.myion->nion; ++ii)
         coutput << Ifmt(5) << ii+1 << " "
                 << Lfmt(2) << mysolid.myion->symbol(ii) << "  ( "
                 << Ffmt(10,5) << fion[3*ii]   << " "
                 << Ffmt(10,5) << fion[3*ii+1] << " "
                 << Ffmt(10,5) << fion[3*ii+2] << " )" << std::endl;
      coutput << "  C.O.M.  ( "
              << Ffmt(10,5) << mysolid.myion->com_fion(fion,0) << " "
              << Ffmt(10,5) << mysolid.myion->com_fion(fion,1) << " "
              << Ffmt(10,5) << mysolid.myion->com_fion(fion,2) << " )" << std::endl;;
      coutput << "|F|/nion  = " << std::setprecision(5) << std::fixed
              << std::setw(10) << mysolid.myion->rms_fion(fion) << std::endl
              << "max|Fatom|= " << std::setprecision(5) << std::fixed
              << std::setw(10) << mysolid.myion->max_fion(fion) << "  ("
              << std::setprecision(3) << std::fixed << std::setw(8)
              << mysolid.myion->max_fion(fion)*(27.2116/0.529177)
              << " eV/Angstrom)" << std::endl
              << std::endl;
   }
   DSCAL_PWDFT(nfsize, mrone, fion, one);


   /* initialize lmbfgs */
   nwpw_lmbfgs geom_lmbfgs(3*myion.nion, lmbfgs_size, myion.rion1, fion);

   /* update coords in mysolid */
   mysolid.myion->fixed_step(-trust, fion);
   mysolid.myion->shift();

   Xmax = mysolid.myion->xmax();
   Xrms = mysolid.myion->xrms();
   Gmax = 0.0;
   Grms = 0.0;
   for (ii = 0; ii < (myion.nion); ++ii) {
      gg = fion[3*ii]*fion[3*ii] + fion[3*ii+1]*fion[3*ii+1] + fion[3*ii+2]*fion[3*ii+2];
      Grms += gg;
      g = sqrt(gg);
      if (g > Gmax)
         Gmax = g;
   }
   Grms = sqrt(Grms) / ((double)myion.nion);
   if ((Gmax <= tol_Gmax) && (Grms <= tol_Grms) && (Xrms <= tol_Xrms) && (Xmax <= tol_Xmax))
      done = true;


   if (oprint)
   {
      if (done) {
        coutput << "      ----------------------\n"
                << "      Optimization converged\n"
                << "      ----------------------\n";
      }
      if (it == 0) {
         coutput << std::endl << std::endl;
         coutput << "@ Step             Energy     Delta E     Gmax     Grms     Xrms     Xmax   Walltime\n";
         coutput << "@ ---- ------------------ ----------- -------- -------- -------- -------- ----------\n";
      } else {
         coutput << std::endl << std::endl;
         coutput << "  Step             Energy     Delta E     Gmax     Grms     Xrms     Xmax   Walltime\n";
         coutput << "  ---- ------------------ ----------- -------- -------- -------- -------- ----------\n";
      }
      seconds(&cpustep);
      coutput << "@ " << Ifmt(4) << it << " "
              << Ffmt(18,9) << EV << " "
              << Efmt(11,3) << EV - Eold << " "
              << Ffmt(8,5) << Gmax << " "
              << Ffmt(8,5) << Grms << " "
              << Ffmt(8,5) << Xrms << " "
              << Ffmt(8,5) << Xmax << " "
              << Ffmt(10,1) << cpustep - cpu1 << " "
              << Ifmt(9) << myelectron.counter << std::endl;
      // printf("@ %4d %18.9lf %11.3le %8.5lf %8.5lf %8.5lf %8.5lf %10.1lf
      // %9d\n",it,EV,(EV-Eold),Gmax,Grms,Xrms,Xmax,cpustep-cpu1,myelectron.counter);

      if ((Gmax<=tol_Gmax) && (Grms<=tol_Grms) && (Xrms<=tol_Xrms) && (Xmax<=tol_Xmax))
      {
         coutput << "                                            ok       ok       ok       ok\n";
      }
      else if ((Gmax <= tol_Gmax) || (Grms <= tol_Grms) || (Xrms <= tol_Xrms) || (Xmax <= tol_Xmax))
      {
         std::string oktag("");
         coutput << "                                     ";
         oktag = "  ";
         if (Gmax <= tol_Gmax)
           oktag = "ok";
         coutput << "       " << oktag;
         oktag = "  ";
         if (Grms <= tol_Grms)
           oktag = "ok";
         coutput << "       " << oktag;
         oktag = "  ";
         if (Xrms <= tol_Xrms)
           oktag = "ok";
         coutput << "       " << oktag;
         oktag = "  ";
         if (Xmax <= tol_Xmax)
           oktag = "ok";
         coutput << "       " << oktag << "\n";
      }
   }


   ++it;

   while ((!done) && (it <= maxit))
   {
      if (oprint)
      {
         coutput << "\n\n";
         coutput << " -----------------------------------------------------------------------------------\n";
         coutput << " ----------------------------- Optimization Step "
                 << std::setw(5) << it << " -----------------------------\n";
         coutput << " -----------------------------------------------------------------------------------\n\n";
      }

      /* print out the current geometry */
      if (oprint)
      {
         coutput << " ---------------------------------\n";
         coutput << "  Geometry for Step " << it << "\n";
         coutput << " ---------------------------------\n";
         coutput << mysolid.myion->print_bond_angle_torsions();
         coutput << std::endl << myion.print_constraints(1);
      }

      /*  calculate energy */
      if (oprint)
      {
         coutput << "\n\n\n";
         coutput << " ---------------------------------\n";
         coutput << "  Calculate Energy for Step " << it << "\n";
         coutput << " ---------------------------------\n\n";
      }
      Eold = EV;
      EV = band_cgsd_energy(control, mysolid, true, coutput);

      /* calculate the gradient */
      if (oprint)
      {
         coutput << "\n";
         coutput << " ---------------------------------\n";
         coutput << "  Calculate Gradient for Step " << it << "\n";
         coutput << " ---------------------------------\n\n";
      }
      band_cgsd_energy_gradient(mysolid, fion);

      if (oprint)
      {
         coutput << " ion forces (au):"
                 << "\n";
         for (auto ii=0; ii<mysolid.myion->nion; ++ii)
            coutput << Ifmt(5) << ii+1 << " "
                    << Lfmt(2) << mysolid.myion->symbol(ii) << "  ( "
                    << Ffmt(10,5) << fion[3*ii] << " "
                    << Ffmt(10,5) << fion[3*ii+1] << " "
                    << Ffmt(10,5) << fion[3*ii+2] << " )" << std::endl;
         coutput << "  C.O.M.  ( "
                 << Ffmt(10,5) << mysolid.myion->com_fion(fion,0) << " "
                 << Ffmt(10,5) << mysolid.myion->com_fion(fion,1) << " "
                 << Ffmt(10,5) << mysolid.myion->com_fion(fion,2) << " )\n";
         coutput << "|F|/nion  = " << Ffmt(10,5) << mysolid.myion->rms_fion(fion) << std::endl
                 << "max|Fatom|= " << Ffmt(10,5) << mysolid.myion->max_fion(fion) << "  ("
                 << Ffmt(8,3) << mysolid.myion->max_fion(fion)*(27.2116/0.529177)
                 << " eV/Angstrom)" << std::endl << std::endl;

         //coutput << std::endl << myion.print_constraints(1);
      }
      DSCAL_PWDFT(nfsize, mrone, fion, one);

      /* lmbfgs gradient */
      geom_lmbfgs.lmbfgs(myion.rion1, fion, sion);

      /* update coords in mysolid */
      mysolid.myion->fixed_step(-trust, sion);
      mysolid.myion->shift();

      Xmax = mysolid.myion->xmax();
      Xrms = mysolid.myion->xrms();
      Gmax = 0.0;
      Grms = 0.0;
      for (ii=0; ii<(myion.nion); ++ii) {
         gg = fion[3*ii]*fion[3*ii] + fion[3*ii+1]*fion[3*ii+1] + fion[3*ii+2]*fion[3*ii+2];
         Grms += gg;
         g = sqrt(gg);
         if (g > Gmax)
            Gmax = g;
      }
      Grms = sqrt(Grms) / ((double)myion.nion);
      if ((Gmax<=tol_Gmax) && (Grms<=tol_Grms) && (Xrms<=tol_Xrms) && (Xmax<=tol_Xmax))
         done = true;


      /* print out the current energy */
      if (oprint)
      {
         if (done)
         {
            coutput << "      ----------------------\n"
                    << "      Optimization converged\n"
                    << "      ----------------------\n";
         }

         if (it == 0)
         {
            coutput << std::endl << std::endl;
            coutput << "@ Step             Energy     Delta E     Gmax     Grms     Xrms     Xmax   Walltime\n";
            coutput << "@ ---- ------------------ ----------- -------- -------- -------- -------- ----------\n";
         }
         else
         {
            coutput << std::endl << std::endl;
            coutput << "  Step             Energy     Delta E     Gmax     Grms     Xrms     Xmax   Walltime\n";
            coutput << "  ---- ------------------ ----------- -------- -------- -------- -------- ----------\n";
         }
         seconds(&cpustep);
         coutput << "@ " << Ifmt(4) << it << " "
                 << Ffmt(18,9) << EV << " "
                 << Efmt(11,3) << EV - Eold << " "
                 << Ffmt(8,5) << Gmax << " "
                 << Ffmt(8,5) << Grms << " "
                 << Ffmt(8,5) << Xrms << " "
                 << Ffmt(8,5) << Xmax << " "
                 << Ffmt(10,1) << cpustep - cpu1
                 << " " << Ifmt(9) << myelectron.counter << std::endl;
          // printf("@ %4d %18.9lf %11.3le %8.5lf %8.5lf %8.5lf %8.5lf %10.1lf
          // %9d\n",it,EV,(EV-Eold),Gmax,Grms,Xrms,Xmax,cpustep-cpu1,myelectron.counter);

         if ((Gmax <= tol_Gmax) && (Grms <= tol_Grms) && (Xrms <= tol_Xrms) && (Xmax <= tol_Xmax))
         {
            coutput << "                                            ok       ok       ok       ok\n";
         }
         else if ((Gmax <= tol_Gmax) || (Grms <= tol_Grms) || (Xrms <= tol_Xrms) || (Xmax <= tol_Xmax))
         {
            std::string oktag = "";
            coutput << "                                     ";
            oktag = "  ";
            if (Gmax <= tol_Gmax)
              oktag = "ok";
            coutput << "       " << oktag;
            oktag = "  ";
            if (Grms <= tol_Grms)
              oktag = "ok";
            coutput << "       " << oktag;
            oktag = "  ";
            if (Xrms <= tol_Xrms)
              oktag = "ok";
            coutput << "       " << oktag;
            oktag = "  ";
            if (Xmax <= tol_Xmax)
              oktag = "ok";
            coutput << "       " << oktag << "\n";
         }
      }

      ++it;
   }

   if (oprint)
   {
      coutput << "\n\n";
      coutput << " ---------------------------------\n";
      coutput << "  Final Geometry \n";
      coutput << " ---------------------------------\n";
      coutput << mysolid.myion->print_bond_angle_torsions();
      coutput << std::endl << myion.print_constraints(1);
      coutput << "\n\n";
   }
   if (myparallel.is_master()) seconds(&cpu3);


   // **************************************************************************
   //
   //             Vibrational analysis (finite-difference Hessian)
   //
   // --------- NOTE -----------------------------------------------------------
   // BAND vibrational modes reuse the same force code as geometry optimization.
   // No k-dependent perturbation is needed.  The only requirement is force 
   // consistency.
   //
   // **************************************************************************
   /*if(control.geovib()) {
       compute_band_hessian(..);
       compute_vib_modes(..);
   }
   */



   //*******************************************************************

  
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
      mysolid.writecpsi(control.output_movecs_filename(), coutput);

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

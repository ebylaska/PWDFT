
#include	<iostream>
#include	<cstdio>
#include	<cmath>
#include	<cstdlib>
#include	<string>
#include	<vector>
//
#include 	"iofmt.hpp"
#include        "util_linesearch.hpp"
#include	"Parallel.hpp"
//#include	"control.hpp"
#include	"Control2.hpp"
#include	"Lattice.hpp"
#include	"util_date.hpp"
#include	"PGrid.hpp"
#include	"Pneb.hpp"
#include	"Ion.hpp"
#include	"Ewald.hpp"
#include	"Strfac.hpp"
#include        "nwpw_Nose_Hoover.hpp"
#include	"Kinetic.hpp"
#include	"Coulomb.hpp"
#include	"exchange_correlation.hpp"
#include	"Pseudopotential.hpp"
#include	"Electron.hpp"
#include	"Molecule.hpp"
#include	"inner_loop.hpp"
#include	"psi.hpp"
#include        "nwpw_aimd_running_data.hpp"
//#include	"rtdb.hpp"
#include	"mpi.h"


#include	"psp_library.hpp"
#include	"psp_file_check.hpp"
#include	"nwpw_timing.hpp"
#include        "gdevice.hpp"

#include        "nwpw_lmbfgs.hpp"

#include	"cgsd_energy.hpp"

#include "json.hpp"
using json = nlohmann::json;




namespace pwdft {

/******************************************
 *                                        *
 *              pspw_bomd                 *
 *                                        *
 ******************************************/
int pspw_bomd(MPI_Comm comm_world0, std::string& rtdbstring)
{
   //Parallel myparallel(argc,argv);
   Parallel myparallel(comm_world0);
   //RTDB myrtdb(&myparallel, "eric.db", "old");

   bool verlet,SA;
   int version,nfft[3],ne[2],ispin;
   int i,ii,ia,nn,ngrid[3],matype,nelem,icount;
   char date[26];
   double sum1,sum2,ev,zv;
   double cpu1,cpu2,cpu3,cpu4,cpustep;
   double E[60],deltae,deltac,deltar,viral,unita[9];
   double sa_alpha[2],sa_decay[2],Te_init,Tr_init,Te_new,Tr_new;
   double kb = 3.16679e-6;

   Control2 control(myparallel.np(),rtdbstring);
   int flag =  control.task();

   bool hprint = (myparallel.is_master() && control.print_level("high"));
   bool oprint = (myparallel.is_master() && control.print_level("medium"));
   bool lprint = (myparallel.is_master() && control.print_level("low"));

   /* reset Parallel base_stdio_print = lprint */
   myparallel.base_stdio_print = lprint;


   //double *psi1,*psi2,*Hpsi,*psi_r;
   //double *dn;
   //double *hml,*lmbda,*eig;

   for (ii=0; ii<60; ++ii) E[ii] = 0.0;

   if (myparallel.is_master()) seconds(&cpu1);
   if (oprint)
   {
      std::ios_base::sync_with_stdio();
      std::cout << "          *****************************************************\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *   PWDFT PSPW Born-Oppenheimer molecular dynamics  *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *  [ (Grassmann/Stiefel manifold implementation) ]  *\n";
      std::cout << "          *  [              C++ implementation             ]  *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *              version #7.00   02/27/21             *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *    This code was developed by Eric J. Bylaska,    *\n";
      std::cout << "          *    Abhishek Bagusetty, David H. Bross, ...        *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *****************************************************\n";
      std::cout << "          >>> job started at       " << util_date() << " <<<\n";
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
   psp_file_check(&myparallel,&myion,control);


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
   Pneb mygrid(&myparallel,&mylattice,control,control.ispin(),control.ne_ptr());

   /* initialize gdevice memory */
   gdevice_psi_alloc(mygrid.npack(1),mygrid.neq[0]+mygrid.neq[1]);


   /* setup structure factor */
   Strfac mystrfac(&myion, &mygrid);
   mystrfac.phafac();

   /* initialize operators */
   Kinetic_Operator mykin(&mygrid);
   Coulomb_Operator mycoulomb(&mygrid);

   /* initialize xc */
   XC_Operator      myxc(&mygrid,control);

   /* initialize psp */
   Pseudopotential mypsp(&myion,&mygrid,&mystrfac,control);

   /* initialize electron operators */
   Electron_Operators myelectron(&mygrid,&mykin, &mycoulomb, &myxc, &mypsp);

   // setup ewald
   Ewald myewald(&myparallel,&myion,&mylattice,control,mypsp.zv);
   myewald.phafac();

   /* initialize thermostats */
   double w = 0.01;
   nwpw_Nose_Hoover mynose(myion,(mygrid.ne[0]+mygrid.ne[1]),w,control);

   /* initialize simulated annealing */
   SA       = false;
   Tr_init  = 0.0;
   sa_alpha[1]     = 1.0;
   if (control.SA())
   {
      sa_decay[1] = control.SA_decay(1);
      if (mynose.on())
      {
         SA      = true;
         Tr_init = mynose.Tr;
      }
      else
      {
         double dt   = control.bo_time_step();
         SA          = false;
         sa_alpha[1] = exp(-(dt/sa_decay[1]));
      }
   }


   // initialize Molecule
   Molecule mymolecule(control.input_movecs_filename(),control.input_movecs_initialize(),
                       &mygrid,&myion,&mystrfac,&myewald,&myelectron,&mypsp);

   MPI_Barrier(comm_world0);


   // driver parameters
   int maxit       = control.driver_maxiter();
   double tol_Gmax = control.driver_gmax();
   double tol_Grms = control.driver_grms();
   double tol_Xrms = control.driver_xrms();
   double tol_Xmax = control.driver_xmax();
   double trust    = control.driver_trust();
   int lmbfgs_size = control.driver_lmbfgs_size();


//                 |**************************|
// *****************   summary of input data  **********************
//                 |**************************|

   if (oprint)
   {
      std::cout << "\n";
      std::cout << "          ==============  summary of input  ==================\n";
      std::cout << "\n input psi filename: " << control.input_movecs_filename() << "\n";
      std::cout << "\n";
      std::cout << " number of processors used: " << myparallel.np() << "\n";
      std::cout << " processor grid           : " << myparallel.np_i() << " x" << myparallel.np_j() << "\n";
      if (mygrid.maptype==1) std::cout << " parallel mapping         : 1d-slab"    << "\n";
      if (mygrid.maptype==2) std::cout << " parallel mapping         : 2d-hilbert" << "\n";
      if (mygrid.maptype==3) std::cout << " parallel mapping         : 2d-hcurve" << "\n";
      if (mygrid.isbalanced()) 
         std::cout << " parallel mapping         : balanced" << "\n";
      else
         std::cout << " parallel mapping         : not balanced" << "\n";

      std::cout << "\n options:\n";
      std::cout << "   boundary conditions  = ";
      if (control.version==3) std::cout << "periodic\n";
      if (control.version==4) std::cout << "aperiodic\n";

      std::cout << "   electron spin        = ";
      if (ispin==1)
         std::cout << "restricted\n";
      else
         std::cout << "unrestricted\n";
      std::cout << myxc;
  
      //std::cout << "\n elements involved in the cluster:\n";
      //for (ia=0; ia<myion.nkatm; ++ia)
      //{
      //   printf("    %2d : %4s   core charge: %4.1lf  lmax=%1d\n",
      //           ia+1,myion.atom(ia),mypsp.zv[ia],mypsp.lmax[ia]);
      //   printf("           comment : %s\n",mypsp.comment[ia]);
      //   printf("           pseudopotential type            : %3d\n",mypsp.psp_type[ia]);
      //   printf("           highest angular component       : %3d\n",mypsp.lmax[ia]);
      //   printf("           local potential used            : %3d\n",mypsp.locp[ia]);
      //   printf("           number of non-local projections : %3d\n",mypsp.nprj[ia]);
      //   if (mypsp.semicore[ia])
      //      printf("           semicore corrections included   : %6.3lf (radius) %6.3lf (charge)\n",mypsp.rcore[ia],mypsp.ncore(ia));
      //   printf("           cutoff = ");
      //   for (ii=0; ii<=mypsp.lmax[ia]; ++ii)
      //      printf("%8.3lf",mypsp.rc[ia][ii]);
      //   printf("\n");
      //}
      std::cout << mypsp.print_pspall();

      std::cout << "\n total charge =" << Ffmt(8,3) << control.total_charge() << std::endl;

      std::cout << "\n atom composition:" << "\n";
      for (ia=0; ia<myion.nkatm; ++ia)
         std::cout << "   " << myion.atom(ia) << " : " << myion.natm[ia];
      std::cout << "\n\n initial ion positions (au):" << "\n";
      for (ii=0; ii<myion.nion; ++ii)
         std::cout << Ifmt(4) << ii+1 << " " << myion.symbol(ii) 
                   << "\t( " << Ffmt(10,5) << myion.rion1[3*ii] << " " << Ffmt(10,5) << myion.rion1[3*ii+1] << " " << Ffmt(10,5) << myion.rion1[3*ii+2]
                   << " ) - atomic mass = " << Ffmt(6,3) << myion.amu(ii) << std::endl;
      std::cout << "   G.C.\t( " << Ffmt(10,5) << myion.gc(0) << " " << Ffmt(10,5) << myion.gc(1) << " " << Ffmt(10,5) << myion.gc(2) << " )" << std::endl;
      std::cout << " C.O.M.\t( " << Ffmt(10,5) << myion.com(0) << " " << Ffmt(10,5) << myion.com(1) << " " << Ffmt(10,5) << myion.com(2) << " )" << std::endl;
      std::cout << "\n";
      std::cout <<" number of electrons: spin up =" << Ifmt(6) << mygrid.ne[0] << " (" << Ifmt(4) << mygrid.neq[0]
                << " per task) down =" << Ifmt(6) << mygrid.ne[ispin-1] << " (" << Ifmt(4) << mygrid.neq[ispin-1] << " per task)" << std::endl;

      std::cout << "\n";
      std::cout << " supercell:\n";
      std::cout << "      volume = " << Ffmt(10,2) << mylattice.omega() << std::endl;
      std::cout << "      lattice:    a1 = < " << Ffmt(8,3) << mylattice.unita(0,0) << " " << Ffmt(8,3) << mylattice.unita(1,0) << " " << Ffmt(8,3) << mylattice.unita(2,0) << " >\n";
      std::cout << "                  a2 = < " << Ffmt(8,3) << mylattice.unita(0,1) << " " << Ffmt(8,3) << mylattice.unita(1,1) << " " << Ffmt(8,3) << mylattice.unita(2,1) << " >\n";
      std::cout << "                  a3 = < " << Ffmt(8,3) << mylattice.unita(0,2) << " " << Ffmt(8,3) << mylattice.unita(1,2) << " " << Ffmt(8,3) << mylattice.unita(2,2) << " >\n";
      std::cout << "      reciprocal: b1 = < " << Ffmt(8,3) << mylattice.unitg(0,0) << " " << Ffmt(8,3) << mylattice.unitg(1,0) << " " << Ffmt(8,3) << mylattice.unitg(2,0) << " >\n";
      std::cout << "                  b2 = < " << Ffmt(8,3) << mylattice.unitg(0,1) << " " << Ffmt(8,3) << mylattice.unitg(1,1) << " " << Ffmt(8,3) << mylattice.unitg(2,1) << " >\n";
      std::cout << "                  b3 = < " << Ffmt(8,3) << mylattice.unitg(0,2) << " " << Ffmt(8,3) << mylattice.unitg(1,2) << " " << Ffmt(8,3) << mylattice.unitg(2,2) << " >\n";

      {double aa1,bb1,cc1,alpha1,beta1,gamma1;
       mylattice.abc_abg(&aa1,&bb1,&cc1,&alpha1,&beta1,&gamma1);
       std::cout << "      lattice:    a =    " << Ffmt(8,3) << aa1    << " b =   " << Ffmt(8,3) << bb1   << " c =    " << Ffmt(8,3) << cc1 << std::endl;
       std::cout << "                  alpha =" << Ffmt(8,3) << alpha1 << " beta =" << Ffmt(8,3) << beta1 << " gamma =" << Ffmt(8,3) << gamma1<< std::endl;}
      std::cout << "      density cutoff =" << Ffmt(7,3) << mylattice.ecut()
                << " fft =" << Ifmt(4) << mygrid.nx << " x " << Ifmt(4) << mygrid.ny << " x " << Ifmt(4) << mygrid.nz
                << "  (" << Ifmt(8) << mygrid.npack_all(0) << " waves " << Ifmt(8) << mygrid.npack(0) << " per task)" << std::endl;
      std::cout << "      wavefnc cutoff =" << Ffmt(7,3) << mylattice.wcut()
                << " fft =" << Ifmt(4) << mygrid.nx << " x " << Ifmt(4) << mygrid.ny << " x " << Ifmt(4) << mygrid.nz
                << "  (" << Ifmt(8) << mygrid.npack_all(1) << " waves " << Ifmt(8) << mygrid.npack(1) << " per task)" << std::endl;
      std::cout << "\n";
      std::cout << " Ewald parameters:\n";
      std::cout << "      energy cutoff = " << Ffmt(7,3) << myewald.ecut()
                << " fft =" << Ifmt(4) << myewald.nx() << " x " << Ifmt(4) << myewald.ny() << " x " << Ifmt(4) << myewald.nz()
                << "  (" << Ifmt(8) << myewald.npack_all() << " waves " << Ifmt(8) << myewald.npack() << " per task)" << std::endl;
      std::cout << "      Ewald summation: cut radius = " << Ffmt(7,3) << myewald.rcut() << " and " << Ifmt(3) << myewald.ncut() << std::endl;
      std::cout << "                       Mandelung Wigner-Seitz =" << Ffmt(12,8) << myewald.mandelung()
                << " (alpha =" << Ffmt(12,8) << myewald.rsalpha() << " rs =" << Ffmt(12,8) << myewald.rs() << ")" << std::endl;

      if (flag > 0) 
      {
         std::cout << std::endl;
         std::cout << " technical parameters:\n";
         std::cout << "      fixed step: time step =" << Ffmt(12,2) << control.time_step() << "  ficticious mass =" << Ffmt(12,2) << control.fake_mass() << std::endl;
         std::cout << "      tolerance =" << Efmt(12,3) << control.tolerances(0) << " (energy) "
                                          << Efmt(12,3) << control.tolerances(1) << " (density) "
                                          << Efmt(12,3) << control.tolerances(2) << " (ion)\n";
         std::cout << "      max iterations = " << Ifmt(10) << control.loop(0)*control.loop(1)
                   << " (" << Ifmt(5) << control.loop(0) << " inner " << Ifmt(5) << control.loop(1) << " outer)\n";
         if (control.minimizer()==1) std::cout << "      minimizer = Grassmann conjugate gradient\n";
         if (control.minimizer()==2) std::cout << "      minimizer = Grassmann lmbfgs\n";
         if (control.minimizer()==4) std::cout << "      minimizer = Stiefel conjugate gradient\n";
         if (control.minimizer()==5) std::cout << "      minimizer = scf (potential)\n";
         if (control.minimizer()==7) std::cout << "      minimizer = Stiefel lmbfgs\n";
         if (control.minimizer()==8) std::cout << "      minimizer = scf (density)\n";
         if ((control.minimizer()==5) || (control.minimizer()==8))
         {
            std::cout << std::endl;
            std::cout << " Kohn-Sham scf parameters:\n";
            std::cout << "     Kohn-Sham algorithm  = conjugate gradient\n";
            std::cout << "     SCF algorithm        = simple mixing\n";
            std::cout << "     SCF mixing parameter =    x.xxxx\n";
            std::cout << "     Kohn-Sham iterations = xxxx\n";
            if (control.minimizer()==5) std::cout << "     SCF mixing type      = potential\n";
            if (control.minimizer()==8) std::cout << "     SCF mixing type      = density\n";
            std::cout << "     Kerker damping       =    x.xxxx\n";
         }
      }
      else
      {
         std::cout << std::endl;
         std::cout << " technical parameters:\n";
         std::cout << "      optimization of psi and densities turned off" << std::endl;
      }
      std::cout << std::endl << std::endl << std::endl;

      std::cout << " -----------------------------------------------------------------------------------\n";
      std::cout << " ----------------------- Born-Oppenheimer molecular dynamics -----------------------\n";
      std::cout << " -----------------------------------------------------------------------------------\n\n";
      std::cout << "\n";
      std::cout << " molecular dynamics parameters:\n";
      if (myion.fix_translation) std::cout << "      translation constrained\n";
      if (myion.fix_rotation)    std::cout << "      rotation constrained\n";

      std::cout << "      time step =" << Ffmt(11,2) << control.bo_time_step()
                << " iterations = "    << Ifmt(10) <<  control.bo_steps(0)*control.bo_steps(1)
                << " ("  << Ifmt(5) << control.bo_steps(0) << " inner " << Ifmt(5) << control.bo_steps(1) << " outer)\n";
      if (control.bo_algorithm()==0) std::cout << "      integration algorithm = position Verlet\n";
      if (control.bo_algorithm()==1) std::cout << "      integration algorithm = velocity Verlet\n";
      if (control.bo_algorithm()==2) std::cout << "      integration algorithm = leap frog\n";


      std::cout << std::endl;
      std::cout << " scaling parameters:   " << std::endl;
      std::cout << "      cooling/heating rate   =" << Efmt(12,5) << control.ion_scaling() << " (ion)"    << std::endl;
      std::cout << "      initial kinetic energy =" << Efmt(12,5) << myion.eki0            << " (ion)"    << std::endl;
      std::cout << "                           "    << Efmt(15,5) << myion.ekg             << " (C.O.M.)" << std::endl;
      std::cout << "      after scaling          =" << Efmt(12,5) << myion.eki1            << " (ion)"    << std::endl;
      std::cout << "      increased energy       =" << Efmt(12,5) << myion.eki1-myion.eki0 << " (ion)"    << std::endl;
      std::cout << std::endl;

      if (mynose.on())
         std::cout << mynose.inputprint();
      else
         std::cout << " constant energy simulation" << std::endl;

      if (SA) std::cout << "      SA decay rate = " << Ffmt(10,3) << sa_decay[1] << "(ion)\n";

   }
   if (myparallel.is_master()) seconds(&cpu2);


//*                |***************************|
//******************     call GeoVibminimizer  **********************
//*                |***************************|

   /* intialize the linesearch */
   util_linesearch_init(); 

   /*  calculate energy and gradient */
   double g,gg,Gmax,Grms,Xrms,Xmax;
   double Eold =  0.0;
   double EV   = 0.0;
   double eki  = 0.0;
   double *Emol  = mymolecule.E;

   int nfsize   = 3*myion.nion;
   int one      = 1;
   double mrone = -1.0;

   // allocate temporary memory from stack
   double fion[3*myion.nion];

   bool done = false;

   /*  calculate energy */
   if (oprint) {

      std::cout << "\n\n";
      std::cout << " -----------------------------------------------------------------------------------\n";
      std::cout << " -----------------------------    Initial Geometry     -----------------------------\n";
      std::cout << " -----------------------------------------------------------------------------------\n\n";
      std::cout << " ---------------------------------\n";
      std::cout << "         Initial Geometry         \n";
      std::cout << " ---------------------------------\n";
      std::cout <<  mymolecule.myion->print_bond_angle_torsions();
      std::cout << "\n\n\n";
      std::cout << "         ================ Born-Oppenheimer AIMD iteration ================\n";
      std::cout << "     >>> iteration started at " << util_date() << " <<<\n";
      std::cout << "     iter.          KE+Energy             Energy        KE_Ion   Temperature\n";
      std::cout << "     -----------------------------------------------------------------------\n";
   }

   nwpw_aimd_running_data mymotion_data(control,&myparallel,&mygrid,&myion,
                                        E,mymolecule.hml,mymolecule.psi1,mymolecule.rho1);

   if (control.bo_steps(1)>0)
   {
      bool nose = mynose.on();
      bool verlet   = (control.bo_algorithm()==0);
      bool vverlet  = (control.bo_algorithm()==1);
      bool leapfrog = (control.bo_algorithm()==2);
      double dt  = control.bo_time_step();
      int it_out = control.bo_steps(1);
      int it_in  = control.bo_steps(0);
      double r   = 1.0;
      double ssr = 1.0;
      double eki1 = 0.0;
      icount = 0;

      double fion1[3*myion.nion];

      EV = cgsd_energy(control,mymolecule,false,std::cout);
      cgsd_energy_gradient(mymolecule,fion);

      if (nose) r = (1.0-0.5*dt*mynose.dXr());
      myion.Newton_step(fion,sa_alpha[1]*r);

      eki = myion.eki1;
      if (nose) mynose.Newton_step(0.01,eki);


      E[0] = EV+eki; 
      E[1] = EV; 
      E[2] = 0.0; 
      E[3] = eki; 
      for (auto i=0; i<56; ++i) E[i+4] = Emol[i];

      done = false;
      // outer loop iteration
      while (!done)
      {
         ++icount;

         // inner loop iteration
         for (auto it=0; it<it_in; ++it)
         {
            if (vverlet)
            {
               memcpy(fion1,fion,3*myion.nion*sizeof(double));
               myion.shift21();
            }
            else 
            {
               myion.shift();
               if (nose) mynose.shift();
            }

            //  calculate the energy and gradient 
            EV = cgsd_energy(control,mymolecule,false,std::cout);
            cgsd_energy_gradient(mymolecule,fion); 

            if (nose)
            {
               ssr = mynose.ssr();
               myion.Nose_step(ssr,fion);
               eki = myion.eki1;
               mynose.Verlet_step(0.01,eki);
            }
            else if (vverlet)
            {
               myion.vVerlet_step(fion,fion1);
               myion.vshift();
               myion.Newton_step(fion,sa_alpha[1]*r);
               eki = myion.eki1;
            }
            else
            {
               myion.Verlet_step(fion,sa_alpha[1]);
               eki = myion.eki1;
            }
  
            E[0] = EV+eki; 
            E[1] = EV; 
            E[2] = 0.0; 
            E[3] = eki; 
            for (auto i=0; i<56; ++i) E[i+4] = Emol[i];

         } // end inner loop

         // Update AIMD Running data
         mymotion_data.update_iteration(icount);
         //eave = mymotion_data.eave; evar = mymotion_data.evar;
         //have = mymotion_data.have; hvar = mymotion_data.hvar;
         //qave = mymotion_data.qave; qvar = mymotion_data.qvar;

         // Write out loop energies
         if (oprint)
         {
            if (SA)
            {
               if (mynose.on())
                  std::cout << Ifmt(10) << (icount*it_in)
                            << Efmt(19,10) << EV+eki+mynose.r_energy()
                            << Efmt(19,10) << EV
                            << Efmt(14,5)  << eki
                            << Ffmt(9,1)   << Tr_new << std::endl;
               else
                  std::cout << Ifmt(10) << (icount*it_in)
                            << Efmt(19,10) << EV+eki
                            << Efmt(19,10) << EV
                            << Efmt(14,5)  << eki
                            << Ffmt(9,1)   << Tr_new << std::endl;
            }
            else
            {
               if (mynose.on())
                  std::cout << Ifmt(10)    << (icount*it_in)
                            << Efmt(19,10) << EV+eki+mynose.r_energy()
                            << Efmt(19,10) << EV
                            << Efmt(14,5)  << eki
                            << Ffmt(14,2)  << myion.Temperature() << std::endl;
               else
                  std::cout << Ifmt(10) << (icount*it_in)
                            << Efmt(19,10) << EV+eki
                            << Efmt(19,10) << EV
                            << Efmt(14,5)  << eki
                            << Ffmt(14,2)  << myion.Temperature() << std::endl;
            }
         }


         // Check for running out of time
         if (control.out_of_time())
         {   
            done = true;
            if (oprint) std::cout << "         *** out of time. iteration terminated." << std::endl;
         }
         // Check for Completion
         else if (icount>=it_out)
         {
            done = true;
            if (oprint) std::cout << "         *** arrived at the Maximum iteration.   terminated." << std::endl;
         }
      } // end outer loop 

   } // end bomd iterations


   if (oprint) {
      std::cout << "\n\n";
      std::cout << " ---------------------------------\n";
      std::cout << "  Final Geometry \n";
      std::cout << " ---------------------------------\n";
      std::cout <<  mymolecule.myion->print_bond_angle_torsions();
      std::cout << "\n\n";
   }
   if (myparallel.is_master()) seconds(&cpu3);


//*******************************************************************


   // write energy results to the json
   auto rtdbjson = json::parse(rtdbstring);
   rtdbjson["pspw"]["energy"]   = EV;
   rtdbjson["pspw"]["energies"] = mymolecule.E;
   rtdbjson["pspw"]["eigenvalues"]  = mymolecule.eig_vector();

   // APC analysis
   if (mypsp.myapc->apc_on)
   {
      if (!(mypsp.myapc->v_apc_on))
         mypsp.myapc->gen_APC(mymolecule.dng1,false);

      // set qion
      double qion[myion.nion];
      for (auto ii=0; ii<myion.nion; ++ii)
         qion[ii] = -mypsp.myapc->Qtot_APC(ii) + mypsp.zv[myion.katm[ii]];
      rtdbjson["nwpw"]["apc"]["q"] = std::vector<double>(qion,&qion[myion.nion]);

      if (oprint)
      {
         std::cout <<  mypsp.myapc->print_APC(mypsp.zv);
      }
   }

   // set rtdbjson initialize_wavefunction option to false
   if (rtdbjson["nwpw"]["initialize_wavefunction"].is_boolean()) rtdbjson["nwpw"]["initialize_wavefunction"] = false;

   MPI_Barrier(comm_world0);

   /* write psi */
   if (flag > 0) mymolecule.writepsi(control.output_movecs_filename());

   /* write rtdbjson */
   rtdbstring    = rtdbjson.dump();
   myion.writejsonstr(rtdbstring);


//                 |**************************|
// *****************   report consumed time   **********************
//                 |**************************|
   if (myparallel.is_master()) seconds(&cpu4);
   if (oprint) 
   {
      double t1 = cpu2-cpu1;
      double t2 = cpu3-cpu2;
      double t3 = cpu4-cpu3;
      double t4 = cpu4-cpu1;
      double av = t2/((double ) myelectron.counter);
      //std::cout.setf(ios::scientific);
      std::cout << std::scientific;
      std::cout << "\n";
      std::cout << " -----------------"    << "\n";
      std::cout << " cputime in seconds"   << "\n";
      std::cout << " prologue    : " << t1 << "\n";
      std::cout << " main loop   : " << t2 << "\n";
      std::cout << " epilogue    : " << t3 << "\n";
      std::cout << " total       : " << t4 << "\n";
      std::cout << " cputime/step: " << av << " ( " << myelectron.counter << " evaluations, " << util_linesearch_counter() << " linesearches)\n";
      std::cout << "\n";

      nwpw_timing_print_final(myelectron.counter);

      std::cout << "\n";
      std::cout << " >>> job completed at     " << util_date() << " <<<\n";

   }

   /* deallocate memory */
   gdevice_psi_dealloc();
   //delete [] fion;
   //delete [] sion;


   MPI_Barrier(comm_world0);

   return 0;
}

}


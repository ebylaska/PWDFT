
#include	<iostream>
#include	<cstdio>
#include	<cmath>
#include	<cstdlib>
#include	<string>
#include	<vector>
//
#include        "iofmt.hpp"
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
#include	"Kinetic.hpp"
#include	"Coulomb.hpp"
#include	"exchange_correlation.hpp"
#include	"Pseudopotential.hpp"
#include	"Electron.hpp"
#include	"Molecule.hpp"
#include	"inner_loop.hpp"
#include	"psi.hpp"
#include	"rtdb.hpp"
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
 *            pspw_geovib                 *
 *                                        *
 ******************************************/
int pspw_geovib(MPI_Comm comm_world0, std::string& rtdbstring)
{
   //Parallel myparallel(argc,argv);
   Parallel myparallel(comm_world0);
   //RTDB myrtdb(&myparallel, "eric.db", "old");

   int version,nfft[3],ne[2],ispin;
   int i,ii,ia,nn,ngrid[3],matype,nelem,icount;
   char date[26];
   double sum1,sum2,ev,zv;
   double cpu1,cpu2,cpu3,cpu4,cpustep;
   double E[50],deltae,deltac,deltar,viral,unita[9];

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

   for (ii=0; ii<50; ++ii) E[ii] = 0.0;

   if (myparallel.is_master()) seconds(&cpu1);
   if (oprint)
   {
      std::ios_base::sync_with_stdio();
      std::cout << "          *****************************************************\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *             PWDFT PSPW GeoVib Calculation         *\n";
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

   // initialize Molecule
   Molecule mymolecule(control.input_movecs_filename(),control.input_movecs_initialize(),
                       &mygrid,&myion,&mystrfac,&myewald,&myelectron,&mypsp);

   /* intialize the linesearch */
   util_linesearch_init();

   MPI_Barrier(comm_world0);


   // driver parameters
   int maxit       = control.driver_maxiter();
   double tol_Gmax = control.driver_gmax();
   double tol_Grms = control.driver_grms();
   double tol_Xrms = control.driver_xrms();
   double tol_Xmax = control.driver_xmax();
   double trust    = control.driver_trust();
   int lmbfgs_size = control.driver_lmbfgs_size();

   //GeoVib mygeovib(&molecule,control);

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
      std::cout << " ----------------------------- Geometry Optimization -------------------------------\n";
      std::cout << " -----------------------------------------------------------------------------------\n\n";
      std::cout << " Optimization parameters:  " << Ifmt(5) << maxit << " (maxiter) " 
                << Efmt(12,3) << tol_Gmax << " (gmax) " 
                << Efmt(12,3) << tol_Grms << " (grms) " 
                << Efmt(12,3) << tol_Xrms << " (xrms) " 
                << Efmt(12,3) << tol_Xmax << " (xmax)\n";

      std::cout << "    maximum number of steps       (maxiter) = " << Ifmt(4) << maxit << std::endl;
      std::cout << "    maximum gradient threshold       (Gmax) = " << Efmt(12,6) << tol_Gmax << std::endl;
      std::cout << "    rms gradient threshold           (Grms) = " << Efmt(12,6) << tol_Grms << std::endl;
      std::cout << "    rms cartesian step threshold     (Xrms) = " << Efmt(12,6) << tol_Xrms << std::endl;
      std::cout << "    maximum cartesian step threshold (Xmax) = " << Efmt(12,6) << tol_Xmax << std::endl;
      std::cout << std::endl;
      std::cout << "    fixed trust radius              (Trust) = " << Efmt(12,6) << trust << std::endl;
      std::cout << "    number lmbfgs histories   (lmbfgs_size) = " << Ifmt(4) << lmbfgs_size << std::endl;

   }
   if (myparallel.is_master()) seconds(&cpu2);


//*                |***************************|
//******************     call GeoVibminimizer  **********************
//*                |***************************|

   /*  calculate energy and gradient */
   double g,gg,Gmax,Grms,Xrms,Xmax;
   double Eold =  0.0;
   double EV   = 0.0;

   int nfsize   = 3*myion.nion;
   int one      = 1;
   double mrone = -1.0;

   // allocate temporary memory from stack
   double fion[3*myion.nion];
   double sion[3*myion.nion];

   bool done = false;
   int it = 0;

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
      std::cout << " ---------------------------------\n";
      std::cout << "     Calculate Initial Energy     \n";
      std::cout << " ---------------------------------\n\n";
   }
   EV = cgsd_energy(control,mymolecule,true,std::cout);
   /*  calculate the gradient */
   if (oprint){
      std::cout << "\n";
      std::cout << " ---------------------------------\n";
      std::cout << "    Calculate Initial Gradient    \n";
      std::cout << " ---------------------------------\n\n";
   }
   cgsd_energy_gradient(mymolecule,fion);
   if (oprint) {
      std::cout << " ion forces (au):" << "\n";
      for (auto ii=0; ii<mymolecule.myion->nion; ++ii)
         std::cout << " " << Ifmt(5) << ii+1 << " " << mymolecule.myion->symbol(ii)
                   << "  ( " << Ffmt(12,6) << fion[3*ii] << " " << Ffmt(12,6) << fion[3*ii+1] << " " << Ffmt(12,6) << fion[3*ii+2] << " )\n";
       std::cout << "      C.O.M. ( " << Ffmt(12,6) << mymolecule.myion->com_fion(fion,0) << " "
                                      << Ffmt(12,6) << mymolecule.myion->com_fion(fion,1) << " "
                                      << Ffmt(12,6) << mymolecule.myion->com_fion(fion,2) << " )\n";
       std::cout << "   |F|/nion  = " << std::setprecision(6) << std::fixed << std::setw(12) << mymolecule.myion->rms_fion(fion) << std::endl
                 << "   max|Fatom|= " << std::setprecision(6) << std::fixed << std::setw(12) << mymolecule.myion->max_fion(fion)
                 << "  (" << std::setprecision(3) << std::fixed << std::setw(8) << mymolecule.myion->max_fion(fion)*(27.2116/0.529177)  << " eV/Angstrom)"
                 << std::endl << std::endl;
   }
   DSCAL_PWDFT(nfsize,mrone,fion,one);

   /* initialize lmbfgs */
   nwpw_lmbfgs geom_lmbfgs(3*myion.nion,lmbfgs_size,myion.rion1,fion);


   /* update coords in mymolecule */
   mymolecule.myion->fixed_step(-trust,fion);
   mymolecule.myion->shift();

   Xmax = mymolecule.myion->xmax();
   Xrms = mymolecule.myion->xrms();
   Gmax = 0.0;
   Grms = 0.0;
   for (ii=0; ii<(myion.nion); ++ii)
   {
      gg = fion[3*ii]*fion[3*ii] + fion[3*ii+1]*fion[3*ii+1] + fion[3*ii+2]*fion[3*ii+2];
      Grms += gg;
      g = sqrt(gg); if (g>Gmax) Gmax = g;
   }
   Grms = sqrt(Grms)/((double) myion.nion);
   if ((Gmax<=tol_Gmax)&&(Grms<=tol_Grms)&&(Xrms<=tol_Xrms)&&(Xmax<=tol_Xmax)) done = true;

   if (oprint) {
       if (done)
       {
	  std::cout << "      ----------------------\n"
		    << "      Optimization converged\n"
		    << "      ----------------------\n";
       }
       if (it==0)
       {
          std::cout << std::endl << std::endl; 
          std::cout << "@ Step             Energy     Delta E     Gmax     Grms     Xrms     Xmax   Walltime\n";
          std::cout << "@ ---- ------------------ ----------- -------- -------- -------- -------- ----------\n";
       }
       else
       {
          std::cout << std::endl << std::endl;
          std::cout << "  Step             Energy     Delta E     Gmax     Grms     Xrms     Xmax   Walltime\n";
          std::cout << "  ---- ------------------ ----------- -------- -------- -------- -------- ----------\n";
       }
       seconds(&cpustep);
       std::cout << "@ " << Ifmt(4) << it << " " << Ffmt(18,9) << EV << " " << Efmt(11,3) << EV-Eold
                   << " " << Ffmt(8,5) << Gmax <<  " " << Ffmt(8,5) << Grms
                   << " " << Ffmt(8,5) << Xrms <<  " " << Ffmt(8,5) << Xmax
                   << " " << Ffmt(10,1) << cpustep-cpu1 << " " << Ifmt(9) << myelectron.counter << std::endl;
       //printf("@ %4d %18.9lf %11.3le %8.5lf %8.5lf %8.5lf %8.5lf %10.1lf %9d\n",it,EV,(EV-Eold),Gmax,Grms,Xrms,Xmax,cpustep-cpu1,myelectron.counter);

       if ((Gmax<=tol_Gmax)&&(Grms<=tol_Grms)&&(Xrms<=tol_Xrms)&&(Xmax<=tol_Xmax))
       {
	 std::cout << "                                            ok       ok       ok       ok\n";
       }
       else if ((Gmax<=tol_Gmax)||(Grms<=tol_Grms)||(Xrms<=tol_Xrms)||(Xmax<=tol_Xmax))
       {
	  std::string oktag("");
          std::cout << "                                     ";
          oktag = "  "; if (Gmax<=tol_Gmax) oktag = "ok"; std::cout << "       " << oktag;
          oktag = "  "; if (Grms<=tol_Grms) oktag = "ok"; std::cout << "       " << oktag;
          oktag = "  "; if (Xrms<=tol_Xrms) oktag = "ok"; std::cout << "       " << oktag;
          oktag = "  "; if (Xmax<=tol_Xmax) oktag = "ok"; std::cout << "       " << oktag << "\n";
       }
    }



   ++it;

   while ((!done) && (it<=maxit))
   {
   
      if (oprint) {
         std::cout << "\n\n";
         std::cout << " -----------------------------------------------------------------------------------\n";
         std::cout << " ----------------------------- Optimization Step " << std::setw(5) << it << " -----------------------------\n";
         std::cout << " -----------------------------------------------------------------------------------\n\n";
      }

      /* print out the current geometry */
      if (oprint) {
         std::cout << " ---------------------------------\n";
         std::cout << "  Geometry for Step " <<  it << "\n";
         std::cout << " ---------------------------------\n";
         std::cout <<  mymolecule.myion->print_bond_angle_torsions();
      }

      /*  calculate energy */
      if (oprint) {
         std::cout << "\n\n\n";
         std::cout << " ---------------------------------\n";
         std::cout << "  Calculate Energy for Step " <<  it << "\n";
         std::cout << " ---------------------------------\n\n";
      }
      Eold = EV;
      EV = cgsd_energy(control,mymolecule,true,std::cout);

      /*  calculate the gradient */
      if (oprint) {
         std::cout << "\n";
         std::cout << " ---------------------------------\n";
         std::cout << "  Calculate Gradient for Step " <<  it << "\n";
         std::cout << " ---------------------------------\n\n";
      }

      cgsd_energy_gradient(mymolecule,fion);

      if (oprint) {
	 std::cout << " ion forces (au):" << "\n";
         for (auto ii=0; ii<mymolecule.myion->nion; ++ii)
            std::cout << " " << Ifmt(5) << ii+1 << " " << mymolecule.myion->symbol(ii)
                      << "  ( " << Ffmt(12,6) << fion[3*ii] << " " << Ffmt(12,6) << fion[3*ii+1] << " " << Ffmt(12,6) << fion[3*ii+2] << " )\n";
          std::cout << "      C.O.M. ( " << Ffmt(12,6) << mymolecule.myion->com_fion(fion,0) << " " 
                                         << Ffmt(12,6) << mymolecule.myion->com_fion(fion,1) << " " 
                                         << Ffmt(12,6) << mymolecule.myion->com_fion(fion,2) << " )\n";
          std::cout << "   |F|/nion  = " << std::setprecision(6) << std::fixed << std::setw(12) << mymolecule.myion->rms_fion(fion) << std::endl 
	            << "   max|Fatom|= " << std::setprecision(6) << std::fixed << std::setw(12) << mymolecule.myion->max_fion(fion) 
		    << "  (" << std::setprecision(3) << std::fixed << std::setw(8) << mymolecule.myion->max_fion(fion)*(27.2116/0.529177)  << " eV/Angstrom)"
		    << std::endl << std::endl;
      }
      DSCAL_PWDFT(nfsize,mrone,fion,one);

      /* lmbfgs gradient */
      geom_lmbfgs.lmbfgs(myion.rion1,fion,sion);


      /* update coords in mymolecule */
      mymolecule.myion->fixed_step(-trust,sion);
      mymolecule.myion->shift();


      Xmax = mymolecule.myion->xmax();
      Xrms = mymolecule.myion->xrms();
      Gmax = 0.0;
      Grms = 0.0;
      for (ii=0; ii<(myion.nion); ++ii)
      {
         gg = fion[3*ii]*fion[3*ii] + fion[3*ii+1]*fion[3*ii+1] + fion[3*ii+2]*fion[3*ii+2];
         Grms += gg; 
         g = sqrt(gg); if (g>Gmax) Gmax = g;
      }
      Grms = sqrt(Grms)/((double) myion.nion);
      if ((Gmax<=tol_Gmax)&&(Grms<=tol_Grms)&&(Xrms<=tol_Xrms)&&(Xmax<=tol_Xmax)) done = true;


      /* print out the current energy */
      if (oprint)
      {
         if (done)
         {
	    std::cout << "      ----------------------\n"
		      << "      Optimization converged\n"
		      << "      ----------------------\n";
         }

         if (it==0) 
         {
            std::cout << std::endl << std::endl;
            std::cout << "@ Step             Energy     Delta E     Gmax     Grms     Xrms     Xmax   Walltime\n";
            std::cout << "@ ---- ------------------ ----------- -------- -------- -------- -------- ----------\n";
         } 
         else 
         {
            std::cout << std::endl << std::endl;
            std::cout << "  Step             Energy     Delta E     Gmax     Grms     Xrms     Xmax   Walltime\n";
            std::cout << "  ---- ------------------ ----------- -------- -------- -------- -------- ----------\n";
         }
         seconds(&cpustep);
         std::cout << "@ " << Ifmt(4) << it << " " << Ffmt(18,9) << EV << " " << Efmt(11,3) << EV-Eold
                   << " " << Ffmt(8,5) << Gmax <<  " " << Ffmt(8,5) << Grms
                   << " " << Ffmt(8,5) << Xrms <<  " " << Ffmt(8,5) << Xmax
                   << " " << Ffmt(10,1) << cpustep-cpu1 << " " << Ifmt(9) << myelectron.counter << std::endl;
         //printf("@ %4d %18.9lf %11.3le %8.5lf %8.5lf %8.5lf %8.5lf %10.1lf %9d\n",it,EV,(EV-Eold),Gmax,Grms,Xrms,Xmax,cpustep-cpu1,myelectron.counter);
       
         if ((Gmax<=tol_Gmax)&&(Grms<=tol_Grms)&&(Xrms<=tol_Xrms)&&(Xmax<=tol_Xmax))
         {
	    std::cout << "                                            ok       ok       ok       ok\n";
         }
         else if ((Gmax<=tol_Gmax)||(Grms<=tol_Grms)||(Xrms<=tol_Xrms)||(Xmax<=tol_Xmax))
         {
            std::string oktag = "";
            std::cout << "                                     ";
            oktag = "  "; if (Gmax<=tol_Gmax) oktag = "ok"; std::cout << "       " << oktag;
            oktag = "  "; if (Grms<=tol_Grms) oktag = "ok"; std::cout << "       " << oktag;
            oktag = "  "; if (Xrms<=tol_Xrms) oktag = "ok"; std::cout << "       " << oktag;
            oktag = "  "; if (Xmax<=tol_Xmax) oktag = "ok"; std::cout << "       " << oktag << "\n";
         }
      }

     ++it;
   }

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




#include	<iostream>
#include	<cstdio>
#include	<cmath>
#include	<cstdlib>
#include	<string>
#include	<vector>
using namespace std;

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
int pspw_geovib(MPI_Comm comm_world0, string& rtdbstring)
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

   //double *psi1,*psi2,*Hpsi,*psi_r;
   //double *dn;
   //double *hml,*lmbda,*eig;

   for (ii=0; ii<50; ++ii) E[ii] = 0.0;

   if (myparallel.is_master())
   {
      seconds(&cpu1);
      ios_base::sync_with_stdio();
      cout << "          *****************************************************\n";
      cout << "          *                                                   *\n";
      cout << "          *             PWDFT PSPW GeoVib Calculation         *\n";
      cout << "          *                                                   *\n";
      cout << "          *  [ (Grassmann/Stiefel manifold implementation) ]  *\n";
      cout << "          *  [              C++ implementation             ]  *\n";
      cout << "          *                                                   *\n";
      cout << "          *              version #7.00   02/27/21             *\n";
      cout << "          *                                                   *\n";
      cout << "          *    This code was developed by Eric J. Bylaska,    *\n";
      cout << "          *    Abhishek Bagusetty, David H. Bross, ...        *\n";
      cout << "          *                                                   *\n";
      cout << "          *****************************************************\n";
      cout << "          >>> job started at       " << util_date() << " <<<\n";
   }

   //control_read(myrtdb);
   //control_read(myparallel.np(),rtdbstring);
   Control2 control(myparallel.np(),rtdbstring);
   int flag =  control.task();

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
   Molecule mymolecule(control.input_movecs_filename(),
                       &mygrid,&myion,&mystrfac,&myewald,&myelectron);

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

   if (myparallel.is_master())
   {
      cout << "\n";
      cout << "          ==============  summary of input  ==================\n";
      cout << "\n input psi filename: " << control.input_movecs_filename() << "\n";
      cout << "\n";
      cout << " number of processors used: " << myparallel.np() << "\n";
      cout << " processor grid           : " << myparallel.np_i() << " x" << myparallel.np_j() << "\n";
      if (mygrid.maptype==1) cout << " parallel mapping         : 1d-slab"    << "\n";
      if (mygrid.maptype==2) cout << " parallel mapping         : 2d-hilbert" << "\n";
      if (mygrid.maptype==3) cout << " parallel mapping         : 2d-hcurve" << "\n";
      if (mygrid.isbalanced()) 
         cout << " parallel mapping         : balanced" << "\n";
      else
         cout << " parallel mapping         : not balanced" << "\n";

      cout << "\n options:\n";
      cout << "   boundary conditions  = ";
      if (control.version==3) cout << "periodic\n";
      if (control.version==4) cout << "aperiodic\n";

      cout << "   electron spin        = ";
      if (ispin==1)
         cout << "restricted\n";
      else
         cout << "unrestricted\n";
      cout << myxc;
  
      //cout << "\n elements involved in the cluster:\n";
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
      cout << mypsp.print_pspall();

      printf("\n total charge:%8.3lf\n", control.total_charge());

      cout << "\n atom composition:" << "\n";
      for (ia=0; ia<myion.nkatm; ++ia)
         cout << "   " << myion.atom(ia) << " : " << myion.natm[ia];
      cout << "\n\n initial ion positions (au):" << "\n";
      for (ii=0; ii<myion.nion; ++ii)
         printf("%4d %s\t( %10.5lf %10.5lf %10.5lf ) - atomic mass = %6.3lf\n",ii+1,myion.symbol(ii),
                                               myion.rion1[3*ii],
                                               myion.rion1[3*ii+1],
                                               myion.rion1[3*ii+2],
                                               myion.amu(ii));
      printf("   G.C.\t( %10.5lf %10.5lf %10.5lf )\n", myion.gc(0), myion.gc(1), myion.gc(2));
      printf(" C.O.M.\t( %10.5lf %10.5lf %10.5lf )\n", myion.com(0),myion.com(1),myion.com(2));
      cout << "\n";
      printf(" number of electrons: spin up=%6d (%4d per task) down=%6d (%4d per task)\n",
             mygrid.ne[0],mygrid.neq[0],mygrid.ne[ispin-1],mygrid.neq[ispin-1]);

      cout << "\n";
      cout << " supercell:\n";
      printf("      volume = %10.2lf\n",mylattice.omega());
      printf("      lattice:    a1=< %8.3lf %8.3lf %8.3lf >\n",mylattice.unita(0,0),mylattice.unita(1,0),mylattice.unita(2,0));
      printf("                  a2=< %8.3lf %8.3lf %8.3lf >\n",mylattice.unita(0,1),mylattice.unita(1,1),mylattice.unita(2,1));
      printf("                  a3=< %8.3lf %8.3lf %8.3lf >\n",mylattice.unita(0,2),mylattice.unita(1,2),mylattice.unita(2,2));
      printf("      reciprocal: b1=< %8.3lf %8.3lf %8.3lf >\n",mylattice.unitg(0,0),mylattice.unitg(1,0),mylattice.unitg(2,0));
      printf("                  b2=< %8.3lf %8.3lf %8.3lf >\n",mylattice.unitg(0,1),mylattice.unitg(1,1),mylattice.unitg(2,1));
      printf("                  b3=< %8.3lf %8.3lf %8.3lf >\n",mylattice.unitg(0,2),mylattice.unitg(1,2),mylattice.unitg(2,2));

      {double aa1,bb1,cc1,alpha1,beta1,gamma1;
       mylattice.abc_abg(&aa1,&bb1,&cc1,&alpha1,&beta1,&gamma1);
       printf("      lattice:    a=    %8.3lf b=   %8.3lf c=    %8.3lf\n",aa1,bb1,cc1);
       printf("                  alpha=%8.3lf beta=%8.3lf gamma=%8.3lf\n",alpha1,beta1,gamma1);}

      printf("      density cutoff= %7.3lf fft= %4d x %4d x %4d  (%8d waves %8d per task)\n",
             mylattice.ecut(),mygrid.nx,mygrid.ny,mygrid.nz,mygrid.npack_all(0),mygrid.npack(0));
      printf("      wavefnc cutoff= %7.3lf fft= %4d x %4d x %4d  (%8d waves %8d per task)\n",
             mylattice.wcut(),mygrid.nx,mygrid.ny,mygrid.nz,mygrid.npack_all(1),mygrid.npack(1));

      cout << "\n";
      cout << " ewald parameters:\n";
      printf("      energy cutoff= %7.3lf fft= %4d x %4d x %4d  (%8d waves %8d per task)\n",
             myewald.ecut(),myewald.nx(),myewald.ny(),myewald.nz(),myewald.npack_all(),myewald.npack());
      printf("      Ewald summation: cut radius=  %7.3lf and %3d\n", myewald.rcut(),myewald.ncut());
      printf("                       Mandelung Wigner-Seitz= %12.8lf (alpha=%12.8lf rs=%11.8lf)\n",myewald.mandelung(),myewald.rsalpha(),myewald.rs());

      if (flag > 0) 
      {
         cout << std::endl;
         cout << " technical parameters:\n";
         printf("      time step= %11.2lf  ficticious mass=%11.2lf\n",
                control.time_step(),control.fake_mass());
         printf("      tolerance=%12.3le (energy) %12.3le (density) %12.3le (ion)\n",
                control.tolerances(0),control.tolerances(1),control.tolerances(2));
         //printf("      tolerance=%12.3le (energy) %12.3le (density)\n",
         //       control.tolerances(0),control.tolerances(1));
         printf("      max iterations = %10d (%5d inner %5d outer)\n",
                control.loop(0)*control.loop(1),control.loop(0),control.loop(1));
         if (control.minimizer()==1) cout << "      minimizer = Grassmann conjugate gradient\n";
         if (control.minimizer()==2) cout << "      minimizer = Grassmann lmbfgs\n";
         if (control.minimizer()==4) cout << "      minimizer = Stiefel conjugate gradient\n";
         if (control.minimizer()==5) cout << "      minimizer = scf (potential)\n";
         if (control.minimizer()==7) cout << "      minimizer = Stiefel lmbfgs\n";
         if (control.minimizer()==8) cout << "      minimizer = scf (density)\n";
         if ((control.minimizer()==5) || (control.minimizer()==8))
         {
            cout << std::endl;
            cout << " Kohn-Sham scf parameters:\n";
            cout << "     Kohn-Sham algorithm  = conjugate gradient\n";
            cout << "     SCF algorithm        = simple mixing\n";
            cout << "     SCF mixing parameter =    x.xxxx\n";
            cout << "     Kohn-Sham iterations = xxxx\n";
            if (control.minimizer()==5) cout << "     SCF mixing type      = potential\n";
            if (control.minimizer()==8) cout << "     SCF mixing type      = density\n";
            cout << "     Kerker damping       =    x.xxxx\n";
         }
      }
      else
      {
         cout << std::endl;
         cout << " technical parameters:\n";
         cout << "      optimization of psi and densities turned off" << std::endl;
      }
      cout << std::endl << std::endl << std::endl;

      cout << " -----------------------------------------------------------------------------------\n";
      cout << " ----------------------------- Geometry Optimization -------------------------------\n";
      cout << " -----------------------------------------------------------------------------------\n\n";
      printf(" Optimization parameters:  %4d (maxiter) %12.3le (gmax) %12.3le (grms) %12.3le (xrms) %12.3le (xmax)\n",
                maxit,tol_Gmax,tol_Grms,tol_Xrms,tol_Xmax);

      printf("    maximum number of steps       (maxiter) = %4d\n",maxit); 
      printf("    maximum gradient threshold       (Gmax) = %12.6le\n",tol_Gmax);
      printf("    rms gradient threshold           (Grms) = %12.6le\n",tol_Grms);
      printf("    rms cartesian step threshold     (Xrms) = %12.6le\n",tol_Xrms);
      printf("    maximum cartesian step threshold (Xmax) = %12.6le\n",tol_Xmax);
      printf("\n");
      printf("    fixed trust radius              (Trust) = %12.6le\n",trust);
      printf("    number lmbfgs histories   (lmbfgs_size) = %4d\n",lmbfgs_size); 

      seconds(&cpu2);
   }


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
   if (myparallel.is_master()){
      cout << "\n\n";
      cout << " -----------------------------------------------------------------------------------\n";
      cout << " -----------------------------    Initial Geometry     -----------------------------\n";
      cout << " -----------------------------------------------------------------------------------\n\n";
      cout << " ---------------------------------\n";
      cout << "         Initial Geometry         \n";
      cout << " ---------------------------------\n";
      cout <<  mymolecule.myion->print_bond_angle_torsions();
      cout << "\n\n\n";
      cout << " ---------------------------------\n";
      cout << "     Calculate Initial Energy     \n";
      cout << " ---------------------------------\n\n";
   }
   EV = cgsd_energy(control,mymolecule);
   /*  calculate the gradient */
   if (myparallel.is_master()){
      cout << "\n";
      cout << " ---------------------------------\n";
      cout << "    Calculate Initial Gradient    \n";
      cout << " ---------------------------------\n\n";
   }
   cgsd_energy_gradient(mymolecule,fion);
   if (myparallel.is_master()) {
      cout << " ion forces (au):" << "\n";
      for (auto ii=0; ii<mymolecule.myion->nion; ++ii)
         printf(" %5d %4s  ( %12.6lf %12.6lf %12.6lf )\n",
                         ii+1,mymolecule.myion->symbol(ii),
                         fion[3*ii],
                         fion[3*ii+1],
                         fion[3*ii+2]);
         printf("      C.O.M. ( %12.6lf %12.6lf %12.6lf )\n",
                         mymolecule.myion->com_fion(fion,0),
                         mymolecule.myion->com_fion(fion,1),
                         mymolecule.myion->com_fion(fion,2));
         cout << "   |F|/nion  = " << std::setprecision(6) << fixed << std::setw(12) << mymolecule.myion->rms_fion(fion) << std::endl
              << "   max|Fatom|= " << std::setprecision(6) << fixed << std::setw(12) << mymolecule.myion->max_fion(fion)
              << "  (" << std::setprecision(3) << fixed << std::setw(8) << mymolecule.myion->max_fion(fion)*(27.2116/0.529177)  << " eV/Angstrom)"
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

   if (myparallel.is_master()){
       if (done)
       {
          cout << "      ----------------------\n"
               << "      Optimization converged\n"
               << "      ----------------------\n";
       }
       if (it==0)
       {
          printf("\n\n");
          printf("@ Step             Energy     Delta E     Gmax     Grms     Xrms     Xmax   Walltime\n");
          printf("@ ---- ------------------ ----------- -------- -------- -------- -------- ----------\n");
       }
       else
       {
          printf("\n\n");
          printf("  Step             Energy     Delta E     Gmax     Grms     Xrms     Xmax   Walltime\n");
          printf("  ---- ------------------ ----------- -------- -------- -------- -------- ----------\n");
       }
       seconds(&cpustep);
       printf("@ %4d %18.9lf %11.3le %8.5lf %8.5lf %8.5lf %8.5lf %10.1lf %9d\n",it,EV,(EV-Eold),Gmax,Grms,Xrms,Xmax,cpustep-cpu1,myelectron.counter);

       if ((Gmax<=tol_Gmax)&&(Grms<=tol_Grms)&&(Xrms<=tol_Xrms)&&(Xmax<=tol_Xmax))
       {
          cout << "                                            ok       ok       ok       ok\n";
       }
       else if ((Gmax<=tol_Gmax)||(Grms<=tol_Grms)||(Xrms<=tol_Xrms)||(Xmax<=tol_Xmax))
       {
          std:string oktag = "";
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
   
      if (myparallel.is_master()){
         cout << "\n\n";
         cout << " -----------------------------------------------------------------------------------\n";
         cout << " ----------------------------- Optimization Step " << std::setw(5) << it << " -----------------------------\n";
         cout << " -----------------------------------------------------------------------------------\n\n";
      }

      /* print out the current geometry */
      if (myparallel.is_master()) {
         cout << " ---------------------------------\n";
         cout << "  Geometry for Step " <<  it << "\n";
         cout << " ---------------------------------\n";
         cout <<  mymolecule.myion->print_bond_angle_torsions();
      }

      /*  calculate energy */
      if (myparallel.is_master()){
         cout << "\n\n\n";
         cout << " ---------------------------------\n";
         cout << "  Calculate Energy for Step " <<  it << "\n";
         cout << " ---------------------------------\n\n";
      }
      Eold = EV;
      EV = cgsd_energy(control,mymolecule);

      /*  calculate the gradient */
      if (myparallel.is_master()){
         cout << "\n";
         cout << " ---------------------------------\n";
         cout << "  Calculate Gradient for Step " <<  it << "\n";
         cout << " ---------------------------------\n\n";
      }
      cgsd_energy_gradient(mymolecule,fion);
      if (myparallel.is_master()){
         cout << " ion forces (au):" << "\n";
         for (auto ii=0; ii<mymolecule.myion->nion; ++ii)
            printf(" %5d %4s  ( %12.6lf %12.6lf %12.6lf )\n",
                            ii+1,mymolecule.myion->symbol(ii),
                            fion[3*ii],
                            fion[3*ii+1],
                            fion[3*ii+2]);
            printf("      C.O.M. ( %12.6lf %12.6lf %12.6lf )\n",
                            mymolecule.myion->com_fion(fion,0),
                            mymolecule.myion->com_fion(fion,1),
                            mymolecule.myion->com_fion(fion,2));
            cout << "   |F|/nion  = " << std::setprecision(6) << fixed << std::setw(12) << mymolecule.myion->rms_fion(fion) << std::endl 
                 << "   max|Fatom|= " << std::setprecision(6) << fixed << std::setw(12) << mymolecule.myion->max_fion(fion) 
                 << "  (" << std::setprecision(3) << fixed << std::setw(8) << mymolecule.myion->max_fion(fion)*(27.2116/0.529177)  << " eV/Angstrom)"
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
      if (myparallel.is_master())
      {
         if (done)
         {
            cout << "      ----------------------\n"
                 << "      Optimization converged\n"
                 << "      ----------------------\n";
         }

         if (it==0) 
         {
            printf("\n\n");
            printf("@ Step             Energy     Delta E     Gmax     Grms     Xrms     Xmax   Walltime\n");
            printf("@ ---- ------------------ ----------- -------- -------- -------- -------- ----------\n");
         } 
         else 
         {
            printf("\n\n");
            printf("  Step             Energy     Delta E     Gmax     Grms     Xrms     Xmax   Walltime\n");
            printf("  ---- ------------------ ----------- -------- -------- -------- -------- ----------\n");
         }
         seconds(&cpustep);
         printf("@ %4d %18.9lf %11.3le %8.5lf %8.5lf %8.5lf %8.5lf %10.1lf %9d\n",it,EV,(EV-Eold),Gmax,Grms,Xrms,Xmax,cpustep-cpu1,myelectron.counter);
       
         if ((Gmax<=tol_Gmax)&&(Grms<=tol_Grms)&&(Xrms<=tol_Xrms)&&(Xmax<=tol_Xmax))
         {
            cout << "                                            ok       ok       ok       ok\n";
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

   if (myparallel.is_master()) {
      cout << "\n\n";
      cout << " ---------------------------------\n";
      cout << "  Final Geometry \n";
      cout << " ---------------------------------\n";
      cout <<  mymolecule.myion->print_bond_angle_torsions();
      cout << "\n\n";
   }
   if (myparallel.is_master()) seconds(&cpu3);

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

      if (myparallel.is_master())
         std::cout <<  mypsp.myapc->print_APC(mypsp.zv);
   }

   MPI_Barrier(comm_world0);


//*******************************************************************


   /* write psi */
   if (flag > 0) mymolecule.writepsi(control.output_movecs_filename());


   /* write rtdbjson */
   rtdbstring    = rtdbjson.dump();
   myion.writejsonstr(rtdbstring);

//                 |**************************|
// *****************   report consumed time   **********************
//                 |**************************|
   if (myparallel.is_master()) 
   {
      seconds(&cpu4);
      double t1 = cpu2-cpu1;
      double t2 = cpu3-cpu2;
      double t3 = cpu4-cpu3;
      double t4 = cpu4-cpu1;
      double av = t2/((double ) myelectron.counter);
      //cout.setf(ios::scientific);
      std::cout << std::scientific;
      cout << "\n";
      cout << " -----------------"    << "\n";
      cout << " cputime in seconds"   << "\n";
      cout << " prologue    : " << t1 << "\n";
      cout << " main loop   : " << t2 << "\n";
      cout << " epilogue    : " << t3 << "\n";
      cout << " total       : " << t4 << "\n";
      cout << " cputime/step: " << av << " ( " << myelectron.counter << " evaluations, " << util_linesearch_counter() << " linesearches)\n";
      cout << "\n";

      nwpw_timing_print_final(myelectron.counter);

      cout << "\n";
      cout << " >>> job completed at     " << util_date() << " <<<\n";

   }

   /* deallocate memory */
   gdevice_psi_dealloc();
   //delete [] fion;
   //delete [] sion;


   MPI_Barrier(comm_world0);

   return 0;
}

}



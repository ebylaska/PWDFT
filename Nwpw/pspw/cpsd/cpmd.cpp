
#include	<iostream>
#include	<cstdio>
#include	<cmath>
#include	<cstdlib>
#include	<string>
using namespace std;

#include	"Parallel.hpp"
//#include	"control.hpp"
#include	"Control2.hpp"
#include	"util_date.hpp"
#include	"Lattice.hpp"
#include	"PGrid.hpp"
#include	"Pneb.hpp"
#include	"Ion.hpp"
#include	"Ewald.hpp"
#include	"Strfac.hpp"
#include	"Kinetic.hpp"
#include	"Coulomb.hpp"
#include	"exchange_correlation.hpp"
#include	"Pseudopotential.hpp"
#include	"inner_loop_md.hpp"
#include	"psi.hpp"
#include	"rtdb.hpp"
#include	"mpi.h"

#include "json.hpp"
using json = nlohmann::json;


/******************************************
 *                                        *
 *                 cpmd                   *
 *                                        *
 ******************************************/
int cpmd(MPI_Comm comm_world0, string& rtdbstring)
{
   //Parallel myparallel(argc,argv);
   Parallel myparallel(comm_world0);
   //RTDB myrtdb(&myparallel, "eric.db", "old");

   bool verlet;
   int version,nfft[3],ne[2],ispin;
   int i,ii,ia,nn,ngrid[3],matype,nelem,icount,done;
   char date[26];
   double sum1,sum2,ev;
   double cpu1,cpu2,cpu3,cpu4;
   double E[50],deltae,deltac,deltar,viral,unita[9];
   double *psi0,*psi1,*psi2,*Hpsi,*psi_r;
   double *dn;
   double *hml,*lmbda,*eig;

   for (ii=0; ii<50; ++ii) E[ii] = 0.0;

   if (myparallel.is_master())
   {
      seconds(&cpu1);
      ios_base::sync_with_stdio();
      cout << "          *****************************************************\n";
      cout << "          *                                                   *\n";
      cout << "          *     Car-Parrinello calculation for molecules,     *\n";
      cout << "          *       microclusters, liquids, and materials       *\n";
      cout << "          *                                                   *\n";
      cout << "          *     [     extended Lagrangian molecular   ]       *\n";
      cout << "          *     [        dynamics simulation          ]       *\n";
      cout << "          *     [          C++ implementation         ]       *\n";
      cout << "          *                                                   *\n";
      cout << "          *            version #7.00   03/20/20               *\n";
      cout << "          *                                                   *\n";
      cout << "          *    This code was developed by Eric J. Bylaska     *\n";
      cout << "          *                                                   *\n";
      cout << "          *****************************************************\n";
      cout << "          >>> job started at       " << util_date() << " <<<\n";
   }

   //control_read(myrtdb);
   //control_read(myparallel.np(),rtdbstring);
   Control2 control(myparallel.np(),rtdbstring);
   Lattice mylattice(control);
   //myparallel.init2d(control_np_orbital());
   myparallel.init2d(control.np_orbital(),control.pfft3_qsize());

   /* initialize lattice, parallel grid structure */
   psi_get_header(&myparallel,&version,nfft,unita,&ispin,ne,control.input_movecs_filename());
   Pneb mygrid(&myparallel,&mylattice,control,ispin,ne);

   /* initialize psi0, psi1, and psi2 */
   psi0  = mygrid.g_allocate(1);
   psi1  = mygrid.g_allocate(1);
   psi2  = mygrid.g_allocate(1);
   Hpsi  = mygrid.g_allocate(1);
   psi_r = mygrid.h_allocate();
   dn    = mygrid.r_nalloc(ispin);
   hml   = mygrid.m_allocate(-1,1);
   lmbda = mygrid.m_allocate(-1,1);
   eig   = new double[ne[0]+ne[1]];

   /* read wavefunction */
   psi_read0(&mygrid,&version,nfft,unita,&ispin,ne,psi2,control.input_movecs_filename());

   /* ortho check */
   sum2  = mygrid.gg_traceall(psi2,psi2);
   sum1  = ne[0]+ne[1];
   if (ispin==1) sum1 *= 2;
   if (fabs(sum2-sum1)>1.0e-10)
   {
      if (myparallel.is_master())
         printf("Warning: Gram-Schmidt Being performed on psi2\n");
   }

   /* read wavefunction velocities */
   mygrid.g_zero(psi0);
   if (psi_filefind(&mygrid,control.input_v_movecs_filename()))
      psi_read0(&mygrid,&version,nfft,unita,&ispin,ne,psi0,control.input_v_movecs_filename());


   /* read in ion structure */
   //Ion myion(myrtdb);
   Ion myion(rtdbstring,control);

   /* setup structure factor */
   Strfac mystrfac(&myion, &mygrid);
   mystrfac.phafac();

   /* initialize operators */
   Kinetic_Operator mykin(&mygrid);
   Coulomb_Operator mycoulomb(&mygrid);
   XC_Operator      myxc(&mygrid,control);
   
   Pseudopotential mypsp(&myion,&mygrid,&mystrfac,control);

   /* setup ewald */
   Ewald myewald(&myparallel,&myion,&mylattice,control,mypsp.zv);
   myewald.phafac();



//                 |**************************|
// *****************   summary of input data  **********************
//                 |**************************|

   if (myparallel.is_master())
   {
      cout << "\n\n";
      cout << "          ==============  summary of input  ==================\n";
      cout << "\n input psi filename:  " << control.input_movecs_filename() << "\n";
      cout << " input vpsi filename: " << control.input_v_movecs_filename() << "\n";
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
      cout << "   ion motion           = ";
      if (control.geometry_optimize())
         cout << "yes\n";
      else
         cout << "no\n";
      cout << "   boundary conditions  = ";
      cout << "periodic\n";

      cout << "   electron spin        = ";
      if (ispin==1)
         cout << "restricted\n";
      else
         cout << "unrestricted\n";
      cout << "   exchange-correlation = ";
         cout << "LDA (Vosko et al) parameterization\n";
  
      cout << "\n elements involved in the cluster:\n";
      for (ia=0; ia<myion.nkatm; ++ia)
      {
         printf("    %2d : %4s   core charge: %4.1lf  lmax=%1d\n",
                 ia+1,myion.atom(ia),mypsp.zv[ia],mypsp.lmax[ia]);
         printf("           comment : %s\n",mypsp.comment[ia]);
         printf("           pseudopotential type            : %3d\n",mypsp.psp_type[ia]);
         printf("           highest angular component       : %3d\n",mypsp.lmax[ia]);
         printf("           local potential used            : %3d\n",mypsp.locp[ia]);
         printf("           number of non-local projections : %3d\n",mypsp.nprj[ia]);
         if (mypsp.semicore[ia])
            printf("           semicore corrections included   : %6.3lf (radius) %6.3lf (charge)\n",mypsp.rcore[ia],mypsp.ncore(ia));
         printf("           cutoff = ");
         for (ii=0; ii<=mypsp.lmax[ia]; ++ii)
            printf("%8.3lf",mypsp.rc[ia][ii]);
         printf("\n");
      }

      cout << "\n atom composition:" << "\n";
      for (ia=0; ia<myion.nkatm; ++ia)
         cout << "   " << myion.atom(ia) << " : " << myion.natm[ia];
      cout << "\n\n initial position of ions (au):" << "\n";
      for (ii=0; ii<myion.nion; ++ii)
         printf("%4d %s\t( %10.5lf %10.5lf %10.5lf ) - atomic mass = %6.3lf\n",ii+1,myion.symbol(ii),
                                               myion.rion(0,ii),
                                               myion.rion(1,ii),
                                               myion.rion(2,ii),
                                               myion.amu(ii));
      printf("   G.C.\t( %10.5lf %10.5lf %10.5lf )\n", myion.gc(0), myion.gc(1), myion.gc(2));
      printf(" C.O.M.\t( %10.5lf %10.5lf %10.5lf )\n", myion.com(0),myion.com(1),myion.com(2));

      cout << "\n\n initial velocity of ions (au):" << "\n";
      for (ii=0; ii<myion.nion; ++ii)
         printf("%4d %s\t( %10.5lf %10.5lf %10.5lf )\n",ii+1,myion.symbol(ii),
                                               myion.vion(0,ii),
                                               myion.vion(1,ii),
                                               myion.vion(2,ii));
      printf("   G.C.\t( %10.5lf %10.5lf %10.5lf )\n", myion.vgc(0), myion.vgc(1), myion.vgc(2));
      printf(" C.O.M.\t( %10.5lf %10.5lf %10.5lf )\n", myion.vcom(0),myion.vcom(1),myion.vcom(2));

      cout << "\n";
      printf(" number of electrons: spin up=%6d (%4d per task) down=%6d (%4d per task)\n",
             mygrid.ne[0],mygrid.neq[0],mygrid.ne[ispin-1],mygrid.neq[ispin-1]);

      cout << "\n";
      cout << " supercell:\n";
      printf("      volume : %10.2lf\n",mylattice.omega());
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

      cout << "\n";
      cout << " technical parameters:\n";
      printf("      time step= %11.2lf  ficticious mass=%11.2lf\n",
             control.time_step(),control.fake_mass());
      printf("      tolerance=%12.3le (energy) %12.3le (density) %12.3le (ion)\n",
             control.tolerances(0),control.tolerances(1),control.tolerances(2));
      printf("      max iterations = %10d (%5d inner %5d outer)\n",
             control.loop(0)*control.loop(1),control.loop(0),control.loop(1));
      cout << "\n\n\n";



   }

//                 |**************************|
// *****************   start iterations       **********************
//                 |**************************|
   if (myparallel.is_master()) 
   {
      seconds(&cpu2);
      cout << "         ================ Car-Parrinello iteration ================\n";
      cout << "     >>> iteration started at " << util_date() << " <<<\n";
      cout << "     iter.          KE+Energy             Energy        KE_psi        KE_Ion   Temperature\n";
      cout << "     -------------------------------------------------------------------------------------\n";


   }
   verlet = false;
   inner_loop_md(verlet,control,&mygrid,&myion,
                 &mykin,&mycoulomb,&myxc,
                 &mypsp,&mystrfac,&myewald,
                 psi1,psi2,Hpsi,psi_r,
                 dn,hml,lmbda,
                 E,&deltae,&deltac,&deltar);

   verlet = true;
   done   = 0;
   icount = 0;
   while (!done)
   {
      ++icount;
      inner_loop_md(verlet,control,&mygrid,&myion,
                 &mykin,&mycoulomb,&myxc,
                 &mypsp,&mystrfac,&myewald,
                 psi1,psi2,Hpsi,psi_r,
                 dn,hml,lmbda,
                 E,&deltae,&deltac,&deltar);

      if (myparallel.is_master())
         printf("%10d%19.10le%19.10le%14.5le%14.5le%14.2lf\n",icount*control.loop(0),
                                       E[0],E[1],E[2],E[3],0.0);

      /* check for competion */
      if (icount>=control.loop(1))
      {
         done = 1;
         if (myparallel.is_master())
            cout << "         *** arrived at the Maximum iteration.   terminated\n";
      }
   }
   if (myparallel.is_master()) 
   {
      seconds(&cpu3);
      cout << "     >>> iteration ended at   " << util_date() << " <<<\n";
   }


   /* diagonalize the hamiltonian */
   mygrid.m_diagonalize(hml,eig);


//                  |***************************|
// ****************** report summary of results **********************
//                  |***************************|
   if (myparallel.is_master()) 
   {
      cout << "\n\n";
      cout << "          =============  summary of results  =================\n";
      cout << "\n final position of ions (au):" << "\n";
      for (ii=0; ii<myion.nion; ++ii)
         printf("%4d %s\t( %10.5lf %10.5lf %10.5lf ) - atomic mass = %6.3lf\n",ii+1,myion.symbol(ii),
                                               myion.rion(0,ii),
                                               myion.rion(1,ii),
                                               myion.rion(2,ii),
                                               myion.amu(ii));
      printf("   G.C.\t( %10.5lf %10.5lf %10.5lf )\n", myion.gc(0), myion.gc(1), myion.gc(2));
      printf(" C.O.M.\t( %10.5lf %10.5lf %10.5lf )\n", myion.com(0),myion.com(1),myion.com(2));
      cout << "\n final velocity of ions (au):" << "\n";
      for (ii=0; ii<myion.nion; ++ii)
         printf("%4d %s\t( %10.5lf %10.5lf %10.5lf )\n",ii+1,myion.symbol(ii),
                                               myion.vion(0,ii),
                                               myion.vion(1,ii),
                                               myion.vion(2,ii));
      printf("   G.C.\t( %10.5lf %10.5lf %10.5lf )\n", myion.vgc(0), myion.vgc(1), myion.vgc(2));
      printf(" C.O.M.\t( %10.5lf %10.5lf %10.5lf )\n", myion.vcom(0),myion.vcom(1),myion.vcom(2));
      cout << "\n\n";
      printf(" total     energy    : %19.10le (%15.5le /ion)\n",      E[0],E[0]/myion.nion);
      printf(" total orbital energy: %19.10le (%15.5le /electron)\n", E[1],E[1]/(mygrid.ne[0]+mygrid.ne[1]));
      printf(" hartree energy      : %19.10le (%15.5le /electron)\n", E[2],E[2]/(mygrid.ne[0]+mygrid.ne[1]));
      printf(" exc-corr energy     : %19.10le (%15.5le /electron)\n", E[3],E[3]/(mygrid.ne[0]+mygrid.ne[1]));
      printf(" ion-ion energy      : %19.10le (%15.5le /ion)\n",      E[4],E[4]/myion.nion);
      printf("\n");
      printf(" K.S. kinetic energy : %19.10le (%15.5le /electron)\n",      E[5],E[5]/(mygrid.ne[0]+mygrid.ne[1]));
      printf(" K.S. V_l energy     : %19.10le (%15.5le /electron)\n",      E[6],E[6]/(mygrid.ne[0]+mygrid.ne[1]));
      printf(" K.S. V_nl energy    : %19.10le (%15.5le /electron)\n",      E[7],E[7]/(mygrid.ne[0]+mygrid.ne[1]));
      printf(" K.S. V_Hart energy  : %19.10le (%15.5le /electron)\n",      E[8],E[8]/(mygrid.ne[0]+mygrid.ne[1]));
      printf(" K.S. V_xc energy    : %19.10le (%15.5le /electron)\n",      E[9],E[9]/(mygrid.ne[0]+mygrid.ne[1]));
      viral = (E[9]+E[8]+E[7]+E[6])/E[5];
      printf(" Viral Coefficient   : %19.10le\n",viral);

      printf("\n orbital energies:\n"); 
      nn = ne[0] - ne[1];
      ev = 27.2116;
      for (i=0; i<nn; ++i)
      {
         printf("%18.7le",eig[i]); printf(" ("); printf("%8.3f",eig[i]*ev); printf("eV)\n");
      }
      for (i=0; i<ne[1]; ++i)
      {
         printf("%18.7le",eig[i+nn]); printf(" ("); printf("%8.3lf",eig[i+nn]*ev); printf("eV) ");
         printf("%18.7le",eig[i+(ispin-1)*ne[0]]); printf(" ("); printf("%8.3lf",eig[i+(ispin-1)*ne[0]]*ev); printf("eV)\n");
      }

      cout << "\n output psi filename:  " << control.output_movecs_filename() << "\n";
      cout << " output vpsi filename: " << control.output_v_movecs_filename() << "\n";
   }

   /* write wavefunction and velocity wavefunction */
   psi_write(&mygrid,&version,nfft,unita,&ispin,ne,psi1,control.output_movecs_filename());
   psi_write(&mygrid,&version,nfft,unita,&ispin,ne,psi0,control.output_v_movecs_filename());

   /* deallocate memory */
   mygrid.g_deallocate(psi0);
   mygrid.g_deallocate(psi1);
   mygrid.g_deallocate(psi2);
   mygrid.g_deallocate(Hpsi);
   mygrid.h_deallocate(psi_r);
   mygrid.r_dealloc(dn);
   mygrid.m_deallocate(hml);
   mygrid.m_deallocate(lmbda);
   delete [] eig;

   // write results to the json
   auto rtdbjson =  json::parse(rtdbstring);
   rtdbjson["pspw"]["energy"]   = E[0];
   rtdbjson["pspw"]["energies"] = E;
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
      double av = t2/((double ) control.loop(0)*icount);
      cout.setf(ios::scientific);
      cout << "\n";
      cout << " -----------------"    << "\n";
      cout << " cputime in seconds"   << "\n";
      cout << " prologue    : " << t1 << "\n";
      cout << " main loop   : " << t2 << "\n";
      cout << " epilogue    : " << t3 << "\n";
      cout << " total       : " << t4 << "\n";
      cout << " cputime/step: " << av << "\n";
      cout << "\n";
      cout << " >>> job completed at     " << util_date() << " <<<\n";

   }

   return 0;
}

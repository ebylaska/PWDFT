
#include	<iostream>
#include	<cstdio>
#include	<cmath>
#include	<cstdlib>
#include	<string>
using namespace std;

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
#include        "exchange_correlation.hpp"
#include	"Pseudopotential.hpp"
#include	"inner_loop.hpp"
#include	"psi.hpp"
#include	"rtdb.hpp"
#include	"mpi.h"

#include	"psp_library.hpp"
#include	"psp_file_check.hpp"
#include	"nwpw_timing.hpp"
#include        "gdevice.hpp"

#include "json.hpp"
using json = nlohmann::json;

namespace pwdft {
using namespace pwdft;

/******************************************
 *                                        *
 *                cpsd                    *
 *                                        *
 ******************************************/
int cpsd(MPI_Comm comm_world0, string& rtdbstring)
{
   //Parallel myparallel(argc,argv);
   Parallel myparallel(comm_world0);
   //RTDB myrtdb(&myparallel, "eric.db", "old");

   int version,nfft[3],ne[2],ispin;
   int i,ii,ia,nn,ngrid[3],matype,nelem,icount,done;
   char date[26];
   double sum1,sum2,ev,zv;
   double cpu1,cpu2,cpu3,cpu4;
   double E[60],deltae,deltac,deltar,viral,unita[9],en[2];
   double *psi1,*psi2,*Hpsi,*psi_r;
   double *dn;
   double *hml,*lmbda,*eig;

   for (ii=0; ii<60; ++ii) E[ii] = 0.0;

   if (myparallel.is_master())
   {
      seconds(&cpu1);
      ios_base::sync_with_stdio();
      cout << "          *****************************************************\n";
      cout << "          *                                                   *\n";
      cout << "          *     Car-Parrinello calculation for molecules,     *\n";
      cout << "          *       microclusters, liquids, and materials       *\n";
      cout << "          *                                                   *\n";
      cout << "          *     [     steepest descent minimization   ]       *\n";
      cout << "          *     [          C++ implementation         ]       *\n";
      cout << "          *                                                   *\n";
      cout << "          *            version #7.00   09/20/18               *\n";
      cout << "          *                                                   *\n";
      cout << "          *    This code was developed by Eric J. Bylaska     *\n";
      cout << "          *                                                   *\n";
      cout << "          *****************************************************\n";
      cout << "          >>> job started at       " << util_date() << " <<<\n";
   }

   //control_read(myrtdb);
   //control_read(myparallel.np(),rtdbstring);
   Control2 control(myparallel.np(),rtdbstring);

   /* initialize processor grid structure */
   myparallel.init2d(control.np_orbital(),control.pfft3_qsize());
   MPI_Barrier(comm_world0);

   /* initialize lattice */
   Lattice mylattice(control);

   /* read in ion structure */
   //Ion myion(myrtdb);
   Ion myion(rtdbstring,control);
   MPI_Barrier(comm_world0);

   /* Check for and generate psp files                       */
   /* - this routine also sets the valence charges in myion, */
   /*   and total_ion_charge and ne in control               */
   psp_file_check(&myparallel,&myion,control);
   MPI_Barrier(comm_world0);


   /* debug output - print charge, ispin, and ne */
/*
   if (myparallel.is_master()) 
   { 
       cout << endl;
       cout << "total_ion_charge = " << myion.total_zv() << endl;
       cout << "ispin = " << control.ispin() << endl;
       cout << "ne = " << control.ne(0) << " " << control.ne(1) << endl;
       cout << "ne = " << control.ne_ptr()[0] << " " << control.ne_ptr()[1] << endl;
   }
*/

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

   /* initialize psi1 and psi2 */
   psi1  = mygrid.g_allocate(1);
   psi2  = mygrid.g_allocate(1);
   Hpsi  = mygrid.g_allocate(1);
   psi_r = mygrid.h_allocate();
   dn    = mygrid.r_nalloc(ispin);
   hml   = mygrid.m_allocate(-1,1);
   lmbda = mygrid.m_allocate(-1,1);
   eig   = new double[ne[0]+ne[1]];
   gdevice_psi_alloc(mygrid.npack(1),mygrid.neq[0]+mygrid.neq[1]);

   //psi_read(&mygrid,&version,nfft,unita,&ispin,ne,psi2,control.input_movecs_filename());
   bool newpsi = psi_read(&mygrid,control.input_movecs_filename(),psi2);
   MPI_Barrier(comm_world0);


   /* setup structure factor */
   Strfac mystrfac(&myion, &mygrid);
   mystrfac.phafac();

   /* initialize operators */
   Kinetic_Operator mykin(&mygrid);
   Coulomb_Operator mycoulomb(&mygrid);
   XC_Operator      myxc(&mygrid,control);

   /* initialize psps */
   Pseudopotential mypsp(&myion,&mygrid,&mystrfac,control);

   /* setup ewald */
   Ewald myewald(&myparallel,&myion,&mylattice,control,mypsp.zv);
   myewald.phafac();



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
      cout << "   ion motion           = ";
      if (control.geometry_optimize())
         cout << "yes\n";
      else
         cout << "no\n";
      cout << "   boundary conditions  = ";
      if (control.version==3) cout << "periodic\n";
      if (control.version==4) cout << "aperiodic\n";

      cout << "   electron spin        = ";
      if (ispin==1)
         cout << "restricted\n";
      else
         cout << "unrestricted\n";
      cout << myxc;
      //cout << "   exchange-correlation = ";
      //cout << "LDA (Vosko et al) parameterization\n";
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

      cout << "\n total charge: " << control.total_charge() << "\n";

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
      printf("      volume : %10.2lf\n",mylattice.omega());
      printf("      lattice:    a1=< %8.3lf %8.3lf %8.3lf >\n",mylattice.unita(0,0),mylattice.unita(1,0),mylattice.unita(2,0));
      printf("                  a2=< %8.3lf %8.3lf %8.3lf >\n",mylattice.unita(0,1),mylattice.unita(1,1),mylattice.unita(2,1));
      printf("                  a3=< %8.3lf %8.3lf %8.3lf >\n",mylattice.unita(0,2),mylattice.unita(1,2),mylattice.unita(2,2));
      printf("      reciprocal: b1=< %8.3lf %8.3lf %8.3lf >\n",mylattice.unitg(0,0),mylattice.unitg(1,0),mylattice.unitg(2,0));
      printf("                  b2=< %8.3lf %8.3lf %8.3lf >\n",mylattice.unitg(0,1),mylattice.unitg(1,1),mylattice.unitg(2,1));
      printf("                  b3=< %8.3lf %8.3lf %8.3lf >\n",mylattice.unitg(0,2),mylattice.unitg(1,2),mylattice.unitg(2,2));
      {  double aa1,bb1,cc1,alpha1,beta1,gamma1;
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
      cout << "     ========================== iteration ==========================\n";
      cout << "          >>> iteration started at " << util_date() << "  <<<\n";
      cout << "     iter.             Energy       DeltaE     DeltaPsi     DeltaIon\n";
      cout << "     ---------------------------------------------------------------\n";


   }
   done   = 0;
   icount = 0;
   while (!done)
   {
      ++icount;
      inner_loop(control,&mygrid,&myion,
                 &mykin,&mycoulomb,&myxc,
                 &mypsp,&mystrfac,&myewald,
                 psi1,psi2,Hpsi,psi_r,
                 dn,hml,lmbda,
                 E,&deltae,&deltac,&deltar);

      if (myparallel.is_master())
         printf("%10d%19.10le%13.5le%13.5le%13.5le\n",icount*control.loop(0),
                                       E[0],deltae,deltac,deltar);

      /* check for competion */
      if ((deltae>0.0)&&(icount>1))
      {
         done = 1;
         if (myparallel.is_master())
            cout << "         *** Energy going up. iteration terminated\n";
      }
      else if ((fabs(deltae)<control.tolerances(0)) &&
               (deltac      <control.tolerances(1)) &&
               (deltar      <control.tolerances(2)))
      {
         done = 1;
         if (myparallel.is_master())
            cout << "         *** tolerance ok.    iteration terminated\n";
      }
      else if (icount>=control.loop(1))
      {
         done = 1;
         if (myparallel.is_master())
            cout << "          *** arrived at the Maximum iteration.   terminated ***\n";
      }
   }
   if (myparallel.is_master()) 
   {
      seconds(&cpu3);
      cout << "          >>> iteration ended at   " << util_date() << "  <<<\n";
   }


   /* diagonalize the hamiltonian */
   mygrid.m_diagonalize(hml,eig);


   /* calculate real-space number of electrons, en */
   {
      double omega = mylattice.omega();
      double scal1 = 1.0/((double) ((mygrid.nx)*(mygrid.ny)*(mygrid.nz)));
      double dv = omega*scal1;

      en[0] = dv*mygrid.r_dsum(dn);
      en[1] = en[0];
      if (ispin > 1)
         en[1] =  dv*mygrid.r_dsum(&dn[mygrid.n2ft3d]);
   }

   // calculate APC if v_apc_on==false
   if (mypsp.myapc->apc_on)
      if (!(mypsp.myapc->v_apc_on)) 
         mypsp.myapc->dngen_APC(dn,false);


//                  |***************************|
// ****************** report summary of results **********************
//                  |***************************|
   if (myparallel.is_master()) 
   {
      cout << "\n\n";
      cout << "          =============  summary of results  =================\n";
      cout << "\n final ion positions (au):" << "\n";
      for (ii=0; ii<myion.nion; ++ii)
         printf("%4d %s\t( %10.5lf %10.5lf %10.5lf ) - atomic mass = %6.3lf\n",ii+1,myion.symbol(ii),
                                               myion.rion1[3*ii],
                                               myion.rion1[3*ii+1],
                                               myion.rion1[3*ii+2],
                                               myion.amu(ii));
      printf("   G.C.\t( %10.5lf %10.5lf %10.5lf )\n", myion.gc(0), myion.gc(1), myion.gc(2));
      printf(" C.O.M.\t( %10.5lf %10.5lf %10.5lf )\n", myion.com(0),myion.com(1),myion.com(2));

      //if (mypsp.myapc->v_apc_on)
      //   cout << mypsp.myapc->shortprint_APC();

      cout << "\n\n";
      cout << fixed << " number of electrons: spin up= " << setw(11) << setprecision(5) << en[0]
                                           << "  down= " << setw(11) << setprecision(5) << en[ispin]
                                           << " (real space)";
      cout << "\n\n";
      printf(" total     energy    : %19.10le (%15.5le /ion)\n",      E[0],E[0]/myion.nion);
      printf(" total orbital energy: %19.10le (%15.5le /electron)\n", E[1],E[1]/(mygrid.ne[0]+mygrid.ne[1]));
      printf(" hartree energy      : %19.10le (%15.5le /electron)\n", E[2],E[2]/(mygrid.ne[0]+mygrid.ne[1]));
      printf(" exc-corr energy     : %19.10le (%15.5le /electron)\n", E[3],E[3]/(mygrid.ne[0]+mygrid.ne[1]));
      printf(" ion-ion energy      : %19.10le (%15.5le /ion)\n",      E[4],E[4]/myion.nion);
      if (mypsp.myapc->v_apc_on) 
         printf(" APC energy          : %19.10le (%15.5le /ion)\n",      E[51],E[51]/myion.nion);

      printf("\n");
      printf(" K.S. kinetic energy : %19.10le (%15.5le /electron)\n",      E[5],E[5]/(mygrid.ne[0]+mygrid.ne[1]));
      printf(" K.S. V_l energy     : %19.10le (%15.5le /electron)\n",      E[6],E[6]/(mygrid.ne[0]+mygrid.ne[1]));
      printf(" K.S. V_nl energy    : %19.10le (%15.5le /electron)\n",      E[7],E[7]/(mygrid.ne[0]+mygrid.ne[1]));
      printf(" K.S. V_Hart energy  : %19.10le (%15.5le /electron)\n",      E[8],E[8]/(mygrid.ne[0]+mygrid.ne[1]));
      printf(" K.S. V_xc energy    : %19.10le (%15.5le /electron)\n",      E[9],E[9]/(mygrid.ne[0]+mygrid.ne[1]));
      if (mypsp.myapc->v_apc_on) 
         printf(" K.S. V_APC energy   : %19.10le (%15.5le /ion)\n",      E[52],E[52]/myion.nion);

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
      cout << std::endl;

      // write APC analysis
      if (mypsp.myapc->apc_on)
         std::cout <<  mypsp.myapc->print_APC(mypsp.zv);
   }

   psi_write(&mygrid,&version,nfft,unita,&ispin,ne,psi1,control.output_movecs_filename());
   MPI_Barrier(comm_world0);

   /* deallocate memory */
   mygrid.g_deallocate(psi1);
   mygrid.g_deallocate(psi2);
   mygrid.g_deallocate(Hpsi);
   mygrid.h_deallocate(psi_r);
   mygrid.r_dealloc(dn);
   mygrid.m_deallocate(hml);
   mygrid.m_deallocate(lmbda);
   delete [] eig;
   gdevice_psi_dealloc();

   // write results to the json
   auto rtdbjson =  json::parse(rtdbstring);
   rtdbjson["pspw"]["energy"]   = E[0];
   rtdbjson["pspw"]["energies"] = E;
   if (mypsp.myapc->apc_on)
   {
      double qion[myion.nion];
      for (auto ii=0; ii<myion.nion; ++ii)
         qion[ii] = -mypsp.myapc->Qtot_APC(ii) + mypsp.zv[myion.katm[ii]];
      rtdbjson["nwpw"]["apc"]["q"] = std::vector<double>(qion,&qion[myion.nion]);
   }
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
      //cout.setf(ios::scientific);
      std::cout << std::scientific;
      std::cout << std::endl;
      std::cout << " -----------------"    << std::endl;
      std::cout << " cputime in seconds"   << std::endl;
      std::cout << " prologue    : " << t1 << std::endl;
      std::cout << " main loop   : " << t2 << std::endl;
      std::cout << " epilogue    : " << t3 << std::endl;
      std::cout << " total       : " << t4 << std::endl;
      std::cout << " cputime/step: " << av << std::endl;
      std::cout << std::endl;

      nwpw_timing_print_final(control.loop(0)*icount);

      std::cout << std::endl;
      std::cout << " >>> job completed at     " << util_date() << " <<<" << std::endl;;

   }

   MPI_Barrier(comm_world0);
   return 0;
}

}

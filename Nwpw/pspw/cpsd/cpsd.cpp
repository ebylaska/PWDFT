
#include	<iostream>
#include	<cstdio>
#include	<stdio.h>
#include	<cmath>
#include	<cstdlib>
#include	<string>
using namespace std;

#include	"Parallel.hpp"
#include	"control.hpp"
#include	"util_date.hpp"
#include	"PGrid.hpp"
#include	"Pneb.hpp"
#include	"Ion.hpp"
#include	"Ewald.hpp"
#include	"Strfac.hpp"
#include	"Kinetic.hpp"
#include	"Coulomb.hpp"
#include	"Pseudopotential.hpp"
#include	"inner_loop.hpp"
#include	"psi.hpp"
#include	"rtdb.hpp"
#include	"mpi.h"


//int cpsd(int argc, char *argv[])
int cpsd(MPI_Comm comm_world0, string rtdbstring)
{
   //Parallel myparallel(argc,argv);
   Parallel myparallel(comm_world0);
   RTDB myrtdb(&myparallel, "eric.db", "old");

   int version,nfft[3],ne[2],ispin;
   int i,ii,ia,nn,ngrid[3],matype,nelem,icount,done;
   char date[26];
   double sum1,sum2,ev;
   double cpu1,cpu2,cpu3,cpu4;
   double E[50],deltae,deltac,deltar,viral,unita[9];
   double *psi1,*psi2,*Hpsi,*psi_r;
   double *dn;
   double *hml,*lmbda,*eig;

   for (ii=0; ii<50; ++ii) E[ii] = 0.0;

   if (myparallel.is_master())
   {
      seconds(&cpu1);
      ios_base::sync_with_stdio();
      cout << "          *****************************************************\n";
      cout << "          *                                                   *\n";
      cout << "          *     Car-Parrinello microcluster calculation       *\n";
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

   control_read(myrtdb);
   lattice_init();
   myparallel.init2d(control_np_orbital());

   /* initialize lattice, parallel grid structure */
   psi_get_header(&myparallel,&version,nfft,unita,&ispin,ne);
   Pneb mygrid(&myparallel,ispin,ne);

   /* initialize psi1 and psi2 */
   psi1  = mygrid.g_allocate(1);
   psi2  = mygrid.g_allocate(1);
   Hpsi  = mygrid.g_allocate(1);
   psi_r = mygrid.h_allocate();
   dn    = mygrid.r_nalloc(ispin);
   hml   = mygrid.m_allocate(-1,1);
   lmbda = mygrid.m_allocate(-1,1);
   eig   = new double[ne[0]+ne[1]];

   psi_read(&mygrid,&version,nfft,unita,&ispin,ne,psi2);

   /* ortho check */
   sum2  = mygrid.gg_traceall(psi2,psi2);
   sum1  = ne[0]+ne[1];
   if (ispin==1) sum1 *= 2;
   if (fabs(sum2-sum1)>1.0e-10)
   {
      if (myparallel.is_master())
         printf("Warning: Gram-Schmidt Being performed on psi2\n");
   }



   /* read in ion structure */
   Ion myion(myrtdb);

   /* setup structure factor */
   Strfac mystrfac(&myion, &mygrid);
   mystrfac.phafac();

   /* initialize operators */
   Kinetic_Operator mykin(&mygrid);
   Coulomb_Operator mycoulomb(&mygrid);

   //Coulomb_Operator mycoul(mygrid);
   //XC_Operator      myxc(mygrid);

   
   Pseudopotential mypsp(&myion,&mygrid,&mystrfac);

   /* setup ewald */
   Ewald myewald(&myparallel,&myion, mypsp.zv);
   myewald.phafac();



//                 |**************************|
// *****************   summary of input data  **********************
//                 |**************************|

   if (myparallel.is_master())
   {
      cout << "\n\n\n";
      cout << "          ==============  summary of input  ==================\n";
      cout << "\n";
      cout << " number of processors used: " << myparallel.np() << "\n";
      cout << " processor grid           : " << myparallel.np_i() << " x" << myparallel.np_j() << "\n";
      if (mygrid.maptype==1) cout << " parallel mapping         : slab"    << "\n";
      if (mygrid.maptype==2) cout << " parallel mapping         : hilbert" << "\n";
      if (mygrid.isbalanced()) 
         cout << " parallel mapping         : balanced" << "\n";
      else
         cout << " parallel mapping         : not balanced" << "\n";

      cout << "\n options:\n";
      cout << "   electron spin = ";
      if (ispin==1)
         cout << "restricted\n";
      else
         cout << "unrestricted\n";
      cout << "\n input movecs:" << control_input_movecs_filename() << "\n";
  
      cout << "\n elements involved in the cluster:\n";
      for (ia=0; ia<myion.nkatm; ++ia)
      {
         printf("    %2d : %4s   core charge: %4.1lf  lmax=%1d\n",
                 ia,myion.atom(ia),mypsp.zv[ia],mypsp.lmax[ia]);
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
      cout << "\n\n initial ion positions (au):" << "\n";
      for (ii=0; ii<myion.nion; ++ii)
         printf("%4d %s\t %8.3lf %8.3lf %8.3lf\n",ii+1,myion.symbol(ii),
                                               myion.rion1[3*ii],
                                               myion.rion1[3*ii+1],
                                               myion.rion1[3*ii+2]);

      cout << "\n";
      printf(" number of electrons: spin up=%6d (%4d per task) down=%6d (%4d per task)\n",
             mygrid.ne[0],mygrid.neq[0],mygrid.ne[ispin-1],mygrid.neq[ispin-1]);

      cout << "\n";
      cout << " supercell:\n";
      printf("      volume : %10.2lf\n",lattice_omega());
      printf("      lattice:    a1=< %8.3lf %8.3lf %8.3lf >\n",lattice_unita(0,0),lattice_unita(1,0),lattice_unita(2,0));
      printf("                  a2=< %8.3lf %8.3lf %8.3lf >\n",lattice_unita(0,1),lattice_unita(1,1),lattice_unita(2,1));
      printf("                  a3=< %8.3lf %8.3lf %8.3lf >\n",lattice_unita(0,2),lattice_unita(1,2),lattice_unita(2,2));
      printf("      reciprocal: b1=< %8.3lf %8.3lf %8.3lf >\n",lattice_unitg(0,0),lattice_unitg(1,0),lattice_unitg(2,0));
      printf("                  b2=< %8.3lf %8.3lf %8.3lf >\n",lattice_unitg(0,1),lattice_unitg(1,1),lattice_unitg(2,1));
      printf("                  b3=< %8.3lf %8.3lf %8.3lf >\n",lattice_unitg(0,2),lattice_unitg(1,2),lattice_unitg(2,2));

      printf("      density cutoff= %7.3lf fft= %4d x %4d x %4d  (%8d waves %8d per task)\n",
             lattice_ecut(),mygrid.nx,mygrid.ny,mygrid.nz,mygrid.npack_all(0),mygrid.npack(0));
      printf("      wavefnc cutoff= %7.3lf fft= %4d x %4d x %4d  (%8d waves %8d per task)\n",
             lattice_wcut(),mygrid.nx,mygrid.ny,mygrid.nz,mygrid.npack_all(1),mygrid.npack(1));

      cout << "\n";
      cout << " ewald parameters:\n";
      printf("      energy cutoff= %7.3lf fft= %4d x %4d x %4d  (%8d waves %8d per task)\n",
             myewald.ecut(),myewald.nx(),myewald.ny(),myewald.nz(),myewald.npack_all(),myewald.npack());
      printf("      summation: cut radius=  %7.3lf and %3d   mandelung= %12.8lf\n",
             myewald.rcut(),myewald.ncut(),myewald.mandelung());

      cout << "\n";
      cout << " technical parameters:\n";
      printf("      time step= %10.2lf  ficticious mass=%10.2lf\n",
             control_time_step(),control_fake_mass());
      printf("      tolerance=%11.3le (energy) %11.3le (density) %11.3le (ion)\n",
             control_tolerances(0),control_tolerances(1),control_tolerances(2));
      printf("      max iterations = %10d (%5d inner %5d outer)\n",
             control_loop(0)*control_loop(1),control_loop(0),control_loop(1));
      cout << "\n\n\n";



   }

//                 |**************************|
// *****************   start iterations       **********************
//                 |**************************|
   if (myparallel.is_master()) 
   {
      seconds(&cpu2);
      cout << "          ================ iteration =========================\n\n";
      cout << "          >>> iteration started at " << util_date() << " <<<\n\n";
      cout << "   iter.             Energy       DeltaE     DeltaPsi     DeltaIon\n";
      cout << "   ---------------------------------------------------------------\n";


   }
   done   = 0;
   icount = 0;
   while (!done)
   {
      ++icount;
      inner_loop(&mygrid,&myion,
                 &mykin,&mycoulomb,
                 &mypsp,&mystrfac,&myewald,
                 psi1,psi2,Hpsi,psi_r,
                 dn,hml,lmbda,
                 E,&deltae,&deltac,&deltar);

      if (myparallel.is_master())
         printf("%8d%19.10le%13.5le%13.5le%13.5le\n",icount*control_loop(0),
                                       E[0],deltae,deltar,deltac);

      /* check for competion */
      if ((deltae>0.0)&&(icount>1))
      {
         done = 1;
         cout << "          *** Energy going up. iteration terminated\n";
      }
      else if ((fabs(deltae)<control_tolerances(0)) &&
               (deltac      <control_tolerances(1)) &&
               (deltar      <control_tolerances(2)))
      {
         done = 1;
         if (myparallel.is_master())
            cout << "          *** tolerance ok.    iteration terminated\n";
      }
      else if (icount>=control_loop(1))
      {
         done = 1;
         if (myparallel.is_master())
            cout << "          *** arrived at the Maximum iteration.   terminated\n";
      }
   }
   if (myparallel.is_master()) 
   {
      seconds(&cpu3);
      cout << "\n          >>> iteration ended at   " << util_date() << " <<<\n";
   }


   /* diagonalize the hamiltonian */
   mygrid.m_diagonalize(hml,eig);


//                  |***************************|
// ****************** report summary of results **********************
//                  |***************************|
   if (myparallel.is_master()) 
   {
      cout << "\n\n\n";
      cout << "          =============  summary of results  =================\n";
      cout << "\n final ion positions (au):" << "\n";
      for (ii=0; ii<myion.nion; ++ii)
         printf("%4d %s\t %8.3lf %8.3lf %8.3lf\n",ii+1,myion.symbol(ii),
                                               myion.rion1[3*ii],
                                               myion.rion1[3*ii+1],
                                               myion.rion1[3*ii+2]);
      cout << "\n\n";
      printf(" total     energy    : %19.10le (%15.5le /ion)\n",      E[0],E[0]/myion.nion);
      printf(" total orbtial energy: %19.10le (%15.5le /electron)\n", E[1],E[1]/(mygrid.ne[0]+mygrid.ne[1]));
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

   }

   /* deallocate memory */
   mygrid.g_deallocate(psi1);
   mygrid.g_deallocate(psi2);
   mygrid.g_deallocate(Hpsi);
   mygrid.h_deallocate(psi_r);
   mygrid.r_dealloc(dn);
   mygrid.m_deallocate(hml);
   mygrid.m_deallocate(lmbda);
   delete [] eig;

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
      double av = t2/((double ) control_loop(0)*icount);
      cout.setf(ios::scientific);
      cout << "\n";
      cout << "-----------------"    << "\n";
      cout << "cputime in seconds"   << "\n";
      cout << "prologue    : " << t1 << "\n";
      cout << "main loop   : " << t2 << "\n";
      cout << "epilogue    : " << t3 << "\n";
      cout << "total       : " << t4 << "\n";
      cout << "cputime/step: " << av << "\n";
      cout << "\n";
      cout << "          >>> job completed at     " << util_date() << " <<<\n";

   }

   return 0;
}

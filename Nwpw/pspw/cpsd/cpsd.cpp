
#include	<iostream>
#include	<cstdio>
#include	<cmath>
#include	<cstdlib>
#include	<string>

#include        "iofmt.hpp"
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
//#include	"rtdb.hpp"
#include	"mpi.h"

#include	"psp_library.hpp"
#include	"psp_file_check.hpp"
#include	"nwpw_timing.hpp"
#include        "gdevice.hpp"

#include "json.hpp"
using json = nlohmann::json;


namespace pwdft {


/******************************************
 *                                        *
 *                cpsd                    *
 *                                        *
 ******************************************/
int cpsd(MPI_Comm comm_world0, std::string& rtdbstring)
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

   Control2 control(myparallel.np(),rtdbstring);

   bool hprint = (myparallel.is_master() && control.print_level("high"));
   bool oprint = (myparallel.is_master() && control.print_level("medium"));
   bool lprint = (myparallel.is_master() && control.print_level("low"));

   /* reset Parallel base_stdio_print = lprint */
   myparallel.base_stdio_print = lprint;


   for (ii=0; ii<60; ++ii) E[ii] = 0.0;

   if (myparallel.is_master()) seconds(&cpu1);
   if (oprint)
   {
      std::ios_base::sync_with_stdio();
      std::cout << "          *****************************************************\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *     Car-Parrinello calculation for molecules,     *\n";
      std::cout << "          *       microclusters, liquids, and materials       *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *     [     steepest descent minimization   ]       *\n";
      std::cout << "          *     [          C++ implementation         ]       *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *            version #7.00   09/20/18               *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *    This code was developed by Eric J. Bylaska     *\n";
      std::cout << "          *                                                   *\n";
      std::cout << "          *****************************************************\n";
      std::cout << "          >>> job started at       " << util_date() << " <<<\n";
   }

   //control_read(myrtdb);
   //control_read(myparallel.np(),rtdbstring);

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
   psp_file_check(&myparallel,&myion,control,std::cout);
   MPI_Barrier(comm_world0);


   /* debug output - print charge, ispin, and ne */
/*
   if (myparallel.is_master()) 
   { 
       std::cout << endl;
       std::cout << "total_ion_charge = " << myion.total_zv() << endl;
       std::cout << "ispin = " << control.ispin() << endl;
       std::cout << "ne = " << control.ne(0) << " " << control.ne(1) << endl;
       std::cout << "ne = " << control.ne_ptr()[0] << " " << control.ne_ptr()[1] << endl;
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
   bool newpsi = psi_read(&mygrid,control.input_movecs_filename(),
                                  control.input_movecs_initialize(),psi2,std::cout);
   MPI_Barrier(comm_world0);


   /* setup structure factor */
   Strfac mystrfac(&myion, &mygrid);
   mystrfac.phafac();

   /* initialize operators */
   Kinetic_Operator mykin(&mygrid);
   Coulomb_Operator mycoulomb(&mygrid);
   XC_Operator      myxc(&mygrid,control);

   /* initialize psps */
   Pseudopotential mypsp(&myion,&mygrid,&mystrfac,control,std::cout);

   /* setup ewald */
   Ewald myewald(&myparallel,&myion,&mylattice,control,mypsp.zv);
   myewald.phafac();



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
      std::cout << "   ion motion           = ";
      if (control.geometry_optimize())
         std::cout << "yes\n";
      else
         std::cout << "no\n";
      std::cout << "   boundary conditions  = ";
      if (control.version==3) std::cout << "periodic\n";
      if (control.version==4) std::cout << "aperiodic\n";

      std::cout << "   electron spin        = ";
      if (ispin==1)
         std::cout << "restricted\n";
      else
         std::cout << "unrestricted\n";
      std::cout << myxc;
      //std::cout << "   exchange-correlation = ";
      //std::cout << "LDA (Vosko et al) parameterization\n";
  
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

      std::cout << "\n atom composition:" << std::endl;
      for (ia=0; ia<myion.nkatm; ++ia)
         std::cout << "   " << myion.atom(ia) << " : " << myion.natm[ia];
      std::cout << "\n\n initial ion positions (au):" << std::endl;
      for (ii=0; ii<myion.nion; ++ii)
         std::cout << Ifmt(4) << ii+1 << " " << myion.symbol(ii)
                   << "\t( " << Ffmt(10,5) << myion.rion1[3*ii] << " " << Ffmt(10,5) << myion.rion1[3*ii+1] << " " << Ffmt(10,5) << myion.rion1[3*ii+2]
                   << " ) - atomic mass = " << Ffmt(6,3) << myion.amu(ii) << std::endl;
      std::cout << "   G.C.\t( " << Ffmt(10,5) << myion.gc(0) << " " << Ffmt(10,5) << myion.gc(1) << " " << Ffmt(10,5) << myion.gc(2) << " )" << std::endl;
      std::cout << " C.O.M.\t( " << Ffmt(10,5) << myion.com(0) << " " << Ffmt(10,5) << myion.com(1) << " " << Ffmt(10,5) << myion.com(2) << " )" << std::endl;
      std::cout << std::endl;
      std::cout <<" number of electrons: spin up =" << Ifmt(6) << mygrid.ne[0] << " (" << Ifmt(4) << mygrid.neq[0]
                << " per task) down =" << Ifmt(6) << mygrid.ne[ispin-1] << " (" << Ifmt(4) << mygrid.neq[ispin-1] << " per task)" << std::endl;

      std::cout << std::endl;
      std::cout << " supercell:" << std::endl;
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
                << " fft= " << Ifmt(4) << myewald.nx() << " x " << Ifmt(4) << myewald.ny() << " x " << Ifmt(4) << myewald.nz()
                << "  (" << Ifmt(8) << myewald.npack_all() << " waves " << Ifmt(8) << myewald.npack() << " per task)" << std::endl;
      std::cout << "      Ewald summation: cut radius = " << Ffmt(7,3) << myewald.rcut() << " and " << Ifmt(3) << myewald.ncut() << std::endl;
      std::cout << "                       Mandelung Wigner-Seitz =" << Ffmt(12,8) << myewald.mandelung()
                << " (alpha =" << Ffmt(12,8) << myewald.rsalpha() << " rs =" << Ffmt(12,8) << myewald.rs() << ")" << std::endl;

      std::cout << std::endl;
      std::cout << " technical parameters:\n";
      std::cout << "      fixed step: time step =" << Ffmt(12,2) << control.time_step() << "  ficticious mass =" << Ffmt(12,2) << control.fake_mass() << std::endl;
      std::cout << "      tolerance =" << Efmt(12,3) << control.tolerances(0) << " (energy) "
                                       << Efmt(12,3) << control.tolerances(1) << " (density) "
                                       << Efmt(12,3) << control.tolerances(2) << " (ion)\n";
      std::cout << "      max iterations = " << Ifmt(10) << control.loop(0)*control.loop(1)
                << " (" << Ifmt(5) << control.loop(0) << " inner " << Ifmt(5) << control.loop(1) << " outer)\n";

   }

//                 |**************************|
// *****************   start iterations       **********************
//                 |**************************|
   if (myparallel.is_master()) seconds(&cpu2);
   if (oprint) 
   {
      std::cout << "     ========================== iteration ==========================\n";
      std::cout << "          >>> iteration started at " << util_date() << "  <<<\n";
      std::cout << "     iter.             Energy       DeltaE     DeltaPsi     DeltaIon\n";
      std::cout << "     ---------------------------------------------------------------\n";
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

      if (oprint)
         std::cout << Ifmt(10) << icount*control.loop(0) << Efmt(19,10) << E[0] 
                   << Efmt(13,5) << deltae << Efmt(13,5) << deltac << Efmt(13,5) << deltar << std::endl; 
         //printf("%10d%19.10le%13.5le%13.5le%13.5le\n",icount*control.loop(0),
         //                              E[0],deltae,deltac,deltar);

      /* check for competion */
      if ((deltae>0.0)&&(icount>1))
      {
         done = 1;
         if (oprint)
            std::cout << "         *** Energy going up. iteration terminated\n";
      }
      else if ((std::fabs(deltae)<control.tolerances(0)) &&
               (deltac      <control.tolerances(1)) &&
               (deltar      <control.tolerances(2)))
      {
         done = 1;
         if (oprint)
            std::cout << "         *** tolerance ok.    iteration terminated\n";
      }
      else if (icount>=control.loop(1))
      {
         done = 1;
         if (oprint)
            std::cout << "          *** arrived at the Maximum iteration.   terminated ***\n";
      }
   }
   if (myparallel.is_master()) seconds(&cpu3);
   if (oprint) 
   {
      std::cout << "          >>> iteration ended at   " << util_date() << "  <<<\n";
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
   if (oprint) 
   {
      std::cout << "\n\n";
      std::cout << "          =============  summary of results  =================\n";
      std::cout << "\n final ion positions (au):" << "\n";
      for (ii=0; ii<myion.nion; ++ii)
         std::cout << Ifmt(4) << ii+1 << " " << myion.symbol(ii) 
                   << "\t( " 
                   << Ffmt(10,5) << myion.rion1[3*ii]   << " " 
                   << Ffmt(10,5) << myion.rion1[3*ii+1] << " " 
                   << Ffmt(10,5) << myion.rion1[3*ii+2] << " ) - atomic mass = " << Ffmt(6,3) << myion.amu(ii) << std::endl; 
      std::cout << "   G.C.\t( " << Ffmt(10,5) << myion.gc(0) << " " 
                                 << Ffmt(10,5) << myion.gc(1) << " " 
                                 << Ffmt(10,5) << myion.gc(2) << " )" <<  std::endl; 
      std::cout << " C.O.M.\t( " << Ffmt(10,5) << myion.com(0) << " " 
                                 << Ffmt(10,5) << myion.com(1) << " " 
                                 << Ffmt(10,5) << myion.com(2) << " )" <<  std::endl; 

      //if (mypsp.myapc->v_apc_on)
      //   std::cout << mypsp.myapc->shortprint_APC();

      std::cout << "\n\n";
      std::cout << std::fixed << " number of electrons: spin up= " << std::setw(11) << std::setprecision(5) << en[0]
		<< "  down= " << std::setw(11) << std::setprecision(5) << en[ispin]
		<< " (real space)";
      std::cout << std::endl << std::endl;
      std::cout << " total     energy    : " << Efmt(19,10) << E[0] << " (" << Efmt(15,5) << E[0]/myion.nion << " /ion)" << std::endl;  
      std::cout << " total orbital energy: " << Efmt(19,10) << E[1] << " (" << Efmt(15,5) << E[1]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      std::cout << " hartree energy      : " << Efmt(19,10) << E[2] << " (" << Efmt(15,5) << E[2]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl; 
      std::cout << " exc-corr energy     : " << Efmt(19,10) << E[3] << " (" << Efmt(15,5) << E[1]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      if (mypsp.myapc->v_apc_on) 
         std::cout << " APC energy          : " << Efmt(19,10) << E[51] << " (" << Efmt(15,5)  << E[51]/myion.nion << " /ion)" << std::endl;
      std::cout << " ion-ion energy      : " << Efmt(19,10) << E[4] << " (" << Efmt(15,5) << E[4]/myion.nion << " /ion)" << std::endl;

      std::cout << std::endl;
      std::cout << " K.S. kinetic energy : " << Efmt(19,10) << E[5] << " (" << Efmt(15,5) << E[5]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      std::cout << " K.S. V_l energy     : " << Efmt(19,10) << E[6] << " (" << Efmt(15,5) << E[6]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      std::cout << " K.S. V_nl_energy    : " << Efmt(19,10) << E[7] << " (" << Efmt(15,5) << E[7]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      std::cout << " K.S. V_Hart energy  : " << Efmt(19,10) << E[8] << " (" << Efmt(15,5) << E[8]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;
      std::cout << " K.S. V_xc energy    : " << Efmt(19,10) << E[9] << " (" << Efmt(15,5) << E[9]/(mygrid.ne[0]+mygrid.ne[1]) << " /electron)" << std::endl;

      if (mypsp.myapc->v_apc_on) 
         std::cout << " K.S. V_APC energy   : " << Efmt(19,10) << E[52] << " (" << Efmt(15,5) << E[52]/myion.nion << " /ion)" << std::endl;

      viral = (E[9]+E[8]+E[7]+E[6])/E[5];
      std::cout << " Viral Coefficient   : " << Efmt(19,10) << viral << std::endl;

      std::cout << "\n orbital energies:\n"; 
      nn = ne[0] - ne[1];
      ev = 27.2116;
      for (i=0; i<nn; ++i)
      {
         std::cout << Efmt(18,7) << eig[i] << " (" << Ffmt(8,3) << eig[i]*ev << "eV)" << std::endl;
      }
      for (i=0; i<ne[1]; ++i)
      {
         std::cout << Efmt(18,7) << eig[i+nn] << eig[i]    << " (" << Ffmt(8,3) << eig[i+nn]*ev << eig[i]*ev << "eV) " 
                   << Efmt(18,7) << eig[i+(ispin-1)*ne[0]] << " (" << Ffmt(8,3) << eig[i+(ispin-1)*ne[0]]*ev << "eV)"  << std::endl;
         //printf("%18.7le",eig[i+nn]); printf(" ("); printf("%8.3lf",eig[i+nn]*ev); printf("eV) ");
         //printf("%18.7le",eig[i+(ispin-1)*ne[0]]); printf(" ("); printf("%8.3lf",eig[i+(ispin-1)*ne[0]]*ev); printf("eV)\n");
      }
      std::cout << std::endl;

      // write APC analysis
      if (mypsp.myapc->apc_on)
         std::cout <<  mypsp.myapc->print_APC(mypsp.zv);
   }

   psi_write(&mygrid,&version,nfft,unita,&ispin,ne,psi1,control.output_movecs_filename(),std::cout);
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

   // set rtdbjson initialize_wavefunction option to false
   if (rtdbjson["nwpw"]["initialize_wavefunction"].is_boolean()) rtdbjson["nwpw"]["initialize_wavefunction"] = false;

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
      double av = t2/((double ) control.loop(0)*icount);
      //std::cout.setf(ios::scientific);
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

      nwpw_timing_print_final(control.loop(0)*icount,std::cout);

      std::cout << std::endl;
      std::cout << " >>> job completed at     " << util_date() << " <<<" << std::endl;;

   }

   MPI_Barrier(comm_world0);
   return 0;
}

}

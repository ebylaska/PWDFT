#include        <iostream>
#include        <cstdio>
#include        <cmath>
#include        <cstdlib>
#include        <string>

#include        "Parallel.hpp"
#include        "Control2.hpp"
#include        "util_date.hpp"
#include        "Ion.hpp"
#include        "Pneb.hpp"
#include        "Molecule.hpp"
#include	"Geodesic.hpp"

#include	"cgsd.hpp"

namespace pwdft {

/******************************************
 *                                        *
 *            cgsd_noit_energy            *
 *                                        *
 ******************************************/
double cgsd_noit_energy(Molecule& mymolecule) 
{
   Parallel *parall = mymolecule.mygrid->d3db::parall;

   /* generate phase factors and local psp and semicore density */
   mymolecule.phafacs_vl_potential_semicore();

   double total_energy  = mymolecule.gen_all_energies();

   /* report summary of results */
   if (parall->is_master())
   {
      std::cout << "     ==================  optimization turned off  ===================\n" << std::endl;
      std::cout << std::endl;
      std::cout << mymolecule;
   }

   return total_energy;
}



/******************************************
 *                                        *
 *              cgsd_energy               *
 *                                        *
 ******************************************/
double cgsd_energy(Control2& control, Molecule& mymolecule)
{
   Parallel *parall = mymolecule.mygrid->d3db::parall;
   Pneb *mygrid = mymolecule.mygrid;
   Ion  *myion  = mymolecule.myion;

   bool   stalled=false;
   int    ne[2],ispin,nion;
   double E[60],total_energy,deltae,deltae_old,deltac;

   int it_in   = control.loop(0);
   int it_out  = control.loop(1);
   double tole = control.tolerances(0);
   double tolc = control.tolerances(1);

   double dt = control.time_step();
   double dte = dt /sqrt(control.fake_mass());

   int minimizer = control.minimizer();

   for (auto ii=0; ii<60; ++ii) E[ii] = 0.0;

   if (parall->is_master())
   {
      if (minimizer==1) std::cout << "     ============ Grassmann conjugate gradient iteration ============" << std::endl;
      if (minimizer==2) std::cout << "     ================== Grassmann lmbfgs iteration ==================" << std::endl;
      if (minimizer==4) std::cout << "     ============= Stiefel conjugate gradient iteration =============" << std::endl;
      if (minimizer==5) std::cout << "     ============= Kohn-Sham scf iteration (potential) ==============" << std::endl;
      if (minimizer==7) std::cout << "     =================== Stiefel lmbfgs iteration ===================" << std::endl;
      if (minimizer==8) std::cout << "     ============== Kohn-Sham scf iteration (density) ===============" << std::endl;

      std::cout << "          >>> iteration started at " << util_date() << "  <<<\n";
      std::cout << "     iter.                 Energy          DeltaE        DeltaRho\n";
      std::cout << "     ----------------------------------------------------------------\n";
      //printf("%10d%25.12le%16.6le%16.6le\n",1000,99.99, 1.33434340e-4, 2.33434211e-6);


   }

   //if (minimizer > 1) pspw_Grsm_list_start()
   if ((minimizer==5) || (minimizer==8)) it_out = 1;

   Geodesic mygeodesic(minimizer,&mymolecule);


   /* generate phase factors and local psp and semicore density */
   mymolecule.phafacs_vl_potential_semicore();

   if (mymolecule.newpsi) 
   {
      for (int it=0; it<it_in; ++it) mymolecule.sd_update(dte);
      if (parall->is_master()) std::cout << "        - " << it_in << " steepest descent iterations performed" << std::endl;
   }

   //std::cout << "cgsd_energy: minimizer = " << minimizer << std::endl;
   deltae = -1.0e-03;
   int bfgscount = 0;
   int icount = 0;
   bool converged = false;
   while ((icount < it_out) && (!converged))
   {
      ++icount;
      if (stalled)
      {
         for (int it=0; it<it_in; ++it) mymolecule.sd_update(dte);
         if (parall->is_master()) std::cout << "        - " << it_in << " steepest descent iterations performed" << std::endl;
         bfgscount = 0;
      }

      deltae_old = deltae;
      if (minimizer==1)
      {
        ++bfgscount; total_energy = cgsd_cgminimize(mymolecule,mygeodesic,E,&deltae,&deltac,bfgscount,it_in,tole,tolc);
      }
      else if (minimizer==2)
      {
        ++bfgscount; // call bfgsminimize(E,deltae,deltac,bfgscount,it_in)
      }
      else if (minimizer==4)
      {
        ++bfgscount; // call cgminimize2(E,deltae,deltac,bfgscount,it_in)
      }
      else if (minimizer==7)
      {
        ++bfgscount; // call bfgsminimize2(E,deltae,deltac,bfgscount,it_in)
      }

      if (parall->is_master())
         printf("%10d%25.12le%16.6le%16.6le\n",icount*it_in,total_energy,deltae, deltac);

      if ((fabs(deltae)>fabs(deltae_old)) || (fabs(deltae)>1.0e-2) || (deltae>0.0))
         stalled = true;
      else
         stalled = false;

      converged = (fabs(deltae)<tole) && (deltac<tolc);


   }
   if (parall->is_master())
   {
      if (converged)  std::cout <<  "     *** tolerance ok. iteration terminated" << std::endl;
      if (!converged) std::cout <<  "     *** arrived at the Maximum iteration.  terminated" << std::endl;
      std::cout << "          >>> iteration ended at   " << util_date() << "  <<<\n";
   }



   /* report summary of results */
   //total_energy  = mymolecule.gen_all_energies();
   if (parall->is_master())
   {
      std::cout << std::endl;
      std::cout << mymolecule;
   }

   return total_energy;
}



/******************************************
 *                                        *
 *           cgsd_energy_gradient         *
 *                                        *
 ******************************************/
void cgsd_energy_gradient(Molecule& mymolecule, double *grad_ion)
{
   mymolecule.psi_1local_force(grad_ion);
   mymolecule.psi_1nonlocal_force(grad_ion);
   mymolecule.ewald_fion(grad_ion);
   mymolecule.semicore_force(grad_ion);
   mymolecule.psi_1apc_force(grad_ion);

}

}


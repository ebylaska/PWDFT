#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "Control2.hpp"
#include "Geodesic12.hpp"
#include "Ion.hpp"
#include "Molecule.hpp"
#include "Parallel.hpp"
#include "Pneb.hpp"
#include "iofmt.hpp"
#include "pspw_lmbfgs.hpp"
#include "pspw_lmbfgs2.hpp"
#include "util_date.hpp"
#include "nwpw_scf_mixing.hpp"

#include "cgsd.hpp"

#include "iofmt.hpp"

namespace pwdft {

/******************************************
 *                                        *
 *            cgsd_steepest_update        *
 *                                        *
 ******************************************/
void static cgsd_steepest_update(Molecule &mymolecule, const bool oprint, const bool nolagrange, 
                                 const int it_in, const double dte, std::ostream &coutput)
{
   if (nolagrange)
   {
      for (int it=0; it<it_in; ++it)
         mymolecule.sd_update2(dte);
      if (oprint)
         coutput << "        - " << it_in << " steepest descent2 iterations performed" << std::endl;
   }
   else
   {
      for (int it=0; it<it_in; ++it)
         mymolecule.sd_update(dte);
      if (oprint) 
         coutput << "        - " << it_in << " steepest descent iterations performed" << std::endl;
   }
}

/******************************************
 *                                        *
 *            cgsd_noit_energy            *
 *                                        *
 ******************************************/
double cgsd_noit_energy(Molecule &mymolecule, bool doprint,
                        std::ostream &coutput) {
  Parallel *parall = mymolecule.mygrid->d3db::parall;

  /* generate phase factors and local psp and semicore density */
  mymolecule.phafacs_vl_potential_semicore();

  double total_energy = mymolecule.gen_all_energies();

  /* report summary of results */
  if (parall->base_stdio_print && doprint) 
  {
     coutput << "     ==================  optimization turned off  ===================" << std::endl << std::endl;
     coutput << std::endl;
     coutput << mymolecule;
  }

  return total_energy;
}

/******************************************
 *                                        *
 *              cgsd_energy               *
 *                                        *
 ******************************************/
double cgsd_energy(Control2 &control, Molecule &mymolecule, bool doprint, std::ostream &coutput) 
{
   Parallel *parall = mymolecule.mygrid->d3db::parall;
   Pneb *mygrid = mymolecule.mygrid;
   Ion *myion = mymolecule.myion;
 
   bool stalled = false;
   int ne[2],ispin,nion;
   double E[70],total_energy,deltae,deltae_old,deltac;
 
   bool nolagrange = control.nolagrange();
   int it_in   = control.loop(0);
   int it_out  = control.loop(1);
   double tole = control.tolerances(0);
   double tolc = control.tolerances(1);
   double scf_alpha = control.scf_alpha();
   double scf_beta  = control.scf_beta();
   double kerker_g0 = control.kerker_g0();
   int diis_histories = control.diis_histories();
   int scf_algorithm = control.scf_algorithm();

   double dt = control.time_step();
   double dte = dt/sqrt(control.fake_mass());
 
   int minimizer   = control.minimizer();
   int lmbfgs_size = control.lmbfgs_size();
 
   bool hprint = (parall->is_master() && control.print_level("high") && doprint);
   bool oprint = (parall->is_master() && control.print_level("medium") && doprint);
   bool lprint = (parall->is_master() && control.print_level("low") && doprint);
 
   for (auto ii=0; ii<70; ++ii)
      E[ii] = 0.0;
 
   if (oprint) 
   {
      if (minimizer == 1) coutput << "     =========== Grassmann conjugate gradient iteration ===========" << std::endl;
      if (minimizer == 2) coutput << "     ================= Grassmann lmbfgs iteration =================" << std::endl;
      if (minimizer == 4) coutput << "     ============ Stiefel conjugate gradient iteration ============" << std::endl;
      if (minimizer == 5) coutput << "     ============ Kohn-Sham scf iteration (potential) =============" << std::endl;
      if (minimizer == 7) coutput << "     ================== Stiefel lmbfgs iteration ==================" << std::endl;
      if (minimizer == 8) coutput << "     ============= Kohn-Sham scf iteration (density) ==============" << std::endl;
      if (minimizer == 9) coutput << "     =========== Grassman cg (Stich linesearch) iteration =========" << std::endl;
      if (minimizer == 10) coutput<< "     ============ Grassman lmbfgs (Stich l.s.) iteration ==========" << std::endl;
     
      coutput << "          >>> iteration started at " << util_date() << "  <<<" << std::endl;;
      coutput << "     iter.                   Energy          DeltaE        DeltaRho" << std::endl;
      coutput << "     --------------------------------------------------------------" << std::endl;
      // printf("%10d%25.12le%16.6le%16.6le\n",1000,99.99, 1.33434340e-4, 2.33434211e-6);
   }
 
   // if (minimizer > 1) pspw_Grsm_list_start()
   //if ((minimizer == 5) || (minimizer == 8))
   //   it_out = 1;
 
   Geodesic12 mygeodesic12(minimizer, &mymolecule, control);
 
   /* generate phase factors and local psp and semicore density */
   mymolecule.phafacs_vl_potential_semicore();
 
   // std::cout << "cgsd_energy: minimizer = " << minimizer << std::endl;
   deltae = -1.0e-03;
   int bfgscount = 0;
   int icount = 0;
   bool converged = false;
 
   if (minimizer == 1) {
      if (mymolecule.newpsi) 
      {
         //int ispin = mymolecule.mygrid->ispin;
         //double *vall =  mygrid->r_nalloc(ispin);
         //mymolecule.gen_vall();
         //mymolecule.get_vall(vall);
         //mymolecule.psi_minimize(vall, coutput);
         int it_in0 = 15;
         cgsd_steepest_update(mymolecule,oprint,nolagrange,it_in0,dte,coutput);
      }
      while ((icount < it_out) && (!converged)) 
      {
         ++icount;
         if (stalled) 
         {
            cgsd_steepest_update(mymolecule,oprint,nolagrange,it_in,dte,coutput);

            bfgscount = 0;
         }
         deltae_old = deltae;
         total_energy = cgsd_cgminimize(mymolecule,mygeodesic12.mygeodesic1,E,&deltae,
                                        &deltac,bfgscount,it_in,tole,tolc);
         ++bfgscount;
         if (oprint)
           coutput << Ifmt(10) << icount*it_in 
                   << Efmt(25,12) << total_energy
                   << Efmt(16,6) << deltae 
                   << Efmt(16,6) << deltac << std::endl;
         if ((std::fabs(deltae) > fabs(deltae_old)) ||
             (std::fabs(deltae) > 1.0e-2) || (deltae > 0.0))
            stalled = true;
         else
            stalled = false;
         converged = (std::fabs(deltae) < tole) && (deltac < tolc);
      }
 
   } else if (minimizer == 2) {
      if (mymolecule.newpsi) {
         int it_in0 = 15;
         cgsd_steepest_update(mymolecule,oprint,nolagrange,it_in0,dte,coutput);
      }
      pspw_lmbfgs psi_lmbfgs(mygeodesic12.mygeodesic1, lmbfgs_size);
      while ((icount < it_out) && (!converged)) {
         ++icount;
         if (stalled) {
            cgsd_steepest_update(mymolecule,oprint,nolagrange,it_in,dte,coutput);
            bfgscount = 0;
         }
         deltae_old = deltae;
         total_energy = cgsd_bfgsminimize(mymolecule, mygeodesic12.mygeodesic1, psi_lmbfgs, E,
                               &deltae, &deltac, bfgscount, it_in, tole, tolc);
         ++bfgscount;
         if (oprint)
            coutput << Ifmt(10) << icount * it_in 
                    << Efmt(25,12) << total_energy
                    << Efmt(16,6) << deltae 
                    << Efmt(16,6) << deltac << std::endl;
         if ((std::fabs(deltae) > fabs(deltae_old)) ||
             (std::fabs(deltae) > 1.0e-2) || (deltae > 0.0))
            stalled = true;
         else
            stalled = false;
         converged = (std::fabs(deltae) < tole) && (deltac < tolc);
      }
   } else if (minimizer == 4) {
      if (mymolecule.newpsi) {
        int it_in0 = 15;
        for (int it = 0; it < it_in0; ++it)
          mymolecule.sd_update_sic(dte);
        if (oprint)
          coutput << "        - " << it_in0
                  << " steepest descent iterations performed" << std::endl;
      }
      while ((icount < it_out) && (!converged)) {
        ++icount;
        if (stalled) {
          for (int it = 0; it < it_in; ++it)
            mymolecule.sd_update_sic(dte);
          if (oprint)
            coutput << "        - " << it_in
                      << " steepest descent iterations performed" << std::endl;
          bfgscount = 0;
        }
        deltae_old = deltae;
        total_energy =
            cgsd_cgminimize2(mymolecule, mygeodesic12.mygeodesic2, E, &deltae,
                             &deltac, bfgscount, it_in, tole, tolc);
        ++bfgscount;
        if (oprint)
          coutput << Ifmt(10) << icount*it_in 
                  << Efmt(25,12) << total_energy
                  << Efmt(16,6) << deltae 
                  << Efmt(16,6) << deltac << std::endl;
        if ((std::fabs(deltae) > fabs(deltae_old)) ||
            (std::fabs(deltae) > 1.0e-2) || (deltae > 0.0))
          stalled = true;
        else
          stalled = false;
        converged = (std::fabs(deltae) < tole) && (deltac < tolc);
      }
   } else if (minimizer == 5) {
      if (mymolecule.newpsi) {
         int it_in0 = 15;
         cgsd_steepest_update(mymolecule,oprint,nolagrange,it_in0,dte,coutput);
      }
      while ((icount < it_out) && (!converged)) 
      {
         ++icount;
         if (stalled) 
         {
            cgsd_steepest_update(mymolecule,oprint,nolagrange,it_in,dte,coutput);

            bfgscount = 0;
         }
         deltae_old = deltae;
         total_energy = cgsd_bybminimize(mymolecule,mygeodesic12.mygeodesic1,E,&deltae,
                                         &deltac,bfgscount,it_in, 
                                         scf_algorithm,scf_alpha,kerker_g0,diis_histories,
                                         tole,tolc);
        ++bfgscount;
        if (oprint)
          coutput << Ifmt(10) << icount*it_in
                  << Efmt(25,12) << total_energy
                  << Efmt(16,6) << deltae
                  << Efmt(16,6) << deltac << std::endl;
        if ((std::fabs(deltae) > fabs(deltae_old)) ||
            (std::fabs(deltae) > 1.0e-2) || (deltae > 0.0))
           stalled = true;
        else
           stalled = false;
        stalled = false;
        converged = (std::fabs(deltae) < tole) && (deltac < tolc);
      }

   } else if (minimizer == 7) {
      if (mymolecule.newpsi) {
         int it_in0 = 15;
         for (int it=0; it<it_in0; ++it)
            mymolecule.sd_update_sic(dte);
         if (oprint)
            coutput << "        - " << it_in0 << " steepest descent iterations performed" << std::endl;
      }
      pspw_lmbfgs2 psi_lmbfgs2(mygeodesic12.mygeodesic2, lmbfgs_size);
      while ((icount < it_out) && (!converged)) {
         ++icount;
         if (stalled) {
            for (int it=0; it<it_in; ++it)
               mymolecule.sd_update_sic(dte);
            if (oprint)
               coutput << "        - " << it_in << " steepest descent iterations performed" << std::endl;
            bfgscount = 0;
         }
         deltae_old = deltae;
         total_energy = cgsd_bfgsminimize2(mymolecule,mygeodesic12.mygeodesic2,psi_lmbfgs2,
                                           E,&deltae,&deltac,bfgscount,it_in,tole,tolc);
         ++bfgscount;
         if (oprint)
            coutput << Ifmt(10) << icount * it_in 
                    << Efmt(25,12) << total_energy
                    << Efmt(16,6) << deltae 
                    << Efmt(16,6) << deltac << std::endl;
         if ((std::fabs(deltae) > fabs(deltae_old)) ||
             (std::fabs(deltae) > 1.0e-2) || (deltae > 0.0))
            stalled = true;
         else
            stalled = false;
         converged = (std::fabs(deltae) < tole) && (deltac < tolc);
      }
   } else if (minimizer == 8) {

      if (mymolecule.newpsi) {

         //total_energy = cgsd_bybminimize0(mymolecule,mygeodesic12.mygeodesic2,psi_lmbfgs2,
         //                                 E,&deltae,&deltac,bfgscount,it_in,tole,tolc);
         int it_in0 = 15;
         cgsd_steepest_update(mymolecule,oprint,nolagrange,it_in0,dte,coutput);
         //double total_energy0 = mymolecule.energy();
         //total_energy0 = cgsd_bybminimize2(mymolecule,mygeodesic12.mygeodesic1,E,&deltae,
         //                                 &deltac,0,5,1,tole,tolc);
      }
      
      double total_energy0 = mymolecule.energy();
      double x,sumxx = 0.0;
      int n2ft3d = mygrid->n2ft3d;
      int ispin = mygrid->ispin;
      double omega = mygrid->lattice->omega();
      double scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
      double dv = omega * scal1;
      for (int i=0; i < n2ft3d; ++i) 
      {
         x = (mymolecule.rho1[i]);
         x += (mymolecule.rho1[i+(ispin-1)*n2ft3d]);
         sumxx += x;
      }
      double sum0  = mygrid->d3db::parall->SumAll(1, sumxx) * dv;
      //std::cout << "total energy0=" << total_energy0 << std::endl;
      //std::cout << "total sum0=" << sum0 << std::endl;

      double scf_error;

      nwpw_scf_mixing scfmix(mygrid,kerker_g0,
                             scf_algorithm,scf_alpha,scf_beta,diis_histories,
                             mygrid->ispin,mygrid->n2ft3d,mymolecule.rho1);

      while ((icount < (it_out*it_in)) && (!converged))
      {
         ++icount;
         if ((stalled) && (icount%it_in==0))
         {
            cgsd_steepest_update(mymolecule,oprint,nolagrange,it_in,dte,coutput);
            scfmix.reset_mix(mymolecule.rho1);

            bfgscount = 0;
            stalled = false;
         }
         int ks_it_in = control.ks_maxit_orb();
         int ks_it_out =  control.ks_maxit_orbs();
         deltae_old = deltae;

         int nsize = mygrid->ispin*mygrid->n2ft3d;
         std::memcpy(mymolecule.rho2,mymolecule.rho1,nsize*sizeof(double));

         total_energy = cgsd_bybminimize2(mymolecule,mygeodesic12.mygeodesic1,E,&deltae,
                                         &deltac,bfgscount,ks_it_in,ks_it_out,
                                         tole,tolc);
         double total_energy2 = mymolecule.energy();
         scfmix.mix(mymolecule.rho1,mymolecule.rho1,deltae,&scf_error);
         deltac = scf_error;
         deltae = total_energy2 - total_energy0;
         total_energy0 = total_energy2;
         ++bfgscount;

         converged = (std::fabs(deltae) < tole) && (deltac < tolc);

         if ((oprint) && ((icount%it_in==0) || converged))
         {
            coutput << Ifmt(10)    << icount
                    << Efmt(25,12) << total_energy
                    << Efmt(16,6)  << deltae 
                    << Efmt(16,6)  << deltac << std::endl;
         }

         if ((std::fabs(deltae) > fabs(deltae_old)) ||
             ((deltac > 1.0e-3) && (std::fabs(deltae)<1.0e-10)) ||
             (std::fabs(deltae) > 1.0e-2) || (deltae > 0.0))
             stalled = true;

      }    
   } else if (minimizer == 9) {
   } else if (minimizer == 10) {
   }

 
   if (oprint) {
      if (converged)  coutput << "     *** tolerance ok. iteration terminated" << std::endl;
      if (!converged) coutput << "     *** arrived at the Maximum iteration.  terminated" << std::endl;
      coutput << "          >>> iteration ended at   " << util_date() << "  <<<" << std::endl;
   }
 
   /* report summary of results */
   // total_energy  = mymolecule.gen_all_energies();
   if (oprint) {
      coutput << std::endl;
      coutput << mymolecule;
   }
 
   return total_energy;
}

/******************************************
 *                                        *
 *           cgsd_energy_gradient         *
 *                                        *
 ******************************************/
void cgsd_energy_gradient(Molecule &mymolecule, double *grad_ion) 
{

   mymolecule.psi_1local_force(grad_ion);
   mymolecule.psi_1nonlocal_force(grad_ion);
   mymolecule.semicore_force(grad_ion);
 
   mymolecule.ion_fion(grad_ion);
   mymolecule.psi_1apc_force(grad_ion);
   mymolecule.psi_1apc_force(grad_ion);
}

} // namespace pwdft

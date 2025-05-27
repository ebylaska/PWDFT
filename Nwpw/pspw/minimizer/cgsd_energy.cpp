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
/**
 * @brief Main driver for energy minimization.
 *
 * This function coordinates the energy minimization process by selecting and executing
 * various minimization strategies (e.g., conjugate gradient, L-BFGS, steepest descent, 
 * and band-by-band methods) to iteratively reduce the total energy of the system.
 *
 * The algorithm monitors convergence through the change in energy (ΔE) and density residuals (Δρ),
 * and can automatically switch minimization strategies (for example, from a Stiefel conjugate gradient
 * update (minimizer=4) to a band-by-band update (minimizer=8)) when stalling is detected.
 *
 * Depending on the control parameters, the function chooses among several minimizer modes:
 *  - Minimizer 1, 2, 9, 10: Grassmann-based updates (with symmetrized Hamiltonian)
 *  - Minimizer 4, 7: Stiefel-based conjugate gradient or L-BFGS updates (used when fractional 
 *                    occupations are present, where orbital details matter)
 *  - Minimizer 5, 8: Band-by-band (density-based) updates that can be activated if full-space 
 *                    updates stall.
 *
 * The function updates the state of the provided molecule object (including orbitals, densities,
 * and SCF-related parameters), performs orbital rotations and reorthonormalization as necessary,
 * and returns the converged total energy.
 *
 * @param control    The Control2 object containing simulation parameters such as iteration limits,
 *                   convergence tolerances, SCF mixing parameters, and minimizer selection flags.
 * @param mymolecule The Molecule object representing the current state of the system (orbitals,
 *                   densities, grid information, etc.). This object is updated during minimization.
 * @param doprint    A boolean flag indicating whether detailed progress and diagnostic output
 *                   should be printed.
 * @param coutput    An output stream (std::ostream) used for logging progress, diagnostic information,
 *                   and final results.
 *
 * @return The converged total energy of the system.
 *
 * @note Depending on the convergence behavior, the algorithm may switch automatically from
 *       minimizer=4 (Stiefel conjugate gradient) to minimizer=8 (band-by-band update) to overcome
 *       stalling.
 *
 * @author Eric J. Bylaska
 * @date 2/2/3025 
 */
double cgsd_energy(Control2 &control, Molecule &mymolecule, bool doprint, std::ostream &coutput) 
{
   Parallel *parall = mymolecule.mygrid->d3db::parall;
   Pneb *mygrid = mymolecule.mygrid;
   Ion *myion = mymolecule.myion;
 
   bool stalled = false;
   int ne[2],ispin,nion;
   double E[70],total_energy,deltae,deltae_old,deltac;
 
   bool nolagrange = control.nolagrange() ||  control.fractional();
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

   bool extra_rotate = control.scf_extra_rotate();
 
   for (auto ii=0; ii<70; ++ii)
      E[ii] = 0.0;
 
   if (oprint) 
   {
      if (minimizer == 1) coutput << "     =========== Grassmann conjugate gradient iteration ===========" << std::endl;
      if (minimizer == 2) coutput << "     ================= Grassmann lmbfgs iteration =================" << std::endl;
      if (minimizer == 3) coutput << "     ======== Kohn-Sham  scf iteration (Grassman) iteration =======" << std::endl;
      if (minimizer == 4) coutput << "     ============ Stiefel conjugate gradient iteration ============" << std::endl;
      if (minimizer == 5) coutput << "     ============ Kohn-Sham scf iteration (potential) =============" << std::endl;
      if (minimizer == 6) coutput << "     ========== Kohn-Sham scf iteration (lmbfgs) iteration ========" << std::endl;
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
 
   if (minimizer == 1) 
   {
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
   } 

   else if (minimizer == 2) 
   {
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
   } 

   else if (minimizer == 3) 
   {
      if (mymolecule.newpsi) {
         int it_in0 = 15;
         for (int it=0; it<it_in0; ++it)
            mymolecule.sd_update(dte);
         if (oprint) coutput << "        - " << it_in0 << " steepest descent iterations performed" << std::endl;
      }

      // Initial SCF setup
      // generate density and rotate orbitals 
      // Generate initial density, potential and then orbital diagonalization
      double total_energy0 = mymolecule.energy(); // Run Hψ = Eψ and compute E[0]
      mymolecule.gen_hml();                       // Generate ⟨ψ|H|ψ⟩ (stored in hml)
      mymolecule.diagonalize();      // Diagonalize H matrix (sets eig)
      mymolecule.rotate1to2();       // Rotate ψ₁ to ψ₂ using eigenvectors
      mymolecule.swap_psi1_psi2();   // Swap ψ₂ → ψ₁ (start clean state)

      // Normalize total density (diagnostic)
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

      // Setup SCF loop
      deltae = -1.0e-03;
      int bfgscount = 0;
      int icount = 0;
      int ks_it_in  = control.ks_maxit_orb();
      int ks_it_out = control.ks_maxit_orbs();
      double scf_error = 0.0;
      double *vout = mymolecule.rho1;
      double *vnew = mymolecule.rho2;

      nwpw_scf_mixing scfmix(mygrid,kerker_g0,
                             scf_algorithm,scf_alpha,scf_beta,diis_histories,
                             mygrid->ispin,mygrid->n2ft3d,vout);

      while ((icount < (it_out*it_in)) && (!converged))
      {
         ++icount;
         if (stalled)
         {
            for (int it=0; it<it_in; ++it)
               mymolecule.sd_update(dte);
            if (oprint)
               coutput << "        - " << it_in << " steepest descent iterations performed" << std::endl;
            bfgscount = 0;
         }
         deltae_old = deltae;

         // minimize ks orbitals it_in steps with fixed ks potential
         double total_energy_fixedV = cgsd_cgksminimize(mymolecule,mygeodesic12.mygeodesic1,E,&deltae,
                                                        &deltac,bfgscount,it_in,tole,tolc);

         // diagonalize psi wrt current ks potential
         mymolecule.gen_hml();
         mymolecule.diagonalize();
         mymolecule.rotate1to2();
         mymolecule.swap_psi1_psi2();

         // Generate updated density from current ψ and occupations, generate ks potential 
         // and then calculate energy, after this rho1, dng1, rho1_all, 
         // and ks potentials and hpsi have been updated
         total_energy = mymolecule.energy();

         // [Insert fractional occupation update here if needed]
         // if (mysolid.fractional) update_occupations(...);


         //std::cout << "total_energy=" << total_energy << " " << total_energy2 
         //          << " " << total_energy2 - total_energy << std::endl;

         //define fractional occupation here
         scfmix.mix(vout,vnew,deltae,&scf_error);
         std::memcpy(vout,vnew,ispin*n2ft3d*sizeof(double));

        // density now updated, and the ks potentials have to be updated here!?

         deltac = scf_error;

         converged = (std::fabs(deltae) < tole) && (deltac < tolc);
         //deltac = mysolid.rho_error();
         deltae = total_energy - total_energy0;
         total_energy0 = total_energy;
         ++bfgscount;

         converged = (std::fabs(deltae) < tole) && (deltac < tolc);

         if ((oprint) && ((icount%it_in==0) || converged))
         {
            coutput << Ifmt(10)    << icount
                    << Efmt(25,12) << total_energy
                    << Efmt(16,6)  << deltae
                    << Efmt(16,6)  << deltac
                    << Efmt(16,6)  << total_energy-total_energy_fixedV << std::endl;
         }

         // Finalize SCF step with updated potentials
         mymolecule.gen_scf_potentials_from_rho1();

      }
   }

   else if (minimizer == 4) 
   {
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
   }

   else if (minimizer == 5) 
   {
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
   } 

   else if (minimizer == 7) 
   {
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

   } 

   else if (minimizer == 8) 
   {
      if (mymolecule.newpsi) {

         //total_energy = cgsd_bybminimize0(mymolecule,mygeodesic12.mygeodesic2,psi_lmbfgs2,
         //                                 E,&deltae,&deltac,bfgscount,it_in,tole,tolc);
         int it_in0 = 15;
         cgsd_steepest_update(mymolecule,oprint,nolagrange,it_in0,dte,coutput);
         //double total_energy0 = mymolecule.energy();
         //total_energy0 = cgsd_bybminimize2(mymolecule,mygeodesic12.mygeodesic1,E,&deltae,
         //                                 &deltac,0,5,1,tole,tolc);
      }
      
      // rotate orbitals 
      double total_energy0 = mymolecule.energy();
      mymolecule.gen_hml();
      mymolecule.diagonalize();
      mymolecule.rotate1to2();
      mymolecule.swap_psi1_psi2();

      //define fractional occupation here

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

      int ks_it_in  = control.ks_maxit_orb();
      int ks_it_out = control.ks_maxit_orbs();
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

            // rotate orbitals
            if (extra_rotate) 
            {
               mymolecule.gen_hml();
               mymolecule.diagonalize();
               mymolecule.rotate1to2();
               mymolecule.swap_psi1_psi2();
            }

            //define fractional occupation here

            scfmix.reset_mix(mymolecule.rho1);

            bfgscount = 0;
            stalled = false;
         }
         deltae_old    = deltae;

         int nsize = mygrid->ispin*mygrid->n2ft3d;
         std::memcpy(mymolecule.rho2,mymolecule.rho1,nsize*sizeof(double));

         total_energy = cgsd_bybminimize2(mymolecule,mygeodesic12.mygeodesic1,E,&deltae,
                                          &deltac,bfgscount,ks_it_in,ks_it_out,
                                          tole,tolc);
         // rotate orbitals
         if (extra_rotate) 
         {
            total_energy = mymolecule.energy();
            mymolecule.gen_hml();
            mymolecule.diagonalize();
            mymolecule.rotate1to2();
            mymolecule.swap_psi1_psi2();
         }

         //std::cout << "total_energy=" << total_energy << " " << total_energy2 
         //          << " " << total_energy2 - total_energy << std::endl;

         //define fractional occupation here

         scfmix.mix(mymolecule.rho1,mymolecule.rho1,deltae,&scf_error);
         deltac = scf_error;
         deltae = total_energy - total_energy0;
         total_energy0 = total_energy;
         ++bfgscount;

         converged = (std::fabs(deltae) < tole) && (deltac < tolc);

         if ((oprint) && ((icount%it_in==0) || converged))
         {
            coutput << Ifmt(10)    << icount
                    << Efmt(25,12) << total_energy
                    << Efmt(16,6)  << deltae 
                    << Efmt(16,6)  << deltac << std::endl;
         }

         stalled = ((std::fabs(deltae) > fabs(deltae_old)) ||
                   ((deltac > 1.0e-3) && (std::fabs(deltae)<1.0e-10)) ||
                   (std::fabs(deltae) > 1.0e-2) || (deltae > 0.0));

        // if ((std::fabs(deltae) > fabs(deltae_old)) ||
         //    ((deltac > 1.0e-3) && (std::fabs(deltae)<1.0e-10)) ||
         //    (std::fabs(deltae) > 1.0e-2) || (deltae > 0.0))
         //    stalled = true;
      }    

   } 

   else if (minimizer == 9) 
   {

   } 

   else if (minimizer == 10) 
   {
   }

 
   if (oprint) {
      if (converged)  coutput << "     *** tolerance ok. iteration terminated" << std::endl;
      if (!converged) coutput << "     *** arrived at the Maximum iteration.  terminated" << std::endl;
      coutput << "          >>> iteration ended at   " << util_date() << "  <<<" << std::endl;
   }
 
   /* report summary of results */
   total_energy  = mymolecule.gen_all_energies();
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

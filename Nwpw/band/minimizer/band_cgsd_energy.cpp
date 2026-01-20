#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "Control2.hpp"
#include "band_Geodesic12.hpp"
#include "Ion.hpp"
#include "Solid.hpp"
#include "Parallel.hpp"
#include "Cneb.hpp"
#include "iofmt.hpp"
#include "band_lmbfgs.hpp"
#include "band_lmbfgs2.hpp"
#include "util_date.hpp"
#include "nwpw_cscf_mixing.hpp"

#include "band_cgsd.hpp"

#include "iofmt.hpp"

namespace pwdft {

/******************************************
 *                                        *
 *            band_cgsd_noit_energy       *
 *                                        *
 ******************************************/
double band_cgsd_noit_energy(Solid &mysolid, bool doprint, std::ostream &coutput) 
{
   Parallel *parall = mysolid.mygrid->c3db::parall;
 
   /* generate phase factors and local psp and semicore density */
   mysolid.phafacs_vl_potential_semicore();
 
   double total_energy = mysolid.gen_all_energies();
 
   /* report summary of results */
   if (parall->base_stdio_print && doprint) 
   {
      coutput << "     ==================  optimization turned off  ===================" << std::endl << std::endl;
      coutput << std::endl;
      coutput << mysolid;
   }
   coutput << mysolid.print_filled(parall->base_stdio_print && doprint);
 
   return total_energy;
}

/******************************************
 *                                        *
 *           band_cgsd_energy             *
 *                                        *
 ******************************************/
double band_cgsd_energy(Control2 &control, Solid &mysolid, bool doprint, std::ostream &coutput) 
{
   Parallel *parall = mysolid.mygrid->c3db::parall;
   Cneb *mygrid = mysolid.mygrid;
   Ion *myion = mysolid.myion;
 
   bool stalled = false;
   int ne[2],ispin,nion;
   double E[80],total_energy,deltae,deltae_old,deltac,deltasmear,deltasmear_old;
 
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

   const double occ_rmsd_abs_thresh   = 1.0e-3;
   const double occ_rmsd_delta_thresh = 5.0e-5;

   const double occ_rmsd_thresh   = 1.0e-3;   // already implicit
   const double occ_delta_thresh  = 5.0e-4;   // oscillation detector

   const int    occ_update_stride   = 5;     // update occupations every 5 outer steps
   const double drho_for_occ_update = 5e-3;  // or when density change is already small



 
   for (auto ii=0; ii<80; ++ii)
      E[ii] = 0.0;
 
   if (oprint) 
   {
      if (minimizer == 1) coutput << "     ======= bundled Grassmann conjugate gradient iteration =======" << std::endl;
      if (minimizer == 2) coutput << "     =========+=== bundled Grassmann lmbfgs iteration =============" << std::endl;
      if (minimizer == 3) coutput << "     ======== Kohn-Sham scf iteration (Grassmann) iteration =======" << std::endl;
      if (minimizer == 4) coutput << "     ============ Stiefel conjugate gradient iteration ============" << std::endl;
      if (minimizer == 5) coutput << "     ============ Kohn-Sham scf iteration (potential) =============" << std::endl;
      if (minimizer == 6) coutput << "     ========== Kohn-Sham scf iteration (lmbfgs) iteration ========" << std::endl;
      if (minimizer == 7) coutput << "     ================== Stiefel lmbfgs iteration ==================" << std::endl;
      if (minimizer == 8) coutput << "     ============= Kohn-Sham scf iteration (density) ==============" << std::endl;
     
      coutput << "          >>> iteration started at " << util_date() << "  <<<" << std::endl;;
      if ((mysolid.fractional) &&  (!mysolid.fractional_frozen))
      {
         if (minimizer==3)
         {
            coutput << "     iter.              free energy                   energy    smear energy           EF_up    delta free_E    delta energy   delta smear_E       delta rho" << std::endl;
            coutput << "     -------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
         }
         else
         {
            coutput << "     iter.              free energy                   energy    smear energy           EF_up    delta free_E    delta energy   delta smear_E       delta rho" << std::endl;
            coutput << "     -------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
         }
      }
      else
      {
         if (minimizer==3)
         {
            coutput << "     iter.                   energy    delta energy       delta rho       delta scf" << std::endl;
            coutput << "     ------------------------------------------------------------------------------" << std::endl;
         }
         else
         {
            coutput << "     iter.                   energy    delta energy       delta rho" << std::endl;
            coutput << "     --------------------------------------------------------------" << std::endl;
         }
      }
      // printf("%10d%25.12le%16.6le%16.6le\n",1000,99.99, 1.33434340e-4, 2.33434211e-6);

   }
 
   // if (minimizer > 1) band_Grsm_list_start()
   //if ((minimizer == 5) || (minimizer == 8))
   //   it_out = 1;
 
   band_Geodesic12 mygeodesic12(minimizer, &mysolid, control);
 
   /* generate phase factors and local psp and semicore density */
   mysolid.phafacs_vl_potential_semicore();
 
   // std::cout << "band_cgsd_energy: minimizer = " << minimizer << std::endl;
   deltae     = -1.0e-03;
   deltasmear = 1.0e-03;
   int bfgscount = 0;
   int icount = 0;
   bool converged = false;
 
   if (minimizer == 1) {
      if (mysolid.newpsi) 
      {
         int it_in0 = 15;
         for (int it=0; it<it_in0; ++it)
            mysolid.sd_update(dte);
         if (oprint) coutput << "        - " << it_in0 << " steepest descent iterations performed" << std::endl;
      }
      double ee=mysolid.energy();
      while ((icount < it_out) && (!converged)) 
      {
         ++icount;
         if (stalled) 
         {
            for (int it=0; it<it_in; ++it)
               mysolid.sd_update(dte);
            if (oprint)
               coutput << "        - " << it_in << " steepest descent iterations performed" << std::endl;
            bfgscount = 0;
         }
         deltae_old = deltae;
         total_energy = band_cgsd_cgminimize(mysolid,mygeodesic12.mygeodesic1,E,&deltae,
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
      if (mysolid.newpsi) {
         int it_in0 = 15;
         for (int it=0; it<it_in0; ++it)
            mysolid.sd_update(dte);
         if (oprint) coutput << "        - " << it_in0 << " steepest descent iterations performed" << std::endl;
      }
      band_lmbfgs psi_lmbfgs(mygeodesic12.mygeodesic1, lmbfgs_size);
      while ((icount < it_out) && (!converged)) {
         ++icount;
         if (stalled) {
            for (int it = 0; it < it_in; ++it)
               mysolid.sd_update(dte);
            if (oprint) coutput << "        - " << it_in << " steepest descent iterations performed" << std::endl;
            bfgscount = 0;
         }
         deltae_old = deltae;
         total_energy = band_cgsd_bfgsminimize(mysolid, mygeodesic12.mygeodesic1, psi_lmbfgs, E,
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


   } else if (minimizer == 3) {
      if (mysolid.newpsi)
      {
         int it_in0 = 15;
         for (int it=0; it<it_in0; ++it)
            mysolid.sd_update(dte);
         if (oprint) coutput << "        - " << it_in0 << " steepest descent iterations performed" << std::endl;
      }



      // Initial SCF setup
      // generate density and rotate orbitals 
      // Generate initial density, potential and then orbital diagonalization
      //double total_energy0 = mysolid.energy(); // Run Hψ = Eψ and compute E[0]
      double total_energy0 = mysolid.gen_all_energies();

      mysolid.gen_hml();                       // Generate ⟨ψ|H|ψ⟩ (stored in hml)
      mysolid.diagonalize();      // Diagonalize H matrix (sets eig)
      mysolid.rotate1to2();       // Rotate ψ₁ to ψ₂ using eigenvectors
      mysolid.swap_psi1_psi2();   // Swap ψ₂ → ψ₁ (start clean state)

      // Normalize total density (diagnostic)
      double x,sumxx = 0.0;
      int nfft3d = mygrid->nfft3d;
      int ispin = mygrid->ispin;
      double omega = mygrid->lattice->omega();
      double scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
      double dv = omega * scal1;
      for (int i=0; i < nfft3d; ++i)
      {
         x = (mysolid.rho1[i]);
         x += (mysolid.rho1[i+(ispin-1)*nfft3d]);
         sumxx += x;
      }
      double sum0  = mygrid->c3db::parall->SumAll(1, sumxx) * dv;

      //std::cout << "total energy0=" << total_energy0 << std::endl;
      //std::cout << "total sum0=" << sum0 << std::endl;

      // Setup SCF loop
      //deltae = -1.0e-03;
      int bfgscount = 0;
      icount = 0;
      int ks_it_in  = control.ks_maxit_orb();
      int ks_it_out = control.ks_maxit_orbs();
      double scf_error = 0.0;
      double *vout = mysolid.rho1;
      double *vnew = mysolid.rho2;

      nwpw_cscf_mixing scfmix(mygrid,kerker_g0,
                             scf_algorithm,scf_alpha,scf_beta,diis_histories,
                             mygrid->ispin,mygrid->nfft3d,vout);


      while ((icount < (it_out*it_in)) && (!converged))
      {
         //++icount;
         icount += it_in;;
         if (stalled) 
         {
            for (int it=0; it<it_in; ++it)
               mysolid.sd_update(dte);
            if (oprint)
               coutput << "        - " << it_in << " steepest descent iterations performed" << std::endl;
            bfgscount = 0;
         }
         deltae_old = deltae;
         deltasmear_old = deltasmear;


         // minimize ks orbitals it_in steps with fixed ks potential
         double smear0 = mysolid.E[28];

         double total_energy_fixedV = band_cgsd_cgksminimize(mysolid,mygeodesic12.mygeodesic1,E,&deltae,
                                                             &deltac,bfgscount,it_in,tole,tolc);

         // diagonalize psi wrt current ks potential
         mysolid.gen_hml();
         mysolid.diagonalize();
         mysolid.rotate1to2();
         mysolid.swap_psi1_psi2();

         // Generate updated density from current ψ and occupations, generate ks potential 
         // and then calculate energy, after this rho1, dng1, rho1_all, 
         // and ks potentials and hpsi have been updated
         //total_energy = mysolid.energy();
         //if (icount == (it_out*it_in))
         //   mysolid.occupation_update = false;
           
         // fractional occupation are updated in gen_all_energies
         total_energy = mysolid.gen_all_energies();

         // compute deltas BEFORE mixing
         //deltae     = total_energy - total_energy0;
         //deltasmear = mysolid.E[28] - smear0;

         double alpha_current = scfmix.get_alpha();

         // adaptive SCF damping based on occupation stability 
         if ((mysolid.fractional) &&  (!mysolid.fractional_frozen))
         {
          
            //double occ_rmsd  = mysolid.rmsd_occupation;
            //double docc_rmsd = std::abs(mysolid.rmsd_occupation -
            //                            mysolid.rmsd_occupation_prev);
          
            bool occ_unstable = (mysolid.rmsd_occupation > occ_rmsd_thresh);
            bool occ_oscillating = (std::abs(mysolid.rmsd_occupation - mysolid.rmsd_occupation_prev) > occ_delta_thresh);

            // If occupations are unstable, damp density mixing
            //if (mysolid.rmsd_occupation > 1.0e-4)
            //if ( (occ_rmsd  > occ_rmsd_abs_thresh) ||
            //     (docc_rmsd > occ_rmsd_delta_thresh) )
            if (occ_unstable)
            {
               // hard reset
               alpha_current = std::max(0.005, 0.5 * alpha_current);
               scfmix.set_alpha(alpha_current);
               //scfmix.reset_history(mysolid.rho1);
               //scfmix.reset_history(mysolid.rho1);
            }
            else if (occ_oscillating)
            {
               // soft response: damp, but keep history
               alpha_current = std::max(0.01, 0.8 * alpha_current);
               scfmix.set_alpha(alpha_current);
            }
            else
            {
               // occupations stable → restore nominal SCF alpha
               //alpha_use = scf_alpha;
               alpha_current = std::min(scf_alpha, alpha_current * 1.2);
               scfmix.set_alpha(alpha_current);
            }

         }

         // define new KS potential here
         //Perform SCF mixing
         scfmix.mix(vout,vnew,deltae,&scf_error);

         // after scfmix.mix(vout,vnew,...)
         double dr2_local = 0.0;
         for (int i=0; i<ispin*nfft3d; ++i) 
         {
            double d = vnew[i] - vout[i];
            dr2_local += d*d;
         }
         double dr2 = mygrid->c3db::parall->SumAll(1, dr2_local);
         double drho = std::sqrt(dr2 / (double)(ispin*nfft3d));

         std::memcpy(vout,vnew,ispin*nfft3d*sizeof(double));
         
        // density now updated, and the ks potentials have to be updated here!?

         //deltac = scf_error;
         deltac = drho;


         //deltac = mysolid.rho_error();
         deltae = total_energy - total_energy0;
         deltasmear = mysolid.E[28] - smear0;
         total_energy0 = total_energy;
         ++bfgscount;

         bool occ_stable = (mysolid.rmsd_occupation < occ_rmsd_thresh);
         bool occ_reasonable = true;
         converged = (std::fabs(deltae) < tole) && (deltac < tolc);

         if (mysolid.fractional && !mysolid.fractional_frozen)
         {
            //occ_reasonable = (mysolid.rmsd_occupation < 5 * occ_rmsd_abs_thresh);
            occ_reasonable = (mysolid.rmsd_occupation < 5 * occ_rmsd_thresh);
            converged &= occ_stable;
         }



         if ((oprint) && ((icount%it_in==0) || converged))
         {
            if ((mysolid.fractional) &&  (!mysolid.fractional_frozen))
            {
               coutput << Ifmt(10)    << icount
                       << Efmt(25,12) << total_energy
                       << Efmt(25,12)  << total_energy - mysolid.E[28]
                       << Efmt(16,6)  << mysolid.E[28]
                       << Efmt(16,6)  << mysolid.smearfermi[0]
                       << Efmt(16,6)  << deltae
                       << Efmt(16,6)  << deltae-deltasmear
                       << Efmt(16,6)  << deltasmear 
                       << Efmt(16,6)  << deltac <<  std::endl;
               coutput << "                 "
                       << "SCF α=" << Efmt(8,3) << alpha_current
                       << "  SCF error=" << Efmt(8,3) << scf_error
                       << "  occ_rmsd="  << Efmt(10,3) << mysolid.rmsd_occupation
                       << "  Δocc_rmsd=" << Efmt(10,3) << (mysolid.rmsd_occupation - mysolid.rmsd_occupation_prev);
               if (!occ_reasonable && !converged)
                  coutput << "  warning: occupations still fluctuating";
               coutput << std::endl;
            }
            else
               coutput << Ifmt(10)    << icount
                       << Efmt(25,12) << total_energy
                       << Efmt(16,6)  << deltae
                       << Efmt(16,6)  << deltac 
                       << Efmt(16,6)  << total_energy-total_energy_fixedV << std::endl;
         }

         bool tolerance_occupation_update =  drho < drho_for_occ_update;
         bool periodic_occupation_update = (bfgscount % occ_update_stride == 0);

         const bool occ_update_allowed = mysolid.initial_occupation_update;
         const bool occ_update_due     = periodic_occupation_update;
         const bool occ_update_safe    = tolerance_occupation_update;

         //mysolid.occupation_update = occ_update_allowed && occ_stable &&  (occ_update_safe || occ_update_due);
         mysolid.occupation_update = occ_update_allowed &&  (occ_update_safe || occ_update_due);

   
         // Finalize SCF step with updated potentials
         mysolid.gen_scf_potentials_from_rho1();
      }



   } else if (minimizer == 4) {
      if (mysolid.newpsi) {
        int it_in0 = 15;
        for (int it = 0; it < it_in0; ++it)
          mysolid.sd_update_sic(dte);
        if (oprint)
          coutput << "        - " << it_in0
                  << " steepest descent iterations performed" << std::endl;
      }
      while ((icount < it_out) && (!converged)) {
        ++icount;
        if (stalled) {
          for (int it = 0; it < it_in; ++it)
            mysolid.sd_update_sic(dte);
          if (oprint)
            coutput << "        - " << it_in
                      << " steepest descent iterations performed" << std::endl;
          bfgscount = 0;
        }
        deltae_old = deltae;
      //  total_energy =
      //      band_cgsd_cgminimize2(mysolid, mygeodesic12.mygeodesic2, E, &deltae,
      //                            &deltac, bfgscount, it_in, tole, tolc);
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
   } else if (minimizer == 6) {
   } else if (minimizer == 7) {
      if (mysolid.newpsi) {
         int it_in0 = 15;
         for (int it=0; it<it_in0; ++it)
            mysolid.sd_update_sic(dte);
         if (oprint)
            coutput << "        - " << it_in0 << " steepest descent iterations performed" << std::endl;
      }
      band_lmbfgs2 psi_lmbfgs2(mygeodesic12.mygeodesic2, lmbfgs_size);
      while ((icount < it_out) && (!converged)) {
         ++icount;
         if (stalled) {
            for (int it=0; it<it_in; ++it)
               mysolid.sd_update_sic(dte);
            if (oprint)
               coutput << "        - " << it_in << " steepest descent iterations performed" << std::endl;
            bfgscount = 0;
         }
         deltae_old = deltae;
       //  total_energy = band_cgsd_bfgsminimize2(mysolid,mygeodesic12.mygeodesic2,psi_lmbfgs2,
       //                                         E,&deltae,&deltac,bfgscount,it_in,tole,tolc);
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

      if (mysolid.newpsi)
      {
         int it_in0 = 15;
         for (int it=0; it<it_in0; ++it)
            mysolid.sd_update(dte);
         if (oprint) coutput << "        - " << it_in0 << " steepest descent iterations performed" << std::endl;
      }

      // Initial SCF setup
      // generate density and rotate orbitals 
      // Generate initial density, potential and then orbital diagonalization
      double total_energy0 = mysolid.energy(); // Run Hψ = Eψ and compute E[0]
      mysolid.gen_hml();                       // Generate ⟨ψ|H|ψ⟩ (stored in hml)
      mysolid.diagonalize();      // Diagonalize H matrix (sets eig)
      mysolid.rotate1to2();       // Rotate ψ₁ to ψ₂ using eigenvectors
      mysolid.swap_psi1_psi2();   // Swap ψ₂ → ψ₁ (start clean state)


      // Fractional occupation: initialize occupation if needed
      // [Insert your occupation initialization routine here if using smeared/fractional occs]

      // Normalize total density (diagnostic)
      double x,sumxx = 0.0;
      int nfft3d = mygrid->nfft3d;
      int ispin = mygrid->ispin;
      double omega = mygrid->lattice->omega();
      double scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
      double dv = omega * scal1;
      for (int i=0; i < nfft3d; ++i)
      {
         x = (mysolid.rho1[i]);
         x += (mysolid.rho1[i+(ispin-1)*nfft3d]);
         sumxx += x;
      }
      double sum0  = mygrid->c3db::parall->SumAll(1, sumxx) * dv;
      
      //std::cout << "total energy0=" << total_energy0 << std::endl;
      //std::cout << "total sum0=" << sum0 << std::endl;

      // Setup SCF loop
      int ks_it_in  = control.ks_maxit_orb();
      int ks_it_out = control.ks_maxit_orbs();
      double scf_error = 0.0;
      double *vout = mysolid.rho1;
      double *vnew = mysolid.rho2;

      nwpw_cscf_mixing scfmix(mygrid,kerker_g0,
                              scf_algorithm,scf_alpha,scf_beta,diis_histories,
                              mygrid->ispin,mygrid->nfft3d,vout);

      while ((icount < (it_out*it_in)) && (!converged))
      {
         ++icount;

         // Restart SCF if stalled
         if ((stalled) && (icount%it_in==0))
         {
            for (int it=0; it<it_in; ++it)
               mysolid.sd_update(dte);

            // rotate orbitals
            if (extra_rotate)
            {
               mysolid.gen_hml();
               mysolid.diagonalize();
               mysolid.rotate1to2();
               mysolid.swap_psi1_psi2();
            }

            // [Optional] Re-initialize fractional occupation
            // define fractional occupation here

            scfmix.reset_mix(vout);

            bfgscount = 0;
            stalled = false;
         }

         deltae_old     = deltae;
         deltasmear_old = deltasmear;

         // Copy current density to old
         std::memcpy(vnew,vout,ispin*nfft3d*sizeof(double));

         // Minimize energy wrt orbitals
         double smear0 = mysolid.E[28];
         //if (icount == (it_out*it_in))
         //   mysolid.occupation_update = false;

         total_energy = band_cgsd_bybminimize2(mysolid,mygeodesic12.mygeodesic1,E,&deltae,
                                          &deltac,bfgscount,ks_it_in,ks_it_out,
                                          tole,tolc);

         //double mu = mysolid.smearfermi[0];
         // Optional orbital rotation post-minimization
         if (extra_rotate)
         {
            total_energy = mysolid.energy();
            mysolid.gen_hml();
            mysolid.diagonalize();
            mysolid.rotate1to2();
            mysolid.swap_psi1_psi2();
         }

         //  Generate updated density from current ψ
         mysolid.gen_rho1(); // updatating rho1==vout

         // [Insert fractional occupation update here if needed]
         // if (mysolid.fractional) update_occupations(...);


         //std::cout << "total_energy=" << total_energy << " " << total_energy2 
         //          << " " << total_energy2 - total_energy << std::endl;

         //define fractional occupation here
         scfmix.mix(vout,vnew,deltae,&scf_error);
         std::memcpy(vout,vnew,ispin*nfft3d*sizeof(double));

         deltac = scf_error;

         //deltac = mysolid.rho_error();
         deltae = total_energy - total_energy0;
         deltasmear = mysolid.E[28] - smear0;
         total_energy0 = total_energy;
         ++bfgscount;

         converged = (std::fabs(deltae) < tole) && (deltac < tolc);

         if ((oprint) && ((icount%it_in==0) || converged))
         {
            if ((mysolid.fractional) &&  (!mysolid.fractional_frozen))
               coutput << Ifmt(10)    << icount
                       << Efmt(25,12) << total_energy
                       << Efmt(25,12)  << total_energy - mysolid.E[28]
                       << Efmt(16,6)  << mysolid.E[28]
                       << Efmt(16,6)  << mysolid.smearfermi[0]
                       << Efmt(16,6)  << deltae
                       << Efmt(16,6)  << deltae-deltasmear
                       << Efmt(16,6)  << deltasmear 
                       << Efmt(16,6)  << deltac << std::endl;
            else
               coutput << Ifmt(10)    << icount
                       << Efmt(25,12) << total_energy
                       << Efmt(16,6)  << deltae
                       << Efmt(16,6)  << deltac << std::endl;

         }

         // Finalize SCF step with updated potentials
         mysolid.gen_scf_potentials_from_rho1();

      }
   }

   else if (minimizer == 9)
   {
   }

   else if (minimizer == 10)
   {
   }

 
   if (oprint) 
   {
      if (converged)  coutput << "     *** tolerance ok. iteration terminated" << std::endl;
      if (!converged) coutput << "     *** arrived at the Maximum iteration.  terminated" << std::endl;
      coutput << "          >>> iteration ended at   " << util_date() << "  <<<" << std::endl;
   }
 
   /* report summary of results */
   total_energy  = mysolid.gen_all_energies();
   if (oprint) 
   {
      coutput << mysolid;
      coutput << std::endl;
   }
   coutput << mysolid.print_filled(oprint);;
   //coutput << mysolid(oprint);
 
   return total_energy;
}

/******************************************
 *                                        *
 *         band_cgsd_energy_gradient      *
 *                                        *
 ******************************************/
void band_cgsd_energy_gradient(Solid &mysolid, double *grad_ion) 
{
   Parallel *parall = mysolid.mygrid->c3db::parall;
   bool oprint = parall->base_stdio_print;
   oprint = false;

   mysolid.psi_1local_force(grad_ion);

   if (oprint)
   {
      std::cout << std::endl << " local Ion Forces (au):" << std::endl;
      for (auto ii=0; ii< mysolid.myion->nion; ++ii)
         std::cout << Ifmt(4) << ii + 1 << " " << mysolid.myion->symbol(ii) << "\t( "
                   << Ffmt(10,5) << grad_ion[3 *ii] << " "
                   << Ffmt(10,5) << grad_ion[3*ii+1] << " "
                   << Ffmt(10,5) << grad_ion[3*ii+2] << " )\n";
      std::cout << std::endl << std::endl;
   }

   mysolid.psi_1nonlocal_force(grad_ion);

   if (oprint)
   {
      std::cout << std::endl << " nonlocal Ion Forces (au):" << std::endl;
      for (auto ii=0; ii< mysolid.myion->nion; ++ii)
         std::cout << Ifmt(4) << ii + 1 << " " << mysolid.myion->symbol(ii) << "\t( "
                   << Ffmt(10,5) << grad_ion[3 *ii] << " "
                   << Ffmt(10,5) << grad_ion[3*ii+1] << " "
                   << Ffmt(10,5) << grad_ion[3*ii+2] << " )\n";
      std::cout << std::endl << std::endl;
   }
   mysolid.semicore_force(grad_ion);
   if (oprint)
   {
      std::cout << std::endl << " semicore Ion Forces (au):" << std::endl;
      for (auto ii=0; ii< mysolid.myion->nion; ++ii)
         std::cout << Ifmt(4) << ii + 1 << " " << mysolid.myion->symbol(ii) << "\t( "
                   << Ffmt(10,5) << grad_ion[3 *ii] << " "
                   << Ffmt(10,5) << grad_ion[3*ii+1] << " "
                   << Ffmt(10,5) << grad_ion[3*ii+2] << " )\n";
      std::cout << std::endl << std::endl;
   }
 
   mysolid.ion_fion(grad_ion);
   if (oprint)
   {
      std::cout << std::endl << " ion-ion Ion Forces (au):" << std::endl;
      for (auto ii=0; ii< mysolid.myion->nion; ++ii)
         std::cout << Ifmt(4) << ii + 1 << " " << mysolid.myion->symbol(ii) << "\t( "
                   << Ffmt(10,5) << grad_ion[3 *ii] << " "
                   << Ffmt(10,5) << grad_ion[3*ii+1] << " "
                   << Ffmt(10,5) << grad_ion[3*ii+2] << " )\n";
      std::cout << std::endl << std::endl;
   }
}

} // namespace pwdft

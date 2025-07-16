#ifndef _SOLID_HPP_
#define _SOLID_HPP_

#pragma once

#include <iomanip>
#include <iostream>
#include <vector>

#include "Control2.hpp"
#include "cElectron.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "Cneb.hpp"
#include "CPseudopotential.hpp"
#include "CStrfac.hpp"
#include "cpsi.hpp"

namespace pwdft {

#define ionstream(A, B, C)                                                     \
  (A) << std::scientific << std::setw(19) << std::setprecision(10) << (B)      \
      << std::setw(0) << " (" << std::setw(15) << std::setprecision(5) << (C)  \
      << std::setw(0) << " /ion)"                                              \
      << "\n"
#define elcstream(A, B, C)                                                     \
  (A) << std::scientific << std::setw(19) << std::setprecision(10) << (B)      \
      << std::setw(0) << " (" << std::setw(15) << std::setprecision(5) << (C)  \
      << std::setw(0) << " /electron)"                                         \
      << "\n"

#define eig1stream(A, B)                                                       \
  std::scientific << std::setw(18) << std::setprecision(7) << (A)              \
                  << std::setw(0) << " (" << std::fixed << std::setw(8)        \
                  << std::setprecision(3) << (B) << std::setw(0) << "eV)\n"
#define eig2stream(A, B, C, D)                                                 \
  std::scientific << std::setw(18) << std::setprecision(7) << (A)              \
                  << std::setw(0) << " (" << std::fixed << std::setw(8)        \
                  << std::setprecision(3) << (B) << std::setw(0) << "eV) "     \
                  << std::scientific << std::setw(18) << std::setprecision(7)  \
                  << (C) << std::setw(0) << " (" << std::fixed << std::setw(8) \
                  << std::setprecision(3) << (D) << std::setw(0) << " eV)\n"

/*
E[0] =
E[1] =
E[2] =
E[3] =
E[4] =
E[5] =
E[6] =
*/


//
// Solid.hpp
//
// The `Solid` class defines the quantum system in a plane-wave DFT simulation.
// It is the central data wall in the Walls and Mirrors design: it owns the evolving
// state of the system (wavefunctions, densities, eigenvalues, occupations, Hamiltonians)
// and supports iterative updates via SCF and orbital optimization.
//
// All physical transformations of the system (e.g., diagonalization, density generation,
// steepest descent, Hamiltonian updates) occur through this class. It mirrors ψ, ρ, and
// occupation vectors across SCF iterations and calls into operator walls like
// `cElectron_Operators` to perform Hamiltonian applications and energy calculations.
//
// This class should remain the single source of truth for the system state.
//


class Solid {

  double omega, scal2, scal1, dv;

  int ispin, ne[2], neall, nbrillq, nbrillouin, n2ft3d, nfft3d, shift1, shift2, mshift;
  int nextra[2]     = {0,0};
  int ne_excited[2] = {0,0};
  int nfft[3];
  int version = 5;

public:
   Cneb *mygrid;
   Ion *myion;
   CStrfac *mystrfac;
   Ewald *myewald;
   cElectron_Operators *myelectron;
   CPseudopotential *mypsp;

 
   // --- Primary quantum state (mirrored between ψ1 and ψ2) ---
   double *psi1,*rho1,*rho1_all,*dng1; // Current orbitals, densities, and fftd densities
   double *psi2,*rho2,*rho2_all,*dng2;

   // --- Hamiltonian matrix and eigenvalues (mirrored with ψ) ---
   // hml  stores the matrix ⟨ψ1|H|ψ1⟩ before diagonalization.
   // In rotate1to2(), it is used to rotate orbitals via:
   //     ψ₂ = hml† ψ₁   (if hml holds eigenvector matrix, V)
   //
   // Note: hml may be overwritten post-diagonalization.
   //       Use hml as a temporary workspace unless separate storage is required.
   double *hml  = nullptr;  // ⟨ψ1|H|ψ1⟩
   double *eig  = nullptr;  // KS eigenvalues of ψ1

   double *hml2 = nullptr;  // ⟨ψ2|H|ψ2⟩ if ψ2 is used in SD or CG updates
   double *eig2 = nullptr;  // KS eigenvalues of ψ2


   // --- Occupation vectors (mirrored for fractional SCF) ---
   double *eig_prev = nullptr;
   double *occ1 = nullptr;
   double *occ2 = nullptr;

   // --- Constraint and transformation matrices (not mirrored) ---
   double *lmbda = nullptr;  // Lagrange multipliers: enforces ⟨ψ|ψ⟩ = I during SD/CG updates
   double *rotmat = nullptr; // Rotation matrix (unitary): diagonalizes hml to rotate ψ1 → ψ2

   // --- Excited (virtual) orbitals generated post-SCF (not used in minimization) ---
   // psi1_excited holds virtual (unoccupied) orbitals generated post-SCF.
   // These are used for diagnostics, spectroscopy, and excited-state analyses.
   // psi2_excited is not currently used — excited states are not evolved or minimized.
   // If excited-state optimization (e.g., TDDFT or constrained DFT) is implemented,
   // psi2_excited can be added as a mirror of psi1_excited.
   double *psi1_excited = nullptr;
   double *psi2_excited = nullptr;  // Unused for now — reserved for future capability

   double *hml_excited = nullptr;   // ⟨ψ_virtual|H|ψ_virtual⟩
   double *eig_excited = nullptr;   // KS eigenvalues for virtual orbitals

   int multiplicity;
   double total_charge;

   // psi smearing block
   bool fractional=false;
   bool fractional_frozen=false;
   int smearoccupation, smeartype;
   double smearfermi[2], smearcorrection, smearkT;
   double fractional_alpha, fractional_alpha_min, fractional_alpha_max, fractional_beta, fractional_gamma, fractional_rmsd_threshold;
   int fractional_it=0;
   bool occupation_update = false;;
 
   double E[80],en[2],ep,sp,tole;
 
   bool newpsi = false;
 
   /* Constructors */
   Solid(char *,bool,Cneb *,Ion *,CStrfac *,Ewald *,cElectron_Operators *,CPseudopotential *,Control2 &, std::ostream &);
 
   /* destructor */
   ~Solid() 
   {
      // --- Deallocate wavefunction/density memory ---
      if (psi1)     mygrid->g_deallocate(psi1);
      if (rho1)     mygrid->r_pack_deallocate(rho1);
      if (rho1_all) mygrid->r_pack_deallocate(rho1_all);
      if (dng1)     mygrid->c_pack_deallocate(dng1);
 
      if (psi2)     mygrid->g_deallocate(psi2);
      if (rho2)     mygrid->r_pack_deallocate(rho2);
      if (rho2_all) mygrid->r_pack_deallocate(rho2_all);
      if (dng2)     mygrid->c_pack_deallocate(dng2);
 
      // --- Hamiltonian and eigenvalue memory ---
      if (hml)      mygrid->w_deallocate(hml);
      if (eig)      delete[] eig;
      if (eig_prev) delete[] eig_prev;
 
      // --- Fractional occupation ---
      // 🔥 THESE ARE new[] ALLOCATED, NOT GRID-ALLOCATED
      if (occ1)     delete[] occ1;
      if (occ2)     delete[] occ2;
 
      if (lmbda)    mygrid->w_deallocate(lmbda);
 
      // --- Excited state orbitals ---
      if (psi1_excited) mygrid->g_deallocate(psi1_excited);
      if (psi2_excited) mygrid->g_deallocate(psi2_excited);
      if (hml_excited)  mygrid->w_deallocate(hml_excited);
      if (eig_excited)  delete[] eig_excited;
 
      // --- Optional / commented out ---
      // if (hml2)   mygrid->w_deallocate(hml2);
      // if (eig2)   delete[] eig2;
   }


   void ecpsi_initialize(char *,bool, const int *, std::ostream &);
   void ecpsi_rename(const char *, std::ostream &);
   void ecpsi_finalize(char *, std::ostream &);
   void ecpsi_minimize(double *, std::ostream &);
   void compute_Horb_for_cg(const int, double *, double *, double *);
   double ecpsi_KS_update_virtual(const int, const int, const int, const int, const double, const double, double *, double *, double *, double *, double *,  std::ostream &);

   void ecpsi_linesearch_update(const int, double, double, double *, double *, double *, double *);
   void ecpsi_sort_virtual(const int, double *, double *);

 
   /* write psi solid */
   void writecpsi(char *output_filename, std::ostream &coutput) {
      //cpsi_write(mygrid,&version,nfft,mygrid->lattice->unita_ptr(),&ispin,ne,&nbrillouin,psi1,output_filename,coutput);
      cpsi_write(mygrid, &version, nfft, mygrid->lattice->unita_ptr(),
                 &ispin, ne, &nbrillouin, psi1, &smearoccupation, occ1,
                 output_filename, coutput);
   }

   void writecpsi_excited(char *output_filename, std::ostream &coutput) {
      cpsi_write(mygrid,&version,nfft,mygrid->lattice->unita_ptr(),&ispin,ne,&nbrillouin,psi1_excited,output_filename,coutput);
   }

   double cpsi_KS_update(const int, const double, const double , double *, const int, const int *, const int,
                         double *,  double *, std::ostream &);

   double cpsi_KS_update_orb(const int, const int, const int, const int, const double,
                             const double, double *, double *, double *, double *, std::ostream &);

   void cpsi_linesearch_update(const int, double, double, double *, double *, double *, double *);

   void replace_excited_psi1(Control2 &, std::ostream &);

 
   /* solid energy */
   double energy() {
      myelectron->run(psi1, rho1, dng1, rho1_all, occ1);
      E[0] = (myelectron->energy(psi1, rho1, dng1, rho1_all, occ1) + myewald->energy());
      return E[0];
   }

   double energy0() {
      myelectron->run0(psi1);
      E[0] = (myelectron->energy(psi1, rho1, dng1, rho1_all, occ1) + myewald->energy());
      return E[0];
   }
 
   double psi2_energy() {
      myelectron->run(psi2, rho2, dng2, rho2_all);
      E[0] = (myelectron->energy(psi2, rho2, dng2, rho2_all) + myewald->energy());
      return E[0];
   }

   double psi2_energy0() {
      myelectron->run0(psi2);
      E[0] = (myelectron->energy(psi2, rho2, dng2, rho2_all) + myewald->energy());
      return E[0];
   }
 
   /* solid energy and eigenvalues */
   double energy_eigenvalues() {
      myelectron->run(psi1, rho1, dng1, rho1_all);
      E[0] = (myelectron->energy(psi1, rho1, dng1, rho1_all) + myewald->energy());
      
      /* generate eigenvalues */
      myelectron->gen_hml(psi1, hml);
      //mygrid->m_diagonalize(hml, eig);
      mygrid->w_diagonalize(hml, eig);
      
      return E[0];
   }
 
   /* solid energy and eigenvalues and other energies and en */
   double gen_all_energies() 
   {
      myelectron->run(psi1, rho1, dng1, rho1_all,occ1);
      myelectron->gen_energies_en(psi1, rho1, dng1, rho1_all, E, en,occ1);
      
      /*  ion-ion energy */
      E[4] = myewald->energy();
      E[0] += E[4];

      /* get contraints energies */
      if (myion->has_ion_bond_constraints()) {
         E[70] = myion->energy_ion_bond_constraints();
         E[0] += E[70];
      }
      if (myion->has_ion_bondings_constraints()) {
         E[71] = myion->energy_ion_bondings_constraints();
         E[0] +=  E[71];
      }

      
      /* generate eigenvalues */
      myelectron->gen_hml(psi1, hml);
      //mygrid->m_diagonalize(hml, eig);
      mygrid->w_diagonalize(hml, eig);
      
      if ((fractional) && (!fractional_frozen))
      {
         //std::cout << "Define occupations, eig=" << eig[0] << " " << eig[1] << " " << eig[2] << " " << eig[3] << " " 
         //                                        << eig[4] << " " << eig[5] << " " << eig[6] << " " << eig[7] << std::endl;
         //std::cout << "Define occupations, fractional_it=" << fractional_it << std::endl;
         if (fractional_it == 0)  // Initialize eig_prev for the first iteration
            std::memcpy(eig_prev,eig,nbrillq*(ne[0]+ne[1])*sizeof(double));
         else             // Smooth eigenvalues in subsequent iterations
            for (size_t i=0; i<nbrillq*(ne[0]+ne[1]); ++i)
               eig_prev[i] = (1.0-fractional_gamma)*eig_prev[i] + fractional_gamma*eig[i];


         if ((fractional_it>0) && (smeartype>=0) && (occupation_update))
         {
            // Define occupations based on smoothed eigenvalues
            double smearcorrection_old = smearcorrection;
            mygrid->m_0define_occupation(-1.0, false,
                                      multiplicity,
                                      myion->total_zv(),total_charge,
                                      eig_prev,hml,occ2,
                                      smeartype,smearkT,smearfermi,&smearcorrection);

            // RMSD occupation computation
            double rmsd_occupation = 0.0;
            for (size_t i=0; i<nbrillq*(ne[0]+ne[1]); ++i)
            {
                double delta_occ = occ2[i] - occ1[i];
                rmsd_occupation += delta_occ * delta_occ;
            }
            rmsd_occupation = std::sqrt(rmsd_occupation / nbrillq*(ne[0]+ne[1]));

            // Adaptive alpha adjustment
            if (rmsd_occupation < fractional_rmsd_threshold)  // Converging well
               fractional_alpha = std::min(fractional_alpha_max, fractional_alpha * (1.0 + fractional_beta));
            else  // Oscillations or divergence
               fractional_alpha = std::max(fractional_alpha_min, fractional_alpha * (1.0 - fractional_beta));




            // Update occupations
            for (auto i=0; i<nbrillq*(ne[0]+ne[1]); ++i)
               occ1[i] = (1.0-fractional_alpha)*occ1[i] + fractional_alpha*occ2[i];
            //std::memcpy(occ2,occ1,(ne[0]+ne[1])*sizeof(double));

            // Debugging output (optional)
           /*  std::cout << " Iteration: " << fractional_it
                      << ", RMSD: " << rmsd_occupation
                      << ", Alpha: " << fractional_alpha
                      << ", Smear Correction: " << smearcorrection
                      << ", Delta Smear Correction: " << smearcorrection - smearcorrection_old << std::endl;;
           */


         }
         else
         {
            smearfermi[0] = 0.0;
            smearfermi[1] = 0.0;
            smearcorrection = 0.0;
            for (auto nbq=0; nbq<nbrillq; ++nbq)
            {
               int ishift = nbq*(ne[0]+ne[1]);
               double *eigk  = eig + ishift;
               double *occ1k = occ1 + ishift;

               smearfermi[0]   +=  mygrid->define_smearfermi(ne[0],eigk,occ1k);
               smearcorrection +=  mygrid->add_smearcorrection(smeartype,ne[0],eigk,occ1k,smearfermi[0],smearkT);
               if (ispin==1)
               {
                  smearcorrection *= 2.0;
               }
               else
               {
                  smearfermi[1]    +=  mygrid->define_smearfermi(ne[1],eigk+ne[0],occ1k+ne[0]);
                  smearcorrection +=  mygrid->add_smearcorrection(smeartype,ne[1],eigk+ne[0],occ1k+ne[0],smearfermi[0],smearkT);
               }
            }

         }
         E[28] = smearcorrection;
         //E[0]  +=  E[28];
         fractional_it++;
      }




      /* generate dipole */
      //mypsp->mydipole->gen_dipole(rho1);
      
      return E[0];
   }
 
   /* various solid energies */
   double eorbit() { return myelectron->eorbit(psi1); }
   double psi2_eorbit() { return myelectron->eorbit(psi2); }
   double ehartree() { return myelectron->ehartree(dng1); }
   double exc() { return myelectron->exc(rho1_all); }
   double pxc() { return myelectron->pxc(rho1); }
   double eke() { return myelectron->eke(psi1); }
   double vl_ave() { return myelectron->vl_ave(dng1); }
   double vnl_ave() { return myelectron->vnl_ave(psi1); }
   double eion() {
     double ee = 0.0;
     ee = myewald->energy();
     return ee;
   }
 
   /* solid - generate current hamiltonian */
   void gen_hml() { myelectron->gen_hml(psi1, hml); }

   void gen_vall() { myelectron->gen_vall(); }
   void get_vall(double *vall_out) { myelectron->get_vall(vall_out); }
   void set_vall(const double *vall_in) { myelectron->set_vall(vall_in); }
   void gen_rho1() { myelectron->genrho(psi1,rho1,occ1); }
   void gen_densities1() { myelectron->gen_psi_r(psi1);
                           myelectron->gen_densities(rho1,dng1,rho1_all,occ1); }
   void gen_scf_potentials_from_rho1() { myelectron->scf_update_from_dn(rho1,dng1,rho1_all);}

 
   /* solid - diagonalize the current hamiltonian */
   void diagonalize() { mygrid->w_diagonalize(hml, eig); }
   void rotate1to2() 
   { 
      int nshift0 = 2*(mygrid->neq[0]+mygrid->neq[1])*mygrid->CGrid::npack1_max();
      int mshift0 = 2*(mygrid->ne[0]*mygrid->ne[0] + mygrid->ne[1]*mygrid->ne[1]);
      double rone[2]  = {1.0,0.0};
      double rzero[2] = {0.0,0.0};
      for (auto nbq=0; nbq<(mygrid->nbrillq); ++nbq)
      {
         double *psi1k = psi1 + nbq*nshift0;
         double *psi2k = psi2 + nbq*nshift0;
         double *hmlk = hml   + nbq*mshift0;
         mygrid->fwf_Multiply(-1,psi1k,hmlk,rone,psi2k,rzero);
      }
   }




 
   /* solid - call phafacs and gen_vl_potential and semicore */
   void phafacs_vl_potential_semicore() {
      mystrfac->phafac();
      myewald->phafac();
      myelectron->gen_vl_potential();
      myelectron->semicore_density_update();
   }
 
   /* apply psi2 = psi1 - dte*Hpsi1 + lmbda*psi1*/
   void sd_update(double dte) {
 
      /* apply psi2 = psi1 + dte*Hpsi1 */
      myelectron->run(psi1, rho1, dng1, rho1_all);
      
      // myelectron->add_dteHpsi((-dte),psi1,psi2);
      myelectron->add_dteHpsi((dte), psi1, psi2);
      
      /* lagrange multiplier - Expensive */
      mygrid->ggw_lambda(dte, psi1, psi2, lmbda);
      
      /* pointer swap of psi2 and psi1 */
      double *t2 = psi2;
      psi2 = psi1;
      psi1 = t2;
   }
 
   /* apply psi2 = psi1 - dte*Hpsi1 + lmbda*psi1*/
   void sd_update_sic(double dte) {
 
      /* apply psi2 = psi1 + dte*Hpsi1 */
      myelectron->run(psi1, rho1, dng1, rho1_all);
      
      // myelectron->add_dteHpsi((-dte),psi1,psi2);
      myelectron->add_dteHpsi((dte), psi1, psi2);
      
      /* lagrange multiplier - Expensive */
      //mygrid->ggm_lambda_sic(dte, psi1, psi2, lmbda);
      
      /* pointer swap of psi2 and psi1 */
      double *t2 = psi2;
      psi2 = psi1;
      psi1 = t2;
   }

   /***********************************************
    *                                             *
    *              _psi_1_get_Tgradient           *
    *                                             *
    ***********************************************/
   /**
    * Calculates the total energy of a system and computes the gradient of the kinetic energy term
    * with respect to the wavefunction `psi1`.
    *
    * This function evaluates the electronic properties of the system by executing a series of
    * computational steps that embody both electronic energy evaluation and kinetic gradient determination:
    * 
    * Steps involved:
    * 1. Executes initial electronic calculations on the wavefunction (`psi1`) using the `run` method,
    *    which updates the charge density (`rho1`), its derivatives (`dng1`), and the all-electron
    *    charge density (`rho1_all`).
    *    - Note: During this step, the density and Kohn-Sham (KS) potentials are dynamically updated.
    * 
    * 2. Computes the total energy by summing the electronic energy obtained from the `myelectron` object 
    *    and the Ewald energy calculated through the `myewald` object, which accounts for long-range 
    *    electrostatic interactions.
    * 
    * 3. Constructs the Hamiltonian matrix (`hml`) for the system using the current wavefunction configuration, 
    *    which encapsulates contributions from updated potentials.
    * 
    * 4. Determines the gradient of the kinetic energy term relative to `psi1` using the Hamiltonian matrix, 
    *    and deposits these calculated values into the provided gradient vector `G1`.
    *
    * @param G1 A pointer to an array where the calculated kinetic energy gradient values will be stored.
    * @return Returns the total energy of the system, encompassing contributions from both electronic
    *         and Ewald energies.
    */
   double psi_1get_Tgradient(double *G1) 
   {
      double total_energy;
      myelectron->run(psi1, rho1, dng1, rho1_all);
      total_energy = myelectron->energy(psi1, rho1, dng1, rho1_all) + myewald->energy();
      myelectron->gen_hml(psi1, hml);
      myelectron->get_Tgradient(psi1, hml, G1);
     
      return total_energy;
   }

   /***********************************************
    *                                             *
    *              psi_1_get_Tgradient0           *
    *                                             *
    ***********************************************/
   /**
    * Calculates the total energy of a system and computes the gradient of the kinetic energy term
    * with respect to the wavefunction psi1 using a specific initialization routine.
    *
    * This function is similar to psi_1get_Tgradient but utilizes a modified or preliminary setup 
    * step found in the `myelectron->run0` method. This method performs initial computations 
    * without updating the density and the Kohn-Sham (KS) potentials.
    *
    * Detailed Steps:
    * 1. Executes initial electronic calculations on the wavefunction (`psi1`) using the `run0` method.
    *    - Note: The density and KS potentials remain unchanged during this step.
    * 2. Computes the total energy by summing the electronic energy (via `myelectron`) with 
    *    the Ewald energy (`myewald`), reflecting long-range electrostatic interactions.
    * 3. Constructs the Hamiltonian matrix (`hml`) for the system based on the current state of the 
    *    wavefunction.
    * 4. Calculates the gradient of the kinetic energy term relative to `psi1` and populates 
    *    the provided gradient vector `G1` with these values.
    *
    * @param G1 A pointer to an array where the calculated T-gradient values will be stored.
    * @return Returns the total energy of the system, which includes contributions from both electronic 
    *         and Ewald energies.
    *
    * Note: Employing `run0` suggests a specialized context or initial conditions, requiring 
    *       understanding of distinct application scenarios within broader computational workflows.
    */
   double psi_1get_Tgradient0(double *G1) 
   {
      double total_energy;
      myelectron->run0(psi1);
      total_energy = myelectron->energy(psi1, rho1, dng1, rho1_all) + myewald->energy();
      myelectron->gen_hml(psi1, hml);
      myelectron->get_Tgradient(psi1, hml, G1);
     
      return total_energy;
   }


 
   double psi_1get_TSgradient(double *G1) 
   {
      double total_energy;
      myelectron->run(psi1, rho1, dng1, rho1_all);
      total_energy = myelectron->energy(psi1, rho1, dng1, rho1_all) + myewald->energy();
      myelectron->gen_hmlt(psi1, hml);
      myelectron->get_Tgradient(psi1, hml, G1);
      
      return total_energy;
   }
 
   void swap_psi1_psi2() {
      /* pointer swap of psi2 and psi1 */
      double *t2 = psi2;
      psi2 = psi1;
      psi1 = t2;
   }
 
   /* calculates the difference squared  error between rho1 and rho2 */
   double rho_error() 
   {
      double x;
      double sumxx = 0.0;
      for (int i=0; i<nfft3d; ++i) 
      {
         x = (rho2[i] - rho1[i]);
         x += (rho2[i+(ispin-1)*nfft3d] - rho1[i+(ispin-1)*nfft3d]);
         sumxx += x * x;
      }
      return mygrid->c3db::parall->SumAll(1,sumxx) * dv;
   }
 
   std::vector<double> eig_vector() {
      int localsize = nbrillq * (ne[0] + ne[1]);
      int globalsize = mygrid->c3db::parall->np_k() * localsize;
      int totalsize = nbrillouin*(ne[0] + ne[1]);

      std::vector<double> eig_global;
      if (mygrid->c3db::parall->taskid_k() == MASTER)
         eig_global.resize(globalsize);

      mygrid->c3db::parall->Vector_GatherAll(3, localsize, eig, eig_global.data(), MASTER);

      if (mygrid->c3db::parall->taskid_k() == MASTER)
         return std::vector<double>(eig_global.begin(), eig_global.begin() + totalsize);
      else
         return {};
   }
 
   void psi_1local_force(double *grad_ion) {
      myelectron->vl_force(dng1, grad_ion);
   }
 
   void psi_1nonlocal_force(double *grad_ion) {
      myelectron->vnl_force(psi1, grad_ion);
   }
   void semicore_force(double *grad_ion) {
      myelectron->semicore_force(grad_ion);
   }
   void ewald_fion(double *grad_ion) { myewald->force(grad_ion); }
   void ion_fion(double *grad_ion) {
      if (myelectron->is_periodic())
         myewald->force(grad_ion);

      myion->add_contraint_force(grad_ion);
   }

   std::string print_virtual() 
   {
      std::stringstream stream;
 
      std::ios init(NULL);
      init.copyfmt(stream);
      std::string eoln = "\n";
      int nexall = ne_excited[0] + ne_excited[1];

      for (auto nbq=0; nbq<nbrillouin; ++nbq)
      {
         stream << eoln;
         stream << " Brillouin zone point: " << nbq + 1 << eoln << eoln;
         stream << " virtual orbital energies:" <<  eoln;
         int nn = ne_excited[0] - ne_excited[1];
         double ev = 27.2116;
      
         // Print the first set of excited states in reverse order without symmetry considerations
         for (int i = ne_excited[0]-1; i>=ne_excited[1]; --i) 
            stream << eig1stream(eig_excited[i + nbq*nexall], eig_excited[i + nbq*nexall]*ev);
      
         // Print the second set of excited states in reverse order without symmetry considerations
         for (int i = ne_excited[1]-1; i>=0; --i) 
         {
            stream << eig2stream(
                        eig_excited[i + nn + nbq*nexall], eig_excited[i + nn + nbq*nexall]*ev,
                        eig_excited[i + (ispin - 1) * ne_excited[0] + nbq*nexall],
                        eig_excited[i + (ispin - 1) * ne_excited[0] + nbq*nexall]*ev);
         }
         stream << eoln;
      }

      return stream.str();
   }


   std::string print_filled(const bool oprint) 
   {
      std::stringstream stream;
 
      std::ios init(NULL);
      init.copyfmt(stream);
      std::string eoln = "\n";

      for (auto nb=0; nb<nbrillouin; ++nb)
      {
         int nbq = mygrid->ktoindex(nb);
         int pk = mygrid->ktop(nb);

         int n = ne[0] + ne[1];
         double tmpeig[n];
         double tmpocc[n];
         std::memset(tmpeig,0,n*sizeof(double));
         std::memset(tmpocc,0,n*sizeof(double));
         if (pk==mygrid->c3db::parall->taskid_k())
         {
            std::memcpy(tmpeig,eig+nbq*n,n*sizeof(double));
            if (occ1) std::memcpy(tmpocc,occ1+nbq*n,n*sizeof(double));
         }
         mygrid->c3db::parall->Vector_SumAll(3,n,tmpeig);
         mygrid->c3db::parall->Vector_SumAll(3,n,tmpocc);
         //std::memcpy(tmpeig,eig+nbq*n,n*sizeof(double));

         if (oprint)
         {
            stream << eoln;
            stream << eoln;
            stream << mygrid->mybrillouin->print_zone_point(nb);
            stream << eoln;
            stream << " orbital energies:" << eoln;
         }
         int nn = ne[0] - ne[1];
         double ev = 27.2116;

         //for (int i=0; i<nn; ++i)
         //   os << eig1stream(tmpeig[mysolid.ne[0]-1-i], tmpeig[mysolid.ne[0]-1-i] * ev);

         if (oprint)
         {
            for (int i=0; i<nn; ++i)
            {
               //os << eig1stream(mymolecule.eig[i], mymolecule.eig[i] * ev);
               if (fractional)
               {
                  if ((tmpocc[i] < 1.e-3) && (tmpocc[i]>1.0e-12))
                     stream << Efmt(18,7) << tmpeig[i] << " (" << Ffmt(8,3) << tmpeig[i] * ev << "eV) occ="
                        << Efmt(9,3) << tmpocc[i] << std::endl;
                  else
                     stream << Efmt(18,7) << tmpeig[i] << " (" << Ffmt(8,3) << tmpeig[i] * ev << "eV) occ="
                        << Ffmt(5,3) << tmpocc[i] << std::endl;
               }
               else
                  stream << Efmt(18,7) << tmpeig[i] << " (" << Ffmt(8,3) << tmpeig[i] * ev << "eV)" << std::endl;
            }
            for (int i=0; i<ne[1]; ++i)
            {
                if (fractional)
                {
                   if ((tmpocc[i+nn] < 1.e-3) && (tmpocc[i+nn]>1.0e-12))
                      stream << Efmt(18,7) << eig[i+nn]    << " ("
                         << Ffmt(8,3)  << eig[i+nn]*ev << "eV)  occ="
                         << Efmt(9,3)  << occ1[i+nn]   << " ";
                   else
                      stream << Efmt(18,7) << eig[i+nn]    << " ("
                         << Ffmt(8,3)  << eig[i+nn]*ev << "eV)  occ="
                         << Ffmt(5,3)  << occ1[i+nn]   << " "
                         << Efmt(18,7) << eig[i+(ispin-1)*ne[0]]    << " ("
                         << Ffmt(8,3)  << eig[i+(ispin-1)*ne[0]]*ev << "eV) occ="
                         << Ffmt(5,3)  << occ1[i+(ispin-1)*ne[0]] << std::endl;
      
                   if ((tmpocc[i+(ispin-1)*ne[0]] < 1.e-3) && (tmpocc[i+(ispin-1)*ne[0]]>1.0e-12))
                      stream << Efmt(18,7) << tmpeig[i+(ispin-1)*ne[0]]    << " ("
                         << Ffmt(8,3)  << tmpeig[i+(ispin-1)*ne[0]]*ev << "eV) occ="
                         << Efmt(9,3)  << tmpocc[i + (ispin-1)*ne[0]] << std::endl;
                   else
                      stream << Efmt(18,7) << tmpeig[i+(ispin-1)*ne[0]]    << " ("
                         << Ffmt(8,3)  << tmpeig[i+(ispin-1)*ne[0]]*ev << "eV) occ="
                         << Ffmt(5,3)  << tmpocc[i + (ispin-1)*ne[0]] << std::endl;
                }
                else
                   stream << Efmt(18,7) << tmpeig[i+nn] << " ("
                      << Ffmt(8,3)  << tmpeig[i + nn] * ev << "eV) "
                      << Efmt(18,7) << tmpeig[i+(ispin-1)*ne[0]] << " ("
                      << Ffmt(8,3)  << tmpeig[i+(ispin-1)*ne[0]]*ev << "eV)" << std::endl;
      
            }
            stream << eoln;
         }

      }
      return stream.str();

   }


 

 
   friend std::ostream &operator<<(std::ostream &os, const Solid &mysolid) 
   {
      /* using old style c++ formatting */
      std::ios init(NULL);
      init.copyfmt(os);
      std::string eoln = "\n";
      os << "     =============  energy results (Solid object)  =============" << eoln;
      os << eoln << eoln;
     
      os << std::fixed << " number of electrons: spin up= " << std::setw(11)
         << std::setprecision(5) << mysolid.en[0]
         << "  down= " << std::setw(11) << std::setprecision(5)
         << mysolid.en[mysolid.ispin - 1] << " (real space)";
      os << eoln << eoln;
      os << eoln;
      os << ionstream(" total     energy    : ", mysolid.E[0],mysolid.E[0]/mysolid.myion->nion);
      os << elcstream(" total orbital energy: ", mysolid.E[1],mysolid.E[1]/mysolid.neall);
      os << elcstream(" hartree energy      : ", mysolid.E[2],mysolid.E[2]/mysolid.neall);
      os << elcstream(" exc-corr energy     : ", mysolid.E[3],mysolid.E[3]/mysolid.neall);

      os << ionstream(" ion-ion energy      : ", mysolid.E[4], mysolid.E[4]/mysolid.myion->nion);
     
      os << eoln;
      os << elcstream(" kinetic (planewave) : ", mysolid.E[5],mysolid.E[5]/mysolid.neall);
      os << elcstream(" V_local (planewave) : ", mysolid.E[6],mysolid.E[6]/mysolid.neall);
      os << elcstream(" V_nl    (planewave) : ", mysolid.E[7],mysolid.E[7]/mysolid.neall);
      os << elcstream(" V_Coul  (planewave) : ", mysolid.E[8],mysolid.E[8]/mysolid.neall);
      os << elcstream(" V_xc    (planewave) : ", mysolid.E[9],mysolid.E[9]/mysolid.neall);

      //if (mysolid.myelectron->is_v_apc_on())
      //   os << ionstream(" K.S. V_APC energy   : ",mysolid.E[52],mysolid.E[52]/mysolid.myion->nion);
      os << " Viral Coefficient   : " << std::setw(19) << std::setprecision(10)
         << (mysolid.E[9]+mysolid.E[8]+mysolid.E[7]+mysolid.E[6])/mysolid.E[5] << std::endl;

      if (mysolid.myion->has_ion_constraints())
      {
         os << std::endl;
         if (mysolid.myion->has_ion_bond_constraints())
            os << " spring bond         : " << Efmt(19,10) << mysolid.E[70] << " ("
                                             << Efmt(15,5)  << mysolid.E[70]/mysolid.myion->nion << " /ion)" << std::endl;
         if (mysolid.myion->has_ion_bondings_constraints())
            os << " spring bondings     : " << Efmt(19,10) << mysolid.E[71] << " ("
                                            << Efmt(15,5)  << mysolid.E[71]/mysolid.myion->nion << " /ion)" << std::endl;
      }

      /*
      for (auto nb=0; nb<mysolid.nbrillouin; ++nb)
      {
         int nbq = mysolid.mygrid->ktoindex(nb);
         int pk = mysolid.mygrid->ktop(nb);
       
         int n = mysolid.ne[0] + mysolid.ne[1];
         double tmpeig[n];
         double tmpocc[n];
         std::memset(tmpeig,0,n*sizeof(double));
         std::memset(tmpocc,0,n*sizeof(double));
         if (pk==mysolid.mygrid->c3db::parall->taskid_k())
         {
            std::memcpy(tmpeig,mysolid.eig+nbq*n,n*sizeof(double));
            if ( mysolid.occ1) std::memcpy(tmpocc,mysolid.occ1+nbq*n,n*sizeof(double));
         }
         //mysolid.mygrid->c3db::parall->Vector_SumAll(3,n,tmpeig);
         //mysolid.mygrid->c3db::parall->Vector_SumAll(3,n,tmpocc);
         //std::memcpy(tmpeig,mysolid.eig+nbq*n,n*sizeof(double));

         os << eoln;
         os << eoln;
         os << mysolid.mygrid->mybrillouin->print_zone_point(nb);
         os << eoln;
         os << " orbital energies:" << eoln;
         int nn = mysolid.ne[0] - mysolid.ne[1];
         double ev = 27.2116;

         //for (int i=0; i<nn; ++i)
         //   os << eig1stream(tmpeig[mysolid.ne[0]-1-i], tmpeig[mysolid.ne[0]-1-i] * ev);

         for (int i=0; i<nn; ++i)
         {
            //os << eig1stream(mymolecule.eig[i], mymolecule.eig[i] * ev);
            if (mysolid.fractional)
            {
               if ((tmpocc[i] < 1.e-3) && (tmpocc[i]>1.0e-12))
                  os << Efmt(18,7) << tmpeig[i] << " (" << Ffmt(8,3) << tmpeig[i] * ev << "eV) occ="
                     << Efmt(9,3) << tmpocc[i] << std::endl;
               else
                  os << Efmt(18,7) << tmpeig[i] << " (" << Ffmt(8,3) << tmpeig[i] * ev << "eV) occ="
                     << Ffmt(5,3) << tmpocc[i] << std::endl;
            }
            else
               os << Efmt(18,7) << tmpeig[i] << " (" << Ffmt(8,3) << tmpeig[i] * ev << "eV)" << std::endl;
         }
         for (int i=0; i<mysolid.ne[1]; ++i)
         {
             if (mysolid.fractional)
             {
                if ((tmpocc[i+nn] < 1.e-3) && (tmpocc[i+nn]>1.0e-12))
                   os << Efmt(18,7) << mysolid.eig[i+nn]    << " ("
                      << Ffmt(8,3)  << mysolid.eig[i+nn]*ev << "eV)  occ="
                      << Efmt(9,3)  << mysolid.occ1[i+nn]   << " ";
                else
                   os << Efmt(18,7) << mysolid.eig[i+nn]    << " ("
                      << Ffmt(8,3)  << mysolid.eig[i+nn]*ev << "eV)  occ="
                      << Ffmt(5,3)  << mysolid.occ1[i+nn]   << " "
                      << Efmt(18,7) << mysolid.eig[i+(mysolid.ispin-1)*mysolid.ne[0]]    << " ("
                      << Ffmt(8,3)  << mysolid.eig[i+(mysolid.ispin-1)*mysolid.ne[0]]*ev << "eV) occ="
                      << Ffmt(5,3)  << mysolid.occ1[i+(mysolid.ispin-1)*mysolid.ne[0]] << std::endl;
      
                if ((tmpocc[i+(mysolid.ispin-1)*mysolid.ne[0]] < 1.e-3) && (tmpocc[i+(mysolid.ispin-1)*mysolid.ne[0]]>1.0e-12))
                   os << Efmt(18,7) << tmpeig[i+(mysolid.ispin-1)*mysolid.ne[0]]    << " ("
                      << Ffmt(8,3)  << tmpeig[i+(mysolid.ispin-1)*mysolid.ne[0]]*ev << "eV) occ="
                      << Efmt(9,3)  << tmpocc[i + (mysolid.ispin-1)*mysolid.ne[0]] << std::endl;
                else
                   os << Efmt(18,7) << tmpeig[i+(mysolid.ispin-1)*mysolid.ne[0]]    << " ("
                      << Ffmt(8,3)  << tmpeig[i+(mysolid.ispin-1)*mysolid.ne[0]]*ev << "eV) occ="
                      << Ffmt(5,3)  << tmpocc[i + (mysolid.ispin-1)*mysolid.ne[0]] << std::endl;
             }
             else
                os << Efmt(18,7) << tmpeig[i+nn] << " ("
                   << Ffmt(8,3)  << tmpeig[i + nn] * ev << "eV) "
                   << Efmt(18,7) << tmpeig[i+(mysolid.ispin-1)*mysolid.ne[0]] << " ("
                   << Ffmt(8,3)  << tmpeig[i+(mysolid.ispin-1)*mysolid.ne[0]]*ev << "eV)" << std::endl;
      
         }
         os << eoln;

      //      os << eig2stream(tmpeig[i+nn], 
      //                       tmpeig[i+nn]*ev,
      //                       tmpeig[i+(mysolid.ispin-1)*mysolid.ne[0]],
      //                       tmpeig[i+(mysolid.ispin-1)*mysolid.ne[0]]*ev);
      //   os << eoln;
      }
    */
     
      // write dipoles
      //os << mysolid.mypsp->mydipole->shortprint_dipole();
     
      os.copyfmt(init);
     
      return os;
   }

   // Add a flag to force reinitialization of the wavefunction
   bool force_reinit_flag = false;
   void force_reinit_wavefunction() { force_reinit_flag = true; }
   bool should_force_reinit_wavefunction() const { return force_reinit_flag; }
   void clear_force_reinit_wavefunction() { force_reinit_flag = false; }

   int get_total_electrons() const { return ne[0] + ne[1]; }
};

} // namespace pwdft

#endif

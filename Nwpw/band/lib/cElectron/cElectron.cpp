
#include "cCoulomb.hpp"
#include "Ion.hpp"
#include "cKinetic.hpp"
#include "Cneb.hpp"
#include "CPseudopotential.hpp"
#include "CStrfac.hpp"
#include "blas.h"
#include "cExchange_Correlation.hpp"

#include "cpsi_H.hpp"

#include "cElectron.hpp"
#define mytaskid 1

namespace pwdft {

/********************************************
 *                                          *
 * cElectron_Operators::Electron_Operators  *
 *                                          *
 ********************************************/
cElectron_Operators::cElectron_Operators(Cneb *mygrid0, cKinetic_Operator *myke0,
                                         cCoulomb_Operator *mycoulomb0,
                                         cXC_Operator *myxc0,
                                         CPseudopotential *mypsp0) 
{
   mygrid = mygrid0;
   myke = myke0;
   mycoulomb = mycoulomb0;
   mypsp = mypsp0;
   myxc = myxc0;
 
   ispin = mygrid->ispin;
   neall = mygrid->neq[0] + mygrid->neq[1];
 
   /* allocate memory */
   Hpsi  = mygrid->g_allocate_nbrillq_all();
   psi_r = mygrid->h_allocate_nbrillq_all();
   xcp = mygrid->r_nalloc(ispin);
   xce = mygrid->r_nalloc(ispin);
 
   x = mygrid->c_alloc();
   rho = mygrid->c_alloc();
   vall = mygrid->c_alloc();
   vl = mygrid->c_pack_allocate(0);

   vc = mygrid->c_pack_allocate(0);
   vcall = mygrid->c_pack_allocate(0);

   hmltmp =  mygrid->w_allocate_nbrillq_all();
 
   omega = mygrid->lattice->omega();
   scal1 = 1.0/((double)((mygrid->nx)*(mygrid->ny)*(mygrid->nz)));
   scal2 = 1.0/omega;
   dv = omega*scal1;
 
   nfft3d = (mygrid->nfft3d);
   shift1 = 2*(mygrid->npack1_max());
   npack1 = shift1;
   shift2 = (mygrid->n2ft3d);
}

/********************************************
 *                                          *
 *      cElectron_Operators::gen_psi_r      *
 *                                          *
 ********************************************/
/**
 * @brief Transform the orbitals ψ from reciprocal space to real space.
 *
 * Converts the orbital ψ(G) in reciprocal space to ψ(r) in real space using
 * a backward FFT. The result is stored in the member variable psi_r.
 *
 * @param psi   Input orbitals in reciprocal space (Fourier coefficients).
 *
 * Notes:
 * - psi_r is typically used for evaluating real-space potentials or densities.
 * - Assumes mygrid is properly initialized with FFT parameters.
 */
void cElectron_Operators::gen_psi_r(double *psi) 
{
   /* convert psi(G) to psi(r) */
   mygrid->gh_fftb(psi,psi_r);
}

/********************************************
 *                                          *
 *     cElectron_Operators::gen_density     *
 *                                          *
 ********************************************/
void cElectron_Operators::gen_density(double *dn, double *occ) 
{
   /* generate dn */
   if (occ)
   {
      // Compute density including occupation numbers
      mygrid->hr_aSumSqr_occ(scal2,occ,psi_r,dn);
   }
   else
   {
     // Compute density without occupation numbers
      mygrid->hr_aSumSqr(scal2, psi_r, dn);
   }
}

/********************************************
 *                                          *
 *     cElectron_Operators::gen_densities   *
 *                                          *
 ********************************************/
/**
 * @brief Computes all density-related quantities from the current real-space wavefunction.
 *
 * This function performs three key steps:
 * 1. Computes the spin-resolved electron density `dn` from the real-space orbitals `psi_r`,
 *    optionally applying fractional occupations `occ`.
 * 2. Computes the total electron density `rho`, transforms it to reciprocal space, and stores
 *    the result in packed G-space format as `dng`, which is used for solving the Coulomb potential.
 * 3. Constructs `dnall`, the full spin-resolved real-space density used for exchange-correlation
 *    and semicore corrections. If a semicore pseudopotential is present, its density contribution
 *    is added accordingly.
 *
 * @param[out] dn     Spin-resolved real-space electron density (ρ↑, ρ↓).
 * @param[out] dng    Packed G-space density used in the Coulomb solver.
 * @param[out] dnall  Full spin-resolved real-space density (with semicore correction if applicable).
 * @param[in]  occ    Optional fractional occupation vector. If null, assumes full occupancy.
 *
 * @note This function assumes that the real-space representation `psi_r` has already been generated.
 */
void cElectron_Operators::gen_densities(double *dn, double *dng, double *dnall, double *occ) 
{
   /* generate dn */
   if (occ)
       mygrid->hr_aSumSqr_occ(scal2,occ,psi_r,dn);
   else
      mygrid->hr_aSumSqr(scal2,psi_r,dn);
     
 
   /* generate rho and dng */
   double *tmp = x;
   mygrid->rrc_Sum(dn, dn+(ispin-1)*nfft3d, rho);
   mygrid->cc_SMul(scal1, rho, tmp);
   // mygrid->rc_fft3d(tmp);
   mygrid->rc_pfft3f(0, tmp);
   mygrid->c_pack(0, tmp);
   mygrid->cc_pack_copy(0, tmp, dng);
 
   /* generate dnall - used for semicore corrections */
   if (mypsp->has_semicore()) 
   {
      for (int ms = 0; ms < ispin; ++ms)
         mygrid->rrr_SMulAdd(0.5, mypsp->semicore_density, dn+ms*nfft3d, dnall+ms*nfft3d);
   } 
   else 
   {
      for (int ms = 0; ms < ispin; ++ms)
         mygrid->rr_copy(dn+ms*nfft3d, dnall+ms*nfft3d);
   }
}

/*********************************************
 *                                           *
 *  cElectron_Operators::dn_to_dng_dnall     *
 *                                           *
 *********************************************/
/**
 * @brief Converts real-space spin-resolved densities into packed G-space
 *        and spin-summed forms required for SCF potential construction.
 *
 * This routine performs the following:
 * 1. Sums spin densities to create the total charge density ρ.
 * 2. Scales and FFTs the result into G-space, storing it in `dng`.
 * 3. Generates `dnall`, the modified real-space density used for
 *    exchange-correlation and semicore correction logic.
 *
 * If semicore pseudopotentials are active, `dnall` includes a scaled
 * addition of the semicore reference density.
 *
 * @param[in]  dn     Real-space spin-resolved electron densities (ρ↑, ρ↓).
 * @param[out] dng    Packed G-space total density (for Coulomb solver).
 * @param[out] dnall  Real-space total density (for XC and semicore corrections).
 */
void cElectron_Operators::dn_to_dng_dnall(double *dn, double *dng, double *dnall)

{

   /* generate rho and dng */
   double *tmp = x;
   mygrid->rrc_Sum(dn, dn+(ispin-1)*nfft3d, rho);
   mygrid->cc_SMul(scal1, rho, tmp);
   // mygrid->rc_fft3d(tmp);
   mygrid->rc_pfft3f(0, tmp);
   mygrid->c_pack(0, tmp);
   mygrid->cc_pack_copy(0, tmp, dng);
 
   /* generate dnall - used for semicore corrections */
   if (mypsp->has_semicore()) 
   {
      for (int ms = 0; ms < ispin; ++ms)
         mygrid->rrr_SMulAdd(0.5, mypsp->semicore_density, dn+ms*nfft3d, dnall+ms*nfft3d);
   } 
   else 
   {
      for (int ms = 0; ms < ispin; ++ms)
         mygrid->rr_copy(dn+ms*nfft3d, dnall+ms*nfft3d);
   }
}




/*********************************************
 *                                           *
 *  cElectron_Operators::gen_scf_potentials  *
 *                                           *
 *********************************************/
/**
 * @brief Computes self-consistent field (SCF) potentials from the density.
 *
 * Given real-space (`dn`, `dnall`) and packed G-space (`dng`) density components,
 * this function evaluates the following SCF contributions:
 *   - The Hartree (Coulomb) potential `vc`, stored in both `vc` and `vcall`.
 *   - The exchange-correlation potentials `xcp` (potential) and `xce` (energy).
 *
 * These quantities are used to construct the total effective potential
 * in later steps (e.g., `gen_vall()` or `gen_Hpsi_k()`).
 *
 * @param[in]  dn     Real-space spin-resolved density (used for p_xc).
 * @param[in]  dng    Packed G-space density (used for Hartree).
 * @param[in]  dnall  Spin-summed density (used for exchange-correlation).
 */
void cElectron_Operators::gen_scf_potentials(double *dn, double *dng, double *dnall)
{
   /* generate coulomb potential */
   mycoulomb->vcoulomb(dng, vc);
   mygrid->cc_pack_copy(0,vc,vcall);

 
   // generate exchange-correlation potential */
   myxc->v_exc_all(ispin, dnall, xcp, xce);
   // v_exc(ispin,shift2,dnall,xcp,xce,x);
 
}

/*********************************************
 *                                           *
 *  cElectron_Operators::scf_update_from_dn  *
 *                                           *
 *********************************************/
/**
 * @brief Updates SCF potentials given an existing density.
 *
 * This function is used when the electron density `dn` is already known,
 * and you want to recompute:
 *   - The Fourier-space (packed) density `dng` for Hartree potential.
 *   - The full spin-summed real-space density `dnall` for XC and semicore.
 *   - The SCF potentials (Hartree and XC) based on these densities.
 *
 * It avoids the need to regenerate `dn` from wavefunctions (ψ),
 * and is useful in orbital minimization algorithms (e.g., steepest descent).
 *
 * @param[in]  dn     Real-space spin densities (per spin channel).
 * @param[out] dng    Packed G-space density used for Coulomb potential.
 * @param[out] dnall  Spin-summed real-space density (with semicore correction if needed).
 */
void cElectron_Operators::scf_update_from_dn(double *dn, double *dng, double *dnall)
{
   dn_to_dng_dnall(dn, dng, dnall);
   gen_scf_potentials(dn, dng, dnall);
}

/********************************************
 *                                          *
 *    cElectron_Operators::gen_vl_potential *
 *                                          *
 ********************************************/
void cElectron_Operators::gen_vl_potential() 
{
   double dng0[1], fion0[1];

   /* generate local psp */
   mypsp->v_local(vl, 0, dng0, fion0);
}

/********************************************
 *                                          *
 *    cElectron_Operators::gen_vall         *
 *                                          *
 ********************************************/
/**
 * @brief Generates the total potential array `vall` in real space.
 *
 * Combines k-space contributions from the local potential (`vl`) and 
 * Coulomb potential (`vc`), scales by `scal2`, and transforms to real space
 * using inverse FFT. Adds exchange-correlation potential (`xcp`) in real space.
 *
 * The resulting `vall` array contains the real-space total effective potential,
 * and must be preallocated with size at least 2 × nfft3d (complex doubles).
 *
 * @note Requires `vall` to be allocated as a complex array.
 * 
 * @pre `vl`, `vc`, `xcp`, and `vall` must be properly initialized and sized.
 * 
 * @warning Incorrect allocation of `vall` (e.g., as real instead of complex)
 *          will result in heap-buffer-overflow errors.
 */
void cElectron_Operators::gen_vall()
{
   // add up k-space potentials, vall = scal2*vl + vc  ****
   mygrid->cc_pack_SMul(0, scal2, vl, vall);
   mygrid->cc_pack_Sum2(0,vc,vall);

   // fourier transform k-space potentials ****
   mygrid->c_unpack(0,vall);
   mygrid->cr_fft3d(vall);


   // add xcp to vall 
   if (ispin==2) mygrid->rcc_Sum(xcp+2*nfft3d,vall, vall+2*nfft3d);
   mygrid->rcc_Sum(xcp,vall,vall);
}

/********************************************
 *                                          *
 *      cElectron_Operators::get_vall       *
 *                                          *
 ********************************************/
void cElectron_Operators::get_vall(double *vall_out)
{
   mygrid->cc_copy(vall,vall_out);
   if (ispin==2) mygrid->cc_copy(vall+2*nfft3d,vall_out+2*nfft3d);
}

/********************************************
 *                                          *
 *      cElectron_Operators::set_vall       *
 *                                          *
 ********************************************/
void cElectron_Operators::set_vall(const double *vall_in)
{
   mygrid->cc_copy(vall_in,vall);
   if (ispin==2) mygrid->cc_copy(vall_in+2*nfft3d,vall+2*nfft3d);
}



/************************************************
 *                                              *
 * cElectron_Operators::semicore_density_update *
 *                                              *
 ************************************************/
void cElectron_Operators::semicore_density_update() 
{
   if (mypsp->has_semicore())
      mypsp->semicore_density_update();
}


/********************************************
 *                                          *
 *      cElectron_Operators::gen_Hpsi_k     *
 *                                          *
 ********************************************/
/* vl, vcall,xcp and psi_r  have been set */
/**
 * @brief Generate the action of the Kohn–Sham Hamiltonian on  orbitals
 *        in reciprocal space.
 *
 * Computes Hψ = (T + V_eff)ψ for every Kohn–Sham orbital ψ in reciprocal
 * space. The result is stored in the member variable Hpsi.
 *
 * @param psi   Input orbitals in reciprocal space (Fourier coefficients).
 * @param occ   Occupation number (used to control Fermi-level dependent terms).
 *
 * Notes:
 * - The Hamiltonian includes kinetic, local potential, and exchange-correlation terms.
 * - Non-local pseudopotential components are handled through cpsi_H().
 * - Hpsi is negated after construction to match minimization convention.
 * - The 'move' and 'fion0' arguments control force terms and are typically
 *   inactive in this context.
 */
void cElectron_Operators::gen_Hpsi_k(double *psi, double *occ) 
{
   bool move = false;
   double fion0[1];
 
   mygrid->g_zero(Hpsi);
 
   /* get Hpsi */
   cpsi_H(mygrid,myke,mypsp,psi,psi_r,vl,vcall,xcp,Hpsi,move,fion0,occ);
 
   mygrid->g_Scale(-1.0,Hpsi);
}


/********************************************
 *                                          *
 *      cElectron_Operators::gen_hml        *
 *                                          *
 ********************************************/
void cElectron_Operators::gen_hml(double *psi, double *hml) 
{
   mygrid->ggw_sym_Multiply(psi,Hpsi,hml);
}

/********************************************
 *                                          *
 *      cElectron_Operators::gen_hmlt       *
 *                                          *
 ********************************************/
void cElectron_Operators::gen_hmlt(double *psi, double *hmlt) 
{
   mygrid->ggw_Multiply(psi, Hpsi, hmlt);
}

/********************************************
 *                                          *
 *     cElectron_Operators::get_Tgradient   *
 *                                          *
 ********************************************/
void cElectron_Operators::get_Tgradient(double *psi, double *hml, double *THpsi) 
{
   double rone[2]  = {1.0,0.0};
   double rzero[2] = {0.0,0.0};
   int npack1 =   mygrid->CGrid::npack1_max();
   int npack2 = 2*mygrid->CGrid::npack1_max();
   int shift2 = (mygrid->neq[0]+mygrid->neq[1])*npack2;
   int shift1 = 2*(mygrid->ne[0]*mygrid->ne[0]+mygrid->ne[1]*mygrid->ne[1]);

  for (auto nbq=0; nbq<(mygrid->nbrillq); ++nbq)
  {
      double *psik = psi + nbq*shift2;
      double *THpsik = THpsi + nbq*shift2;
      double *hmlk = hml + nbq*shift1;
      mygrid->fwf_Multiply(-1,psik,hmlk,rone,THpsik,rzero);
      //mygrid->fwf_Multiply(-1,psi,hml,rone,THpsi,rzero);
   // mygrid->g_Scale(-1.0,THpsi);
   }
   mygrid->gg_Minus2(Hpsi,THpsi);
}


/********************************************
 *                                          *
 *     cElectron_Operators::gen_Tangent     *
 *                                          *
 ********************************************/
void cElectron_Operators::gen_Tangent(double *psi, double *hml, double *THpsi) 
{
   double rone[2]  = {1.0,0.0};
   double rmone[2]  = {-1.0,0.0};
   int npack1 =   mygrid->CGrid::npack1_max();
   int npack2 = 2*mygrid->CGrid::npack1_max();
   int shift2 = (mygrid->neq[0]+mygrid->neq[1])*npack2;
   int shift1 = 2*(mygrid->ne[0]*mygrid->ne[0]+mygrid->ne[1]*mygrid->ne[1]);

   for (auto nbq=0; nbq<(mygrid->nbrillq); ++nbq)
   {
      double *psik = psi + nbq*shift2;
      double *THpsik = THpsi + nbq*shift2;
      double *hmlk = hml + nbq*shift1;
      mygrid->fwf_Multiply(-1,psik,hmlk,rone,THpsik,rmone);
   }
}


/********************************************
 *                                          *
 *     cElectron_Operators::get_Gradient    *
 *                                          *
 ********************************************/
 /**
  * @brief Copies the result of the Hamiltonian operator Hψ to the provided output buffer.
  *
  * This function retrieves the latest result of the Hamiltonian applied to ψ (stored in `Hpsi`)
  * and copies it into the user-supplied array `THpsi`.
  *
  * @param[out] THpsi Destination array for the Hamiltonian-applied wavefunction (Hψ).
  *
  * @note This is typically used in gradient calculations where THpsi is compared to ψ for descent updates.
  */
void cElectron_Operators::get_Gradient(double *THpsi) 
{
   mygrid->gg_copy(Hpsi, THpsi);
}

/********************************************
 *                                          *
 *      cElectron_Operators::genrho         *
 *                                          *
 ********************************************/
 /**
  * @brief Generates the real-space electron density from a given wavefunction.
  *
  * This function performs two steps:
  * 1. Converts the complex-valued ψ (in reciprocal or packed form) into a real-space representation.
  * 2. Computes the electron density from the resulting ψ₍ᵣ₎, optionally using occupation numbers.
  *
  * @param[in]  psi Complex orbital wavefunction ψ.
  * @param[out] dn  Output real-space electron density array.
  * @param[in]  occ Optional occupation array (used for fractional occupations).
  *
  * @note This function is commonly called during SCF steps to update the density 
  *       from the current orbital guess.
  */
void cElectron_Operators::genrho(double *psi, double *dn, double *occ) 
{
   this->gen_psi_r(psi);
   this->gen_density(dn,occ);
}

/********************************************
 *                                          *
 *      cElectron_Operators::run            *
 *                                          *
 ********************************************/
/**
 * @brief Executes a full SCF operator application cycle for the current wavefunction.
 *
 * This routine:
 * 1. Converts the complex ψ into real-space form.
 * 2. Computes all relevant electron densities (ρ, ∇ρ, ρ_all).
 * 3. Constructs the SCF potentials (Hartree, XC, pseudopotential, etc.).
 * 4. Applies the Hamiltonian to ψ in reciprocal space (generates Hψ).
 *
 * @param[in]  psi     Complex orbital wavefunction ψ.
 * @param[out] dn      Real-space electron density ρ(𝐫).
 * @param[out] dng     Gradient-related or G-space modified density (for XC, etc.).
 * @param[out] dnall   Optional full spin-summed total density (used for post-SCF diagnostics).
 * @param[in]  occ     Optional occupation vector (used for fractional or smearing occupations).
 *
 * @note This is the primary interface used in SCF cycles to update all potentials 
 *       and generate Hψ = H[ρ]ψ.
 */
void cElectron_Operators::run(double *psi, double *dn, double *dng, double *dnall, double *occ) 
{
   ++counter;
   this->gen_psi_r(psi);
   // this->gen_density(dn);
   this->gen_densities(dn, dng, dnall,occ);
   this->gen_scf_potentials(dn, dng, dnall);
   this->gen_Hpsi_k(psi);
}

/********************************************
 *                                          *
 *      cElectron_Operators::run0           *
 *                                          *
 ********************************************/
/* densities and potentials have been set */
void cElectron_Operators::run0(double *psi)
{
   ++counter;
   this->gen_psi_r(psi);
   // this->gen_density(dn);
   //this->gen_densities(dn, dng, dnall,occ);
   //this->gen_scf_potentials(dn, dng, dnall);
   this->gen_Hpsi_k(psi);
}


/********************************************
 *                                          *
 *       cElectron_Operators::vl_ave        *
 *                                          *
 ********************************************/
/**
 * @brief Computes the average of the local potential over a given density.
 *
 * Calculates the expectation value ⟨ρ|V_local⟩ by integrating the
 * local potential `vl` against a given electron density `dng`.
 *
 * @param[in] dng  Gradient-corrected or packed density array.
 * @return         The average value ⟨ρ|V_local⟩, used for energy diagnostics.
 *
 * @note This is typically used in evaluating energy contributions from the local pseudopotential.
 */
double cElectron_Operators::vl_ave(double *dng) 
{
   return mygrid->cc_pack_dot(0, dng, vl);
}


/********************************************
 *                                          *
 *      cElectron_Operators::vnl_ave        *
 *                                          *
 ********************************************/
/**
 * @brief Computes the average nonlocal pseudopotential energy.
 *
 * Evaluates the expectation value ⟨ψ|V_nonlocal|ψ⟩ using the current wavefunctions `psi`
 * and (optionally) their occupations `occ`. This term represents the nonlocal contribution
 * from the pseudopotential, which typically includes angular-momentum projectors.
 *
 * @param[in] psi   Pointer to the wavefunction array.
 * @param[in] occ   (Optional) Pointer to the occupation array; if null, assumes full occupancy.
 * @return          Nonlocal pseudopotential energy contribution ⟨ψ|V_nonlocal|ψ⟩.
 *
 * @note This is evaluated via projectors stored in `CPseudopotential`.
 */
double cElectron_Operators::vnl_ave(double *psi, double *occ) 
{
   return mypsp->e_nonlocal(psi,occ);
}

/********************************************
 *                                          *
 *      cElectron_Operators::eorbit         *
 *                                          *
 ********************************************/
double cElectron_Operators::eorbit(double *psi, double *occ) 
{

   mygrid->ggw_sym_Multiply(psi,Hpsi,hmltmp);


   // mygrid->m_scal(-1.0,hmltmp);
   //double eorbit0 = mygrid->w_trace(hmltmp);
   double eorbit0 = occ ? mygrid->w_trace_occ(hmltmp,occ) : mygrid->w_trace(hmltmp);
   if (ispin==1)
      eorbit0 = eorbit0 + eorbit0;

   return eorbit0;
}

/********************************************
 *                                          *
 *      cElectron_Operators::ehartree       *
 *                                          *
 ********************************************/
double cElectron_Operators::ehartree(double *dng) 
{
   return mycoulomb->ecoulomb(dng);
}


/********************************************
 *                                          *
 *         cElectron_Operators::exc         *
 *                                          *
 ********************************************/
double cElectron_Operators::exc(double *dnall) 
{
   double excsum = mygrid->rr_dot(dnall, xce);
   if (ispin==1) 
   {
      excsum *= 2.0;
   } 
   else 
   {
      excsum += mygrid->rr_dot(dnall+nfft3d, xce);
   }
   excsum *= dv;
 
   return excsum;
}


/********************************************
 *                                          *
 *         cElectron_Operators::pxc         *
 *                                          *
 ********************************************/
double cElectron_Operators::pxc(double *dn) 
{
   double pxcsum = mygrid->rr_dot(dn, xcp);
   if (ispin==1) 
   {
      pxcsum *= 2.0;
   } 
   else 
   {
      pxcsum += mygrid->rr_dot(dn+nfft3d, xcp+nfft3d);
   }
   pxcsum *= dv;
 
   return pxcsum;
}

/********************************************
 *                                          *
 *        cElectron_Operators::eke          *
 *                                          *
 ********************************************/
double cElectron_Operators::eke(double *psi, double *occ) 
{
   return occ ? myke->ke_ave(psi,occ) : myke->ke_ave(psi);
}

/********************************************
 *                                          *
 *        cElectron_Operators::energy       *
 *                                          *
 ********************************************/
double cElectron_Operators::energy(double *psi, double *dn, double *dng, double *dnall, double *occ) 
{
   double total_energy, eorbit0, ehartr0, exc0, pxc0;
 
   /* total energy calculation */
   mygrid->ggw_sym_Multiply(psi, Hpsi, hmltmp);
   // mygrid->m_scal(-1.0,hmltmp);
   //eorbit0 = mygrid->w_trace(hmltmp);
   eorbit0 = occ ? mygrid->w_trace_occ(hmltmp,occ) : mygrid->w_trace(hmltmp);
   if (ispin == 1)
      eorbit0 = eorbit0 + eorbit0;
 

   ehartr0 = mycoulomb->ecoulomb(dng);
 
   exc0 = mygrid->rr_dot(dnall, xce);
   pxc0 = mygrid->rr_dot(dn, xcp);
   if (ispin == 1) 
   {
      exc0 = exc0 + exc0;
      pxc0 = pxc0 + pxc0;
   } 
   else 
   {
      exc0 += mygrid->rr_dot(dnall+nfft3d, xce);
      pxc0 += mygrid->rr_dot(dn+nfft3d, xcp+nfft3d);
   }
   exc0 *= dv;
   pxc0 *= dv;
 
   total_energy = eorbit0 + exc0 - ehartr0 - pxc0;
 
   return total_energy;
}

/********************************************
 *                                          *
 *   cElectron_Operators::gen_energies_en   *
 *                                          *
 ********************************************/
void cElectron_Operators::gen_energies_en(double *psi, double *dn, double *dng,
                                          double *dnall, double *E, double *en, double *occ) 
{
   double total_energy, eorbit0, ehartr0, exc0, pxc0;
   
 
   /* total energy calculation */
   mygrid->ggw_sym_Multiply(psi, Hpsi, hmltmp);
   // mygrid->m_scal(-1.0,hmltmp);
   //eorbit0 = mygrid->w_trace(hmltmp);
   eorbit0 = occ ? mygrid->w_trace_occ(hmltmp,occ) : mygrid->w_trace(hmltmp);
   if (ispin==1) eorbit0 = eorbit0 + eorbit0;
 
   ehartr0 = mycoulomb->ecoulomb(dng);
 
   exc0 = mygrid->rr_dot(dnall,xce);
   pxc0 = mygrid->rr_dot(dn,xcp);
   if (ispin==1) 
   {
      exc0 = exc0 + exc0;
      pxc0 = pxc0 + pxc0;
   } 
   else 
   {
      exc0 += mygrid->rr_dot(dnall+nfft3d, xce);
      pxc0 += mygrid->rr_dot(dn+nfft3d, xcp+nfft3d);
   }
   exc0 *= dv;
   pxc0 *= dv;
 
   total_energy = eorbit0 + exc0 - ehartr0 - pxc0;
 
   /* set the E[] energies */
   E[0] = total_energy;
   E[1] = eorbit0;
   E[2] = ehartr0;
   E[3] = exc0;
   E[4] = 0.0;
 
   E[5] = occ ? myke->ke_ave(psi,occ) :  myke->ke_ave(psi);
   E[6] = this->vl_ave(dng);
   E[7] = mypsp->e_nonlocal(psi,occ);
   E[8] = 2 * ehartr0;
   E[9] = pxc0;
 
   en[0] = dv * mygrid->r_dsum(dn);
   en[1] = en[0];
   if (ispin>1) en[1] = dv*mygrid->r_dsum(dn+nfft3d);
}



/********************************************
 *                                          *
 *    cElectron_Operators::add_dteHpsi      *
 *                                          *
 ********************************************/
void cElectron_Operators::add_dteHpsi(double dte, double *psi1, double *psi2) 
{
   /* do a steepest descent step */
   mygrid->gg_SMul(dte, Hpsi, psi2);
   mygrid->gg_Sum2(psi1, psi2);
}


/*********************************************
 *                                           *
 *    cElectron_Operators::semicore_force    *
 *                                           *
 *********** *********************************/
void cElectron_Operators::semicore_force(double *fion) 
{
   if (mypsp->has_semicore())
      mypsp->semicore_xc_fion(xcp, fion);
}


/********************************************
 *                                          *
 *     cElectron_Operators::vl_force        *
 *                                          *
 ********************************************/
void cElectron_Operators::vl_force(double *dng, double *fion) 
{
   mypsp->f_local(dng, fion);
}


/********************************************
 *                                          *
 *      Electron_Operators::vnl_force       *
 *                                          *
 ********************************************/
void cElectron_Operators::vnl_force(double *psi, double *fion) 
{
   mypsp->f_nonlocal_fion(psi, fion);
}

} // namespace pwdft

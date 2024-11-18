
#include "Coulomb12.hpp"
#include "Ion.hpp"
#include "Kinetic.hpp"
#include "Pneb.hpp"
#include "Pseudopotential.hpp"
#include "Strfac.hpp"
#include "blas.h"
#include "exchange_correlation.hpp"
//#include        "v_exc.hpp"

#include "psi_H.hpp"

#include "Electron.hpp"
#define mytaskid 1

namespace pwdft {

/********************************************
 *                                          *
 *  Electron_Operators::Electron_Operators  *
 *                                          *
 ********************************************/
Electron_Operators::Electron_Operators(Pneb *mygrid0, Kinetic_Operator *myke0,
                                       Coulomb12_Operator *mycoulomb120,
                                       XC_Operator *myxc0,
                                       Pseudopotential *mypsp0) 
{
   mygrid = mygrid0;
   myke = myke0;
   mycoulomb12 = mycoulomb120;
   mypsp = mypsp0;
   myxc = myxc0;
   periodic = mycoulomb12->has_coulomb1;
   aperiodic = mycoulomb12->has_coulomb2;
 
   ispin = mygrid->ispin;
   neall = mygrid->neq[0] + mygrid->neq[1];
 
   /* allocate memory */
   Hpsi = mygrid->g_allocate(1);
   psi_r = mygrid->h_allocate();
   xcp = mygrid->r_nalloc(ispin);
   xce = mygrid->r_nalloc(ispin);
 
   x = mygrid->r_alloc();
   rho = mygrid->r_alloc();
   vall = mygrid->r_nalloc(ispin);
   vl = mygrid->c_pack_allocate(0);
   if (periodic)
   {
      vc = mygrid->c_pack_allocate(0);
      vcall = mygrid->c_pack_allocate(0);
      if (mycoulomb12->dielectric_on())
         vdielec = mygrid->c_pack_allocate(0);
   }
   if (aperiodic)
   {
      vc = mygrid->r_alloc();
      vcall= mygrid->r_alloc();
      if (mycoulomb12->dielectric_on())
         vdielec = mygrid->r_alloc();
   }

   if (aperiodic)
      vlr_l = mygrid->r_alloc();
 
   hmltmp = mygrid->m_allocate(-1, 1);
 
   omega = mygrid->lattice->omega();
   scal1 = 1.0/((double)((mygrid->nx)*(mygrid->ny)*(mygrid->nz)));
   scal2 = 1.0/omega;
   dv = omega*scal1;
 
   n2ft3d = (mygrid->n2ft3d);
   shift1 = 2 * (mygrid->npack(1));
   npack1 = shift1;
   shift2 = (mygrid->n2ft3d);
}

/********************************************
 *                                          *
 *      Electron_Operators::gen_psi_r       *
 *                                          *
 ********************************************/
void Electron_Operators::gen_psi_r(double *psi) 
{
   /* convert psi(G) to psi(r) */
   mygrid->gh_fftb(psi,psi_r);
 
  /* 
    int indx1 = 0;
    int indx2 = 0;
    for (int i=0; i<neall; ++i)
    {
       double sum1 = mygrid->cc_pack_dot(1,psi+indx1,psi+indx1);
       std::cout << "sum1 = " << sum1 << std::endl; 

       mygrid->cc_pack_copy(1,psi+indx1,psi_r+indx2);
       double sum2 = mygrid->cc_pack_dot(1,psi_r+indx2,psi_r+indx2);
       std::cout << "sum2 = " << sum2 << std::endl; 

       mygrid->c_unpack(1,psi_r+indx2);

       double sum3 = mygrid->cc_dot(psi_r+indx2,psi_r+indx2);
       std::cout << "sum3 = " << sum3 << std::endl; 

       mygrid->cr_fft3d(psi_r+indx2);
       //mygrid->cr_pfft3b(1,psi_r + indx2);

       double sum4 = mygrid->rr_dot(psi_r+indx2,psi_r+indx2);
       std::cout << "sum4 = " << sum4*scal2*dv << " dv*scal2=" << dv*scal2 <<  std::endl; 

       indx1 += shift1;
       indx2 += shift2;
    }
    */
}

/********************************************
 *                                          *
 *      Electron_Operators::gen_density     *
 *                                          *
 ********************************************/
void Electron_Operators::gen_density(double *dn) {
  /* generate dn */
  mygrid->hr_aSumSqr(scal2, psi_r, dn);
}

/********************************************
 *                                          *
 *      Electron_Operators::gen_densities   *
 *                                          *
 ********************************************/
void Electron_Operators::gen_densities(double *dn, double *dng, double *dnall) {
   /* generate dn */
   mygrid->hr_aSumSqr(scal2, psi_r, dn);
 
   /* generate rho and dng */
   double *tmp = x;
   mygrid->rrr_Sum(dn, dn+(ispin-1)*n2ft3d, rho);
   mygrid->rr_SMul(scal1, rho, tmp);
   // mygrid->rc_fft3d(tmp);
   mygrid->rc_pfft3f(0, tmp);
   mygrid->c_pack(0, tmp);
   mygrid->cc_pack_copy(0, tmp, dng);
 
   /* generate dnall - used for semicore corrections */
   if (mypsp->has_semicore()) {
      for (int ms = 0; ms < ispin; ++ms)
         mygrid->rrr_SMulAdd(0.5, mypsp->semicore_density, dn+ms*n2ft3d, dnall+ms*n2ft3d);
   } else {
      for (int ms = 0; ms < ispin; ++ms)
         mygrid->rr_copy(dn+ms*n2ft3d, dnall+ms*n2ft3d);
   }
}

/********************************************
 *                                          *
 *  Electron_Operators::gen_scf_potentials  *
 *                                          *
 ********************************************/
void Electron_Operators::gen_scf_potentials(double *dn, double *dng, double *dnall)
{
   /* generate coulomb potential */
   if (periodic)
   {
      mycoulomb12->mycoulomb1->vcoulomb(dng, vc);
      mygrid->cc_pack_copy(0,vc,vcall);
   }

   if (aperiodic)
   {
      mygrid->rrr_Sum(dn, &dn[(ispin - 1) * n2ft3d], rho);
      mycoulomb12->mycoulomb2->vcoulomb(rho, vc);

      std::memcpy(vcall,vc,n2ft3d*sizeof(double));
      if (mycoulomb12->dielectric_on())
      {
         double fion_tmp[3];
         mycoulomb12->v_dielectric_aperiodic(rho,dng,vc,vdielec,false,fion_tmp);
         mygrid->rr_Sum(vdielec,vcall);
      }
   }
 
   // generate exchange-correlation potential */
   myxc->v_exc_all(ispin, dnall, xcp, xce);
   // v_exc(ispin,shift2,dnall,xcp,xce,x);
 
   // generate apc potential */
   if (mypsp->myapc->v_apc_on)
   {
      double fion0[1];
      gen_vl_potential();
      mypsp->myapc->V_APC(dng, mypsp->zv, vl, false, fion0);
   }
}

/********************************************
 *                                          *
 *    Electron_Operators::gen_vl_potential  *
 *                                          *
 ********************************************/
void Electron_Operators::gen_vl_potential() {
   double dng0[1], fion0[1];

   /* generate local psp */
   mypsp->v_local(vl, 0, dng0, fion0);

   /* generate long range local psp */
   if (aperiodic)
      mypsp->v_lr_local(vlr_l);
}

/********************************************
 *                                          *
 *      Electron_Operators::gen_vall        *
 *                                          *
 ********************************************/
/**
 * @brief Generates the potential `vall` by combining various potential components.
 *
 * This function calculates the potential `vall` for the electron operators, depending
 * on whether the system is periodic or aperiodic. The function handles the summation
 * of k-space potentials, performs Fourier transforms, and optionally adds contributions
 * from an external electric field. If the system has spin (`ispin == 2`), additional
 * exchange-correlation potentials are added.
 *
 * @details
 * - **Periodic Case**:
 *   - Adds k-space potentials (`vall = scal2 * vl + vc`).
 *   - Performs a Fourier transform on `vall`.
 *   - Adds external electric field (`v_field`) if `efield_on` is true.
 * - **Aperiodic Case**:
 *   - Adds k-space potentials (`vall = scal2 * vsr_l`).
 *   - Performs a Fourier transform on `vall`.
 *   - Adds long-range and core potentials (`vlr_l + vc`) to `vall`.
 *   - Adds external electric field (`v_field`) if `efield_on` is true.
 * - Adds exchange-correlation potential `xcp` to `vall` and handles spin case if applicable.
 *
 * @note The function assumes that the grid operations and potentials are properly initialized
 *       before calling this function.
 *
 * @param None
 * @return void
 */
void Electron_Operators::gen_vall() 
{

   // periodic
   if (periodic)
   {
      // add up k-space potentials, vall = scal2*vl + vc  ****
      mygrid->cc_pack_SMul(0, scal2, vl, vall);
      mygrid->cc_pack_Sum2(0,vc,vall);

      // fourier transform k-space potentials ****
      mygrid->c_unpack(0,vall);
      mygrid->cr_fft3d(vall);
      mygrid->r_zero_ends(vall);

      // add v_field to vall */
      if (mypsp->myefield->efield_on)
         mygrid->rr_Sum(mypsp->myefield->v_field,vall);
   }
   // aperiodic
   else
   {
      /* add up k-space potentials, vall = scal2*v_l */
      mygrid->cc_pack_SMul(0, scal2, vl, vall);
      mygrid->c_unpack(0, vall);
      mygrid->cr_fft3d(vall);
      mygrid->r_zero_ends(vall);

      /* add vall += vlr_l + vc */
      mygrid->rrr_Sum2Add(vlr_l, vc, vall);

      /* add v_field to vall */
      if (mypsp->myefield->efield_on)
         mygrid->rr_Sum(mypsp->myefield->v_field, vall);
   }

   // add xcp to vall 
   if (ispin==2) mygrid->rrr_Sum(vall,xcp+n2ft3d,vall+n2ft3d);
   mygrid->rr_Sum(xcp,vall);
}

/********************************************
 *                                          *
 *      Electron_Operators::get_vall        *
 *                                          *
 ********************************************/
/**
 * @brief Copies the calculated potential `vall` to the output array.
 *
 * This function retrieves the potential `vall` and copies it to the provided
 * output array `vall_out`. If the system has spin (`ispin == 2`), it also copies
 * the spin-dependent part of the potential to the corresponding location in
 * the output array.
 *
 * @param vall_out A pointer to an array where the potential `vall` will be copied.
 *                 The array must be appropriately sized to hold the potential data.
 * 
 * @return void
 */
void Electron_Operators::get_vall(double *vall_out) 
{
   mygrid->rr_copy(vall,vall_out);
   if (ispin==2) mygrid->rr_copy(vall+n2ft3d,vall_out+n2ft3d);
}

/********************************************
 *                                          *
 *      Electron_Operators::set_vall        *
 *                                          *
 ********************************************/
/**
 * @brief Sets the potential `vall` from the input array.
 *
 * This function sets the potential `vall` by copying data from the provided
 * input array `vall_in`. If the system has spin (`ispin == 2`), it also copies
 * the spin-dependent part of the potential from the corresponding location in
 * the input array to the internal `vall` array.
 *
 * @param vall_in A pointer to an array containing the potential data to be set.
 *                The array must be appropriately sized to hold the potential data.
 * 
 * @return void
 */
void Electron_Operators::set_vall(const double *vall_in) 
{
   mygrid->rr_copy(vall_in,vall);
   if (ispin==2) mygrid->rr_copy(vall_in+n2ft3d,vall+n2ft3d);
}

/***********************************************
 *                                             *
 * Electron_Operators::semicore_density_update *
 *                                             *
 ***********************************************/
void Electron_Operators::semicore_density_update() {
   if (mypsp->has_semicore())
      mypsp->semicore_density_update();
}

/********************************************
 *                                          *
 *      Electron_Operators::gen_Hpsi_k      *
 *                                          *
 ********************************************/
void Electron_Operators::gen_Hpsi_k(double *psi) {
   bool move = false;
   double fion0[1];
 
   mygrid->g_zero(Hpsi);
 
   /* get Hpsi */
   if (periodic)
      psi_H(mygrid,myke,mypsp,psi,psi_r,vl,vcall,xcp,Hpsi,move,fion0);
   if (aperiodic)
      psi_Hv4(mygrid,myke,mypsp,psi,psi_r,vl,vlr_l,vcall,xcp,Hpsi,move,fion0);
 
   mygrid->g_Scale(-1.0,Hpsi);
}

/********************************************
 *                                          *
 *      Electron_Operators::gen_hml         *
 *                                          *
 ********************************************/
void Electron_Operators::gen_hml(double *psi, double *hml) {
   mygrid->ggm_sym_Multiply(psi,Hpsi,hml);
}

/********************************************
 *                                          *
 *      Electron_Operators::gen_hmlt        *
 *                                          *
 ********************************************/
void Electron_Operators::gen_hmlt(double *psi, double *hmlt) {
   mygrid->ggm_Multiply(psi, Hpsi, hmlt);
}

/********************************************
 *                                          *
 *      Electron_Operators::get_Tgradient   *
 *                                          *
 ********************************************/
void Electron_Operators::get_Tgradient(double *psi, double *hml, double *THpsi) 
{
   mygrid->fmf_Multiply(-1,psi,hml,1.0,THpsi,0.0);
   mygrid->gg_Minus2(Hpsi,THpsi);
   // mygrid->g_Scale(-1.0,THpsi);
}

/********************************************
 *                                          *
 *      Electron_Operators::gen_Tangent     *
 *                                          *
 ********************************************/
void Electron_Operators::gen_Tangent(double *psi,double *hml,double *THpsi) 
{
   mygrid->fmf_Multiply(-1,psi,hml,1.0,THpsi,-1.0);
}

/********************************************
 *                                          *
 *      Electron_Operators::get_Gradient   *
 *                                          *
 ********************************************/
void Electron_Operators::get_Gradient(double *THpsi) 
{
   mygrid->gg_copy(Hpsi, THpsi);
}

/********************************************
 *                                          *
 *       Electron_Operators::genrho         *
 *                                          *
 ********************************************/
void Electron_Operators::genrho(double *psi, double *dn) 
{
   this->gen_psi_r(psi);
   this->gen_density(dn);
}

/********************************************
 *                                          *
 *       Electron_Operators::run            *
 *                                          *
 ********************************************/
void Electron_Operators::run(double *psi, double *dn, double *dng, double *dnall) 
{
   ++counter;
   this->gen_psi_r(psi);
   // this->gen_density(dn);
   this->gen_densities(dn, dng, dnall);
   this->gen_scf_potentials(dn, dng, dnall);
   this->gen_Hpsi_k(psi);
}

/********************************************
 *                                          *
 *       Electron_Operators::vl_ave         *
 *                                          *
 ********************************************/
double Electron_Operators::vl_ave(double *dng) 
{
   return mygrid->cc_pack_dot(0, dng, vl);
}

/********************************************
 *                                          *
 *       Electron_Operators::vlr_ave        *
 *                                          *
 ********************************************/
double Electron_Operators::vlr_ave(double *dn) 
{
   return (mygrid->rr_dot(dn,vlr_l)+mygrid->rr_dot(&dn[(ispin-1)*n2ft3d],vlr_l))*dv;
}

/********************************************
 *                                          *
 *       Electron_Operators::vnl_ave        *
 *                                          *
 ********************************************/
double Electron_Operators::vnl_ave(double *psi) 
{
   return mypsp->e_nonlocal(psi);
}

/********************************************
 *                                          *
 *       Electron_Operators::eorbit         *
 *                                          *
 ********************************************/
double Electron_Operators::eorbit(double *psi) 
{
   //if (mygrid->d3db::parall->is_master())
   //   std::cout << "Eorbit into ggm_sym_Multiply" << std::endl;

   mygrid->ggm_sym_Multiply(psi,Hpsi,hmltmp);

   //if (mygrid->d3db::parall->is_master())
   //   std::cout << "OUT Eorbit into ggm_sym_Multiply" << std::endl;

   // mygrid->m_scal(-1.0,hmltmp);
   double eorbit0 = mygrid->m_trace(hmltmp);
   if (ispin==1)
      eorbit0 = eorbit0 + eorbit0;

 
   return eorbit0;
}

/********************************************
 *                                          *
 *       Electron_Operators::ehartree       *
 *                                          *
 ********************************************/
double Electron_Operators::ehartree(double *dng) 
{
   return mycoulomb12->mycoulomb1->ecoulomb(dng);
}

/********************************************
 *                                          *
 *       Electron_Operators::ehartree2      *
 *                                          *
 ********************************************/
double Electron_Operators::ehartree2(double *dn) 
{
   return 0.5*(mygrid->rr_dot(dn,vc)+mygrid->rr_dot(&dn[(ispin-1)*n2ft3d],vc))*dv;
}

/********************************************
 *                                          *
 *          Electron_Operators::exc         *
 *                                          *
 ********************************************/
double Electron_Operators::exc(double *dnall) 
{
   double excsum = mygrid->rr_dot(dnall, xce);
   if (ispin==1) {
      excsum *= 2.0;
   } else {
      excsum += mygrid->rr_dot(&dnall[n2ft3d], xce);
   }
   excsum *= dv;
 
   return excsum;
}

/********************************************
 *                                          *
 *          Electron_Operators::pxc         *
 *                                          *
 ********************************************/
double Electron_Operators::pxc(double *dn) 
{
   double pxcsum = mygrid->rr_dot(dn, xcp);
   if (ispin==1) {
      pxcsum *= 2.0;
   } else {
      pxcsum += mygrid->rr_dot(&dn[n2ft3d], &xcp[n2ft3d]);
   }
   pxcsum *= dv;
 
   return pxcsum;
}

/********************************************
 *                                          *
 *         Electron_Operators::eke          *
 *                                          *
 ********************************************/
double Electron_Operators::eke(double *psi) { return myke->ke_ave(psi); }

/********************************************
 *                                          *
 *        Electron_Operators::energy        *
 *                                          *
 ********************************************/
double Electron_Operators::energy(double *psi, double *dn, double *dng, double *dnall) 
{
   double total_energy, eorbit0, ehartr0, exc0, pxc0;
 
   /* total energy calculation */
   mygrid->ggm_sym_Multiply(psi, Hpsi, hmltmp);
   // mygrid->m_scal(-1.0,hmltmp);
   eorbit0 = mygrid->m_trace(hmltmp);
   if (ispin == 1)
      eorbit0 = eorbit0 + eorbit0;
 
   if (periodic)
      ehartr0 = mycoulomb12->mycoulomb1->ecoulomb(dng);
   else
      ehartr0 = 0.5*(mygrid->rr_dot(dn,vc)+mygrid->rr_dot(&dn[(ispin-1)*n2ft3d],vc))*dv;
 
   exc0 = mygrid->rr_dot(dnall, xce);
   pxc0 = mygrid->rr_dot(dn, xcp);
   if (ispin == 1) {
     exc0 = exc0 + exc0;
     pxc0 = pxc0 + pxc0;
   } else {
     exc0 += mygrid->rr_dot(&dnall[n2ft3d], xce);
     pxc0 += mygrid->rr_dot(&dn[n2ft3d], &xcp[n2ft3d]);
   }
   exc0 *= dv;
   pxc0 *= dv;
 
   total_energy = eorbit0 + exc0 - ehartr0 - pxc0;
 
   if (mypsp->myapc->v_apc_on) {
     double eapc = mypsp->myapc->Eapc;
     double papc = mypsp->myapc->Papc;
     total_energy += (eapc - papc);
   }

   /* get dielectric energies */
   if (mycoulomb12->dielectric_on())
   {
      double edielec =  mycoulomb12->edielec;
      double pdielec = mycoulomb12->pdielec;
      total_energy += edielec - pdielec;
   }

 
   return total_energy;
}

/********************************************
 *                                          *
 *   Electron_Operators::gen_energies_en    *
 *                                          *
 ********************************************/
void Electron_Operators::gen_energies_en(double *psi, double *dn, double *dng,
                                         double *dnall, double *E, double *en) 
{
   double total_energy, eorbit0, ehartr0, exc0, pxc0;
 
   /* total energy calculation */
   mygrid->ggm_sym_Multiply(psi, Hpsi, hmltmp);
   // mygrid->m_scal(-1.0,hmltmp);
   eorbit0 = mygrid->m_trace(hmltmp);
   if (ispin == 1) eorbit0 = eorbit0 + eorbit0;
 
   if (periodic) ehartr0 = mycoulomb12->mycoulomb1->ecoulomb(dng);
   if (aperiodic) ehartr0 = 0.5*(mygrid->rr_dot(dn,vc)+mygrid->rr_dot(&dn[(ispin-1)*n2ft3d],vc))*dv;
 
   exc0 = mygrid->rr_dot(dnall,xce);
   pxc0 = mygrid->rr_dot(dn,xcp);
   if (ispin == 1) {
      exc0 = exc0 + exc0;
      pxc0 = pxc0 + pxc0;
   } else {
      exc0 += mygrid->rr_dot(&dnall[n2ft3d], xce);
      pxc0 += mygrid->rr_dot(&dn[n2ft3d], &xcp[n2ft3d]);
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
 
   E[5] = myke->ke_ave(psi);
   E[6] = this->vl_ave(dng);
   if (aperiodic)
     E[6] += this->vlr_ave(dn);
   E[7] = mypsp->e_nonlocal(psi);
   E[8] = 2 * ehartr0;
   E[9] = pxc0;
 
   /* get APC energies */
   if (mypsp->myapc->v_apc_on) {
     E[51] = mypsp->myapc->Eapc;
     E[52] = mypsp->myapc->Papc;
     E[0] = E[0] + E[51] - E[52];
   }

   /* get dielectric energies */
   if (mycoulomb12->dielectric_on())
   {
      E[61] = mycoulomb12->edielec;
      E[62] = mycoulomb12->pdielec;
      E[0] = E[0] + E[61] - E[62];
   }
 
   /* get Efield energies */
   if (mypsp->myefield->efield_on) {
     if (mypsp->myefield->efield_type == 0) {
       E[48] = 0.0;
       E[49] = 0.0;
       E[0] += E[48] - E[49];
     } else {
       E[48] = dv * mygrid->rr_dot(rho, mypsp->myefield->v_field);
       E[49] = mypsp->myefield->efield_ion_energy();
       E[0] += E[49];
     }
   }
 
   en[0] = dv * mygrid->r_dsum(dn);
   en[1] = en[0];
   if (ispin > 1) en[1] = dv*mygrid->r_dsum(&dn[n2ft3d]);
}

/********************************************
 *                                          *
 *     Electron_Operators::add_dteHpsi      *
 *                                          *
 ********************************************/
void Electron_Operators::add_dteHpsi(double dte, double *psi1, double *psi2) 
{
   /* do a steepest descent step */
   mygrid->gg_SMul(dte, Hpsi, psi2);
   mygrid->gg_Sum2(psi1, psi2);
}

/********************************************
 *                                          *
 *    Electron_Operators::semicore_force    *
 *                                          *
 ********************************************/
void Electron_Operators::semicore_force(double *fion) 
{
   if (mypsp->has_semicore())
      mypsp->semicore_xc_fion(xcp, fion);
}

/********************************************
 *                                          *
 *      Electron_Operators::vl_force        *
 *                                          *
 ********************************************/
void Electron_Operators::vl_force(double *dng, double *fion) 
{
   mypsp->f_local(dng, fion);
}

/********************************************
 *                                          *
 *      Electron_Operators::vlr_force       *
 *                                          *
 ********************************************/
void Electron_Operators::vlr_force(double *dn, double *fion) 
{
   if (aperiodic) 
   {
      mygrid->rrr_Sum(dn, &dn[(ispin-1)*n2ft3d],rho);
      mypsp->grad_v_lr_local(rho,fion);
   }
}

/********************************************
 *                                          *
 *      Electron_Operators::apc_force       *
 *                                          *
 ********************************************/
void Electron_Operators::apc_force(double *dng, double *fion) 
{
   // generate apc force */
   if (mypsp->myapc->v_apc_on) {
      mypsp->myapc->f_APC(dng, mypsp->zv, fion);
   }
}

/********************************************
 *                                          *
 *      Electron_Operators::dielectric_force       *
 *                                          *
 ********************************************/
void Electron_Operators::dielectric_force(double *fion)
{
   // generate apc force */
   if (mycoulomb12->dielectric_on())
   {
      mycoulomb12->dielectric_fion(fion);
   }
}


/********************************************
 *                                          *
 *      Electron_Operators::vnl_force       *
 *                                          *
 ********************************************/
void Electron_Operators::vnl_force(double *psi, double *fion) 
{
   mypsp->f_nonlocal_fion(psi, fion);
}

} // namespace pwdft

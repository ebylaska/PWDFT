
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
void cElectron_Operators::gen_density(double *dn) 
{
   /* generate dn */
   mygrid->hr_aSumSqr(scal2, psi_r, dn);
}

/********************************************
 *                                          *
 *     cElectron_Operators::gen_densities   *
 *                                          *
 ********************************************/
void cElectron_Operators::gen_densities(double *dn, double *dng, double *dnall) 
{
   /* generate dn */
   mygrid->hr_aSumSqr(scal2, psi_r, dn);
 
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
void cElectron_Operators::gen_scf_potentials(double *dn, double *dng, double *dnall)
{
   /* generate coulomb potential */
   mycoulomb->vcoulomb(dng, vc);
   mygrid->cc_pack_copy(0,vc,vcall);

 
   // generate exchange-correlation potential */
   myxc->v_exc_all(ispin, dnall, xcp, xce);
   // v_exc(ispin,shift2,dnall,xcp,xce,x);
 
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
void cElectron_Operators::gen_Hpsi_k(double *psi) 
{
   bool move = false;
   double fion0[1];
 
   mygrid->g_zero(Hpsi);
 
   /* get Hpsi */
   cpsi_H(mygrid,myke,mypsp,psi,psi_r,vl,vcall,xcp,Hpsi,move,fion0);
 
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
void cElectron_Operators::get_Gradient(double *THpsi) 
{
   mygrid->gg_copy(Hpsi, THpsi);
}

/********************************************
 *                                          *
 *      cElectron_Operators::genrho         *
 *                                          *
 ********************************************/
void cElectron_Operators::genrho(double *psi, double *dn) 
{
   this->gen_psi_r(psi);
   this->gen_density(dn);
}

/********************************************
 *                                          *
 *      cElectron_Operators::run            *
 *                                          *
 ********************************************/
void cElectron_Operators::run(double *psi, double *dn, double *dng, double *dnall) 
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
 *       cElectron_Operators::vl_ave        *
 *                                          *
 ********************************************/
double cElectron_Operators::vl_ave(double *dng) 
{
   return mygrid->cc_pack_dot(0, dng, vl);
}


/********************************************
 *                                          *
 *      cElectron_Operators::vnl_ave        *
 *                                          *
 ********************************************/
double cElectron_Operators::vnl_ave(double *psi) 
{
   return mypsp->e_nonlocal(psi);
}

/********************************************
 *                                          *
 *      cElectron_Operators::eorbit         *
 *                                          *
 ********************************************/
double cElectron_Operators::eorbit(double *psi) 
{

   mygrid->ggw_sym_Multiply(psi,Hpsi,hmltmp);


   // mygrid->m_scal(-1.0,hmltmp);
   double eorbit0 = mygrid->w_trace(hmltmp);
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
double cElectron_Operators::eke(double *psi) { return myke->ke_ave(psi); }

/********************************************
 *                                          *
 *        cElectron_Operators::energy       *
 *                                          *
 ********************************************/
double cElectron_Operators::energy(double *psi, double *dn, double *dng, double *dnall) 
{
   double total_energy, eorbit0, ehartr0, exc0, pxc0;
 
   /* total energy calculation */
   mygrid->ggw_sym_Multiply(psi, Hpsi, hmltmp);
   // mygrid->m_scal(-1.0,hmltmp);
   eorbit0 = mygrid->w_trace(hmltmp);
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
                                          double *dnall, double *E, double *en) 
{
   double total_energy, eorbit0, ehartr0, exc0, pxc0;
 
   /* total energy calculation */
   mygrid->ggw_sym_Multiply(psi, Hpsi, hmltmp);
   // mygrid->m_scal(-1.0,hmltmp);
   eorbit0 = mygrid->w_trace(hmltmp);
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
 
   E[5] = myke->ke_ave(psi);
   E[6] = this->vl_ave(dng);
   E[7] = mypsp->e_nonlocal(psi);
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

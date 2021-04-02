
#include	"blas.h"
#include	"Pneb.hpp"
#include	"Ion.hpp"
#include	"Kinetic.hpp"
#include	"Coulomb.hpp"
#include        "Strfac.hpp"
#include	"Pseudopotential.hpp"
#include        "v_exc.hpp"

#include	"psi_H.hpp"

#include	"Electron.hpp"

/********************************************
 *                                          *
 *  Electron_Operators::Electron_Operators  *
 *                                          *
 ********************************************/
Electron_Operators::Electron_Operators(Pneb *mygrid0, Kinetic_Operator *myke0, Coulomb_Operator *mycoulomb0,  Pseudopotential *mypsp0)
{
   mygrid  = mygrid0;
   myke      = myke0;
   mycoulomb = mycoulomb0;
   mypsp     = mypsp0;

   ispin = mygrid->ispin;
   neall = mygrid->neq[0] + mygrid->neq[1];

   /* allocate memory */
   Hpsi  = mygrid->g_allocate(1);
   psi_r = mygrid->h_allocate();
   xcp   = mygrid->r_nalloc(ispin);
   xce   = mygrid->r_nalloc(ispin);

   x   = mygrid->r_alloc();
   vall= mygrid->r_alloc();
   vl  = mygrid->c_pack_allocate(0);
   vc  = mygrid->c_pack_allocate(0);

   hmltmp = mygrid->m_allocate(-1,1);

   omega = mygrid->lattice->omega();
   scal1 = 1.0/((double) ((mygrid->nx)*(mygrid->ny)*(mygrid->nz)));
   scal2 = 1.0/omega;
   dv = omega*scal1;

   n2ft3d = (mygrid->n2ft3d);
   shift1 = 2*(mygrid->npack(1));
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
   //mygrid->gh_fftb(psi1,psi_r);
   int indx1 = 0;
   int indx2 = 0;
   for (int i=0; i<neall; ++i)
   {
      mygrid->cc_pack_copy(1,&psi[indx1],&psi_r[indx2]);
      mygrid->c_unpack(1,&psi_r[indx2]);
      mygrid->cr_fft3d(&psi_r[indx2]);
      indx1 += shift1;
      indx2 += shift2;
   }

}

/********************************************
 *                                          *
 *      Electron_Operators::gen_density     *
 *                                          *
 ********************************************/
void Electron_Operators::gen_density(double *dn)
{
   /* generate dn */
   mygrid->hr_aSumSqr(scal2,psi_r,dn);
}

/********************************************
 *                                          *
 *      Electron_Operators::gen_densities   *
 *                                          *
 ********************************************/
void Electron_Operators::gen_densities(double *dn, double *dng, double *dnall)
{

   /* generate dn */
   mygrid->hr_aSumSqr(scal2,psi_r,dn);

   /* generate dng */
   double *tmp = x;
   mygrid->rrr_Sum(dn,&dn[(ispin-1)*n2ft3d],tmp);
   mygrid->r_SMul(scal1,tmp);
   mygrid->rc_fft3d(tmp);
   mygrid->c_pack(0,tmp);
   mygrid->cc_pack_copy(0,tmp,dng);

   /* generate dnall - used for semicore corrections */
   if (mypsp->has_semicore())
   {
      for (int ms=0; ms<ispin; ++ms)
         mygrid->rrr_SMulAdd(0.5,mypsp->semicore_density,&dn[ms*n2ft3d],&dnall[ms*n2ft3d]);
   }
   else
   {
       for (int ms=0; ms<ispin; ++ms)
          mygrid->rr_copy(&dn[ms*n2ft3d],&dnall[ms*n2ft3d]);
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
   mycoulomb->vcoulomb(dng,vc);

   /* generate exchange-correlation potential */
   v_exc(ispin,shift2,dnall,xcp,xce,x);

}

/********************************************
 *                                          *
 *    Electron_Operators::gen_vl_potential  *
 *                                          *
 ********************************************/
void Electron_Operators::gen_vl_potential()
{
   double dng0[1], fion0[1];

   /* generate local psp*/
   mypsp->v_local(vl,0,dng0,fion0);
}

/***********************************************
 *                                             *
 * Electron_Operators::semicore_density_update *
 *                                             *
 ***********************************************/
void Electron_Operators::semicore_density_update()
{
   if (mypsp->has_semicore())
      mypsp->semicore_density_update();
}

/********************************************
 *                                          *
 *      Electron_Operators::gen_Hpsi_k      *
 *                                          *
 ********************************************/
void Electron_Operators::gen_Hpsi_k(double *psi)
{
   bool move=false;
   double fion0[1];

   mygrid->g_zero(Hpsi);

   psi_H(mygrid,myke,mypsp,
         psi, psi_r, vl, vc, xcp,
         Hpsi, move,fion0);

   mygrid->g_Scale(-1.0,Hpsi);
}

/********************************************
 *                                          *
 *      Electron_Operators::gen_hml         *
 *                                          *
 ********************************************/
void Electron_Operators::gen_hml(double *psi, double *hml) 
{
   mygrid->ggm_sym_Multiply(psi,Hpsi,hml);
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
   //mygrid->g_Scale(-1.0,THpsi);
}

/********************************************
 *                                          *
 *      Electron_Operators::gen_Tangent     *
 *                                          *
 ********************************************/
void Electron_Operators::gen_Tangent(double *psi, double *hml, double *THpsi)
{
   mygrid->fmf_Multiply(1,psi,hml,1.0,THpsi,-1.0);
}

/********************************************
 *                                          *
 *      Electron_Operators::get_Gradient   *
 *                                          *
 ********************************************/
void Electron_Operators::get_Gradient(double *THpsi)
{
   mygrid->gg_copy(Hpsi,THpsi);
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
   
   this->gen_psi_r(psi);
   //this->gen_density(dn);
   this->gen_densities(dn,dng,dnall);
   this->gen_scf_potentials(dn,dng,dnall);
   this->gen_Hpsi_k(psi);
}


/********************************************
 *                                          *
 *       Electron_Operators::vl_ave         *
 *                                          *
 ********************************************/
double Electron_Operators::vl_ave(double *dng)
{
   return mygrid->cc_pack_dot(0,dng,vl);
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
   mygrid->ggm_sym_Multiply(psi,Hpsi,hmltmp);
   //mygrid->m_scal(-1.0,hmltmp);
   double eorbit0 = mygrid->m_trace(hmltmp);
   if (ispin==1) eorbit0 = eorbit0+eorbit0;

   return eorbit0;
}

/********************************************
 *                                          *
 *       Electron_Operators::ehartree       *
 *                                          *
 ********************************************/
double Electron_Operators::ehartree(double *dng)
{
   return mycoulomb->ecoulomb(dng);
}
  

/********************************************
 *                                          *
 *          Electron_Operators::exc         *
 *                                          *
 ********************************************/
double Electron_Operators::exc(double *dnall)
{
   double excsum     = mygrid->rr_dot(dnall,xce);
   if (ispin==1)
   {
      excsum = excsum+excsum;
   }
   else
   {
      excsum += mygrid->rr_dot(&dnall[n2ft3d],xce);
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
   double pxcsum = mygrid->rr_dot(dn,xcp);
   if (ispin==1)
   {
      pxcsum = pxcsum+pxcsum;
   }
   else
   {
      pxcsum += mygrid->rr_dot(&dn[n2ft3d],&xcp[n2ft3d]);
   }
   pxcsum *= dv;

   return pxcsum;
}

/********************************************
 *                                          *
 *         Electron_Operators::eke          *
 *                                          *
 ********************************************/
double Electron_Operators::eke(double *psi)
{
   return myke->ke_ave(psi);
}


/********************************************
 *                                          *
 *        Electron_Operators::energy        *
 *                                          *
 ********************************************/
double Electron_Operators::energy(double *psi, double *dn, double *dng, double *dnall)
{
   double total_energy, eorbit0, ehartr0, exc0, pxc0;

   /* total energy calculation */
   mygrid->ggm_sym_Multiply(psi,Hpsi,hmltmp);
   //mygrid->m_scal(-1.0,hmltmp);
   eorbit0  = mygrid->m_trace(hmltmp);
   if (ispin==1) eorbit0 = eorbit0+eorbit0;

   ehartr0 = mycoulomb->ecoulomb(dng);
   exc0    = mygrid->rr_dot(dnall,xce);
   pxc0    = mygrid->rr_dot(dn,xcp);
   if (ispin==1)
   {
      exc0 = exc0+exc0;
      pxc0 = pxc0+pxc0;
   }
   else
   {
      exc0 += mygrid->rr_dot(&dnall[n2ft3d],xce);
      pxc0 += mygrid->rr_dot(&dn[n2ft3d],&xcp[n2ft3d]);
   }
   exc0 *= dv;
   pxc0 *= dv;

   total_energy = eorbit0 + exc0 - ehartr0 - pxc0;
  
   return total_energy;
}


/********************************************
 *                                          *
 *   Electron_Operators::gen_energies_en    *
 *                                          *
 ********************************************/
void Electron_Operators::gen_energies_en(double *psi, double *dn, double *dng, double *dnall,
                                           double *E, double *en)
{
   double total_energy, eorbit0, ehartr0, exc0, pxc0;

   /* total energy calculation */
   mygrid->ggm_sym_Multiply(psi,Hpsi,hmltmp);
   //mygrid->m_scal(-1.0,hmltmp);
   eorbit0  = mygrid->m_trace(hmltmp);
   if (ispin==1) eorbit0 = eorbit0+eorbit0;

   ehartr0 = mycoulomb->ecoulomb(dng);
   exc0    = mygrid->rr_dot(dnall,xce);
   pxc0    = mygrid->rr_dot(dn,xcp);
   if (ispin==1)
   {
      exc0 = exc0+exc0;
      pxc0 = pxc0+pxc0;
   }
   else
   {
      exc0 += mygrid->rr_dot(&dnall[n2ft3d],xce);
      pxc0 += mygrid->rr_dot(&dn[n2ft3d],&xcp[n2ft3d]);
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

   E[5] =  myke->ke_ave(psi);
   E[6] = mygrid->cc_pack_dot(0,dng,vl);
   E[7] = mypsp->e_nonlocal(psi);
   E[8] = 2*ehartr0;
   E[9] = pxc0;


   en[0] = dv*mygrid->r_dsum(dn);
   en[1] = en[0];
   if (ispin > 1)
      en[1] =  dv*mygrid->r_dsum(&dn[n2ft3d]);

}


/********************************************
 *                                          *
 *     Electron_Operators::add_dteHpsi      *
 *                                          *
 ********************************************/
void Electron_Operators::add_dteHpsi(double dte, double *psi1, double *psi2)
{
   /* do a steepest descent step */
   mygrid->gg_SMul(dte,Hpsi,psi2);
   mygrid->gg_Sum2(psi1,psi2);
}





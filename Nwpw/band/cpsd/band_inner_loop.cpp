

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "Control2.hpp"
#include "cCoulomb.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "cKinetic.hpp"
#include "CGrid.hpp"
#include "Parallel.hpp"
#include "CPseudopotential.hpp"
#include "cExchange_Correlation.hpp"
#include "band_inner_loop.hpp"
#include "iofmt.hpp"
#include "cpsi_H.hpp"

//#include "nwpw_dplot.hpp"

namespace pwdft {

void band_inner_loop(Control2 &control, Cneb *mygrid, Ion *myion,
                     cKinetic_Operator *myke, cCoulomb_Operator *mycoulomb,
                     cXC_Operator *myxc, CPseudopotential *mypsp, CStrfac *mystrfac,
                     Ewald *myewald, double *psi1, double *psi2, double *Hpsi,
                     double *psi_r, double *dn, double *hml, double *lmbda,
                     double E[], double *deltae, double *deltac, double *deltar) 
{
   int it_in, k, ms;
   int indx1, indx2;
   int one = 1;
   double dc;
   double eorbit, eion, econstraint, exc, ehartr, pxc;
   double eke, elocal, enlocal, dt, dte, Eold;
   double *vl,*vc,*xcp,*xce,*dnall,*x,*dng,*rho,*tmp,*vcall,*vpsi,*sumi;
   double *fion;
   bool move = control.geometry_optimize();
   double omega = mygrid->lattice->omega();
 
   int ispin = mygrid->ispin;
   int neall = mygrid->neq[0] + mygrid->neq[1];
   int shift1 = 2*mygrid->npack1_max(); 
   int shift2 = (mygrid->n2ft3d);
   int nfft3d = (mygrid->nfft3d);
   int n2ft3d = (mygrid->n2ft3d);

   double scal1 = 1.0/((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
   double scal2 = 1.0/omega;
   double dv = omega*scal1;
 
   dt = control.time_step();
   dte = dt/sqrt(control.fake_mass());
   it_in = control.loop(0);
   //allocate temporary memory 
   rho = mygrid->c_alloc();
   tmp = mygrid->c_alloc();

   xcp = mygrid->r_nalloc(ispin);
   xce = mygrid->r_nalloc(ispin);
   dnall = mygrid->r_nalloc(ispin);
   x     = mygrid->r_alloc();

   dng   = mygrid->c_pack_allocate(0);
   vl    = mygrid->c_pack_allocate(0);
   vc    = mygrid->c_pack_allocate(0);
   vcall = mygrid->c_pack_allocate(0);
 
   vpsi = x;
 
   //fion = new double[3 * (myion->nion)]();
   fion = myion->fion1;
 
   // generate local psp
   mypsp->v_local(vl,false,dng,fion);

 
   // myewald->phafac();
 
   //|-\____|\/-----\/\/->    Start Parallel Section    <-\/\/-----\/|____/-|
 
   for (auto it=0; it<it_in; ++it) 
   {
      mygrid->g_zero(Hpsi);
      mygrid->gg_copy(psi2, psi1);
     
      if (move)
      {
         myion->shift();
         mystrfac->phafac();
         mystrfac->phafac_k();
         myewald->phafac();
         // for (int ii=0; ii<(3*myion->nion); ++ii) fion[ii] = 0.0;
      }
     
      // convert psi(G) to psi(r) - Expensive 
      mygrid->gh_fftb(psi1,psi_r);

      // generate dn
      mygrid->hr_aSumSqr(scal2,psi_r,dn);

      // generate dng 
      mygrid->rrc_Sum(dn,dn+(ispin-1)*nfft3d,rho);

      mygrid->rc_pfft3f(0,rho);
      //mygrid->rc_fft3d(rho);

      mygrid->cc_SMul(scal1, rho, tmp);

      mygrid->c_pack(0,tmp);
      mygrid->cc_pack_copy(0,tmp,dng);
      //mygrid->c_pack_SMul(0,scal1, dng);

      // generate dnall - used for semicore corrections
      if (mypsp->has_semicore()) 
      {
         if ((move) || (it == 0))
            mypsp->semicore_density_update();
         for (ms=0; ms<ispin; ++ms)
            mygrid->rrr_SMulAdd(0.5,mypsp->semicore_density,dn+ms*nfft3d,dnall+ms*nfft3d);
      } 
      else 
      {
         for (ms=0; ms<ispin; ++ms)
            mygrid->rr_copy(dn+ms*nfft3d, dnall+ms*nfft3d);
      }
     
      // generate local potentials
      if (move) 
      {
         mypsp->v_local(vl,move,dng,fion);
      }
     
      // apply k-space operators
      // myke->ke(psi1,Hpsi);
     
      // apply non-local PSP  - Expensive
      // mypsp->v_nonlocal_fion(psi1,Hpsi,move,fion);
     
      // generate coulomb potential 
      mycoulomb->vcoulomb(dng,vc);
      mygrid->cc_pack_copy(0,vc,vcall);

      // generate exchange-correlation potential 
      std::memset(xcp,0,ispin*nfft3d*sizeof(double));
      std::memset(xce,0,ispin*nfft3d*sizeof(double));
      myxc->v_exc_all(ispin,dnall,xcp,xce);
     
      // get Hpsi
      cpsi_H(mygrid,myke,mypsp,psi1,psi_r,vl,vcall,xcp,Hpsi,move,fion);

      // do a steepest descent step
      mygrid->gg_SMul(dte,Hpsi,psi2);
      mygrid->gg_Sum2(psi1,psi2);

      if (move)
      {
         // get the ion-ion force 
         myewald->force(fion);
        
         // get the semicore force - needs to be checked 
         if (mypsp->has_semicore())
            mypsp->semicore_xc_fion(xcp, fion);
        
        
         // steepest descent step 
         myion->add_contraint_force(fion);

         // steepest descent step 
         myion->optimize_step(fion);
      }

      // lagrange multiplier - Expensive 
      mygrid->ggw_lambda(dte, psi1, psi2, lmbda);
   }

   //|-\____|\/-----\/\/->    End Parallel Section    <-\/\/-----\/|____/-|

   // total energy calculation 
   mygrid->ggw_sym_Multiply(psi1, Hpsi, hml);

   mygrid->w_scal(-1.0, hml);
   eorbit = mygrid->w_trace(hml);
   if (ispin==1) eorbit = eorbit + eorbit;

   eorbit = mygrid->c3db::parall->SumAll(3,eorbit);

   // hartree energy and ion-ion energy 
   ehartr = mycoulomb->ecoulomb(dng);
   eion = myewald->energy();
   
 
   // xc energy
   exc = mygrid->rr_dot(dnall, xce);
   pxc = mygrid->rr_dot(dn, xcp);
   if (ispin == 1) {
      exc = exc + exc;
      pxc = pxc + pxc;
   } else {
      exc += mygrid->rr_dot(dnall+nfft3d, xce);
      pxc += mygrid->rr_dot(dn+nfft3d, xcp+nfft3d);
   }
   exc *= dv;
   pxc *= dv;
 
   // average Kohn-Sham kineticl energy 
   eke = myke->ke_ave(psi1);
 
   // average Kohn-Sham local psp energy
   elocal = mygrid->cc_pack_dot(0, dng, vl);
 
 
   // average Kohn-Sham v_nonlocal energy 
   enlocal = mypsp->e_nonlocal(psi1);
 
   Eold = E[0];
   E[0] = eorbit + eion + exc - ehartr - pxc;
   E[1] = eorbit;
   E[2] = ehartr;
   E[3] = exc;
   E[4] = eion;
   E[5] = eke;
   E[6] = elocal;
   E[7] = enlocal;
   E[8] = 2 * ehartr;
   E[9] = pxc;

   // get contraints energies
   if (myion->has_ion_bond_constraints()) 
   {
      E[70] = myion->energy_ion_bond_constraints();
      E[0] = E[0] + E[70];
   }
   if (myion->has_ion_bondings_constraints()) 
   {
      E[71] = myion->energy_ion_bondings_constraints();
      E[0] = E[0] + E[71];
   }

 
   // set convergence variables
   *deltae = (E[0] - Eold) / (dt * control.loop(0));
 
   // deltac 
   dc = 0.0;
   mygrid->ggg_Minus(psi2, psi1, Hpsi);
   for (auto nbq=0; nbq<mygrid->nbrillq; ++nbq)
   {
      int nbq1 = nbq+1;
      sumi = new double[neall]();
      for (auto i=0; i<neall; ++i)
         sumi[i] = mygrid->cc_pack_idot(nbq1, 
                                        Hpsi+(i+nbq*neall)*shift1, 
                                        Hpsi+(i+nbq*neall)*shift1);
      
      mygrid->c3db::parall->Vector_SumAll(1, neall, sumi);

      for (auto i = 0; i < neall; ++i)
         if (sumi[i] > dc)
            dc = sumi[i];
      delete[] sumi;
   }
   dc = mygrid->c3db::parall->MaxAll(3, dc);
   *deltac = dc / dte;

 
   // deltar 
   *deltar = 0.0;
   if (move) {
     double sum;
     for (auto ii = 0; ii < (myion->nion); ++ii) {
       sum = sqrt(fion[3*ii]  *fion[3*ii] +
                  fion[3*ii+1]*fion[3*ii+1] +
                  fion[3*ii+2]*fion[3*ii+2]);
       if (sum > (*deltar))
         *deltar = sum;
     }
   }
 
   //delete[] fion;
 
   mygrid->c_dealloc(tmp);
   mygrid->r_dealloc(xcp);
   mygrid->r_dealloc(xce);
   mygrid->r_dealloc(dnall);
   mygrid->r_dealloc(x);
   mygrid->c_dealloc(rho);
   mygrid->c_pack_deallocate(dng);
   mygrid->c_pack_deallocate(vl);

   mygrid->c_pack_deallocate(vc);
   mygrid->c_pack_deallocate(vcall);

}
} // namespace pwdft

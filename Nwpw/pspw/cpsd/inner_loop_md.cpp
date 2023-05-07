
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

#include "Control2.hpp"
#include "Coulomb12.hpp"
#include "Ewald.hpp"
#include "Ion.hpp"
#include "Kinetic.hpp"
#include "PGrid.hpp"
#include "Parallel.hpp"
#include "Pseudopotential.hpp"
#include "exchange_correlation.hpp"
#include "nwpw_Nose_Hoover.hpp"
#include "psi_H.hpp"
//#include	"v_exc.hpp"
#include "inner_loop_md.hpp"

namespace pwdft {

void inner_loop_md(const bool verlet, double *sa_alpha, Control2 &control,
                   Pneb *mygrid, Ion *myion, nwpw_Nose_Hoover *mynose,
                   Kinetic_Operator *myke, Coulomb12_Operator *mycoulomb12,
                   XC_Operator *myxc, Pseudopotential *mypsp, Strfac *mystrfac,
                   Ewald *myewald, double *psi0, double *psi1, double *psi2,
                   double *Hpsi, double *psi_r, double *dn, double *hml,
                   double *lmbda, const int it_in, double E[]) 
{
   int n2ft3d,neall,ispin;
   int shift1, shift2, indx1, indx2;
   int one = 1;
   double scal1,scal2,dv,dc;
   double eorbit,eion,exc,ehartr,pxc;
   double eke,eki,elocal,enlocal,dt,dte,fmass,Eold;
   double *vl,*vlr_l,*vc,*vcall,*xcp,*xce,*dnall,*x,*dng,*rho,*tmp,*vpsi,*sumi;
   double *vdielec;
   double *fion;
   bool move = true;
   bool nose = mynose->on();
   bool periodic  = (control.version==3);
   bool aperiodic = (control.version==4);
   double sse = 1.0;
   double ssr = 1.0;
   double omega = mygrid->lattice->omega();
 
   ispin = mygrid->ispin;
   neall = mygrid->neq[0] + mygrid->neq[1];
   shift1 = 2 * (mygrid->npack(1));
   shift2 = (mygrid->n2ft3d);
   n2ft3d = (mygrid->n2ft3d);
   scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
   // scal2 = 1.0/lattice_omega();
   // dv = lattice_omega()*scal1;
   scal2 = 1.0 / omega;
   dv = omega * scal1;
 
   double sa1 = 1.0;
   double sa2 = 1.0;
   if (!nose) {
     sa1 = 1.0 / (2.0 - sa_alpha[0]);
     sa2 = sa_alpha[0] / (2.0 - sa_alpha[0]);
   }
 
   dt = control.time_step();
   fmass = control.fake_mass();
   dte = dt * dt / fmass;
   if (!verlet) dte = 0.5 * dte;
 
   /* allocate temporary memory */
   rho = mygrid->r_alloc();
   tmp = mygrid->r_alloc();
   xcp = mygrid->r_nalloc(ispin);
   xce = mygrid->r_nalloc(ispin);
   dnall = mygrid->r_nalloc(ispin);
   x = mygrid->r_alloc();
   dng = mygrid->c_pack_allocate(0);
   vl = mygrid->c_pack_allocate(0);
   if (periodic)
   {
      vc = mygrid->c_pack_allocate(0);
      vcall = mygrid->c_pack_allocate(0);
      if (mycoulomb12->dielectric_on())
         vdielec = mygrid->c_pack_allocate(0);
   }
   else if (aperiodic) 
   {
      vc    = mygrid->r_alloc();
      vcall = mygrid->r_alloc();
      vlr_l = mygrid->r_alloc();
      if (mycoulomb12->dielectric_on())
         vdielec = mygrid->r_alloc();
   }
   vpsi = x;
 
   // new double[3*(myion->nion)]();
   fion = myion->fion1;
 
   /* generate local psp*/
   // mypsp->v_local(vl,0,dng,fion);
 
   // myewald->phafac();
 
   //|-\____|\/-----\/\/->    Start Parallel Section    <-\/\/-----\/|____/-|
 
   for (auto it=0; it<it_in; ++it) 
   {
      /* shift wavefuntion */
      mygrid->g_zero(Hpsi);
      mygrid->gg_copy(psi1,psi0);
      mygrid->gg_copy(psi2,psi1);
      // mygrid->gh_fftb(psi1,psi_r);
      
      /* skip shift ion if newton step */
      if (verlet) myion->shift();
      if (nose && verlet) mynose->shift();
      
      mystrfac->phafac();
      myewald->phafac();
      
      indx1 = 0;
      indx2 = 0;
      for (auto i=0; i<neall; ++i) 
      {
         mygrid->cc_pack_copy(1,&psi1[indx1],&psi_r[indx2]);
         mygrid->c_unpack(1,&psi_r[indx2]);
         mygrid->cr_fft3d(&psi_r[indx2]);
         indx1 += shift1;
         indx2 += shift2;
      }
      
      /* generate dn */
      mygrid->hr_aSumSqr(scal2, psi_r, dn);
      
      /* generate dng */
      mygrid->rrr_Sum(dn,&dn[(ispin-1)*n2ft3d],rho);
      mygrid->rr_SMul(scal1,rho,tmp);
      mygrid->rc_fft3d(tmp);
      mygrid->c_pack(0, tmp);
      mygrid->cc_pack_copy(0, tmp, dng);
      
      /* generate dnall - used for semicore corrections */
      if (mypsp->has_semicore()) 
      {
         mypsp->semicore_density_update();
         for (auto ms=0; ms<ispin; ++ms)
            mygrid->rrr_SMulAdd(0.5,mypsp->semicore_density,&dn[ms*n2ft3d],&dnall[ms*n2ft3d]);
      } 
      else 
      {
         for (auto ms=0; ms<ispin; ++ms)
            mygrid->rr_copy(&dn[ms*n2ft3d], &dnall[ms*n2ft3d]);
      }
      
      /* generate local potential */
      if (move) mypsp->v_local(vl,move,dng,fion);
      if (mypsp->myapc->v_apc_on) mypsp->myapc->V_APC(dng,mypsp->zv,vl,move,fion);
      
      /* long-range psp for charge systems */
      if (aperiodic) {
         mypsp->v_lr_local(vlr_l);
         if (move)
           mypsp->grad_v_lr_local(rho,fion);
      }
      
      /* apply k-space operators */
      // myke->ke(psi1,Hpsi);
      // mypsp->v_nonlocal_fion(psi1,Hpsi,move,fion);
      
      /* generate coulomb potential */
      if (periodic)
      {
         mycoulomb12->mycoulomb1->vcoulomb(dng, vc);
         mygrid->cc_pack_copy(0,vc,vcall);
      }
      else if (aperiodic)
      {
         mycoulomb12->mycoulomb2->vcoulomb(rho,vc);
         std::memcpy(vcall,vc,n2ft3d*sizeof(double));
         if (mycoulomb12->dielectric_on())
         {
            mycoulomb12->v_dielectric_aperiodic(rho,dng,vc,vdielec,move,fion);
            mygrid->rr_Sum(vdielec,vcall);
         }
      } 
      
      /* generate exchange-correlation potential */
      myxc->v_exc_all(ispin, dnall, xcp, xce);
      
      /* get Hpsi */
      if (periodic)
         psi_H(mygrid,myke,mypsp,psi1,psi_r,vl,vcall,xcp,Hpsi,move,fion);
      else if (aperiodic)
         psi_Hv4(mygrid,myke,mypsp,psi1,psi_r,vl,vlr_l,vcall,xcp,Hpsi,move,fion);
      
      
      /* get the ion-ion force */
      if (periodic)
         myewald->force(fion); /* get the ewald force */
      else if (aperiodic)
         myion->ion_ion_force(fion); /* get the ion-ion force */
      
      /* get the semicore force - needs to be checked */
      if (mypsp->has_semicore()) mypsp->semicore_xc_fion(xcp,fion);
      
      /* get forces from external Efield */
      if (mypsp->myefield->efield_on)
         mypsp->myefield->efield_ion_fion(fion);
      
      /* car-parrinello Verlet step */
      if (verlet) 
      {
         /* constant temperature */
         if (nose) 
         {
            sse = mynose->sse();
            ssr = mynose->ssr();
           
            mygrid->gg_SMul(0.5 * dte, Hpsi, psi2);
            mygrid->gg_daxpy(-1.0, psi0, psi2);
            mygrid->gg_daxpy(1.0, psi1, psi2);
           
            mygrid->g_Scale(2.0 * sse, psi2);
            mygrid->gg_daxpy(1.0, psi0, psi2);
           
            // myion->Verlet_step(ssr,fion);
            myion->Nose_step(ssr, fion);
         }
         /* constant energy */
         else 
         {
            mygrid->gg_SMul(dte * sa1, Hpsi, psi2);
            mygrid->gg_daxpy(-1.0 * sa2, psi0, psi2);
            mygrid->gg_daxpy(2.0 * sa1, psi1, psi2);
           
            myion->Verlet_step(fion, sa_alpha[1]);
         }
      }
 
      /* car-parrinello Newton step */
      else 
      {
         double r = 1.0;
         double s = 1.0;
         /* constant temperature */
         if (nose) {
           r = (1.0-0.5*dt*mynose->dXr());
           s = (1.0-0.5*dt*mynose->dXe());
         }
        
         mygrid->gg_SMul(dte,Hpsi,psi2);
         mygrid->gg_daxpy(s*dt*sa_alpha[0],psi0,psi2);
         mygrid->gg_Sum2(psi1,psi2);
        
         myion->Newton_step(fion,sa_alpha[1]*r);
      }
 
      /* lagrange multiplier */
      double dte0 = dte;
      if (nose && verlet) dte0 *= sse;
      mygrid->ggm_lambda(dte0, psi1, psi2, lmbda);
     
      /* update thermostats */
      if (nose) 
      {
         if (verlet) {
            double nesum = 1.0*(mygrid->ne[0]+mygrid->ne[ispin-1]);
            double kefac = 0.5*fmass/(dt*dt);
            eke = kefac*(nesum-mygrid->gg_traceall(psi2,psi0));
            eki = myion->eki1;
            mynose->Verlet_step(eke, eki);
         } else {
            eke = fmass * mygrid->gg_traceall(psi0, psi0);
            // eki = myion->ke();
            eki = myion->eki1;
            mynose->Newton_step(eke, eki);
         }
      }
 
   } /* it, innerloop */
 
   //|-\____|\/-----\/\/->    End Parallel Section    <-\/\/-----\/|____/-|
 
   /****************************/
   /* total energy calculation */
   /****************************/
 
   /* if newton then skip energy calculations */
   if (verlet) 
   {
      mygrid->ggm_sym_Multiply(psi1, Hpsi, hml);
      mygrid->m_scal(-1.0,hml);
     
      eorbit = mygrid->m_trace(hml);
      if (ispin==1) eorbit = eorbit + eorbit;
     
      /* hartree energy and ion-ion energy */
      if (periodic) {
         ehartr = mycoulomb12->mycoulomb1->ecoulomb(dng);
         eion = myewald->energy();
      } else if (aperiodic) {
         ehartr = 0.5*mygrid->rr_dot(rho,vc)*dv;
         eion = myion->ion_ion_energy();
      }



      exc = mygrid->rr_dot(dnall,xce);
      pxc = mygrid->rr_dot(dn,xcp);
      if (ispin == 1) 
      {
         exc = exc + exc;
         pxc = pxc + pxc;
      } 
      else 
      {
         exc += mygrid->rr_dot(&dnall[n2ft3d], xce);
         pxc += mygrid->rr_dot(&dn[n2ft3d], &xcp[n2ft3d]);
      }
      exc *= dv;
      pxc *= dv;
     
      elocal = mygrid->cc_pack_dot(0, dng, vl);
     
      /* add in long range part here*/
      if (aperiodic)
         elocal += dv * mygrid->rr_dot(rho, vlr_l);
     
      /* add in other real-space fields here*/
      if (mypsp->myefield->efield_on)
         elocal += dv * mygrid->rr_dot(rho, mypsp->myefield->v_field);
     
      mygrid->g_zero(Hpsi);
      mypsp->v_nonlocal(psi1, Hpsi);
      enlocal = -mygrid->gg_traceall(psi1, Hpsi);
     
      /* set wavefunction velocity and kinetic enegy of psi */
      double h = 1.0 / (2.0 * dt);
      mygrid->g_Scale(-h, psi0);
      mygrid->gg_daxpy(h, psi2, psi0);
      eke = fmass * mygrid->gg_traceall(psi0, psi0);
     
      Eold = E[0];
      E[1] = eorbit + eion + exc - ehartr - pxc;
      E[2] = eke;
      E[3] = myion->ke();
      // E[3] = myion->eki1;
     
      E[4] = eorbit;
      E[5] = ehartr;
      E[6] = exc;
      E[7] = eion;
      E[8] = elocal;
      E[9] = enlocal;
      E[10] = 2 * ehartr;
      E[11] = pxc;
     
      /* get APC energies */
      if (mypsp->myapc->v_apc_on) 
      {
         E[51] = mypsp->myapc->Eapc;
         E[52] = mypsp->myapc->Papc;
         E[1] = E[1] + E[51] - E[52];
      }
     
      /* get dielectric energies */
      if (mycoulomb12->dielectric_on())
      {
         E[61] = mycoulomb12->edielec;
         E[62] = mycoulomb12->pdielec;
         E[0] = E[0] + E[61] - E[62];
      }
     
      /* get Efield energies */
      if (mypsp->myefield->efield_on) 
      {
         if (mypsp->myefield->efield_type==0) 
         {
            E[48] = 0.0;
            E[49] = 0.0;
            E[1] += E[48] - E[49];
         } 
         else 
         {
            E[48] = dv * mygrid->rr_dot(rho, mypsp->myefield->v_field);
            E[49] = mypsp->myefield->efield_ion_energy();
            E[1] += E[49];
         }
      }
     
      /* Running Sums - Energy and Energy**2 sum */
      E[24] += E[1];
      E[25] += E[1] * E[1];
      E[26] += (E[1] + E[2] + E[3]);
      E[27] += (E[1] + E[2] + E[3]) * (E[1] + E[2] + E[3]);
     
      /* Nose thermostat energies */
      if (nose) {
        E[8] = mynose->e_energy();
        E[9] = mynose->r_energy();
        E[0] = E[1] + E[2] + E[3] + E[8] + E[9];
      }
      /* add kinetic energies */
      else 
      {
        E[0] = E[1] + E[2] + E[3];
      }
   }
 
   mygrid->r_dealloc(tmp);
   mygrid->r_dealloc(xcp);
   mygrid->r_dealloc(xce);
   mygrid->r_dealloc(dnall);
   mygrid->r_dealloc(x);
   mygrid->r_dealloc(rho);
   mygrid->c_pack_deallocate(dng);
   mygrid->c_pack_deallocate(vl);
   if (periodic)
   {
      mygrid->c_pack_deallocate(vc);
      mygrid->c_pack_deallocate(vcall);
   }
   else if (aperiodic) 
   {
      mygrid->r_dealloc(vc);
      mygrid->r_dealloc(vcall);
      mygrid->r_dealloc(vlr_l);
   }
   if (mycoulomb12->dielectric_on())
      mygrid->r_dealloc(vdielec);
}
} // namespace pwdft

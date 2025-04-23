

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
#include "inner_loop.hpp"
#include "iofmt.hpp"
#include "psi_H.hpp"
#include "v_exc.hpp"

#include "nwpw_dplot.hpp"

namespace pwdft {

void inner_loop(Control2 &control, Pneb *mygrid, Ion *myion,
                Kinetic_Operator *myke, Coulomb12_Operator *mycoulomb12,
                XC_Operator *myxc, Pseudopotential *mypsp, Strfac *mystrfac,
                Ewald *myewald, double *psi1, double *psi2, double *Hpsi,
                double *psi_r, double *dn, double *hml, double *lmbda,
                double E[], double *deltae, double *deltac, double *deltar,
                bool fractional, double *occ1, double *occ2) 
{
   int it, it_in, i, n2ft3d, neall, ispin, k, ms;
   int shift1, shift2, indx1, indx2;
   int one = 1;
   double scal1, scal2, dv, dc;
   double eorbit, eion, econstraint, exc, ehartr, pxc;
   double eke, elocal, enlocal, dt, dte, Eold;
   double *vl,*vlr_l,*vc,*xcp,*xce,*dnall,*x,*dng,*rho,*tmp,*vcall,*vfield,*vpsi,*sumi;
   double *vdielec;
   double *fion;
   bool periodic = (control.version == 3);
   bool aperiodic = (control.version == 4);
   bool move = control.geometry_optimize();
   double omega = mygrid->lattice->omega();
 
   ispin = mygrid->ispin;
   neall = mygrid->neq[0] + mygrid->neq[1];
   shift1 = 2*(mygrid->npack(1));
   shift2 = (mygrid->n2ft3d);
   n2ft3d = (mygrid->n2ft3d);
   scal1 = 1.0/((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
   // scal2 = 1.0/lattice_omega();
   // dv = lattice_omega()*scal1;
   scal2 = 1.0/omega;
   dv = omega*scal1;
 
   dt = control.time_step();
   dte = dt/sqrt(control.fake_mass());
   it_in = control.loop(0);
 
   /* allocate temporary memory */
   rho = mygrid->r_alloc();
   tmp = mygrid->r_alloc();
   xcp = mygrid->r_nalloc(ispin);
   xce = mygrid->r_nalloc(ispin);
   dnall = mygrid->r_nalloc(ispin);
   x = mygrid->r_alloc();
   vfield = mygrid->r_alloc();
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
 
   //fion = new double[3 * (myion->nion)]();
   fion = myion->fion1;
 
   /* generate local psp*/
   mypsp->v_local(vl,false,dng,fion);
 
   // myewald->phafac();
 
   //|-\____|\/-----\/\/->    Start Parallel Section    <-\/\/-----\/|____/-|
 
   for (it = 0; it < it_in; ++it) 
   {
      mygrid->g_zero(Hpsi);
      mygrid->gg_copy(psi2, psi1);

      if (fractional)
         std::memcpy(occ1,occ2,(mygrid->ne[0]+mygrid->ne[1])*sizeof(double));
     
      if (move)
      {
         myion->shift();
         mystrfac->phafac();
         myewald->phafac();
         // for (int ii=0; ii<(3*myion->nion); ++ii) fion[ii] = 0.0;
      }
     
      /* convert psi(G) to psi(r) - Expensive */
      mygrid->gh_fftb(psi1,psi_r);
      /*
      indx1 = 0;
      indx2 = 0;
      for (i = 0; i < neall; ++i) 
      {
         mygrid->cc_pack_copy(1,psi1+indx1,psi_r+indx2);
         mygrid->c_unpack(1,psi_r+indx2);
         mygrid->cr_pfft3b(1,psi_r+indx2);
         indx1 += shift1;
         indx2 += shift2;
      }
      */
     
      /* generate dn */
      if (fractional)
         mygrid->hr_aSumSqr_occ(scal2,occ1,psi_r,dn);
      else
         mygrid->hr_aSumSqr(scal2,psi_r,dn);
     
      /* generate dng */
      mygrid->rrr_Sum(dn,dn+(ispin-1)*n2ft3d,rho);
      mygrid->rr_SMul(scal1, rho, tmp);
      mygrid->rc_pfft3f(0,tmp);
      //mygrid->rc_fft3d(tmp);
      mygrid->c_pack(0,tmp);
      mygrid->cc_pack_copy(0,tmp,dng);

      /* generate dnall - used for semicore corrections */
      if (mypsp->has_semicore()) 
      {
         if ((move) || (it == 0))
            mypsp->semicore_density_update();
         for (ms=0; ms<ispin; ++ms)
            mygrid->rrr_SMulAdd(0.5,mypsp->semicore_density,dn+ms*n2ft3d,dnall+ms*n2ft3d);
      } 
      else 
      {
         for (ms=0; ms<ispin; ++ms)
            mygrid->rr_copy(dn+ms*n2ft3d, dnall+ms*n2ft3d);
      }
     
      /* generate local potentials */
      if ((move) || (mypsp->myapc->v_apc_on)) 
      {
         mypsp->v_local(vl,move,dng,fion);
         if (mypsp->myapc->v_apc_on)
            mypsp->myapc->V_APC(dng,mypsp->zv,vl,move,fion);
      }
     
      /* long-range psp for charge systems */
      if (aperiodic) 
      {
         mypsp->v_lr_local(vlr_l);
         if (move)
            mypsp->grad_v_lr_local(rho, fion);
      }
     
      /* apply k-space operators */
      // myke->ke(psi1,Hpsi);
     
      /* apply non-local PSP  - Expensive */
      // mypsp->v_nonlocal_fion(psi1,Hpsi,move,fion);
     
      /* generate coulomb potential */
      if (periodic)
      {
         mycoulomb12->mycoulomb1->vcoulomb(dng,vc);
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
      std::memset(xcp,0,ispin*n2ft3d*sizeof(double));
      std::memset(xce,0,ispin*n2ft3d*sizeof(double));
      myxc->v_exc_all(ispin,dnall,xcp,xce);
     
     
      /* get Hpsi */
      if (periodic)
      {
         psi_H(mygrid,myke,mypsp,psi1,psi_r,vl,vcall,xcp,Hpsi,move,fion,occ1);
      }
      else if (aperiodic)
      {
         psi_Hv4(mygrid,myke,mypsp,psi1,psi_r,vl,vlr_l,vcall,xcp,Hpsi,move,fion,occ1);
      }
     
      /* do a steepest descent step */
      mygrid->gg_SMul(dte,Hpsi,psi2);
      mygrid->gg_Sum2(psi1,psi2);
     
      if (move)
      {
         /* get the ion-ion force */
         if (periodic)
            myewald->force(fion);
         else if (aperiodic)
            myion->ion_ion_force(fion);
        
         /* get the semicore force - needs to be checked */
         if (mypsp->has_semicore())
            mypsp->semicore_xc_fion(xcp, fion);
        
         /* get forces from external Efield */
         if (mypsp->myefield->efield_on)
            mypsp->myefield->efield_ion_fion(fion);
        
         /* steepest descent step */
         myion->add_contraint_force(fion);

         /* steepest descent step */
         myion->optimize_step(fion);
      }
     
      /* lagrange multiplier - Expensive */
      if (fractional)
         //mygrid->ggm_occ_lambda(dte, psi1, psi2, occ1, lmbda);
         mygrid->g_ortho(-1,psi2);
      else
         mygrid->ggm_lambda(dte, psi1, psi2, lmbda);
   }

 
   //|-\____|\/-----\/\/->    End Parallel Section    <-\/\/-----\/|____/-|
 
   /* total energy calculation */
   mygrid->ggm_sym_Multiply(psi1, Hpsi, hml);
   mygrid->m_scal(-1.0, hml);
   if (fractional)
      eorbit = mygrid->m_trace_occ(hml,occ1);
   else
      eorbit = mygrid->m_trace(hml);

   if (ispin == 1)
     eorbit = eorbit + eorbit;
 
   /* hartree energy and ion-ion energy */
   if (periodic) {
     ehartr = mycoulomb12->mycoulomb1->ecoulomb(dng);
     eion = myewald->energy();
   } else if (aperiodic) {
     ehartr = 0.5 * mygrid->rr_dot(rho, vc) * dv;
     eion = myion->ion_ion_energy();
   }
 
   /* xc energy */
   exc = mygrid->rr_dot(dnall, xce);
   pxc = mygrid->rr_dot(dn, xcp);
   if (ispin == 1) {
      exc = exc + exc;
      pxc = pxc + pxc;
   } else {
      exc += mygrid->rr_dot(dnall+n2ft3d, xce);
      pxc += mygrid->rr_dot(dn+n2ft3d, xcp+n2ft3d);
   }
   exc *= dv;
   pxc *= dv;
 
   /* average Kohn-Sham kineticl energy */
   eke = myke->ke_ave(psi1,occ1);
 
   /* average Kohn-Sham local psp energy */
   elocal = mygrid->cc_pack_dot(0, dng, vl);
 
   /* add in long range part here*/
   if (aperiodic)
      elocal += dv * mygrid->rr_dot(rho, vlr_l);
 
   /* add in other real-space fields here*/
   if (mypsp->myefield->efield_on)
     elocal += dv* mygrid->rr_dot(rho, mypsp->myefield->v_field);
 
   /* average Kohn-Sham v_nonlocal energy */
   mygrid->g_zero(Hpsi);
   mypsp->v_nonlocal(psi1, Hpsi);
   enlocal = -mygrid->gg_traceall_occ(psi1, Hpsi,occ1);
 
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

 
   /* get APC energies */
   if (mypsp->myapc->v_apc_on) 
   {
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
      if (mypsp->myefield->efield_type == 0) 
      {
         E[48] = 0.0;
         E[49] = 0.0;
         E[0] = E[0] + E[48] - E[49];
      } 
      else
      {
         E[48] = dv * mygrid->rr_dot(rho, mypsp->myefield->v_field);
         E[49] = mypsp->myefield->efield_ion_energy();
         E[0] = E[0] + E[49];
      }
   }

   /* get contraints energies */
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

 
   /* set convergence variables */
   *deltae = (E[0] - Eold) / (dt * control.loop(0));
 
   /* deltac */
   sumi = new double[neall]();
   mygrid->ggg_Minus(psi2, psi1, Hpsi);
   for (i=0; i<neall; ++i)
      sumi[i] = mygrid->cc_pack_idot(1, Hpsi+i*shift1, Hpsi+i*shift1 );
   mygrid->d3db::parall->Vector_SumAll(1, neall, sumi);
   dc = 0.0;
   for (i = 0; i < neall; ++i)
     if (sumi[i] > dc)
       dc = sumi[i];
   dc = mygrid->d3db::parall->MaxAll(2, dc);
   *deltac = dc / dte;
   delete[] sumi;
 
   /* deltar */
   *deltar = 0.0;
   if (move) {
     double sum;
     for (auto ii = 0; ii < (myion->nion); ++ii) {
       sum = sqrt(fion[3 * ii] * fion[3 * ii] +
                  fion[3 * ii + 1] * fion[3 * ii + 1] +
                  fion[3 * ii + 2] * fion[3 * ii + 2]);
       if (sum > (*deltar))
         *deltar = sum;
     }
   }
 
   // Plotting Fattebert dielectric function
    if (false)
    {
       double *sw = mygrid->r_alloc();
       double *sw2 = mygrid->r_alloc();
       // double *r_grid = mygrid->r_nalloc(3);
       // double *grx    = r_grid;
       // double *gry    = r_grid+n2ft3d;
       // double *grz    = r_grid+2*n2ft3d;
       double *rho_ion = mygrid->r_alloc();
       double *dng_ion = mygrid->r_alloc();
       double *v_ion = mygrid->r_alloc();
 
       // mygrid->generate_r_sym_grid(r_grid);
 
       // for (auto k=0; k<n2ft3d; ++k)
       //    sw[k] = util_switching_function(-5.0,1.2,r_grid[3*k+2]);
 
       /* initialize nwpw_dplot */
       double epsmin = 999.0;
       nwpw_dplot mydplot(myion, mygrid, control);
 
       mycoulomb12->generate_dng_ion(dng_ion);
       std::memcpy(rho_ion, dng_ion, n2ft3d * sizeof(double));
       mygrid->c_unpack(0, rho_ion);
       //mygrid->cr_fft3d(rho_ion);
       mygrid->cr_pfft3b(0,rho_ion);
       mygrid->r_zero_ends(rho_ion);
       mycoulomb12->mycoulomb2->vcoulomb(rho_ion, v_ion);
 
       mycoulomb12->v_dielectric2_aperiodic(rho, dng, rho_ion, dng_ion, vc, v_ion,
                                            false, sw, sw2, &mydplot);
 
       mycoulomb12->v_dielectric2_aperiodic(rho, dng, rho_ion, dng_ion, vc, v_ion,
                                            true, sw, sw2, &mydplot);
 
      mycoulomb12->v_dielectric_aperiodic(rho,dng,vc,sw,false,fion);
      std::cout << "Edielec = " << mycoulomb12->edielec << std::endl;
      std::cout << "Pdielec = " << mycoulomb12->pdielec << std::endl;
 
      mycoulomb12->v_dielectric_aperiodic(rho,dng,vc,sw,false,fion);
      std::cout << "Edielec = " << mycoulomb12->edielec << std::endl;
      std::cout << "Pdielec = " << mycoulomb12->pdielec << std::endl;
 
     mygrid->r_dealloc(v_ion);
     mygrid->r_dealloc(rho_ion);
     mygrid->r_dealloc(dng_ion);
     // mygrid->r_dealloc(r_grid);
     mygrid->r_dealloc(sw);
     mygrid->r_dealloc(sw2);
    }
 
   //delete[] fion;
 
   mygrid->r_dealloc(vfield);
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

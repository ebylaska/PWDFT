

#include        <iostream>
#include        <cstdio>
#include        <cmath>
#include        <cstdlib>
using namespace std;

#include        "Parallel.hpp"
#include        "Control2.hpp"
#include        "PGrid.hpp"
#include        "Ion.hpp"
#include        "Ewald.hpp"
#include        "Kinetic.hpp"
#include        "Coulomb.hpp"
#include        "exchange_correlation.hpp"
#include        "Pseudopotential.hpp"
//#include	"v_exc.hpp"
#include	"inner_loop_md.hpp"


namespace pwdft {
using namespace pwdft;

void inner_loop_md(bool verlet, Control2& control, Pneb *mygrid, Ion *myion, 
                Kinetic_Operator *myke, 
                Coulomb_Operator *mycoulomb, 
                XC_Operator      *myxc, 
                Pseudopotential *mypsp, Strfac *mystrfac, Ewald *myewald,
                double *psi0, double *psi1, double *psi2, double *Hpsi, double *psi_r,
                double *dn, double *hml,double *lmbda,
                double E[])
{
   int it,it_in,i,n2ft3d,neall,ispin,k,ms;
   int shift1,shift2,indx1,indx2;
   int one=1;
   double scal1,scal2,dv,dc;
   double eorbit,eion,exc,ehartr,pxc;
   double eke,elocal,enlocal,dt,dte,fmass,Eold;
   double *vl,*vc,*xcp,*xce,*dnall,*x,*dng,*rho,*tmp,*vall,*vpsi,*sumi;
   double *fion;
   bool move = true;
   double omega = mygrid->lattice->omega();

   ispin = mygrid->ispin;
   neall = mygrid->neq[0] + mygrid->neq[1];
   shift1 = 2*(mygrid->npack(1));
   shift2 = (mygrid->n2ft3d);
   n2ft3d = (mygrid->n2ft3d);
   scal1 = 1.0/((double) ((mygrid->nx)*(mygrid->ny)*(mygrid->nz)));
   //scal2 = 1.0/lattice_omega();
   //dv = lattice_omega()*scal1;
   scal2 = 1.0/omega;
   dv = omega*scal1;

   dt = control.time_step();
   fmass = control.fake_mass();
   dte = dt*dt/fmass;

   it_in = control.loop(0);

   /* allocate temporary memory */
   rho   = mygrid->r_alloc();
   tmp   = mygrid->r_alloc();
   xcp   = mygrid->r_nalloc(ispin);
   xce   = mygrid->r_nalloc(ispin);
   dnall = mygrid->r_nalloc(ispin);
   x   = mygrid->r_alloc();
   vall= mygrid->r_alloc();
   dng = mygrid->c_pack_allocate(0);
   vl  = mygrid->c_pack_allocate(0);
   vc  = mygrid->c_pack_allocate(0);
   vpsi=x;

   fion = new double[3*(myion->nion)]();

   /* generate local psp*/
   mypsp->v_local(vl,0,dng,fion);
   

   //myewald->phafac();

   for (it=0; it<it_in; ++it)
   {
      /* shift wavefuntion */
      mygrid->g_zero(Hpsi);
      mygrid->gg_copy(psi1,psi0);
      mygrid->gg_copy(psi2,psi1);
      //mygrid->gh_fftb(psi1,psi_r);

      /* skip shift ion if newton step */
      if (verlet) myion->shift();

      mystrfac->phafac();
      myewald->phafac();

      indx1 = 0;
      indx2 = 0;
      for (i=0; i<neall; ++i) 
      {
         mygrid->cc_pack_copy(1,&psi1[indx1],&psi_r[indx2]);
         mygrid->c_unpack(1,&psi_r[indx2]);
         mygrid->cr_fft3d(&psi_r[indx2]);
         indx1 += shift1;
         indx2 += shift2;
      }

      /* generate dn */
      mygrid->hr_aSumSqr(scal2,psi_r,dn);

      /* generate dng */
      mygrid->rrr_Sum(dn,&dn[(ispin-1)*n2ft3d],tmp);
      mygrid->r_SMul(scal1,tmp);
      mygrid->rc_fft3d(tmp);
      mygrid->c_pack(0,tmp);
      mygrid->cc_pack_copy(0,tmp,dng);

      /* generate dnall - used for semicore corrections */
      if (mypsp->has_semicore())
      {
         if ((move)||(it==0)) mypsp->semicore_density_update();
         for (ms=0; ms<ispin; ++ms)
            mygrid->rrr_SMulAdd(0.5,mypsp->semicore_density,&dn[ms*n2ft3d],&dnall[ms*n2ft3d]);
      }
      else
      {
          for (ms=0; ms<ispin; ++ms)
             mygrid->rr_copy(&dn[ms*n2ft3d],&dnall[ms*n2ft3d]);
      }

      /* generate local potential */
      if (move) mypsp->v_local(vl,move,dng,fion);
      if (mypsp->myapc->v_apc_on) mypsp->myapc->V_APC(dng,mypsp->zv,vl,move,fion);

      /* apply k-space operators */
      myke->ke(psi1,Hpsi);
      mypsp->v_nonlocal_fion(psi1,Hpsi,move,fion);

      /* generate coulomb potential */
      mycoulomb->vcoulomb(dng,vc);

      /* generate exchange-correlation potential */
      myxc->v_exc_all(ispin,dnall,xcp,xce);

      /* apply r-space operators */
     mygrid->cc_SMul(0,scal2,vl,vall);
     mygrid->cc_Sum2(0,vc,vall);
     mygrid->c_unpack(0,vall);
     mygrid->cr_fft3d(vall);
     indx1 = 0;
     indx2 = 0;
     for (ms=0; ms<ispin; ++ms) 
     { 
         mygrid->rrr_Sum(vall,&xcp[ms*n2ft3d],tmp);
         for (i=0; i<(mygrid->neq[ms]); ++i)
         {
            mygrid->rrr_Mul(tmp,&psi_r[indx2],vpsi);
            mygrid->rc_fft3d(vpsi);
            mygrid->c_pack(1,vpsi);
            mygrid->cc_daxpy(1,(-scal1),vpsi,&Hpsi[indx1]);

            indx1 += shift1;
            indx2 += shift2;
         }
     }

     myewald->force(fion);

     /* car-parrinello Verlet step */
     if (verlet) 
     {
        mygrid->gg_SMul(dte,Hpsi,psi2);
        mygrid->gg_daxpy(-1.0,psi0,psi2);
        mygrid->gg_daxpy( 2.0,psi1,psi2);

        myion->Verlet_step(fion);
     }

     /* car-parrinello Newton step */
     else
     {
        mygrid->gg_SMul(dte,Hpsi,psi2);
        mygrid->gg_daxpy(dt,psi0,psi2);
        mygrid->gg_Sum2(psi1,psi2);

        myion->Newton_step(fion);
     }

     /* lagrange multiplier */
     mygrid->ggm_lambda(dte,psi1,psi2,lmbda);
   }

   /* total energy calculation */
   mygrid->ggm_sym_Multiply(psi1,Hpsi,hml);
   mygrid->m_scal(-1.0,hml);

   /* if newton then skip energy calculations */
   if (verlet)
   {
      eorbit  = mygrid->m_trace(hml);
      if (ispin==1) eorbit = eorbit+eorbit;

      ehartr  = mycoulomb->ecoulomb(dng);
      exc     = mygrid->rr_dot(dnall,xce);
      pxc     = mygrid->rr_dot(dn,xcp);
      if (ispin==1)
      {
         exc = exc+exc;
         pxc = pxc+pxc;
      }
      else
      {
         exc += mygrid->rr_dot(&dnall[n2ft3d],xce);
         pxc += mygrid->rr_dot(&dn[n2ft3d],&xcp[n2ft3d]);
      }
      exc *= dv;
      pxc *= dv;

      eion    = myewald->energy();
      elocal  = mygrid->cc_pack_dot(0,dng,vl);
      mygrid->g_zero(Hpsi);
      mypsp->v_nonlocal(psi1,Hpsi);
      enlocal = -mygrid->gg_traceall(psi1,Hpsi);

      /* velocity and kinetic enegy of psi */
      double h = 1.0/(2.8*dt);
      mygrid->g_Scale(-h,psi0);
      mygrid->gg_daxpy(h,psi2,psi0);
      eke = fmass*mygrid->gg_traceall(psi0,psi0);


      Eold = E[0];
      E[1] = eorbit + eion + exc - ehartr - pxc;
      E[2] = eke;
      E[3] = myion->ke();

      E[4] = eorbit;
      E[5] = ehartr;
      E[6] = exc;
      E[7] = eion;

      E[8] = elocal;
      E[9] = enlocal;
      E[10] = 2*ehartr;
      E[11] = pxc;

      if (mypsp->myapc->v_apc_on)
      {
         E[51] = mypsp->myapc->Eapc;
         E[52] = mypsp->myapc->Papc;
         E[1]  = E[1] + E[51] - E[52];
      }


      /* add kinetic energies */
      E[0] = E[1] + E[2] + E[3];

   }

   /* deallocate local heap data */
   delete [] fion;

   mygrid->r_dealloc(tmp);
   mygrid->r_dealloc(xcp);
   mygrid->r_dealloc(xce);
   mygrid->r_dealloc(dnall);
   mygrid->r_dealloc(x);
   mygrid->r_dealloc(vall);
   mygrid->r_dealloc(rho);
   mygrid->c_pack_deallocate(dng);
   mygrid->c_pack_deallocate(vl);
   mygrid->c_pack_deallocate(vc);
}
}

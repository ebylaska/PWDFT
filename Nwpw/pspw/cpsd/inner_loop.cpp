

#include        <iostream>
#include        <cstdio>
#include        <cmath>
#include        <cstdlib>


#include        "Parallel.hpp"
#include        "Control2.hpp"
#include        "PGrid.hpp"
#include        "Ion.hpp"
#include        "Ewald.hpp"
#include        "Kinetic.hpp"
#include        "Coulomb12.hpp"
#include        "exchange_correlation.hpp"
#include        "Pseudopotential.hpp"
#include        "psi_H.hpp"
#include	"v_exc.hpp"
#include	"inner_loop.hpp"
#include        "iofmt.hpp"

namespace pwdft {



void inner_loop(Control2& control, Pneb *mygrid, Ion *myion, 
                Kinetic_Operator *myke, 
                Coulomb12_Operator *mycoulomb12, 
                XC_Operator *myxc, 
                Pseudopotential *mypsp, Strfac *mystrfac, Ewald *myewald,
                double *psi1, double *psi2, double *Hpsi, double *psi_r,
                double *dn, double *hml,double *lmbda,
                double E[], double *deltae, double *deltac, double *deltar)
{
   int it,it_in,i,n2ft3d,neall,ispin,k,ms;
   int shift1,shift2,indx1,indx2;
   int one=1;
   double scal1,scal2,dv,dc;
   double eorbit,eion,exc,ehartr,pxc;
   double eke,elocal,enlocal,dt,dte,Eold;
   double *vl,*vlr_l,*vc,*xcp,*xce,*dnall,*x,*dng,*rho,*tmp,*vall,*vfield,*vpsi,*sumi;
   double *fion;
   bool periodic  = (control.version==3);
   bool aperiodic = (control.version==4);
   bool move = control.geometry_optimize();
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
   dte = dt /sqrt(control.fake_mass());
   it_in = control.loop(0);

   /* allocate temporary memory */
   rho   = mygrid->r_alloc();
   tmp   = mygrid->r_alloc();
   xcp   = mygrid->r_nalloc(ispin);
   xce   = mygrid->r_nalloc(ispin);
   dnall = mygrid->r_nalloc(ispin);
   x     = mygrid->r_alloc();
   vall  = mygrid->r_alloc();
   vfield= mygrid->r_alloc();
   dng = mygrid->c_pack_allocate(0);
   vl  = mygrid->c_pack_allocate(0);
   if (periodic)
   {
      vc  = mygrid->c_pack_allocate(0);
   }
   else if (aperiodic)
   {
      vc    = mygrid->r_alloc();
      vlr_l = mygrid->r_alloc();
   }

   vpsi=x;

   fion = new double[3*(myion->nion)]();

   /* generate local psp*/
   mypsp->v_local(vl,0,dng,fion);
   


   //myewald->phafac();

   for (it=0; it<it_in; ++it)
   {
      mygrid->g_zero(Hpsi);
      mygrid->gg_copy(psi2,psi1);
      //mygrid->gh_fftb(psi1,psi_r);

      if (move)
      {
         myion->shift();
         mystrfac->phafac();
         myewald->phafac();
         //for (int ii=0; ii<(3*myion->nion); ++ii) fion[ii] = 0.0;
      }

      /* convert psi(G) to psi(r) - Expensive */
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
      mygrid->rrr_Sum(dn,&dn[(ispin-1)*n2ft3d],rho);
      mygrid->rr_SMul(scal1,rho,tmp);
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
         if (move) mypsp->grad_v_lr_local(rho,fion);
      }

      /* apply k-space operators */
      //myke->ke(psi1,Hpsi);

      /* apply non-local PSP  - Expensive */
      //mypsp->v_nonlocal_fion(psi1,Hpsi,move,fion);

      /* generate coulomb potential */
      if (periodic)  
         mycoulomb12->mycoulomb1->vcoulomb(dng,vc);
      else if (aperiodic) 
         mycoulomb12->mycoulomb2->vcoulomb(rho,vc);

      /* generate exchange-correlation potential */
      std::memset(xcp,0,ispin*n2ft3d);
      std::memset(xce,0,ispin*n2ft3d);
      myxc->v_exc_all(ispin,dnall,xcp,xce);


      //v_exc(ispin,shift2,dnall,xcp,xce,x);

      //for (k=0; k<n2ft3d; ++k)
      //   std::cout << "k=" << k 
      //             << " xce=" << xce[k] 
      //             << " xcp=" << xcp[k] << " " << xcp[k+n2ft3d] << std::endl;
      //for (k=0; k<n2ft3d; ++k)
      //   std::cout << "k=" << k 
      //             << " dnall=" << dnall[k] 
      //             << " " << dnall[k+n2ft3d]  << std::endl;


      /* get Hpsi */
      if (periodic)
         psi_H(mygrid,myke,mypsp,psi1,psi_r,vl,vc,xcp,Hpsi,move,fion);
      else if (aperiodic) 
         psi_Hv4(mygrid,myke,mypsp,psi1,psi_r,vl,vlr_l,vc,xcp,Hpsi,move,fion);

     /* apply r-space operators  - Expensive*/
     //mygrid->cc_pack_SMul(0,scal2,vl,vall);
     //mygrid->cc_pack_Sum2(0,vc,vall);
     //mygrid->c_unpack(0,vall);
     //mygrid->cr_fft3d(vall);
     //indx1 = 0;
     //indx2 = 0;
     //for (ms=0; ms<ispin; ++ms) 
     //{ 
     //    mygrid->rrr_Sum(vall,&xcp[ms*n2ft3d],tmp);
     //    for (i=0; i<(mygrid->neq[ms]); ++i)
     //    {
     //       mygrid->rrr_Mul(tmp,&psi_r[indx2],vpsi);
     //       mygrid->rc_fft3d(vpsi);
     //       mygrid->c_pack(1,vpsi);
     //       mygrid->cc_pack_daxpy(1,(-scal1),vpsi,&Hpsi[indx1]);
     //
     //       indx1 += shift1;
     //       indx2 += shift2;
     //    }
     //}

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
           mypsp->semicore_xc_fion(xcp,fion);

        /* get forces from external Efield */
        if (mypsp->myefield->efield_on)
           mypsp->myefield->efield_ion_fion(fion);

        /* steepest descent step */
        myion->optimize_step(fion);
     }

     /* lagrange multiplier - Expensive */
     mygrid->ggm_lambda(dte,psi1,psi2,lmbda);
   }

   /* total energy calculation */
   mygrid->ggm_sym_Multiply(psi1,Hpsi,hml);
   mygrid->m_scal(-1.0,hml);
   eorbit  = mygrid->m_trace(hml);
   if (ispin==1) eorbit = eorbit+eorbit;
   /*
   std::cout << "hmlup=" << Ffmt(20,15) << hml[0] << " " << hml[1] << " " << hml[2] << " " << hml[3] << std::endl; 
   std::cout << "      " << hml[4] << " " << hml[5] << " " << hml[6] << " " << hml[7] << std::endl; 
   std::cout << "      " << hml[8] << " " << hml[9] << " " << hml[10] << " " << hml[11] << std::endl; 
   std::cout << "      " << hml[12] << " " << hml[13] << " " << hml[14] << " " << hml[15] << std::endl; 

   std::cout << "hmldn=" << Ffmt(20,15) << hml[16+0] << " " << hml[16+1] << " " << hml[16+2] << std::endl; 
   std::cout << "      " << hml[16+3] << " " << hml[16+4] << " " << hml[16+5] << std::endl; 
   std::cout << "      " << hml[16+6] << " " << hml[16+7] << " " << hml[16+8] << std::endl; 
   */


   /* hartree energy and ion-ion energy */
   if (periodic) 
   {
      ehartr  = mycoulomb12->mycoulomb1->ecoulomb(dng);
      eion    = myewald->energy();
   }
   else if (aperiodic)
   {
       ehartr  = 0.5*mygrid->rr_dot(rho,vc)*dv;
       eion = myion->ion_ion_energy();
   }

   /* xc energy */
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


   /* average Kohn-Sham kineticl energy */
   eke     = myke->ke_ave(psi1);

   /* average Kohn-Sham local psp energy */
   elocal  = mygrid->cc_pack_dot(0,dng,vl);

   /* add in long range part here*/
   if (aperiodic)
      elocal += dv*mygrid->rr_dot(rho,vlr_l);

   /* add in other real-space fields here*/
   if (mypsp->myefield->efield_on)
      elocal += dv*mygrid->rr_dot(rho,mypsp->myefield->v_field);

   /* average Kohn-Sham v_nonlocal energy */
   mygrid->g_zero(Hpsi);
   mypsp->v_nonlocal(psi1,Hpsi);
   enlocal = -mygrid->gg_traceall(psi1,Hpsi);

   Eold = E[0];
   E[0] = eorbit + eion + exc - ehartr - pxc;
   E[1] = eorbit;
   E[2] = ehartr;
   E[3] = exc;
   E[4] = eion;
   E[5] = eke;
   E[6] = elocal;
   E[7] = enlocal;
   E[8] = 2*ehartr;
   E[9] = pxc;

   /* get APC energies */
   if (mypsp->myapc->v_apc_on)
   {
      E[51] = mypsp->myapc->Eapc;
      E[52] = mypsp->myapc->Papc;
      E[0]  = E[0] + E[51] - E[52];
   }


   /* get Efield energies */
   if (mypsp->myefield->efield_on)
   {
      if (mypsp->myefield->efield_type==0)
      {  
         E[48] = 0.0;
         E[49] = 0.0;
         E[0]  = E[0] + E[48] - E[49];
      } 
      else
      {
         E[48] = dv*mygrid->rr_dot(rho,mypsp->myefield->v_field);
         E[49] = mypsp->myefield->efield_ion_energy();
         E[0]  = E[0] + E[49];
      }
   }



   /* set convergence variables */
   *deltae = (E[0]-Eold)/(dt*control.loop(0));

   /* deltac */
   sumi    = new double[neall]();
   mygrid->ggg_Minus(psi2,psi1,Hpsi);
   for (i=0; i<neall; ++i) 
      sumi[i] = mygrid->cc_pack_idot(1,&Hpsi[i*shift1],&Hpsi[i*shift1]);
   mygrid->d3db::parall->Vector_SumAll(1,neall,sumi);
   dc = 0.0;
   for (i=0; i<neall; ++i) 
      if (sumi[i]>dc) dc=sumi[i];
   dc = mygrid->d3db::parall->MaxAll(2,dc);
   *deltac = dc/dte;
   delete [] sumi;

   /* deltar */ 
   *deltar = 0.0;
   if (move)
   {
      double sum;
      for (auto ii=0; ii<(myion->nion); ++ii)
      {
         sum = sqrt( fion[3*ii]  *fion[3*ii] 
                   + fion[3*ii+1]*fion[3*ii+1] 
                   + fion[3*ii+2]*fion[3*ii+2]);
         if (sum>(*deltar)) *deltar = sum;
      }
   }


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
   if (periodic)
      mygrid->c_pack_deallocate(vc);
   else if (aperiodic)
   {
      mygrid->r_dealloc(vc);
      mygrid->r_dealloc(vlr_l);
   }
}
}

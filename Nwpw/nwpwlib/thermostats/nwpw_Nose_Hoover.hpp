#ifndef _NWPW_NOSE_HOOVER_HPP_
#define _NWPW_NOSE_HOOVER_HPP_

// ********************************************************************
// *                                                                  *
// *       nwpw_nose_hoover : implementation of Nose-Hoover dual      *
// *                          thermostats.                            *
// *                                                                  *
// *   The algorithms used in this module are based on the work of    *
// *   P.E. Blochl and M.Parrinello,  J. Chem. Phys.                   *
// *                                                                  *
// ********************************************************************
#include        <iostream>
#include        "blas.h"

#include        "Control2.hpp"
#include        "Ion.hpp"


namespace pwdft {
using namespace pwdft;

class nwpw_Nose_Hoover {


public:
   bool nose_on = false;
   bool nosers  = false;

   int    Mchain=0;
   int    Nchain=0;
   int    nion,nemax;
   double *Xem,*Xe0,*Xe1,*Xe2,*Ee0,*Qe,*Pe;
   double *Xrm,*Xr0,*Xr1,*Xr2,*Er0,*Qr,*Pr;

   double eke0_init,Ne_chain,Pe_init,Pr_init,g_dof;

   double Te,Tr,time_step,fmass,total_mass;

   /* constructor */
   nwpw_Nose_Hoover(Ion&,int,double,Control2&);

   /* destructor */
   ~nwpw_Nose_Hoover() {
      if (nose_on) {
         delete [] Xem;
         delete [] Xe0;
         delete [] Xe1;
         delete [] Xe2;
         delete [] Ee0;
         delete [] Qe;
         delete [] Pe;

         delete [] Xrm;
         delete [] Xr0;
         delete [] Xr1;
         delete [] Xr2;
         delete [] Er0;
         delete [] Qr;
         delete [] Pr;
      }
   }

   /* access functions */
   bool on()      { return nose_on; }
   bool restart() { return nosers; }

   void reset_T(double Te_new, double Tr_new) {

      Te = Te_new;
      Tr = Tr_new;

      double kb = 3.16679e-6;

      /* Set Er0(1) = (1/2)*(g)*k*T), where g=number of degrees of freedom */
      Er0[0] = 0.5*g_dof*kb*Tr;

      /* Set Er0(1:Nchain-1) = 1/2*(k*T) */
      if (Nchain>1) for (auto n=1; n<Nchain; ++n) Er0[n] = 0.5*kb*Tr;

      /* Set Ee0(0) - read from rtdb otherwise use current KE */
      Ee0[0] = 4.0*kb*Te*fmass*(nion/total_mass)*eke0_init;

      /* Set Ee0(1:Mchain-1) = 1/2*(1/betae), where 1/betae = 2*Ee/Ne */
      if (Mchain>1) 
      {
         double betae = Ee0[0]/Ne_chain;
         for (auto m=1; m<Mchain; ++m) Ee0[m] = betae;
      }
   }

   void zero_thermostats() {
      if (nose_on) {
         for (auto m=0; m<Mchain; ++m) Xe2[m] = 0.0;
         for (auto m=0; m<Mchain; ++m) Xe1[m] = 0.0;
         for (auto m=0; m<Mchain; ++m) Xe0[m] = 0.0;
         for (auto m=0; m<Mchain; ++m) Xem[m] = 0.0;
         for (auto m=0; m<Mchain; ++m) Pe[m]  = Pe_init;

         for (auto n=0; n<Nchain; ++n) Xr2[n] = 0.0;
         for (auto n=0; n<Nchain; ++n) Xr1[n] = 0.0;
         for (auto n=0; n<Nchain; ++n) Xr0[n] = 0.0;
         for (auto n=0; n<Nchain; ++n) Xrm[n] = 0.0;
         for (auto n=0; n<Nchain; ++n) Pr[n]  = Pr_init;
      }
   }

   void Newton_Step(double eke, double eki) {
      if (nose_on) {
         double eke_tmp = eke;
         for (auto m=0; m<(Mchain-1); ++m)
         {
            /* integrate thermostat using newton step */
            double FXe = 2.0*(eke_tmp-Ee0[m]);
            double a   = time_step*(1.0 - 0.5*time_step*Xem[m]);
            Xe2[m]     =  Xe1[m] + a*Xem[m] + (0.5*time_step*time_step/Qe[m])*FXe;

            /* define kinetic energy for next link in the chain */
            eke_tmp = Xem[m];
            eke_tmp = 0.5*Qe[m]*(eke_tmp*eke_tmp);
         }

         /* Last link of chain */
         double FXe    = 2.0*(eke_tmp-Ee0[Mchain-1]);
         Xe2[Mchain-1] = Xe1[Mchain-1] + time_step*Xem[Mchain-1] + (0.5*time_step*time_step/Qe[Mchain-1])*FXe;

 
         double ekr_tmp = eki;
         for (auto n=0; n<(Nchain-1); ++n)
         {
            /* integrate thermostat using newton step */
            double FXr = 2.0*(ekr_tmp-Er0[n]);
            double a   = time_step*(1.0 - 0.5*time_step*Xrm[n]);
            Xr2[n]     =  Xr1[n] + a*Xrm[n] + (0.5*time_step*time_step/Qr[n])*FXr;

            /* define kinetic energy for next link in the chain */
            ekr_tmp = Xrm[n];
            ekr_tmp = 0.5*Qr[n]*(ekr_tmp*ekr_tmp);
         }

         /* Last link of chain */
         double FXr    = 2.0*(ekr_tmp-Er0[Nchain-1]);
         Xr2[Nchain-1] = Xr1[Nchain-1] + time_step*Xrm[Nchain-1] + (0.5*time_step*time_step/Qr[Nchain-1])*FXr;
      }
   }

   void Verlet_Step(double eke, double eki) {
      if (nose_on) {
         double eke_tmp = eke;
         for (auto m=0; m<(Mchain-1); ++m)
         {
            /* define dXe/dt = (3*Xe(t) - 4*Xe(t-dt) + Xe(t-2*dt))/(2*dt) */
            double dXe = (3.0*Xe1[m] - 4.0*Xe0[m] +Xem[m])/(2.0*time_step);
            double sse = 1.0/(1.0 +0.5*dXe*time_step);

            /* integrate thermostat using modified verlet */
            double  FXe  = 2.0*(eke_tmp-Ee0[m]);
            Xe2[m] = Xe0[m] + (Xe1[m] - Xe0[m] +  (0.5*time_step*time_step/Qe[m])*FXe)*2.0*sse;

            /* define kinetic energy for next link in the chain */
            eke_tmp = (Xe2[m]-Xe0[m])/(2.0*time_step);
            eke_tmp = 0.5*Qe[m]*(eke_tmp*eke_tmp);
         }

         /* Last link of chain */
         double FXe    = 2.0*(eke_tmp-Ee0[Mchain-1]);
         Xe2[Mchain-1] = 2.0*Xe1[Mchain-1] - Xe0[Mchain-1] + (time_step*time_step/Qe[Mchain-1])*FXe;


         double ekr_tmp = eki;
         for (auto n=0; n<(Nchain-1); ++n)
         {
            /* define dXr/dt = (3*Xr(t) - 4*Xr(t-dt) + Xr(t-2*dt))/(2*dt) */
            double dXr = (3.0*Xr1[n] - 4.0*Xr0[n] +Xrm[n])/(2.0*time_step);
            double ssr = 1.0/(1.0 +0.5*dXr*time_step);

            /* integrate thermostat using modified verlet */
            double  FXr  = 2.0*(ekr_tmp-Er0[n]);
            Xr2[n] = Xr0[n] + (Xr1[n] - Xr0[n] +  (0.5*time_step*time_step/Qr[n])*FXr)*2.0*ssr;

            /* define kinetic energy for next link in the chain */
            ekr_tmp = (Xr2[n]-Xr0[n])/(2.0*time_step);
            ekr_tmp = 0.5*Qr[n]*(ekr_tmp*ekr_tmp);
         }

         /* Last link of chain */
         double FXr    = 2.0*(ekr_tmp-Er0[Mchain-1]);
         Xr2[Mchain-1] = 2.0*Xr1[Mchain-1] - Xr0[Mchain-1] + (time_step*time_step/Qr[Mchain-1])*FXr;
      }
   }

   void shift() {
      if (nose_on) 
      {
         for (auto m=0; m<Mchain; ++m) { Xem[m] = Xe0[m]; Xe0[m] = Xe1[m]; Xe1[m] = Xe2[m]; }
         for (auto n=0; n<Nchain; ++n) { Xrm[n] = Xr0[n]; Xr0[n] = Xr1[n]; Xr1[n] = Xr2[n]; }
      }

   }

   double sse() {
      double dXe = (3.0*Xe1[0] - 4.0*Xe0[0] + Xem[0])/(2.0*time_step);
      return (1.0/(1.0+0.5*dXe*time_step));
   }

   double ssr() {
      double dXr = (3.0*Xr1[0] - 4.0*Xr0[0] + Xrm[0])/(2.0*time_step);
      return (1.0/(1.0+0.5*dXr*time_step));
   }

   double e_energy() { 
      double esum = 0.0;
      for (auto m=0; m<Mchain; ++m)
      {
         double dXe = (3.0*Xe1[m] - 4.0*Xe0[m] + Xem[m])/(2.0*time_step);
         esum += (0.5*Qe[m]*dXe*dXe + 2.0*Ee0[m]*Xe1[m]);
      }
      return esum;
   }

   double r_energy() { 
      double esum = 0.0;
      for (auto n=0; n<Nchain; ++n)
      {
         double dXr  = (3.0*Xr1[n] - 4.0*Xr0[n] + Xrm[n])/(2.0*time_step);
         esum += (0.5*Qr[n]*dXr*dXr + 2.0*Er0[n]*Xr1[n]);
      }
      return esum;
   }

   void writejsonstr(string&);

   std::string inputprint();
   std::string thermostatprint();


};

}

#endif



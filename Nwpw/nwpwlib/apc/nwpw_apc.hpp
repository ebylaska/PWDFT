#ifndef _NWPW_APC_HPP_
#define _NWPW_APC_HPP_

// ********************************************************************
// *                                                                  *
// *       nwpw_apc : used to generate derived atomic point charges   *
// *                  from a plane-wave density.                      *
// *                                                                  *
// *   The algorithms used in this module are based on the work of    *
// *   P.E. Blochl, J. Chem. Phys. vol. 103, page 7422 (1995).        *
// *                                                                  *
// ********************************************************************
#include        <iostream>
#include        "blas.h"

#include        "gdevice.hpp"
#include        "Control2.hpp"
#include        "Ion.hpp"
#include        "Pneb.hpp"
#include        "Strfac.hpp"


namespace pwdft {
using namespace pwdft;

class nwpw_apc {

   Pneb   *mypneb;
   Ion    *myion;
   Strfac *mystrfac;

public:
   bool apc_on;
   int nga,ngs;
   double Gc;

   double *gamma,*A,*Am,*b,*q,*u,*w,*gaus;
   double Eapc,Papc;

   /* constructor */
   nwpw_apc(Ion *, Pneb *, Strfac *, Control2&);

   /* destructor */
   ~nwpw_apc() {
      if (apc_on) {
         delete [] gamma;
         delete [] A;
         delete [] Am;
         delete [] b;
         delete [] q;
         delete [] u;
         delete [] w;
         delete [] gaus;
      }
   }

   void gen_APC(double *, bool);
   void dngen_APC(double *, bool);


};

}

#endif


//pspw_gen_APC(ispin,ne,dng,move)
//pspw_dngen_APC(ispin,ne,dn,move)
//pspw_sumAm_APC(ngs,Am)
//pspw_Amtimesu_APC(ngs,Am,u)
//pspw_Vfac_APC(ngs,Am,u,i)

//pspw_VQ_APC(nion,u,VQ)

//pspw_cosmo_V0_APC(vcosmo)
//pspw_cosmo_V_APC(ispin,ne,dng,vcosmo,ecosmo,pcosmo,move,fion)
//pspw_cosmo_force_APC(ispin,ne,dng,fion)


//pspw_born_V_APC(ispin,ne,dng,vborn,eborn,pborn,move,fion)
//pspw_born_force_APC(ispin,ne,dng,fion)


//pspw_cdft_V_APC(ispin,ne,dng,vcdft,ecdft,pcdft,move,fion)
//pspw_cdft_force_APC(ispin,ne,dng,fion)

//pspw_dQdR_APC(nion,u,fion)
//pspw_generate_dQdR(lb,nga,nion,db,q,u,dA,Am,SumAm,dbtmp,dAtmp)

//pspw_V_APC_on()
//pspw_V_APC(ispin,ne,dng,vapc,Eapc,Papc,move,fion)
//pspw_force_APC(ispin,ne,dng,fion)

//pspw_E_APC(Eapc,Papc)


//pspw_getQ_APC(ii,n)
//pspw_getQtot_APC(ii)

//pspw_shortprint_APC(unit)
//pspw_print_APC(unit)


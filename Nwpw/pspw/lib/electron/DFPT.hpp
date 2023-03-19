#ifndef _DFPT_HPP_
#define _DFPT_HPP_

#pragma once

#include	"Pneb.hpp"
#include	"Kinetic.hpp"
#include	"Coulomb12.hpp"
#include	"exchange_correlation.hpp"
#include        "Strfac.hpp"
#include	"Pseudopotential.hpp"

namespace pwdft {


class	DFPT_Operators {

   Pneb   *mygrid;

   Kinetic_Operator   *myke;
   Coulomb12_Operator *mycoulomb12;
   XC_Operator        *myxc;
   Pseudopotential    *mypsp;

   double *psi0,*psi0_r,*dnall0,*dEpertdpsi0_r;
   double *Hpsi1,*psi1,*b1,*psi1_r;
   double *dn1,*rho1,*dng1,*dnallp,*dnallm;
   double *vc1,*vxc1,*xcpp,*xcpm;
   double *hmltmp,*tmp;

   double omega,scal2,scal1,dv,eps;

   int ispin,neall,n2ft3d,shift1,shift2,npack1;
   bool aperiodic = false;
   bool  periodic = false;

public:
   int counter=0;

   /* Constructors */
   DFPT_Operators(Pneb *, Kinetic_Operator *, Coulomb12_Operator *, XC_Operator *, Pseudopotential *, const double);

   /* destructor */
   ~DFPT_Operators() {
         delete [] Hpsi1;
         delete [] psi1;
         delete [] b1;
         delete [] psi1_r;
         delete [] dn1;
         delete [] rho1;
         delete [] dng1;
         delete [] dnallm;
         delete [] xcpp;
         delete [] xcpm;
         delete [] vc1;
         delete [] vxc1;
         delete [] hmltmp;
         delete [] tmp;
    }


    void start(double *, double *, double *, double *);
    void run();

};

}

#endif

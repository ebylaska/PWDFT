#ifndef _CELECTRON_HPP_
#define _CELECTRON_HPP_

#pragma once

#include "cCoulomb.hpp"
#include "cKinetic.hpp"
#include "Cneb.hpp"
#include "CPseudopotential.hpp"
#include "CStrfac.hpp"
#include "cExchange_Correlation.hpp"

namespace pwdft {

class cElectron_Operators {

   Cneb *mygrid;
 
   cKinetic_Operator *myke;
   cCoulomb_Operator *mycoulomb;
   cXC_Operator *myxc;
   CPseudopotential *mypsp;
 
   double *Hpsi, *psi_r, *vl, *vall, *vc, *xcp, *xce, *x, *rho, *hmltmp;
   double *vcall;
 
   double omega, scal2, scal1, dv;
 
   int ispin, neall, nfft3d, shift1, shift2, npack1;
   bool periodic = true;
 
public:
   int counter = 0;
 
   /* Constructors */
   cElectron_Operators(Cneb *, cKinetic_Operator *, cCoulomb_Operator *,
                       cXC_Operator *, CPseudopotential *);
 
   /* destructor */
   ~cElectron_Operators() 
   {
      delete[] Hpsi;
      delete[] psi_r;
      delete[] vl;
      delete[] x;
      delete[] rho;
      delete[] vall;
      delete[] vc;
      delete[] vcall;
      delete[] xcp;
      delete[] xce;
      delete[] hmltmp;
      //  delete[] vdielec;
   }
 
   void genrho(double *, double *);
   void run(double *, double *, double *, double *);
   void get_Tgradient(double *, double *, double *);
   void get_Gradient(double *);
   void gen_Tangent(double *, double *, double *);
   void gen_hml(double *, double *);
   void gen_hmlt(double *, double *);
   void gen_Hpsi_k(double *);
   void gen_psi_r(double *);
   void gen_density(double *);
   void gen_densities(double *, double *, double *);
   void gen_scf_potentials(double *, double *, double *);
   void gen_vl_potential();
   void semicore_density_update();
   void gen_vall();
   void get_vall(double *);
   void set_vall(const double *);
 
   double vl_ave(double *);
   double vnl_ave(double *);
 
   double energy(double *, double *, double *, double *);
   double eorbit(double *);
   double ehartree(double *);
   double ehartree2(double *);
   double exc(double *);
   double pxc(double *);
   double eke(double *);
 
   void gen_energies_en(double *, double *, double *, double *, double *,
                        double *);
 
   void add_dteHpsi(double, double *, double *);
 
   void vl_force(double *, double *);
   void vnl_force(double *, double *);
   void semicore_force(double *);
 
   bool is_periodic() { return periodic; }

   cKinetic_Operator *get_myke() {return myke;}
};

} // namespace pwdft

#endif

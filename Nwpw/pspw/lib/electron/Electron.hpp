#ifndef _ELECTRON_HPP_
#define _ELECTRON_HPP_

#pragma once

#include "Coulomb12.hpp"
#include "Kinetic.hpp"
#include "Pneb.hpp"
#include "Pseudopotential.hpp"
#include "Strfac.hpp"
#include "exchange_correlation.hpp"

namespace pwdft {

class Electron_Operators {

   Pneb *mygrid;
   Kinetic_Operator *myke;
 
   Coulomb12_Operator *mycoulomb12;
   XC_Operator *myxc;
   Pseudopotential *mypsp;
 
   double *Hpsi, *psi_r, *vl, *vall, *vc, *xcp, *xce, *x, *rho, *hmltmp;
   double *vlr_l; 
   double *vcall, *vdielec;
 
   double omega, scal2, scal1, dv;
 
   int ispin, neall, n2ft3d, shift1, shift2, npack1;
   bool aperiodic = false;
   bool periodic = false;
 
public:

   int counter = 0;
 
   /* Constructors */
   Electron_Operators(Pneb *, Kinetic_Operator *, Coulomb12_Operator *,
                      XC_Operator *, Pseudopotential *);
 
   /* destructor */
   ~Electron_Operators() {
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
      if (aperiodic)
         delete[] vlr_l;
      if (mycoulomb12->dielectric_on())
         delete[] vdielec;
   }
 
   //void genrho(double *, double *);
   //void run(double *, double *, double *, double *);
   void genrho(double *psi, double *dn, double *occ = nullptr);
   void run(double *psi, double *dn, double *dng, double *dnall, double *occ = nullptr);

   void get_Tgradient(double *, double *, double *);
   void get_Gradient(double *);
   void gen_Tangent(double *, double *, double *);
   void gen_hml(double *, double *);
   void gen_hmlt(double *, double *);
   void gen_Hpsi_k(double *psi, double *occ = nullptr);
   void gen_psi_r(double *);
   //void gen_density(double *);
   //void gen_densities(double *, double *, double *);
   void gen_density(double *dn, double *occ = nullptr);
   void gen_densities(double *dn, double *dng, double *dnall, double *occ = nullptr);

   void gen_scf_potentials(double *, double *, double *);
   void gen_vl_potential();
   void semicore_density_update();
   void gen_vall();
   void get_vall(double *);
   void set_vall(const double *);
 
   double vl_ave(double *);
   double vlr_ave(double *);
   double vnl_ave(double *psi, double *occ = nullptr);
 
   double energy(double *psi, double *dn, double *dng, double *dnall, double *occ = nullptr);
   double eorbit(double *psi, double *occ = nullptr);
   double ehartree(double *);
   double ehartree2(double *);
   double exc(double *);
   double pxc(double *);
   double eke(double *psi, double *occ = nullptr);
 
   void gen_energies_en(double *psi, double *dn, double *dng,
                       double *dnall, double *E, double *en, double *occ = nullptr);
 
   void add_dteHpsi(double, double *, double *);
 
   void vl_force(double *, double *);
   void vlr_force(double *, double *);
   void vnl_force(double *, double *);
   void semicore_force(double *);
 
   bool is_dielectric_on() { return mycoulomb12->dielectric_on(); }
   void dielectric_force(double *);
 
   bool is_v_apc_on() { return mypsp->myapc->v_apc_on; }
   void apc_force(double *, double *);

   bool is_aperiodic() { return aperiodic; }
   bool is_periodic() { return periodic; }

   Kinetic_Operator *get_myke() {return myke;}
};

} // namespace pwdft

#endif

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
 
   double *Hpsi=nullptr, 
          *psi_r=nullptr, 
          *vl=nullptr, 
          *vall=nullptr, 
          *vc=nullptr, 
          *xcp=nullptr, 
          *xce=nullptr, 
          *x=nullptr, 
          *rho=nullptr, 
          *hmltmp=nullptr;
   double *vcall=nullptr;
 
   double omega, scal2, scal1, dv;
 
   int ispin, neall, nfft3d, shift1, shift2, npack1;
   bool periodic = true;
 
public:

   // counter is incremented when this->run is called
   int counter = 0;
 
   /* Constructors */
   cElectron_Operators(Cneb *, cKinetic_Operator *, cCoulomb_Operator *,
                       cXC_Operator *, CPseudopotential *);
 
   /* destructor */
   ~cElectron_Operators() 
   {
      if (Hpsi)   delete[] Hpsi;
      if (psi_r)  delete[] psi_r;
      if (vl)     delete[] vl;
      if (x)      delete[] x;
      if (rho)    delete[] rho;
      if (vall)   delete[] vall;
      if (vc)     delete[] vc;
      if (vcall)  delete[] vcall;
      if (xcp)    delete[] xcp;
      if (xce)    delete[] xce;
      if (hmltmp) delete[] hmltmp;
      //  delete[] vdielec;
   }
 
   //void genrho(double *, double *);
   void genrho(double *psi, double *dn, double *occ = nullptr);
   //void run(double *, double *, double *, double *);
   void run(double *psi, double *dn, double *dng, double *dnall, double *occ = nullptr);

   void get_Tgradient(double *, double *, double *);
   void get_Gradient(double *);
   void gen_Tangent(double *, double *, double *);
   void gen_hml(double *, double *);
   void gen_hmlt(double *, double *);
   //void gen_Hpsi_k(double *);
   void gen_Hpsi_k(double *psi, double *occ = nullptr);
   void gen_psi_r(double *);
   //void gen_density(double *);
   //void gen_densities(double *, double *, double *);
   void gen_density(double *dn, double *occ = nullptr);
   void gen_densities(double *dn, double *dng, double *dnall, double *occ = nullptr);
   void dn_to_dng_dnall(double *, double *, double *);
   void scf_update_from_dn(double *, double *, double *);


   void gen_scf_potentials(double *, double *, double *);
   void gen_vl_potential();
   void semicore_density_update();
   void gen_vall();
   void get_vall(double *);
   void set_vall(const double *);
 
   double vl_ave(double *);
   //double vnl_ave(double *);
   double vnl_ave(double *psi, double *occ = nullptr);
 
   //double energy(double *, double *, double *, double *);
   //double eorbit(double *);
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
   void vnl_force(double *, double *);
   void semicore_force(double *);
 
   bool is_periodic() { return periodic; }

   cKinetic_Operator *get_myke() {return myke;}
};

} // namespace pwdft

#endif

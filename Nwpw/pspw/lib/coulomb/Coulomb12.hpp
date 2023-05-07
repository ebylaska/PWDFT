#ifndef _COULOMB12_HPP_
#define _COULOMB12_HPP_

#pragma once

#include "Control2.hpp"
#include "Coulomb.hpp"
#include "Coulomb2.hpp"
#include "Ion.hpp"
#include "Pneb.hpp"
#include "Strfac.hpp"

#include "nwpw_dplot.hpp"

namespace pwdft {

class Coulomb12_Operator {
   Pneb *mypneb;


   // Dielectric variables
   bool has_dielec = false;
   bool relax_dielec  = true;
   bool notset_dielec = true;
   bool cube_dielec  = false;
   double dielec, rhomin, rhomax, beta, rho0;
   double filter_dielec = 0.0;
   double *epsilon, *depsilon, *ddepsilon, *sw, *p;
   double *epsilon_x, *epsilon_y, *epsilon_z, *epsilon_lap, *epsilon_screen;
   double *w_x, *w_y, *w_z;
   double *rho_ind0, *rho_ind1;
   double tole_pol  = 1.0e-7;
   double alpha_pol = 0.41;;
   double rcenter_pol[3] = {0.0,0.0,0.0};
   double rmin_pol = 1.2;
   double rmax_pol = 1.6;
   int model_pol = 0;
   int maxit_pol = 2000;

   bool vdielec0_set = false;
   bool rho_ion_set  = false;
   double *rho_ion,*dng_ion,*v_ion,*vdielec0,*vks0;
   double rcut_ion;

   Ion    *myion;
   Strfac *mystrfac;
   Control2 *tmpcontrol;

public:
   double edielec,pdielec;

   bool has_coulomb1 = false;
   Coulomb_Operator *mycoulomb1;

   bool has_coulomb2 = false;
   Coulomb2_Operator *mycoulomb2;

   /* Constructors */
   Coulomb12_Operator(Pneb *, Control2 &);

   /* destructor */
   ~Coulomb12_Operator() {
      // plot cube files before exiting
      if (cube_dielec) {
         nwpw_dplot mydplot(myion,mypneb,*tmpcontrol);
         mydplot.gcube_write("rho_ion.cube", -1, "SCF dielec function", rho_ion);
         mydplot.gcube_write("rho_ind.cube", -1, "SCF dielec function", rho_ind1);
         mydplot.gcube_write("vdielec.cube", -1, "SCF dielec function", vdielec0);
         mydplot.gcube_write("vks.cube", -1, "SCF dielec function", vks0);
         mydplot.gcube_write("epsilon.cube", -1, "SCF dielec function", epsilon);
         mydplot.gcube_write("epsilon_x.cube", -1, "SCF dielec function",epsilon_x);
         mydplot.gcube_write("epsilon_y.cube", -1, "SCF dielec function",epsilon_y);
         mydplot.gcube_write("epsilon_z.cube", -1, "SCF dielec function",epsilon_z);
      }
      if (has_coulomb1)
         delete mycoulomb1;
      if (has_coulomb2)
         delete mycoulomb2;
      if (has_dielec) {
         mypneb->r_dealloc(epsilon);
         mypneb->r_dealloc(depsilon);
         mypneb->r_dealloc(epsilon_screen);
         mypneb->r_dealloc(epsilon_x);
         mypneb->r_dealloc(epsilon_y);
         mypneb->r_dealloc(epsilon_z);
         mypneb->r_dealloc(sw);
         mypneb->r_dealloc(p);
         mypneb->r_dealloc(w_x);
         mypneb->r_dealloc(w_y);
         mypneb->r_dealloc(w_z);
         mypneb->r_dealloc(rho_ind0);
         mypneb->r_dealloc(rho_ind1);
         mypneb->r_dealloc(rho_ion);
         mypneb->r_dealloc(vdielec0);
         mypneb->r_dealloc(vks0);
      }
   }

   bool dielectric_on() { return has_dielec; }
   void initialize_dielectric(Ion *, Strfac *);

   double v_dielectric(const double *, const double *, const double *,
                       const double *, double *);

   void v_dielectric_aperiodic(const double *, const double *, const double *,
                               double *, bool, double *);

   double v_dielectric2_aperiodic(const double *, const double *, const double *,
                                  const double *, const double *, const double *,
                                  const bool, double *, double *, nwpw_dplot *);

   void dielectric_generate(const double *, const double *);
   void dielectric_fion(double *);
   void generate_dng_ion(double *);
   void dng_ion_vdielec0_fion(const double *, double *);

   std::string shortprint_dielectric();

};

} // namespace pwdft

#endif

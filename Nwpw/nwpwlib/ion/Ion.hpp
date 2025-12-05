#ifndef _IONS_HPP_
#define _IONS_HPP_

#include <cmath>
#include <cstdlib>
#include <cstring>

//#include        "iofmt.hpp"
//#include	"rtdb.hpp"
#include "Control2.hpp"
#include "ion_bond.hpp"
//#include "ion_cbond.hpp"
//#include "ion_angle.hpp"
#include "ion_bondings.hpp"
#include "ion_rcovalent.hpp"

namespace pwdft {

class Ion {

  char *atomarray;

  bool has_constraints = false;
  bool has_bond_constraints = false;
  bool has_bondings_constraints = false;
  ion_bond     *mybond;
  //ion_cbond    *mycbond;
  //ion_angle    *myangle;
  ion_bondings *mybondings;


public:
   int ion_total_charge;
   int nion, nkatm; // number of ions
   int *katm;       // element ID number
   int *natm;
   double *charge, *zv_psp;
   double *mass;
   double *dti;
   double *rion0, *rion1, *rion2; // coordinates of ions
   double *rion_incell0;
   double *vionhalf;              // temp velocities
   double *fion1;                 // forces of ions
   double time_step;

   //dispersion and grimme2
   std::string disp_options;
   bool disp_on = false;
   double *ua_disp;

   bool   is_grimme2 = false;
   int    nion_grimme2;
   int    *indx_grimme2;
   double *rion_grimme2;
 
   /* init_ke variables */
   int ke_count, seed, Tf;
   double ekg, eki0, eki1, ke_total, kg_total, mass_total;
   double kb = 3.16679e-6;
   double g_dof = 1.0;
 
   bool fix_translation = true;
   bool dof_translation = false;
 
   bool fix_rotation = false;
   bool dof_rotation = false;

   bool is_crystal = false;

   //Molecular point group symmetry
   double sym_tolerance = 0.001;
   double inertia_tensor[9], inertia_axes[9], inertia_moments[3];
   std::string rotation_type;
   std::string group_name;
   int group_rank;

 
   /* Constructors */
   // Ion(RTDB&, Control2&);
   Ion(std::string, Control2 &);
 
   /* destructor */
   ~Ion() {
     delete[] katm;
     delete[] natm;
     delete[] atomarray;
     delete[] charge;
     delete[] zv_psp;
     delete[] mass;
     delete[] dti;
     delete[] rion0;
     delete[] rion1;
     delete[] rion2;
     delete[] rion_incell0;
     delete[] vionhalf;
     delete[] fion1;
     delete mybond;
     //delete mycbond;
     //delete myangle;
     delete mybondings;
     if (is_grimme2)
     {
        delete[] indx_grimme2;
        delete[] rion_grimme2;
     }
   }
 
   /* functions */
 
   void shift() {
      for (auto i=0; i<(3*nion); ++i)
         rion0[i] = rion1[i];
      for (auto i=0; i<(3*nion); ++i)
         rion1[i] = rion2[i];
   }
   void shift21() {
      // for (auto i=0; i<(3*nion); ++i) rion1[i] = rion2[i];
      std::memcpy(rion1, rion2, 3 * nion * sizeof(double));
   }
   void vshift() {
      for (auto i = 0; i < (3 * nion); ++i)
         rion0[i] = vionhalf[i];
   }
 
   char *symbol(const int i) { return &atomarray[3*katm[i]]; }
   char *atom(const int ia) { return &atomarray[3*ia]; }
   double amu(const int i) { return mass[i]/1822.89; }
   void writefilename(std::string &);
   void writejsonstr(std::string &);
   double rion(int i, int ii) { return rion1[3*ii+i]; }
   double vion(int i, int ii) { return rion0[3*ii+i]; }
   double fion(int i, int ii) { return fion1[3*ii+i]; }
   double rion_incell(int i, int ii) { return rion_incell0[3*ii+i]; }
   void set_rion_incell(const int ,double *);
 
   int ndof() {
      int dof = 3*nion - 6;
      if (dof_translation) dof += 3;
      if (dof_rotation)    dof += 3;
      if (dof < 1)         dof = 1;
      return dof;
   }
   void remove_com_translation();
   void remove_rotation();
 
   double total_mass() {
      double tmass = 0.0;
      for (auto ii = 0; ii < nion; ++ii)
         tmass += mass[ii];
      return tmass;
   }
   double gc(int i) {
      double ss = 0.0;
      for (auto ii = 0; ii < nion; ++ii)
         ss += rion1[3 * ii + i];
      return ss / ((double)nion);
   }
   double vgc(int i) {
      double ss = 0.0;
      for (auto ii = 0; ii < nion; ++ii)
         ss += rion0[3 * ii + i];
      return ss / ((double)nion);
   }
   double com(int i) {
      double tmass = 0.0;
      double ss = 0.0;
      for (auto ii = 0; ii < nion; ++ii) {
         tmass += mass[ii];
        ss += mass[ii] * rion1[3 * ii + i];
      }
      return ss / tmass;
   }
   double vcom(int i) {
      double tmass = 0.0;
      double ss = 0.0;
      for (auto ii = 0; ii < nion; ++ii) {
         tmass += mass[ii];
         ss += mass[ii] * rion0[3 * ii + i];
      }
      return ss / tmass;
   }
 
   double ke() {
      double eki = 0.0;
      for (auto ii = 0; ii < nion; ++ii) {
         eki += 0.5*mass[ii]*( rion0[3*ii]*rion0[3*ii] 
                             + rion0[3*ii+1]*rion0[3*ii+1] 
                             + rion0[3*ii+2]*rion0[3*ii+2]);
      }
      return eki;
   }
 
   double ke_com() {
     double tmass = 0.0;
     double vgx = 0.0;
     double vgy = 0.0;
     double vgz = 0.0;
     for (auto ii = 0; ii < nion; ++ii) {
       tmass += mass[ii];
       vgx += mass[ii] * rion0[3 * ii];
       vgy += mass[ii] * rion0[3 * ii + 1];
       vgz += mass[ii] * rion0[3 * ii + 2];
     }
     vgx /= tmass;
     vgy /= tmass;
     vgz /= tmass;
 
     return 0.5 * tmass * (vgx * vgx + vgy * vgy + vgz * vgz);
   }
 
   double Temperature() {
     return (2.0 * (ke_total / ((double)ke_count)) / g_dof / kb);
   }
   double com_Temperature() {
     return (2.0 * (kg_total / ((double)ke_count)) / kb);
   }
 
   double ion_ion_energy() 
   {
      double eion=0.0;
      for (auto ii=0; ii<(nion-1); ++ii) 
      {
         double qii = zv_psp[katm[ii]];
         for (auto jj=ii+1; jj<nion; ++jj) 
         {
            double qjj = zv_psp[katm[jj]];
            double x = rion1[3*ii]   - rion1[3*jj];
            double y = rion1[3*ii+1] - rion1[3*jj+1];
            double z = rion1[3*ii+2] - rion1[3*jj+2];
            double r = std::sqrt(x*x + y*y + z*z);
            if (r>1.0e-6) 
               eion += qii*qjj/r;
         }
      }
      return eion;
   }
 
   void ion_ion_force(double *fion) 
   {
      for (auto ii=0; ii<(nion-1); ++ii) 
      {
         double qii = zv_psp[katm[ii]];
         for (auto jj=ii+1; jj<nion; ++jj) 
         {
            double qjj = zv_psp[katm[jj]];
            double x = rion1[3 * ii] - rion1[3 * jj];
            double y = rion1[3 * ii + 1] - rion1[3 * jj + 1];
            double z = rion1[3 * ii + 2] - rion1[3 * jj + 2];
            double r = std::sqrt(x * x + y * y + z * z);
            if (r > 1.0e-6) 
            {
               double v = qii * qjj / (r * r * r);
               fion[3 * ii] += (x * v);
               fion[3 * ii + 1] += (y * v);
               fion[3 * ii + 2] += (z * v);
              
               fion[3 * jj] -= (x * v);
               fion[3 * jj + 1] -= (y * v);
               fion[3 * jj + 2] -= (z * v);
            }
         }
      }
   }
 
   void optimize_step(const double *fion) {
     for (auto ii = 0; ii < nion; ++ii) {
       double alpha = sqrt(dti[ii]);
       rion2[3 * ii] = rion1[3 * ii] + alpha * fion[3 * ii];
       rion2[3 * ii + 1] = rion1[3 * ii + 1] + alpha * fion[3 * ii + 1];
       rion2[3 * ii + 2] = rion1[3 * ii + 2] + alpha * fion[3 * ii + 2];
     }
   }
 
   void Verlet_step(const double *fion, const double alpha) {
     double sa1 = 1.0 / (2.0 - alpha);
     double sa2 = alpha / (2.0 - alpha);
 
     // Generate R(t+delta) = r2
     for (auto ii = 0; ii < nion; ++ii) {
       double scale = sa1 * dti[ii];
       rion2[3*ii]   = 2*sa1*rion1[3*ii]   - sa2*rion0[3*ii]   + scale*fion[3*ii];
       rion2[3*ii+1] = 2*sa1*rion1[3*ii+1] - sa2*rion0[3*ii+1] + scale*fion[3*ii+1];
       rion2[3*ii+2] = 2*sa1*rion1[3*ii+2] - sa2*rion0[3*ii+2] + scale*fion[3*ii+2];
     }
 
     // remove translation
     if (fix_translation) remove_com_translation();
 
     // remove rotation
     if (fix_rotation) remove_rotation();
 
     // update current velocities - place in r0
     double h = 1.0/(2.0*time_step);
     for (auto i=0; i<(3*nion); ++i)
        rion0[i] = h*(rion2[i] - rion0[i]);
 
     // add current kinetic energies to running averages
     eki1 = this->ke();
     ekg = this->ke_com();
     ++ke_count;
     ke_total += eki1;
     kg_total += ekg;
   }
 
   void vVerlet_step(const double fion[], const double fion1[]) {
     double sa = 1.0 / time_step;
     for (auto ii = 0; ii < nion; ++ii) {
       double scale = 0.5 * sa * dti[ii];
       vionhalf[3 * ii] = rion0[3 * ii] + scale * (fion1[3 * ii] + fion[3 * ii]);
       vionhalf[3 * ii + 1] =
           rion0[3 * ii + 1] + scale * (fion1[3 * ii + 1] + fion[3 * ii + 1]);
       vionhalf[3 * ii + 2] =
           rion0[3 * ii + 2] + scale * (fion1[3 * ii + 2] + fion[3 * ii + 2]);
     }
 
     // ***** impose molecular constraints - need to implement rattle here***
   }
 
   void Newton_step(const double *fion, const double alpha) {
     // Generate R(t+delta) = r2
     for (auto ii = 0; ii < nion; ++ii) {
       double scale = 0.5 * dti[ii];
       rion2[3 * ii] = rion1[3 * ii] + alpha * time_step * rion0[3 * ii] +
                       scale * fion[3 * ii];
       rion2[3 * ii + 1] = rion1[3 * ii + 1] +
                           alpha * time_step * rion0[3 * ii + 1] +
                           scale * fion[3 * ii + 1];
       rion2[3 * ii + 2] = rion1[3 * ii + 2] +
                           alpha * time_step * rion0[3 * ii + 2] +
                           scale * fion[3 * ii + 2];
     }
 
     // remove translation
     if (fix_translation)
       remove_com_translation();
 
     // add current kinetic energies to running averages
     eki1 = this->ke();
     ekg = this->ke_com();
     ++ke_count;
     ke_total += eki1;
     kg_total += ekg;
   }
 
   void Nose_step(const double ssr, const double *fion) {
     // Generate R(t+delta) = r2
     double smr = 2.0*ssr-1.0;
     for (auto ii = 0; ii < nion; ++ii) {
        double scale = ssr * dti[ii];
        rion2[3*ii]   = 2*ssr*rion1[3*ii]   - smr*rion0[3*ii]   + scale*fion[3*ii];
        rion2[3*ii+1] = 2*ssr*rion1[3*ii+1] - smr*rion0[3*ii+1] + scale*fion[3*ii+1];
        rion2[3*ii+2] = 2*ssr*rion1[3*ii+2] - smr*rion0[3*ii+2] + scale*fion[3*ii+2];
     }
 
     // remove translation
     if (fix_translation)
       remove_com_translation();
 
     // remove rotation
     if (fix_rotation)
       remove_rotation();
 
     // update current velocities - place in r0
     double h = 1.0 / (2.0 * time_step);
     for (auto i = 0; i < (3 * nion); ++i)
       rion0[i] = h * (rion2[i] - rion0[i]);
 
     // add current kinetic energies to running averages
     eki1 = this->ke();
     ekg = this->ke_com();
     ++ke_count;
     ke_total += eki1;
     kg_total += ekg;
   }
 
   void fixed_step(const double alpha, const double *fion) {
     for (auto i = 0; i < (3 * nion); ++i)
       rion2[i] = rion1[i] + alpha * fion[i];
   }
 
   void set_zv_psp(const int ia, const double zv) { zv_psp[ia] = zv; }
   double total_zv() {
     double tcharge = 0.0;
     for (auto ia = 0; ia < nkatm; ++ia)
       tcharge += zv_psp[ia] * natm[ia];
     return tcharge;
   }
 
   double xrms() {
     double dx, dy, dz;
     double xx = 0.0;
     for (auto ii = 0; ii < nion; ++ii) {
       dx = rion1[3 * ii] - rion0[3 * ii];
       dy = rion1[3 * ii + 1] - rion0[3 * ii + 1];
       dz = rion1[3 * ii + 2] - rion0[3 * ii + 2];
       xx += dx * dx + dy * dy + dz * dz;
     }
     return sqrt(xx) / ((double)nion);
   }
 
   double xmax() {
     double dx, dy, dz, x, xx;
     double yy = 0.0;
     for (auto ii = 0; ii < nion; ++ii) {
       dx = rion1[3 * ii] - rion0[3 * ii];
       dy = rion1[3 * ii + 1] - rion0[3 * ii + 1];
       dz = rion1[3 * ii + 2] - rion0[3 * ii + 2];
       xx = dx * dx + dy * dy + dz * dz;
       x = sqrt(xx);
       if (x > yy)
         yy = x;
     }
     return yy;
   }
 
   std::string print_bond_angle_torsions() {
      return ion_print_bond_angle_torsions(nion, katm, atomarray, rion1);
   }
 
   double com_fion(const double *fion, const int i) {
     double tmass = 0.0;
     double ss = 0.0;
     for (auto ii = 0; ii < nion; ++ii) {
       tmass += mass[ii];
       ss += mass[ii] * fion[3 * ii + i];
     }
     return ss / tmass;
   }
 
   void remove_com_fion(double *fion) {
     double tmass = 0.0;
     double gx = 0.0;
     double gy = 0.0;
     double gz = 0.0;
     for (auto ii = 0; ii < nion; ++ii) {
       tmass += mass[ii];
       gx += mass[ii] * fion[3 * ii];
       gy += mass[ii] * fion[3 * ii + 1];
       gz += mass[ii] * fion[3 * ii + 2];
     }
     gx /= tmass;
     gy /= tmass;
     gz /= tmass;
     for (auto ii = 0; ii < nion; ++ii) {
       fion[3 * ii] -= gx;
       fion[3 * ii + 1] -= gy;
       fion[3 * ii + 2] -= gz;
     }
   }
 
   double max_fion(const double *fion) {
     double x, xx;
     double y = 0.0;
     for (auto ii = 0; ii < nion; ++ii) {
       xx = fion[3 * ii] * fion[3 * ii] + fion[3 * ii + 1] * fion[3 * ii + 1] +
            fion[3 * ii + 2] * fion[3 * ii + 2];
       x = sqrt(xx);
       if (x > y)
         y = x;
     }
     return y;
   }
 
   double rms_fion(const double *fion) {
     double ss = 0.0;
     for (auto ii = 0; ii < nion; ++ii) {
       ss += fion[3*ii]*fion[3*ii] + fion[3*ii+1]*fion[3*ii+1] + fion[3*ii+2]*fion[3*ii+2];
     }
     return sqrt(ss) / ((double)nion);
   }

   // ion contraints access functions
   bool has_ion_constraints() {return has_constraints;}
   bool has_ion_bond_constraints()     {return has_bond_constraints;}
   bool has_ion_bondings_constraints() {return has_bondings_constraints;}
   double energy_ion_bond_constraints()    { return mybond->energy();}
   double energy_ion_bondings_constraints() { return mybondings->energy();}
   std::string print_constraints(const int);

   void add_contraint_force(double *fion) {
      if (has_constraints) {
          double econstraint = 0.0;
         if (has_bond_constraints)     econstraint += mybond->energyfion(fion);
         if (has_bondings_constraints) econstraint += mybondings->energyfion(fion);
      }
   }
   double add_contraint_energyforce(double *fion) {
      double econstraint = 0.0;
      if (has_constraints) {
         if (has_bond_constraints)     econstraint += mybond->energyfion(fion);
         if (has_bondings_constraints) econstraint += mybondings->energyfion(fion);
      }
      return econstraint;
   }



   // ion dispersion access functions
   bool has_ion_disp() {return disp_on;}
   double disp_energy();
   void disp_force(double *);
   void disp_stress(double *);
   
   std::string print_symmetry_group();


};
} // namespace pwdft

#endif

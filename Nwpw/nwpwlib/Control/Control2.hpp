#ifndef _CONTROL2_HPP_
#define _CONTROL2_HPP_

#pragma once

/* Control2.hpp
   Author - Eric Bylaska
*/

#include <algorithm>
#include <map>
#include <string>
#include <vector>
#include <iostream>

namespace pwdft {

class Control2 {

   std::string myrtdbstring, xcstring;
   double punita[9], ptolerances[3], pscaling[2];
   double ptime_step, pfake_mass, pscf_alpha, pscf_beta, pecut, pwcut, prcut;
   double pkerker_g0,pfractional_kT,pfractional_temperature,pfractional_alpha;
   double pfractional_gamma,pfractional_alpha_min,pfractional_alpha_max,pfractional_beta,pfractional_rmsd_threshold,pfractional_rmsd_tolerance;
   double pbo_time_step;
   double ptotal_charge;
   double peprecondition, psprecondition;
   std::vector<double> pfractional_filling;
 
   bool pnolagrange   = false;
   bool puse_grid_cmp = false;
   
   bool pfast_erf = false;
 
   bool pdeltae_check = true;
   bool pis_crystal = false;

   bool ptwodfractional = false;
   bool pfractional = false;
   bool pfractional_frozen = false;
   int pfractional_smeartype;
   int pfractional_orbitals[2];
   int pks_maxit_orb,pks_maxit_orbs;
   int pdiis_histories;
 
   int pbo_steps[2], pbo_algorithm;
   int ploop[2], pngrid[3], pnpsp, pncut, pmapping, pmapping1d, ptile_factor;
  
   int pnp_dimensions[3], pewald_grid[3];
   int pcode, ptask;
   int pispin, pmultiplicity, pne[2], pnexcited[2], ptotal_ion_charge, plmax_multipole;
   int pmove, pfrac_coord, pgram_schmidt;
   int protation, ptranslation, pbalance, pspin_orbit;
   int pmaxit_orb, pmaxit_orbs, pscf_algorithm, pks_algorithm;
   int psymm_number;
   int pqsize;
   int pprint_level = 2;
 
   int pminimizer = 1;
   int plmbfgs_size = 2;
   int pinitial_psi_random_algorithm = 1;
 
   int pdriver_maxiter = 30;
   int pdriver_lmbfgs_size = 9;
   double pdriver_gmax = 0.00045;
   double pdriver_grms = 0.00030;
   double pdriver_xrms = 0.00120;
   double pdriver_xmax = 0.00180;
   double pdriver_trust = 0.3;
 
   bool pgeometry_optimize;
 
   bool pinput_movecs_initialize = false;
   char ppermanent_dir[256];
   char pscratch_dir[256];
   char pinput_movecs_filename[80];
   char poutput_movecs_filename[80];
   char pinput_v_movecs_filename[80];
   char poutput_v_movecs_filename[80];
   char pinput_e_movecs_filename[80];
   char poutput_e_movecs_filename[80];

   bool pefield_on;
   int pefield_type;
   std::vector<double> pefield_vector, pefield_center;
 
   bool papc_on;
   int papc_nga;
   double papc_Gc;
   std::vector<double> papc_gamma, papc_u, papc_q;
 
   bool pborn_on, pborn_relax;
   double pborn_dielec = 78.4;
   double pborn_rcut = 0.1;
   std::vector<double> pborn_bradii, pborn_vradii;
 
   bool pnose_on, pnose_restart;
   int pnose_mchain, pnose_nchain;
   double pnose_pe, pnose_pr, pnose_te, pnose_tr, pnose_ne_chain, pnose_eke0;
   std::vector<double> pnose_xem, pnose_xe0, pnose_xe1, pnose_qe, pnose_xrm,
       pnose_xr0, pnose_xr1, pnose_qr;
 
   bool pgpoisson_on = false;
   bool pgpoisson_relax_dielec = true;
   bool pgpoisson_cube_dielec = false;
   int    pgpoisson_model  = 0;
   int    pgpoisson_maxit  = 2000;
   double pgpoisson_dielec = 78.4;
   double pgpoisson_filter = 0.0;
   double pgpoisson_rho0 = 0.0004;
   double pgpoisson_beta = 1.3;
   double pgpoisson_rhomin = 0.0001;
   double pgpoisson_rhomax = 0.0035;
   double pgpoisson_alpha  = 0.41;
   double pgpoisson_rcut_ion = 1.0;
   double pgpoisson_rmin = 1.0;
   double pgpoisson_rmax = 2.2;
   double pgpoisson_rcenter[3] = {0.0,0.0,0.0};
 
   bool psa_on;
   double psa_decay[2] = {1.0, 1.0};
 
   bool pfei_on;
   bool pcif_on;
   bool pcif_shift_cell = true;
   bool pdipole_on;
   bool pmulliken_on;

   bool pstaged_gpu_fft = true;
   int pfft_container_size = 1;

   // vdw variables - resetable
   std::string poptions_disp;
   bool phas_disp = false;
   bool phas_vdw = false;
   bool pis_vdw2 = false;
   bool pis_grimme2 = false;

   // HFX variables - resetable
   bool pHFX = false;
   bool pHFX_relax = false;
   int    pHFX_screening_type = 0;
   double pHFX_parameter = 0.0;
   double pHFX_screening_radius = 0.0;

   int pio_norbs_max = 100;
   int pio_buffer = true;

   // Brillouin variables 
   int pnbrillouin=0;

   //scf extra scf
   bool pscf_extra_rotate = false;

   //pspspin variables
   bool p_pspspin = false;

   //psputerm variables
   bool p_psputerm = false;
   int  p_pspnuterms = 0;
   

public:
   int version = 3;

   /* constructor */
   Control2(const int, const std::string);
 
   /* public variables */
   std::string psp_library_dir;
   std::map<std::string,std::string> psp_libraries;
 
   std::string xyz_filename,ion_motion_filename,emotion_filename,fei_filename,
               cif_filename,omotion_filename,hmotion_filename,eigmotion_filename,
               dipole_motion_filename;
   std::string permanent_dir_str,scratch_dir_str;
 
   // Access functions
   double * unita_ptr() { return punita; }
   double unita(const int i, const int j) { return punita[i + j * 3]; }
   double unita1d(const int ii) { return punita[ii]; }
   double tolerances(const int i) { return ptolerances[i]; }
 
   double scaling(const int i) { return pscaling[i]; }
   double elc_scaling() { return pscaling[0]; }
   double ion_scaling() {
     if (ptask == 7)
       return pscaling[0];
     else
       return pscaling[1];
   }
 
   double time_step() { return ptime_step; }
   double bo_time_step() { return pbo_time_step; }
   double fake_mass() { return pfake_mass; }
   double ecut() { return pecut; }
   double wcut() { return pwcut; }
   double ewald_rcut() { return prcut; }
   double total_charge() { return ptotal_charge; }
   double Eprecondition() { return peprecondition;}
   double Sprecondition() { return psprecondition;}

   double scf_alpha() { return pscf_alpha; }
   double scf_beta() { return pscf_beta; }
   double kerker_g0() { return pkerker_g0; }
   double fractional_kT() { return pfractional_kT; }
   double fractional_temperature() { return pfractional_temperature; }
   double fractional_alpha() { return pfractional_alpha; }
   double fractional_alpha_min() { return pfractional_alpha_min; }
   double fractional_alpha_max() { return pfractional_alpha_max; }
   double fractional_beta() { return pfractional_beta; }
   double fractional_gamma() { return pfractional_gamma; }
   double fractional_rmsd_threshold() { return pfractional_rmsd_threshold; }
   double fractional_rmsd_tolerance() { return pfractional_rmsd_tolerance; }
 
   int minimizer() { return pminimizer; }
   int lmbfgs_size() { return plmbfgs_size; }
   int ks_algorithm() { return pks_algorithm; }
   int scf_algorithm() { return pscf_algorithm; }
   int fractional_smeartype() { return pfractional_smeartype; }
   int ks_maxit_orb() { return pks_maxit_orb; }
   int ks_maxit_orbs() { return pks_maxit_orbs; }
   int diis_histories() { return pdiis_histories; }
   int task() { return ptask; }
   int np_orbital() { return pnp_dimensions[1]; }
   int np_dimensions(const int i) { return pnp_dimensions[i]; }
   int mapping1d() { return pmapping1d; }
   int mapping() { return pmapping; }
   int tile_factor() { return ptile_factor; }
   int balance() { return pbalance; }
   int ngrid(const int i) { return pngrid[i]; }
   int loop(const int i) { return ploop[i]; }
   int bo_steps(const int i) { return pbo_steps[i]; }
   int bo_algorithm() { return pbo_algorithm; }
   int pfft3_qsize() { return pqsize; }
   int ewald_ngrid(const int i) { return pewald_grid[i]; }
   int ewald_ncut() { return pncut; }
   int multiplicity() { return pmultiplicity; }
   int ispin() { return pispin; }
   int ne(const int i) { return pne[i]; }
   int nexcited(const int i) { return pnexcited[i]; }
   int total_ion_charge() { return ptotal_ion_charge; }
   int lmax_multipole() { return plmax_multipole; }
   std::string xc_name() { return xcstring; }
   int initial_psi_random_algorithm() { return pinitial_psi_random_algorithm; }
   int io_norbs_max() { return pio_norbs_max; }
   bool io_buffer() { return pio_buffer; }
   bool twodfractional() { return ptwodfractional; }
   bool fractional() { return pfractional; }
   bool fractional_frozen() { return pfractional_frozen; }
   int  fractional_orbitals(const int i) { return pfractional_orbitals[i]; }

   bool pspspin()   { return p_pspspin; }
   bool psputerm()  { return p_psputerm; }
   int  pspnuterms() { return p_pspnuterms; }

   bool scf_extra_rotate() { return pscf_extra_rotate; }
 
   int *ne_ptr() { return pne; }

   int nbrillouin() { return pnbrillouin; }
 
   bool deltae_check() { return pdeltae_check; } 
   bool geometry_optimize() { return pgeometry_optimize; }
   bool use_grid_cmp() { return puse_grid_cmp; }
   bool fast_erf() { return pfast_erf; }
   bool is_crystal() { return pis_crystal; }

   bool nolagrange() { return pnolagrange; }
 
   int driver_maxiter() { return pdriver_maxiter; }
   int driver_lmbfgs_size() { return pdriver_lmbfgs_size; }
   double driver_gmax() { return pdriver_gmax; }
   double driver_grms() { return pdriver_grms; }
   double driver_xmax() { return pdriver_xmax; }
   double driver_xrms() { return pdriver_xrms; }
   double driver_trust() { return pdriver_trust; }
 
   bool input_movecs_initialize() { return pinput_movecs_initialize; }
   char *input_movecs_filename() { return pinput_movecs_filename; }
   char *output_movecs_filename() { return poutput_movecs_filename; }
   char *input_v_movecs_filename() { return pinput_v_movecs_filename; }
   char *output_v_movecs_filename() { return poutput_v_movecs_filename; }
   char *input_e_movecs_filename() { return pinput_e_movecs_filename; }
   char *output_e_movecs_filename() { return poutput_e_movecs_filename; }
   char *permanent_dir() { return ppermanent_dir; }
   char *scratch_dir() { return pscratch_dir; }
 
   void add_permanent_dir(char *);
   void add_scratch_dir(char *);
 
   bool out_of_time() { return false; } // needs to be implemented

   bool staged_gpu_fft() { return pstaged_gpu_fft; }
   int  fft_container_size() { return pfft_container_size; }
 
   int print_level(std::string plevel) {
     int doprt = 0;
     if (plevel == "low") {
       doprt = (pprint_level > 0);
 
     } else if (plevel == "medium") {
       doprt = (pprint_level > 1);
 
     } else if (plevel == "high") {
       doprt = (pprint_level > 2);
 
     } else if (plevel == "debug") {
       doprt = (pprint_level > 3);
     }
     return doprt;
   }
 
   void set_total_ion_charge(int icharge) 
   {
      ptotal_ion_charge = icharge;
     
      /* x = total number of electrons */
      int x = (int)(ptotal_ion_charge - ptotal_charge);

     
      /* reassign ispin to agree with total number electrons - odd number of
       * electrons ==> ispin=2 */
      if (((x % 2) != 0) && (pispin == 1))
         pispin = 2;
     
      /* reassign multiplicity to  agree with total number electrons */
      /*  -- odd number of electrons and mult odd */
      if (((x % 2) != 0) && ((pmultiplicity % 2) != 0)) 
      {
         pmultiplicity -= 1;
         while (pmultiplicity > (x + 1))
            pmultiplicity -= 2;
         if (pmultiplicity < 1)
            pmultiplicity = 2;
      }
      /*  -- even number of electrons and mult even */
      if (((x % 2) == 0) && ((pmultiplicity % 2) == 0)) 
      {
         pmultiplicity -= 1;
         while (pmultiplicity > (x + 1))
            pmultiplicity -= 2;
         if (pmultiplicity < 1)
            pmultiplicity = 1;
      }
     
      /* assign ne */
      if (pispin == 1) 
      {
         pne[0] = x / 2;
         pne[1] = 0;
      } 
      else 
      {
         int dx = pmultiplicity - 1;
         pne[0] = (x + dx) / 2;
         pne[1] = (x - dx) / 2;
      }
   }
 
   // Efield
   bool Efield_on() { return pefield_on; }
   int Efield_type() { return pefield_type; }
   double Efield_center(const int i) {
     if (i >= pefield_center.size())
       return 0.0;
     else
       return pefield_center[i];
   }
   double Efield_vector(const int i) {
     if (i >= pefield_vector.size())
       return 0.0;
     else
       return pefield_vector[i];
   }
 
   // APC
   bool APC_on() { return papc_on; }
   int APC_nga() { return papc_nga; }
   double APC_Gc() { return papc_Gc; }
   double APC_gamma(const int i) { return papc_gamma[i]; }
   double APC_q(const int i) {
     if (i >= papc_q.size())
       return 0.0;
     else
       return papc_q[i];
   }
   double APC_u(const int i) {
     if (i >= papc_u.size())
       return 0.0;
     else
       return papc_u[i];
   }
 
   // Born
   bool born_on() { return pborn_on; }
   bool born_relax() { return pborn_relax; }
   double born_bradii(const int i) {
     return ((i >= pborn_bradii.size()) ? 0.0 : pborn_bradii[i]);
   }
   double born_vradii(const int i) {
     return ((i >= pborn_vradii.size()) ? 0.0 : pborn_vradii[i]);
   }
   double born_dielec() { return pborn_dielec; }
   double born_rcut() { return pborn_rcut; }
 
   // Generalized Poisson
   bool gpoisson_on() { return pgpoisson_on; }
   bool gpoisson_relax_dielec() { return pgpoisson_relax_dielec; }
   bool gpoisson_cube_dielec() { return pgpoisson_cube_dielec; }
   int  gpoisson_model() { return pgpoisson_model; }
   int  gpoisson_maxit() { return pgpoisson_maxit; }
   double gpoisson_dielec() { return pgpoisson_dielec; }
   double gpoisson_filter() { return pgpoisson_filter; }
   double gpoisson_rho0() { return pgpoisson_rho0; }
   double gpoisson_beta() { return pgpoisson_beta; }
   double gpoisson_rhomin() { return pgpoisson_rhomin; }
   double gpoisson_rhomax() { return pgpoisson_rhomax; }
   double gpoisson_alpha()    { return pgpoisson_alpha; }
   double gpoisson_rcut_ion() { return pgpoisson_rcut_ion; }
   double gpoisson_rmin() { return pgpoisson_rmin; }
   double gpoisson_rmax() { return pgpoisson_rmax; }
   double gpoisson_rcenter(const int i) { return pgpoisson_rcenter[i]; }
 
 
   // Nose
   bool Nose() { return pnose_on; }
   bool Nose_restart() { return pnose_restart; }
   int Nose_Mchain() { return pnose_mchain; }
   int Nose_Nchain() { return pnose_nchain; }
   double Nose_Te() { return pnose_te; }
   double Nose_Tr() { return pnose_tr; }
   double Nose_Pe() { return pnose_pe; }
   double Nose_Pr() { return pnose_pr; }
   double Nose_eke0() { return pnose_eke0; }
   double Nose_Ne_chain() { return pnose_ne_chain; }
   double Nose_Xem(const int i) {
     return ((i >= pnose_xem.size()) ? 0.0 : pnose_xem[i]);
   }
   double Nose_Xe0(const int i) {
     return ((i >= pnose_xe0.size()) ? 0.0 : pnose_xe0[i]);
   }
   double Nose_Xe1(const int i) {
     return ((i >= pnose_xe1.size()) ? 0.0 : pnose_xe1[i]);
   }
   double Nose_Qe(const int i) {
     return ((i >= pnose_qe.size()) ? 0.0 : pnose_qe[i]);
   }
   double Nose_Xrm(const int i) {
     return ((i >= pnose_xrm.size()) ? 0.0 : pnose_xrm[i]);
   }
   double Nose_Xr0(const int i) {
     return ((i >= pnose_xr0.size()) ? 0.0 : pnose_xr0[i]);
   }
   double Nose_Xr1(const int i) {
     return ((i >= pnose_xr1.size()) ? 0.0 : pnose_xr1[i]);
   }
   double Nose_Qr(const int i) {
     return ((i >= pnose_qr.size()) ? 0.0 : pnose_qr[i]);
   }
 
   // SA - simulated annealing
   bool SA() { return psa_on; }
   double SA_decay(const int i) { return psa_decay[i]; }
 
   // Fei
   bool Fei_on() { return pfei_on; }
 
   // CIF
   bool CIF_on() { return pcif_on; }
   bool CIF_shift_cell() { return pcif_shift_cell; }
 
   // Mulliken
   bool Mulliken_on() { return pmulliken_on; }
 
   // Dipole
   bool dipole_on() { return pdipole_on; }
 
   // GGA value
   int get_gga() {
      std::string myxc_name = this->xc_name();
     
      std::transform(myxc_name.begin(), myxc_name.end(), myxc_name.begin(),
                     ::tolower);
     
      int gga = 0;
     
      if (myxc_name.compare("vosko") == 0)
        gga = 0;
      if (myxc_name.compare("lda") == 0)
        gga = 0;
     
      if (myxc_name.compare("pbe") == 0)
        gga = 10;
      if (myxc_name.compare("pbe96") == 0)
        gga = 10;
      if (myxc_name.compare("blyp") == 0)
        gga = 11;
      if (myxc_name.compare("revpbe") == 0)
        gga = 12;
      if (myxc_name.compare("pbesol") == 0)
        gga = 13;
      if (myxc_name.compare("hser") == 0)
        gga = 14;
      if (myxc_name.compare("b3lypr") == 0)
        gga = 15;
      if (myxc_name.compare("beef") == 0)
        gga = 16;
      if (myxc_name.compare("xbeef-cpbe") == 0)
        gga = 17;
     
      if (myxc_name.compare("pbe0") == 0)
         gga = 110;
      if (myxc_name.compare("blyp0") == 0)
         gga = 111;
      if (myxc_name.compare("revpbe0") == 0)
         gga = 112;
      if (myxc_name.compare("bnl") == 0)
         gga = 113;
      if (myxc_name.compare("hse") == 0)
         gga = 114;
      if (myxc_name.compare("b3lyp") == 0)
         gga = 115;
      if (myxc_name.compare("hartree-fock") == 0)
         gga = 200;
      if (myxc_name.compare("hf") == 0)
         gga = 200;
     
      return gga;
   }

   // vdw
   std::string options_disp() { return poptions_disp; }
   bool has_disp()           { return phas_disp; }
   bool has_vdw()            { return phas_vdw; }
   bool is_vdw2()            { return pis_vdw2; }
   bool is_grimme2()         { return pis_grimme2; }

 
   // cubefiles
   int number_cubefiles();
   int cubetype_cubefiles(const int);
   std::string cubename_cubefiles(const int);
   std::string cubekey_cubefiles(const int);
   double position_tolerance_cubefiles();
   double origin_cubefiles(const int);
   int ncell_cubefiles(const int);

   // bond
   int nhb_bond();
   int    i0_bond(const int);
   int    j0_bond(const int);
   double Kspring0_bond(const int);
   double R0_bond(const int);

   // cbond
   int nhcb_cbond();
   int    i0_cbond(const int);
   int    j0_cbond(const int);
   int    k0_cbond(const int);
   int    l0_cbond(const int);
   double Kspring0_cbond(const int);
   double Rij0_cbond(const int);
   double Rkl0_cbond(const int);

   // angle
   int nha_angle();
   int    i0_angle(const int);
   int    j0_angle(const int);
   int    k0_angle(const int);
   double Kspring0_angle(const int);
   double Theta0_angle(const int);

   // bondings
   int nhc_bondings();
   std::vector<double> coef_bondings(const int);
   std::vector<int> indx_bondings(const int);
   double Kspring0_bondings(const int);
   double gamma0_bondings(const int);

   // smear
   std::vector<double> fractional_filling() {return pfractional_filling;}

   // pspspin
   std::string set_pspspin(int, double *, double *, int *, int *, int *, int *, bool *, bool *);

   // psputerm
   std::string set_psputerm(int, int, int *, double *, double *, bool *);

   // remove virtual from rtdbstring
   //void remove_virtual(){
  // }

};

} // namespace pwdft

#endif

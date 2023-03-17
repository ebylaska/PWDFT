#ifndef _CONTROL2_HPP_
#define _CONTROL2_HPP_
/* Control2.hpp
   Author - Eric Bylaska
*/

#include        <algorithm>
#include	<map>
#include        <string>
#include        <vector>



namespace pwdft {



class Control2 {

   std::string myrtdbstring,xcstring;
   double punita[9],ptolerances[3],pscaling[2];
   double ptime_step,pfake_mass,pks_alpha,pecut,pwcut,prcut;
   double pbo_time_step;
   double ptotal_charge;

   bool   puse_grid_cmp = false;;
   bool   pfast_erf = false;

   int pbo_steps[2],pbo_algorithm;
   int ploop[2],pngrid[3],pnpsp,pncut,pmapping,pmapping1d,ptile_factor;;
   int pnp_dimensions[3],pewald_grid[3];
   int pcode,ptask;
   int pispin,pmultiplicity,pne[2],ptotal_ion_charge,plmax_multipole;
   int pmove,pfrac_coord,pgram_schmidt;
   int protation,ptranslation,pbalance,pspin_orbit;
   int pmaxit_orb,pmaxit_orbs,pscf_algorithm,pks_algorithm;
   int psymm_number;
   int pqsize;
   int pprint_level = 2;

   int pminimizer = 1;
   int plmbfgs_size = 2;
   int  pinitial_psi_random_algorithm = 1;


   int    pdriver_maxiter     = 30;
   int    pdriver_lmbfgs_size = 9;
   double pdriver_gmax  = 0.00045;
   double pdriver_grms  = 0.00030;
   double pdriver_xrms  = 0.00120;
   double pdriver_xmax  = 0.00180;
   double pdriver_trust = 0.3;

   bool pgeometry_optimize;


   bool pinput_movecs_initialize = false;
   char ppermanent_dir[256];
   char pscratch_dir[256];
   char pinput_movecs_filename[80];
   char poutput_movecs_filename[80];
   char pinput_v_movecs_filename[80];
   char poutput_v_movecs_filename[80];

   bool   pefield_on;
   int    pefield_type;
   std::vector<double> pefield_vector,pefield_center;

   bool   papc_on;
   int    papc_nga;
   double papc_Gc;
   std::vector<double> papc_gamma,papc_u,papc_q;

   bool   pborn_on,pborn_relax;
   double pborn_dielec = 78.4;
   double pborn_rcut   =  0.1;
   std::vector<double> pborn_bradii,pborn_vradii;

   bool   pnose_on,pnose_restart;
   int    pnose_mchain,pnose_nchain;
   double pnose_pe,pnose_pr,pnose_te,pnose_tr,pnose_ne_chain,pnose_eke0;
   std::vector<double> pnose_xem,pnose_xe0,pnose_xe1,pnose_qe,
                       pnose_xrm,pnose_xr0,pnose_xr1,pnose_qr;

   bool   psa_on;
   double psa_decay[2] = {1.0,1.0};

   bool   pfei_on;
   bool   pcif_on;
   bool   pdipole_on;
   bool   pmulliken_on;

public:

   int version = 3;

   /* constructor */
   Control2(const int, const std::string);

   /* public variables */
   std::string psp_library_dir;
   std::map<std::string, std::string> psp_libraries;

   std::string xyz_filename,
               ion_motion_filename,
               emotion_filename,
               fei_filename,
               cif_filename,
               omotion_filename,
               hmotion_filename,
               eigmotion_filename,
               dipole_motion_filename;
   std::string permanent_dir_str,scratch_dir_str;


   // Access functions
   double unita(const int i,const int j) { return punita[i+j*3]; }
   double unita1d(const int ii)   { return punita[ii]; }
   double tolerances(const int i) { return ptolerances[i]; }

   double scaling(const int i)    { return pscaling[i]; }
   double elc_scaling()           { return pscaling[0]; }
   double ion_scaling()           { if (ptask==7) return pscaling[0]; else return pscaling[1]; }

   double time_step()             { return ptime_step; }
   double bo_time_step()          { return pbo_time_step; }
   double fake_mass()             { return pfake_mass; }
   double ecut()                  { return pecut; }
   double wcut()                  { return pwcut; }
   double ewald_rcut()            { return prcut; }
   double total_charge()          { return ptotal_charge; }

   int minimizer()             { return pminimizer; }
   int lmbfgs_size()           { return plmbfgs_size; }
   int task()                  { return ptask; }
   int np_orbital()            { return pnp_dimensions[1]; }
   int mapping1d()             { return pmapping1d; }
   int mapping()               { return pmapping; }
   int tile_factor()           { return ptile_factor; }
   int balance()               { return pbalance; }
   int ngrid(const int i)      { return pngrid[i]; }
   int loop(const int i)       { return ploop[i]; }
   int bo_steps(const int i)   { return pbo_steps[i]; }
   int bo_algorithm()          { return pbo_algorithm; }
   int pfft3_qsize()           { return pqsize; }
   int ewald_ngrid(const int i) { return pewald_grid[i]; }
   int ewald_ncut()             { return pncut; }
   int multiplicity()           { return pmultiplicity; }
   int ispin()                  { return pispin; }
   int ne(const int i)          { return pne[i]; }
   int total_ion_charge()       { return ptotal_ion_charge; }
   int lmax_multipole()         { return plmax_multipole; }
   std::string xc_name()	{ return xcstring; }
   int initial_psi_random_algorithm() { return pinitial_psi_random_algorithm; }

   int *ne_ptr()                { return pne; }

   bool geometry_optimize()     { return pgeometry_optimize; }
   bool use_grid_cmp()          { return puse_grid_cmp; }
   bool fast_erf()              { return pfast_erf; }

   int driver_maxiter()    { return pdriver_maxiter; }
   int driver_lmbfgs_size(){ return pdriver_lmbfgs_size; }
   double driver_gmax()    { return pdriver_gmax; }
   double driver_grms()    { return pdriver_grms; }
   double driver_xmax()    { return pdriver_xmax; }
   double driver_xrms()    { return pdriver_xrms; }
   double driver_trust()   { return pdriver_trust; }

   bool   input_movecs_initialize() { return pinput_movecs_initialize; }
   char   *input_movecs_filename() { return pinput_movecs_filename; }
   char   *output_movecs_filename() { return poutput_movecs_filename; }
   char   *input_v_movecs_filename() { return pinput_v_movecs_filename; }
   char   *output_v_movecs_filename() { return poutput_v_movecs_filename; }
   char   *permanent_dir() { return ppermanent_dir;}
   char   *scratch_dir() { return pscratch_dir;}


   void add_permanent_dir(char *);
   void add_scratch_dir(char *);

   bool out_of_time() { return false; } // needs to be implemented

   int  print_level(std::string plevel) {
      int doprt = 0;
      if (plevel=="low") {
         doprt = (pprint_level>0);

      } else if (plevel=="medium") {
         doprt = (pprint_level>1);

      } else if (plevel=="high") {
         doprt =  (pprint_level>2);

      } else if (plevel=="debug") {
         doprt =  (pprint_level>3);
      }
      return doprt;
   }


   void set_total_ion_charge(int icharge)
   {
      ptotal_ion_charge = icharge;

      /* x = total number of electrons */
      int x = (int) (ptotal_ion_charge - ptotal_charge);

      /* reassign ispin to agree with total number electrons - odd number of electrons ==> ispin=2 */
      if (((x%2)!=0) && (pispin==1)) pispin = 2;

      /* reassign multiplicity to  agree with total number electrons */
      /*  -- odd number of electrons and mult odd */
      if (((x%2)!=0) && ((pmultiplicity%2)!=0))
      {
         pmultiplicity -= 1;
         while (pmultiplicity>(x+1)) 
            pmultiplicity -= 2;
         if (pmultiplicity<1)
            pmultiplicity = 2;
      }
      /*  -- even number of electrons and mult even */
      if (((x%2)==0) && ((pmultiplicity%2)==0))
      {
         pmultiplicity -= 1;
         while (pmultiplicity>(x+1)) 
            pmultiplicity -= 2;
         if (pmultiplicity<1)
            pmultiplicity = 1;
      }

      /* assign ne */
      if (pispin==1) 
      {
         pne[0] = x/2;
         pne[1] = 0;
      }
      else
      {
         int dx = pmultiplicity - 1;
         pne[0] = (x+dx)/2;
         pne[1] = (x-dx)/2;
      }
   }

   // Efield
   bool Efield_on()    { return pefield_on; }
   int Efield_type()    { return pefield_type; }
   double Efield_center(const int i){ 
      if (i>=pefield_center.size())
         return 0.0;
      else
         return pefield_center[i];
   }
   double Efield_vector(const int i){ 
      if (i>=pefield_vector.size())
         return 0.0;
      else
         return pefield_vector[i];
   }

   // APC
   bool APC_on()    { return papc_on; }
   int APC_nga()    { return papc_nga; }
   double APC_Gc()  { return papc_Gc; }
   double APC_gamma(const int i)  { return papc_gamma[i]; }
   double APC_q(const int i){ 
      if (i>=papc_q.size())
         return 0.0;
      else
         return papc_q[i];
   }
   double APC_u(const int i)      { 
      if (i>=papc_u.size())
         return 0.0;
      else
         return papc_u[i];
   }

   // Born
   bool born_on()    { return pborn_on; }
   bool born_relax() { return pborn_relax; }
   double born_bradii(const int i) { return ((i>=pborn_bradii.size()) ? 0.0 : pborn_bradii[i]); }
   double born_vradii(const int i) { return ((i>=pborn_vradii.size()) ? 0.0 : pborn_vradii[i]); }
   double born_dielec() { return pborn_dielec; }
   double born_rcut()   { return pborn_rcut; }


   // Nose
   bool Nose() { return pnose_on; }
   bool Nose_restart() { return pnose_restart; }
   int  Nose_Mchain()  { return pnose_mchain; }
   int  Nose_Nchain()  { return pnose_nchain; }
   double Nose_Te()    { return pnose_te; }
   double Nose_Tr()    { return pnose_tr; }
   double Nose_Pe()    { return pnose_pe; }
   double Nose_Pr()    { return pnose_pr; }
   double Nose_eke0()  { return pnose_eke0; }
   double Nose_Ne_chain() { return pnose_ne_chain; }
   double Nose_Xem(const int i) { return ((i>=pnose_xem.size()) ? 0.0 : pnose_xem[i]); }
   double Nose_Xe0(const int i) { return ((i>=pnose_xe0.size()) ? 0.0 : pnose_xe0[i]); }
   double Nose_Xe1(const int i) { return ((i>=pnose_xe1.size()) ? 0.0 : pnose_xe1[i]); }
   double Nose_Qe(const int i)  { return ((i>=pnose_qe.size())  ? 0.0 : pnose_qe[i]);  }
   double Nose_Xrm(const int i) { return ((i>=pnose_xrm.size()) ? 0.0 : pnose_xrm[i]); }
   double Nose_Xr0(const int i) { return ((i>=pnose_xr0.size()) ? 0.0 : pnose_xr0[i]); }
   double Nose_Xr1(const int i) { return ((i>=pnose_xr1.size()) ? 0.0 : pnose_xr1[i]); }
   double Nose_Qr(const int i)  { return ((i>=pnose_qr.size())  ? 0.0 : pnose_qr[i]);  }

   // SA - simulated annealing
   bool SA() { return psa_on; }
   double SA_decay(const int i) { return psa_decay[i]; }

   // Fei
   bool Fei_on() { return pfei_on; }

   // CIF
   bool CIF_on() { return pcif_on; }

   // Mulliken
   bool Mulliken_on() { return pmulliken_on; }

   // Dipole
   bool dipole_on() { return pdipole_on; }



   //GGA value
   int get_gga() {
      std::string myxc_name = this->xc_name();

      std::transform(myxc_name.begin(), myxc_name.end(), myxc_name.begin(), ::tolower);

      int gga = 0;

      if (myxc_name.compare("vosko") == 0) gga=0;
      if (myxc_name.compare("lda")   == 0) gga=0;

      if (myxc_name.compare("pbe")        == 0) gga=10;
      if (myxc_name.compare("pbe96")      == 0) gga=10;
      if (myxc_name.compare("blyp")       == 0) gga=11;
      if (myxc_name.compare("revpbe")     == 0) gga=12;
      if (myxc_name.compare("pbesol")     == 0) gga=13;
      if (myxc_name.compare("hser")       == 0) gga=14;
      if (myxc_name.compare("b3lypr")     == 0) gga=15;
      if (myxc_name.compare("beef")       == 0) gga=16;
      if (myxc_name.compare("xbeef-cpbe") == 0) gga=17;

      if (myxc_name.compare("pbe0")    == 0) gga=110;
      if (myxc_name.compare("blyp0")   == 0) gga=111;
      if (myxc_name.compare("revpbe0") == 0) gga=112;
      if (myxc_name.compare("bnl")     == 0) gga=113;
      if (myxc_name.compare("hse")     == 0) gga=114;
      if (myxc_name.compare("b3lyp")   == 0) gga=115;

      if (myxc_name.compare("hartree-fock") == 0) gga=200;
      if (myxc_name.compare("hf")           == 0) gga=200;

      return gga;
   }

   // cubefiles
   int number_cubefiles();
   int cubetype_cubefiles(const int);
   std::string cubename_cubefiles(const int);
   std::string cubekey_cubefiles(const int);


};

}

#endif

#ifndef _CONTROL2_HPP_
#define _CONTROL2_HPP_
/* Control2.hpp
   Author - Eric Bylaska
*/

#include	<map>
#include        <string>
#include        <vector>

using namespace std;

namespace pwdft {
using namespace pwdft;


class Control2 {

   string myrtdbstring,xcstring;
   double punita[9],ptolerances[3],pscaling[2];
   double ptime_step,pfake_mass,pks_alpha,pecut,pwcut,prcut;
   double pbo_time_step;
   double ptotal_charge;

   int pbo_steps[2],pbo_algorithm;
   int ploop[2],pngrid[3],pnpsp,pncut,pmapping,pmapping1d;
   int pnp_dimensions[3],pewald_grid[3];
   int pcode,ptask;
   int pispin,pmultiplicity,pne[2],ptotal_ion_charge;
   int pmove,pfrac_coord,pSA,pfei,pfei_quench,pgram_schmidt;
   int protation,ptranslation,pbalance,pspin_orbit;
   int pmaxit_orb,pmaxit_orbs,pscf_algorithm,pks_algorithm;
   int psymm_number;
   int pqsize;
   int pprint_level = 2;

   int pminimizer = 1;


   int    pdriver_maxiter     = 30;
   int    pdriver_lmbfgs_size = 9;
   double pdriver_gmax  = 0.00045;
   double pdriver_grms  = 0.00030;
   double pdriver_xrms  = 0.00120;
   double pdriver_xmax  = 0.00180;
   double pdriver_trust = 0.3;

   bool pgeometry_optimize;


   char ppermanent_dir[256];
   char pscratch_dir[256];
   char pinput_movecs_filename[80];
   char poutput_movecs_filename[80];
   char pinput_v_movecs_filename[80];
   char poutput_v_movecs_filename[80];

   bool   papc_on;
   int    papc_nga;
   double papc_Gc;
   std::vector<double> papc_gamma,papc_u,papc_q;


public:

   /* constructor */
   Control2(const int, const string);

   /* public variables */
   string psp_library_dir;
   std::map<std::string, std::string> psp_libraries;

   // Access functions
   double unita(const int i,const int j) { return punita[i+j*3]; }
   double unita1d(const int ii)   { return punita[ii]; }
   double tolerances(const int i) { return ptolerances[i]; }
   double scaling(const int i)    { return pscaling[i]; }
   double time_step()             { return ptime_step; }
   double fake_mass()             { return pfake_mass; }
   double ecut()                  { return pecut; }
   double wcut()                  { return pwcut; }
   double ewald_rcut()            { return prcut; }
   double total_charge()          { return ptotal_charge; }

   int minimizer()             { return pminimizer; }
   int task()                  { return ptask; }
   int np_orbital()            { return pnp_dimensions[1]; }
   int mapping1d()             { return pmapping1d; }
   int mapping()               { return pmapping; }
   int balance()               { return pbalance; }
   int ngrid(const int i)      { return pngrid[i]; }
   int loop(const int i)       { return ploop[i]; }
   int pfft3_qsize()           { return pqsize; }
   int ewald_ngrid(const int i) { return pewald_grid[i]; }
   int ewald_ncut()             { return pncut; }
   int multiplicity()           { return pmultiplicity; }
   int ispin()                  { return pispin; }
   int ne(const int i)          { return pne[i]; }
   int total_ion_charge()       { return ptotal_ion_charge; }
   string xc_name()		{ return xcstring; }

   int *ne_ptr()                { return pne; }

   bool geometry_optimize()       { return pgeometry_optimize; }

   int driver_maxiter()    { return pdriver_maxiter; }
   int driver_lmbfgs_size(){ return pdriver_lmbfgs_size; }
   double driver_gmax()    { return pdriver_gmax; }
   double driver_grms()    { return pdriver_grms; }
   double driver_xmax()    { return pdriver_xmax; }
   double driver_xrms()    { return pdriver_xrms; }
   double driver_trust()   { return pdriver_trust; }

   char   *input_movecs_filename() { return pinput_movecs_filename; }
   char   *output_movecs_filename() { return poutput_movecs_filename; }
   char   *input_v_movecs_filename() { return pinput_v_movecs_filename; }
   char   *output_v_movecs_filename() { return poutput_v_movecs_filename; }
   char   *permanent_dir() { return ppermanent_dir;}
   char   *scratch_dir() { return pscratch_dir;}

   void add_permanent_dir(char *);
   void add_scratch_dir(char *);

   int  print_level(string plevel) {
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

};

}

#endif

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "json.hpp"
//#include	"rtdb.hpp"
#include "Control2.hpp"
#include "Parallel.hpp"

#include "parsestring.hpp"

using json = nlohmann::json;

namespace pwdft {

/***********************************
 *                                 *
 *         factor_count2           *
 *                                 *
 ***********************************/
static int factor_count2(int n, int m) {
  int f = 0;
  int nn = n;
  while ((nn % m) == 0) {
    nn /= m;
    ++f;
  }
  return f;
}
static int dum_ipow(int a, int n) {
  int npa = 1;
  for (auto i = 0; i < n; ++i)
    npa *= a;
  return npa;
}

/***********************************
 *                                 *
 *         control_set_ngrid       *
 *                                 *
 ***********************************/
/* return n so that it is a multiple of 2,3,5,7
 */
static int control_set_ngrid(double x, bool mult2) {
  int nx = (int)floor(x + 0.5);
  if (mult2 && ((nx % 2) != 0))
    ++nx;

  int nf2 = factor_count2(nx, 2);
  int nf3 = factor_count2(nx, 3);
  int nf5 = factor_count2(nx, 5);
  int nf7 = factor_count2(nx, 7);
  int ntest =
      dum_ipow(2, nf2) * dum_ipow(3, nf3) * dum_ipow(5, nf5) * dum_ipow(7, nf7);
  while (nx != ntest) {
    ++nx;
    if (mult2)
      ++nx;
    nf2 = factor_count2(nx, 2);
    nf3 = factor_count2(nx, 3);
    nf5 = factor_count2(nx, 5);
    nf7 = factor_count2(nx, 7);
    ntest = dum_ipow(2, nf2) * dum_ipow(3, nf3) * dum_ipow(5, nf5) *
            dum_ipow(7, nf7);
  }

  return nx;
}

/***********************************
 *                                 *
 *      control_ngrid_default      *
 *                                 *
 ***********************************/
/* defines the default ngrid given unita, ecut and mapping
 */
static void control_ngrid_default(double *unita, double ecut, int mapping,
                                  int *ngrid) {
  double gx, gy, gz, xh, yh, zh, unitg[9];
  double twopi = 8.0 * atan(1.0);
  unitg[0] = unita[4] * unita[8] - unita[5] * unita[7];
  unitg[1] = unita[5] * unita[6] - unita[3] * unita[8];
  unitg[2] = unita[3] * unita[7] - unita[4] * unita[6];
  unitg[3] = unita[7] * unita[2] - unita[8] * unita[1];
  unitg[4] = unita[8] * unita[0] - unita[6] * unita[2];
  unitg[5] = unita[6] * unita[1] - unita[7] * unita[0];
  unitg[6] = unita[1] * unita[5] - unita[2] * unita[4];
  unitg[7] = unita[2] * unita[3] - unita[0] * unita[5];
  unitg[8] = unita[0] * unita[4] - unita[1] * unita[3];
  double volume =
      unita[0] * unitg[0] + unita[1] * unitg[1] + unita[2] * unitg[2];
  for (int i = 0; i < 9; ++i)
    unitg[i] *= twopi / volume;

  gx = unitg[0];
  gy = unitg[1];
  gz = unitg[2];
  xh = sqrt(2.00 * ecut / (gx * gx + gy * gy + gz * gz)) + 0.5;

  gx = unitg[3];
  gy = unitg[4];
  gz = unitg[5];
  yh = sqrt(2.00 * ecut / (gx * gx + gy * gy + gz * gz)) + 0.5;

  gx = unitg[6];
  gy = unitg[7];
  gz = unitg[8];
  zh = sqrt(2.00 * ecut / (gx * gx + gy * gy + gz * gz)) + 0.5;

  ngrid[0] = control_set_ngrid(2.0 * xh, true);
  ngrid[1] = control_set_ngrid(2.0 * yh, true);
  ngrid[2] = control_set_ngrid(2.0 * zh, true);
  if (mapping == 1) {
    if (ngrid[1] > ngrid[2])
      ngrid[2] = ngrid[1];
    else
      ngrid[1] = ngrid[2];
  }
}

/***********************************
 *                                 *
 *         Constructors            *
 *                                 *
 ***********************************/

Control2::Control2(const int np0, const std::string rtdbstring) 
{
   myrtdbstring = rtdbstring;
   json rtdbjson = json::parse(rtdbstring);
 
   int np = np0;
   bool is_cpmd;
 
   ptotal_ion_charge = 0;
   pne[0] = 0;
   pne[1] = 0;
   pnexcited[0] = 0;
   pnexcited[1] = 0;
   pfractional_orbitals[0] = 0;
   pfractional_orbitals[1] = 0;
 
   pfei_on = false;
   pcif_on = false;
   pcif_shift_cell = true;
   pmulliken_on = false;
   pdipole_on = false;
 
   ptask = 0;
   if (!rtdbjson["current_task"].is_null()) {
     if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),
                           "energy"))
       ptask = 1;
     if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),
                           "gradient"))
       ptask = 2;
     if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),
                           "optimize"))
       ptask = 3;
     if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]), "freq"))
       ptask = 4;
     if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),
                           "steepest_descent"))
       ptask = 5;
     if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),
                           "car-parrinello"))
       ptask = 6;
     if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),
                           "born-oppenheimer"))
       ptask = 7;
     if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),
                           "noit_"))
       ptask *= -1;
   }
 
   ptotal_ion_charge = -1.0;
   ptotal_charge = 0.0;
   // if (rtdbjson["charge"].is_number_float()) ptotal_charge =
   // rtdbjson["charge"];
   if (rtdbjson["charge"].is_number_integer() ||
       rtdbjson["charge"].is_number_float())
     ptotal_charge = rtdbjson["charge"];
 
   /* get parallel mappings */
   //pmapping = 1;
   pmapping = 3;
   if (rtdbjson["nwpw"]["mapping"].is_number_integer())
     pmapping = rtdbjson["nwpw"]["mapping"];
 
   /* set mapping1d */
   pmapping1d = 1;
   if (rtdbjson["nwpw"]["mapping1d"].is_number_integer())
     pmapping1d = rtdbjson["nwpw"]["mapping1d"];
 
   /* set ptile_factor for linear algebra gpu pipeline */
   ptile_factor = 1;
   if (rtdbjson["nwpw"]["tile_factor"].is_number_integer())
     ptile_factor = rtdbjson["nwpw"]["tile_factor"];
 
   /* set pinitial_psi_random_algorithm for wavefunction initialization */
   pinitial_psi_random_algorithm = 1;
   if (rtdbjson["nwpw"]["initial_psi_random_algorithm"].is_number_integer())
      pinitial_psi_random_algorithm = rtdbjson["nwpw"]["initial_psi_random_algorithm"];
 
   /* qsize */
   pqsize = 5;
   if (rtdbjson["nwpw"]["pfft3_qsize"].is_number_integer())
     pqsize = rtdbjson["nwpw"]["pfft3_qsize"];
 
   /* get np_dimensions */
   // if (!myrtdb.get("nwpw:np_dimensions",rtdb_int,3,np_dimensions))
   if (rtdbjson["nwpw"]["np_dimensions"].is_null()) 
   {
      pnp_dimensions[0] = np;
      pnp_dimensions[1] = 1;
      pnp_dimensions[2] = 1;
   } 
   else 
   {
      pnp_dimensions[0] = rtdbjson["nwpw"]["np_dimensions"][0];
      pnp_dimensions[1] = rtdbjson["nwpw"]["np_dimensions"][1];
      pnp_dimensions[2] = rtdbjson["nwpw"]["np_dimensions"][2];

      if (pnp_dimensions[2] < 1) pnp_dimensions[2] = 1;
      if (pnp_dimensions[1] < 1) pnp_dimensions[1] = 1;
      if (pnp_dimensions[0] < 1) pnp_dimensions[0] = 1;

      int  nbrill = 1;
      if (!rtdbjson["nwpw"]["brillouin_zone"]["kvectors"].is_null())
           nbrill = rtdbjson["nwpw"]["brillouin_zone"]["kvectors"].size();

      // reset np_dimensions(3) if larger than nbrill 
      if (pnp_dimensions[2]>nbrill) pnp_dimensions[2] = nbrill;

      // reset np_dimensions(2) if it is not a  multiple of np 
      while (np % pnp_dimensions[2] != 0 && pnp_dimensions[2] > 1) {
         pnp_dimensions[2]--;
      }
     
      /* reset np_dimensions[1] if it is not a  multiple of np2 */
      np = np / pnp_dimensions[2];
      while (((np % pnp_dimensions[1]) != 0) && (pnp_dimensions[1] > 1))
         pnp_dimensions[1] = pnp_dimensions[1] - 1;
      pnp_dimensions[0] = np / pnp_dimensions[1];
   }

   if (rtdbjson["nwpw"]["io_norbs_max"].is_number_integer())
      pio_norbs_max = rtdbjson["nwpw"]["io_norbs_max"];
 
   puse_grid_cmp = false;
   if (rtdbjson["nwpw"]["use_grid_cmp"].is_boolean())
      puse_grid_cmp = rtdbjson["nwpw"]["use_grid_cmp"];
 
   pfast_erf = false;
   if (rtdbjson["nwpw"]["fast_erf"].is_boolean())
      pfast_erf = rtdbjson["nwpw"]["fast_erf"];
 
   plmax_multipole = 0;
   if (rtdbjson["nwpw"]["lmax_multipole"].is_number_integer())
      plmax_multipole = rtdbjson["nwpw"]["lmax_multipole"];

   pnolagrange = false;
   if (rtdbjson["nwpw"]["nolagrange"].is_boolean())
      pnolagrange = rtdbjson["nwpw"]["nolagrange"];
 
   pbalance = 1;
   if (rtdbjson["nwpw"]["nobalance"].is_boolean())
     if (rtdbjson["nwpw"]["nobalance"])
       pbalance = 0;
 
   pgeometry_optimize = false;
   if (rtdbjson["nwpw"]["steepest_descent"]["geometry_optimize"].is_boolean())
     pgeometry_optimize =
         rtdbjson["nwpw"]["steepest_descent"]["geometry_optimize"];
 
   psp_library_dir = "";
   if (rtdbjson["psp_library_dir"].is_string())
     psp_library_dir = rtdbjson["psp_library_dir"];
 
   if (!rtdbjson["nwpw"]["pseudopotentials"].is_null())
     for (auto &el : rtdbjson["nwpw"]["pseudopotentials"].items()) {
       psp_libraries[el.key()] = el.value();
     }
 
   pprint_level = 2;
   if (rtdbjson["print"].is_string()) {
     if (mystring_contains(mystring_lowercase(rtdbjson["print"]), "debug"))
       pprint_level = 4;
     if (mystring_contains(mystring_lowercase(rtdbjson["print"]), "high"))
       pprint_level = 3;
     if (mystring_contains(mystring_lowercase(rtdbjson["print"]), "medium"))
       pprint_level = 2;
     if (mystring_contains(mystring_lowercase(rtdbjson["print"]), "low"))
       pprint_level = 1;
     if (mystring_contains(mystring_lowercase(rtdbjson["print"]), "off"))
       pprint_level = 0;
     if (mystring_contains(mystring_lowercase(rtdbjson["print"]), "none"))
       pprint_level = 0;
   }
 
   permanent_dir_str = "";
   scratch_dir_str = "";
   if (rtdbjson["permanent_dir"].is_string())
      permanent_dir_str = rtdbjson["permanent_dir"];
 
   if (rtdbjson["scratch_dir"].is_string())
      scratch_dir_str = rtdbjson["scratch_dir"];

   std::string input_movecs = "eric.movecs";
   std::string output_movecs = "eric.movecs";
   std::string input_v_movecs = "eric.vmovecs";
   std::string output_v_movecs = "eric.vmovecs";
   std::string input_e_movecs = "eric.emovecs";
   std::string output_e_movecs = "eric.emovecs";
   if (rtdbjson["dbname"].is_string()) {
      std::string dbname = rtdbjson["dbname"];
      input_movecs = dbname + ".movecs";
      output_movecs = dbname + ".movecs";
      input_v_movecs = dbname + ".vmovecs";
      output_v_movecs = dbname + ".vmovecs";
      input_e_movecs = dbname + ".emovecs";
      output_e_movecs = dbname + ".emovecs";
   }
 
   // read from nwpw block
   if (rtdbjson["nwpw"]["input_wavefunction_filename"].is_string())
      input_movecs = rtdbjson["nwpw"]["input_wavefunction_filename"];
   if (rtdbjson["nwpw"]["output_wavefunction_filename"].is_string())
      output_movecs = rtdbjson["nwpw"]["output_wavefunction_filename"];
   if (rtdbjson["nwpw"]["input_v_wavefunction_filename"].is_string())
      input_v_movecs = rtdbjson["nwpw"]["input_v_wavefunction_filename"];
   if (rtdbjson["nwpw"]["output_v_wavefunction_filename"].is_string())
      output_v_movecs = rtdbjson["nwpw"]["output_v_wavefunction_filename"];
   if (rtdbjson["nwpw"]["input_e_wavefunction_filename"].is_string())
      input_e_movecs = rtdbjson["nwpw"]["input_e_wavefunction_filename"];
   if (rtdbjson["nwpw"]["output_e_wavefunction_filename"].is_string())
      output_e_movecs = rtdbjson["nwpw"]["output_e_wavefunction_filename"];
 
   // from car-parrinello block
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["input_wavefunction_filename"]
             .is_string())
       input_movecs =
           rtdbjson["nwpw"]["car-parrinello"]["input_wavefunction_filename"];
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["output_wavefunction_filename"]
             .is_string())
       output_movecs =
           rtdbjson["nwpw"]["car-parrinello"]["output_wavefunction_filename"];
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["input_v_wavefunction_filename"]
             .is_string())
       input_v_movecs =
           rtdbjson["nwpw"]["car-parrinello"]["input_v_wavefunction_filename"];
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["output_v_wavefunction_filename"]
             .is_string())
       output_v_movecs =
           rtdbjson["nwpw"]["car-parrinello"]["output_v_wavefunction_filename"];
 
   if (rtdbjson["dbname"].is_string()) {
     std::string dbname = rtdbjson["dbname"];
     xyz_filename = dbname + ".xyz";
     ion_motion_filename = dbname + ".ion_motion";
     emotion_filename = dbname + ".emotion";
     fei_filename = dbname + ".fei";
     cif_filename = dbname + ".cif";
     omotion_filename = dbname + ".omotion";
     hmotion_filename = dbname + ".hmotion";
     fei_filename = dbname + ".fei";
     eigmotion_filename = dbname + ".eigmotion";
     dipole_motion_filename = dbname + ".dipole_motion";
   }
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["xyz_filename"].is_string())
       xyz_filename = rtdbjson["nwpw"]["car-parrinello"]["xyz_filename"];
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["motion_filename"].is_string())
       ion_motion_filename = rtdbjson["nwpw"]["car-parrinello"]["ion_motion_filename"];
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["emotion_filename"].is_string())
       emotion_filename = rtdbjson["nwpw"]["car-parrinello"]["emotion_filename"];
 
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["omotion_filename"].is_string())
       omotion_filename = rtdbjson["nwpw"]["car-parrinello"]["omotion_filename"];
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["hmotion_filename"].is_string())
       hmotion_filename = rtdbjson["nwpw"]["car-parrinello"]["hmotion_filename"];
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["eigmotion_filename"].is_string())
       eigmotion_filename =
           rtdbjson["nwpw"]["car-parrinello"]["eigmotion_filename"];
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["fei_filename"].is_string())
       fei_filename = rtdbjson["nwpw"]["car-parrinello"]["fei_filename"];
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["cif_filename"].is_string())
       cif_filename = rtdbjson["nwpw"]["car-parrinello"]["cif_filename"];
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["dipole_motion_filename"]
             .is_string())
       dipole_motion_filename =
           rtdbjson["nwpw"]["car-parrinello"]["dipole_motion_filename"];
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["fei_on"].is_boolean())
       pfei_on = rtdbjson["nwpw"]["car-parrinello"]["fei_on"];
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["cif_on"].is_boolean())
       pcif_on = rtdbjson["nwpw"]["car-parrinello"]["cif_on"];
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["cif_shift_cell"].is_boolean())
       pcif_shift_cell = rtdbjson["nwpw"]["car-parrinello"]["cif_shift_cell"];

   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["dipole_on"].is_boolean())
       pdipole_on = rtdbjson["nwpw"]["car-parrinello"]["dipole_on"];
   // if (ptask==6) if (rtdbjson["nwpw"]["mulliken"].is_boolean()) pmulliken_on =
   // rtdbjson["nwpw"]["mulliken_on"];
 
   // from steepest_descent block
   if (ptask == 5)
     if (rtdbjson["nwpw"]["steepest_descent"]["input_wavefunction_filename"]
             .is_string())
       input_movecs =
           rtdbjson["nwpw"]["steepest_descent"]["input_wavefunction_filename"];
   if (ptask == 5)
     if (rtdbjson["nwpw"]["steepest_descent"]["output_wavefunction_filename"]
             .is_string())
       output_movecs =
           rtdbjson["nwpw"]["steepest_descent"]["output_wavefunction_filename"];
 
   if (permanent_dir_str.size() > 0) 
   {
      if (input_movecs[0] != '/')    input_movecs    = permanent_dir_str + "/" + input_movecs;
      if (output_movecs[0] != '/')   output_movecs   = permanent_dir_str + "/" + output_movecs;
      if (input_v_movecs[0] != '/')  input_v_movecs  = permanent_dir_str + "/" + input_v_movecs;
      if (output_v_movecs[0] != '/') output_v_movecs = permanent_dir_str + "/" + output_v_movecs;
      if (input_e_movecs[0] != '/')  input_e_movecs  = permanent_dir_str + "/" + input_e_movecs;
      if (output_e_movecs[0] != '/') output_e_movecs = permanent_dir_str + "/" + output_e_movecs;
   }
 
   strcpy(pinput_movecs_filename,const_cast<char *>(input_movecs.data()));
   strcpy(poutput_movecs_filename,const_cast<char *>(output_movecs.data()));
   strcpy(pinput_v_movecs_filename,const_cast<char *>(input_v_movecs.data()));
   strcpy(poutput_v_movecs_filename,const_cast<char *>(output_v_movecs.data()));
   strcpy(pinput_e_movecs_filename,const_cast<char *>(input_e_movecs.data()));
   strcpy(poutput_e_movecs_filename,const_cast<char *>(output_e_movecs.data()));

   strcpy(ppermanent_dir,const_cast<char *>(permanent_dir_str.data()));
   strcpy(pscratch_dir,const_cast<char *>(scratch_dir_str.data()));
   // strcpy(ppsp_library_dir,const_cast<char*>(psp_library_dir_str.data()));
 
   // force wavefunction initializion
   if (rtdbjson["nwpw"]["initialize_wavefunction"].is_boolean())
      pinput_movecs_initialize = rtdbjson["nwpw"]["initialize_wavefunction"];
 
   pfake_mass = 400000.0;
   if (rtdbjson["nwpw"]["fake_mass"].is_number_float())
      pfake_mass = rtdbjson["nwpw"]["fake_mass"];
   if (ptask == 5)
      if (rtdbjson["nwpw"]["steepest_descent"]["fake_mass"].is_number_float())
         pfake_mass = rtdbjson["nwpw"]["steepest_descent"]["fake_mass"];
   if (ptask == 6)
      if (rtdbjson["nwpw"]["car-parrinello"]["fake_mass"].is_number_float())
         pfake_mass = rtdbjson["nwpw"]["car-parrinello"]["fake_mass"];
 
   ptime_step = 5.8;
   if (rtdbjson["nwpw"]["time_step"].is_number_float())
     ptime_step = rtdbjson["nwpw"]["time_step"];
   if (ptask == 5)
     if (rtdbjson["nwpw"]["steepest_descent"]["time_step"].is_number_float())
       ptime_step = rtdbjson["nwpw"]["steepest_descent"]["time_step"];
   if (ptask == 6)
     if (rtdbjson["nwpw"]["car-parrinello"]["time_step"].is_number_float())
       ptime_step = rtdbjson["nwpw"]["car-parrinello"]["time_step"];

   if (rtdbjson["nwpw"]["virtual"][0].is_number_integer())
      pnexcited[0] = rtdbjson["nwpw"]["virtual"][0];
   if (rtdbjson["nwpw"]["virtual"][1].is_number_integer())
      pnexcited[1] = rtdbjson["nwpw"]["virtual"][1];
 
   pscf_algorithm = 0;
   if (rtdbjson["nwpw"]["scf_algorithm"].is_number_integer())
       pscf_algorithm = rtdbjson["nwpw"]["scf_algorithm"];

   if (rtdbjson["nwpw"]["fractional_smeartype"].is_number_integer())
       pfractional_smeartype = rtdbjson["nwpw"]["fractional_smeartype"];


   if (rtdbjson["nwpw"]["fractional_orbitals"][0].is_number_integer())
      pfractional_orbitals[0] = rtdbjson["nwpw"]["fractional_orbitals"][0];
   if (rtdbjson["nwpw"]["fractional_orbitals"][1].is_number_integer())
      pfractional_orbitals[1] = rtdbjson["nwpw"]["fractional_orbitals"][1];


   pks_maxit_orb = 5;
   if (rtdbjson["nwpw"]["ks_maxit_orb"].is_number_integer())
       pks_maxit_orb = rtdbjson["nwpw"]["ks_maxit_orb"];

   pks_maxit_orbs = 1;
   if (rtdbjson["nwpw"]["ks_maxit_orbs"].is_number_integer())
       pks_maxit_orbs = rtdbjson["nwpw"]["ks_maxit_orbs"];

   pdiis_histories = 15;
   if (rtdbjson["nwpw"]["diis_histories"].is_number_integer())
       pdiis_histories = rtdbjson["nwpw"]["diis_histories"];

   pscf_alpha = 0.25;
   if (rtdbjson["nwpw"]["scf_alpha"].is_number_float())
      pscf_alpha = rtdbjson["nwpw"]["scf_alpha"];

   pscf_beta = 0.25;
   if (rtdbjson["nwpw"]["scf_beta"].is_number_float())
      pscf_beta = rtdbjson["nwpw"]["scf_beta"];

   pkerker_g0 = 0.0;
   if (rtdbjson["nwpw"]["kerker_g0"].is_number_float())
      pkerker_g0 = rtdbjson["nwpw"]["kerker_g0"];

   if (rtdbjson["nwpw"]["fractional_kT"].is_number_float())
      pfractional_kT = rtdbjson["nwpw"]["fractional_kT"];
   if (rtdbjson["nwpw"]["fractional_temperature"].is_number_float())
      pfractional_temperature = rtdbjson["nwpw"]["fractional_temperature"];

   pfractional_alpha = 0.5;
   if (rtdbjson["nwpw"]["fractional_alpha"].is_number_float())
      pfractional_alpha = rtdbjson["nwpw"]["fractional_alpha"];

   // Adaptive alpha parameters
   pfractional_alpha_min = 0.1;
   if (rtdbjson["nwpw"]["fractional_alpha_min"].is_number_float())
      pfractional_alpha_min = rtdbjson["nwpw"]["fractional_alpha_min"];

   pfractional_alpha_max = 0.5;
   if (rtdbjson["nwpw"]["fractional_alpha_max"].is_number_float())
      pfractional_alpha_max = rtdbjson["nwpw"]["fractional_alpha_max"];


   pfractional_beta = 0.1;
   if (rtdbjson["nwpw"]["fractional_beta"].is_number_float())
      pfractional_beta = rtdbjson["nwpw"]["fractional_beta"];

   pfractional_gamma = 0.2;
   if (rtdbjson["nwpw"]["fractional_gamma"].is_number_float())
      pfractional_gamma = rtdbjson["nwpw"]["fractional_gamma"];

   pfractional_rmsd_threshold = 1.0e-3;
   if (rtdbjson["nwpw"]["fractional_rmsd_threshold"].is_number_float())
      pfractional_rmsd_threshold = rtdbjson["nwpw"]["fractional_rmsd_threshold"];

   pfractional_rmsd_tolerance = 1.0e-3;
   if (rtdbjson["nwpw"]["fractional_rmsd_tolerance"].is_number_float())
      pfractional_rmsd_tolerance = rtdbjson["nwpw"]["fractional_rmsd_tolerance"];

   if (rtdbjson["nwpw"]["fractional_orbitals"][0].is_number_integer())


   pfractional_filling = {};
   if (!rtdbjson["nwpw"]["fractional_filling"].is_null())
      pfractional_filling = rtdbjson["nwpw"]["fractional_filling"].get<std::vector<double>>();

   pfractional_frozen = false;
   if (rtdbjson["nwpw"]["fractional_frozen"].is_boolean())
      pfractional_frozen = rtdbjson["nwpw"]["fractional_frozen"];

   pfractional = false;
   if (rtdbjson["nwpw"]["fractional"].is_boolean())
      pfractional = rtdbjson["nwpw"]["fractional"];

   ptwodfractional = false;
   if (rtdbjson["nwpw"]["twodfractional"].is_boolean())
      ptwodfractional = rtdbjson["nwpw"]["twodfractional"];

   if (rtdbjson["nwpw"]["scf_extra_rotate"].is_boolean())
      pscf_extra_rotate = rtdbjson["nwpw"]["scf_extra_rotate"];


   ploop[0] = 10;
   ploop[1] = 100;
   if (rtdbjson["nwpw"]["loop"][0].is_number_integer())
      ploop[0] = rtdbjson["nwpw"]["loop"][0];
   if (rtdbjson["nwpw"]["loop"][1].is_number_integer())
      ploop[1] = rtdbjson["nwpw"]["loop"][1];
 
   pbo_time_step = 5.0;
   if (rtdbjson["nwpw"]["bo_time_step"].is_number_float())
      pbo_time_step = rtdbjson["nwpw"]["bo_time_step"];
 
   pbo_steps[0] = 10;
   pbo_steps[1] = 100;
   if (rtdbjson["nwpw"]["bo_steps"][0].is_number_integer())
      pbo_steps[0] = rtdbjson["nwpw"]["bo_steps"][0];
   if (rtdbjson["nwpw"]["bo_steps"][1].is_number_integer())
      pbo_steps[1] = rtdbjson["nwpw"]["bo_steps"][1];
 
   pbo_algorithm = 0;
   if (rtdbjson["nwpw"]["bo_algorithm"].is_number_integer())
     pbo_algorithm = rtdbjson["nwpw"]["bo_algorithm"];
 
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["loop"][0].is_number_integer())
      ploop[0] = rtdbjson["nwpw"]["steepest_descent"]["loop"][0];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["loop"][1].is_number_integer())
      ploop[1] = rtdbjson["nwpw"]["steepest_descent"]["loop"][1];
   if (ptask==6) if (rtdbjson["nwpw"]["car-parrinello"]["loop"][0].is_number_integer())
      ploop[0] = rtdbjson["nwpw"]["car-parrinello"]["loop"][0];
   if (ptask==6) if (rtdbjson["nwpw"]["car-parrinello"]["loop"][1].is_number_integer())
      ploop[1] = rtdbjson["nwpw"]["car-parrinello"]["loop"][1];
 
   pscaling[0] = 1.0;
   pscaling[1] = 1.0;
   if (ptask==6) if (rtdbjson["nwpw"]["car-parrinello"]["scaling"][0].is_number_float())
      pscaling[0] = rtdbjson["nwpw"]["car-parrinello"]["scaling"][0];
   if (ptask==6) if (rtdbjson["nwpw"]["car-parrinello"]["scaling"][0].is_number_integer())
      pscaling[0] = rtdbjson["nwpw"]["car-parrinello"]["scaling"][0];
   if (ptask==6) if (rtdbjson["nwpw"]["car-parrinello"]["scaling"][1].is_number_float())
      pscaling[1] = rtdbjson["nwpw"]["car-parrinello"]["scaling"][1];
   if (ptask==6) if (rtdbjson["nwpw"]["car-parrinello"]["scaling"][1].is_number_integer())
      pscaling[1] = rtdbjson["nwpw"]["car-parrinello"]["scaling"][1];
 
   if (ptask==7)
      if (rtdbjson["nwpw"]["bo_scaling"][0].is_number_float())
         pscaling[0] = rtdbjson["nwpw"]["bo_scaling"][0];
   if (ptask==7)
      if (rtdbjson["nwpw"]["bo_scaling"][0].is_number_integer())
         pscaling[0] = rtdbjson["nwpw"]["bo_scaling"][0];
 
   ptolerances[0] = 1.0e-7;
   ptolerances[1] = 1.0e-7;
   ptolerances[2] = 1.0e-4;
   if (rtdbjson["nwpw"]["tolerances"][0].is_number_float())
      ptolerances[0] = rtdbjson["nwpw"]["tolerances"][0];
   if (rtdbjson["nwpw"]["tolerances"][1].is_number_float())
      ptolerances[1] = rtdbjson["nwpw"]["tolerances"][1];
   if (rtdbjson["nwpw"]["tolerances"][2].is_number_float())
      ptolerances[2] = rtdbjson["nwpw"]["tolerances"][2];
 
   if (rtdbjson["nwpw"]["deltae_check"].is_boolean())
      pdeltae_check = rtdbjson["nwpw"]["deltae_check"];
 
   if (ptask==5) 
   {
      ptolerances[0] = 1.0e-7;
      ptolerances[1] = 1.0e-7;
      ptolerances[2] = 1.0e-4;
   }
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["tolerances"][0].is_number_float())
      ptolerances[0] = rtdbjson["nwpw"]["steepest_descent"]["tolerances"][0];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["tolerances"][1].is_number_float())
      ptolerances[1] = rtdbjson["nwpw"]["steepest_descent"]["tolerances"][1];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["tolerances"][2].is_number_float())
      ptolerances[2] = rtdbjson["nwpw"]["steepest_descent"]["tolerances"][2];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["deltae_check"].is_boolean())
      pdeltae_check = rtdbjson["nwpw"]["steepest_descent"]["deltae_check"];
 
   // pecut=9000.0;
   // pwcut=9000.0;
   // pecut= 60.0; pwcut=30.0;
   pecut = 6000.0;
   pwcut = 3000.0;
   if (rtdbjson["nwpw"]["cutoff"][0].is_number_float())
      pwcut = rtdbjson["nwpw"]["cutoff"][0];
   if (rtdbjson["nwpw"]["cutoff"][1].is_number_float())
      pecut = rtdbjson["nwpw"]["cutoff"][1];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["cutoff"][0].is_number_float())
      pwcut = rtdbjson["nwpw"]["steepest_descent"]["cutoff"][0];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["cutoff"][1].is_number_float())
      pecut = rtdbjson["nwpw"]["steepest_descent"]["cutoff"][1];

   peprecondition = 20.0;
   psprecondition = 200.0;
   if (rtdbjson["nwpw"]["eprecondition"].is_number_float())
      peprecondition = rtdbjson["nwpw"]["eprecondition"];
   if (rtdbjson["nwpw"]["sprecondition"].is_number_float())
      psprecondition = rtdbjson["nwpw"]["sprecondition"];
 
   prcut = 0.0;
   if (rtdbjson["nwpw"]["ewald_rcut"].is_number_float())
      prcut = rtdbjson["nwpw"]["ewald_rcut"];
   pncut = 1;
   if (rtdbjson["nwpw"]["ewald_ncut"].is_number_integer())
      pncut = rtdbjson["nwpw"]["ewald_ncut"];
   pmultiplicity = 1;
   if (rtdbjson["nwpw"]["mult"].is_number_integer())
      pmultiplicity = rtdbjson["nwpw"]["mult"];
   pispin = 1;
   if (rtdbjson["nwpw"]["ispin"].is_number_integer())
      pispin = rtdbjson["nwpw"]["ispin"];
 
   xcstring = "";
   if (rtdbjson["nwpw"]["xc"].is_string())
      xcstring = rtdbjson["nwpw"]["xc"];

   // check for HFX stuff
   std::string myxc_name = xcstring;
   std::transform(myxc_name.begin(), myxc_name.end(), myxc_name.begin(), ::tolower);



 
   punita[0] = 20.0;
   punita[1] = 0.0;
   punita[2] = 0.0;
   punita[3] = 0.0;
   punita[4] = 20.0;
   punita[5] = 0.0;
   punita[6] = 0.0;
   punita[7] = 0.0;
   punita[8] = 20.0;
 
   std::string geomname = "geometry";
   if (rtdbjson["geometry"].is_string())
     geomname = rtdbjson["geometry"];
 
   if (rtdbjson["geometries"][geomname]["is_crystal"].is_boolean())
     pis_crystal = rtdbjson["geometries"][geomname]["is_crystal"];

   if (rtdbjson["geometries"][geomname]["unita"][0].is_number_float())
     punita[0] = rtdbjson["geometries"][geomname]["unita"][0];
   if (rtdbjson["geometries"][geomname]["unita"][1].is_number_float())
     punita[1] = rtdbjson["geometries"][geomname]["unita"][1];
   if (rtdbjson["geometries"][geomname]["unita"][2].is_number_float())
     punita[2] = rtdbjson["geometries"][geomname]["unita"][2];
 
   if (rtdbjson["geometries"][geomname]["unita"][3].is_number_float())
     punita[3] = rtdbjson["geometries"][geomname]["unita"][3];
   if (rtdbjson["geometries"][geomname]["unita"][4].is_number_float())
     punita[4] = rtdbjson["geometries"][geomname]["unita"][4];
   if (rtdbjson["geometries"][geomname]["unita"][5].is_number_float())
     punita[5] = rtdbjson["geometries"][geomname]["unita"][5];
 
   if (rtdbjson["geometries"][geomname]["unita"][6].is_number_float())
     punita[6] = rtdbjson["geometries"][geomname]["unita"][6];
   if (rtdbjson["geometries"][geomname]["unita"][7].is_number_float())
     punita[7] = rtdbjson["geometries"][geomname]["unita"][7];
   if (rtdbjson["geometries"][geomname]["unita"][8].is_number_float())
     punita[8] = rtdbjson["geometries"][geomname]["unita"][8];
 
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][0].is_number_float())
     punita[0] = rtdbjson["nwpw"]["simulation_cell"]["unita"][0];
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][1].is_number_float())
     punita[1] = rtdbjson["nwpw"]["simulation_cell"]["unita"][1];
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][2].is_number_float())
     punita[2] = rtdbjson["nwpw"]["simulation_cell"]["unita"][2];
 
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][3].is_number_float())
     punita[3] = rtdbjson["nwpw"]["simulation_cell"]["unita"][3];
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][4].is_number_float())
     punita[4] = rtdbjson["nwpw"]["simulation_cell"]["unita"][4];
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][5].is_number_float())
     punita[5] = rtdbjson["nwpw"]["simulation_cell"]["unita"][5];
 
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][6].is_number_float())
     punita[6] = rtdbjson["nwpw"]["simulation_cell"]["unita"][6];
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][7].is_number_float())
     punita[7] = rtdbjson["nwpw"]["simulation_cell"]["unita"][7];
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][8].is_number_float())
     punita[8] = rtdbjson["nwpw"]["simulation_cell"]["unita"][8];

 
   pngrid[0] = -1;
   pngrid[1] = -1;
   pngrid[2] = -1;
   if (rtdbjson["nwpw"]["simulation_cell"]["ngrid"][0].is_number_integer())
     pngrid[0] = rtdbjson["nwpw"]["simulation_cell"]["ngrid"][0];
   if (rtdbjson["nwpw"]["simulation_cell"]["ngrid"][1].is_number_integer())
     pngrid[1] = rtdbjson["nwpw"]["simulation_cell"]["ngrid"][1];
   if (rtdbjson["nwpw"]["simulation_cell"]["ngrid"][2].is_number_integer())
     pngrid[2] = rtdbjson["nwpw"]["simulation_cell"]["ngrid"][2];
   if ((pngrid[0] < 0) || (pngrid[1] < 0) || (pngrid[2] < 0)) {
     if (pecut > 5000.0) {
       pngrid[0] = 32;
       pngrid[1] = 32;
       pngrid[2] = 32;
     } else
       control_ngrid_default(punita, pecut, pmapping, pngrid);
   }
 
   pewald_grid[0] = pngrid[0];
   pewald_grid[1] = pngrid[1];
   pewald_grid[2] = pngrid[2];
   if (rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][0].is_number_integer())
     pewald_grid[0] = rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][0];
   if (rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][1].is_number_integer())
     pewald_grid[1] = rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][1];
   if (rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][2].is_number_integer())
     pewald_grid[2] = rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][2];
 
   if (rtdbjson["nwpw"]["simulation_cell"]["boundary_conditions"].is_string()) {
     std::string bc_str = mystring_lowercase(
         rtdbjson["nwpw"]["simulation_cell"]["boundary_conditions"]);
     if (mystring_contains(bc_str, "aperiodic"))
       version = 4;
     if (mystring_contains(bc_str, "free-space"))
       version = 4;
     if (mystring_contains(bc_str, "freespace"))
       version = 4;
     if (mystring_contains(bc_str, "free space"))
       version = 4;
   }
 
   pminimizer = 1;
   if (rtdbjson["nwpw"]["minimizer"].is_number_integer())
     pminimizer = rtdbjson["nwpw"]["minimizer"];
   if (rtdbjson["nwpw"]["lmbfgs_size"].is_number_integer())
     plmbfgs_size = rtdbjson["nwpw"]["lmbfgs_size"];
 
   // Efield data
   pefield_on = false;
   pefield_type = 2;
   if (!rtdbjson["nwpw"]["efield"].is_null()) {
     if (rtdbjson["nwpw"]["efield"]["on"].is_boolean())
       pefield_on = rtdbjson["nwpw"]["efield"]["on"];
     if (rtdbjson["nwpw"]["efield"]["type"].is_number_integer())
       pefield_type = rtdbjson["nwpw"]["efield"]["type"];
 
     if (!rtdbjson["nwpw"]["efield"]["center"].is_null()) {
       size_t nu = rtdbjson["nwpw"]["efield"]["center"].size();
       for (size_t i = 0; i < nu; ++i)
         pefield_center.push_back(rtdbjson["nwpw"]["efield"]["center"][i]);
     }
     if (!rtdbjson["nwpw"]["efield"]["vector"].is_null()) {
       size_t nu = rtdbjson["nwpw"]["efield"]["vector"].size();
       for (size_t i = 0; i < nu; ++i)
         pefield_vector.push_back(rtdbjson["nwpw"]["efield"]["vector"][i]);
     }
   }
 
   // APC data
   papc_on = false;
   papc_nga = 0;
   papc_Gc = 2.5;
   if (!rtdbjson["nwpw"]["apc"].is_null()) {
     if (rtdbjson["nwpw"]["apc"]["on"].is_boolean())
       papc_on = rtdbjson["nwpw"]["apc"]["on"];
     if (rtdbjson["nwpw"]["apc"]["Gc"].is_number_float())
       papc_Gc = rtdbjson["nwpw"]["apc"]["Gc"];
     if (!rtdbjson["nwpw"]["apc"]["gamma"].is_null()) {
       papc_nga = rtdbjson["nwpw"]["apc"]["gamma"].size();
       for (size_t i = 0; i < papc_nga; ++i)
         papc_gamma.push_back(rtdbjson["nwpw"]["apc"]["gamma"][i]);
     }
     if (!rtdbjson["nwpw"]["apc"]["u"].is_null()) {
       size_t nu = rtdbjson["nwpw"]["apc"]["u"].size();
       for (size_t i = 0; i < nu; ++i)
         papc_u.push_back(rtdbjson["nwpw"]["apc"]["u"][i]);
     }
     if (!rtdbjson["nwpw"]["apc"]["q"].is_null()) {
       size_t nq = rtdbjson["nwpw"]["apc"]["q"].size();
       for (size_t i = 0; i < nq; ++i)
         papc_q.push_back(rtdbjson["nwpw"]["apc"]["q"][i]);
     }
   }
 
   // Born data
   pborn_on = false;
   pborn_relax = false;
   if (!rtdbjson["nwpw"]["born"].is_null()) {
     auto bornjson = rtdbjson["nwpw"]["born"];
     if (bornjson["on"].is_boolean())
       pborn_on = bornjson["on"];
     if (bornjson["relax"].is_boolean())
       pborn_relax = bornjson["relax"];
     if (bornjson["dielec"].is_number_float())
       pborn_dielec = bornjson["dielec"];
     if (bornjson["rcut"].is_number_float())
       pborn_rcut = bornjson["rcut"];
 
     if (!bornjson["bradii"].is_null()) {
       size_t nu = bornjson["bradii"].size();
       for (size_t i = 0; i < nu; ++i)
         pborn_bradii.push_back(bornjson["bradii"][i]);
     }
     if (!bornjson["vradii"].is_null()) {
       size_t nu = bornjson["vradii"].size();
       for (size_t i = 0; i < nu; ++i)
         pborn_vradii.push_back(bornjson["vradii"][i]);
     }
   }
   if (pborn_on)
     papc_on = true;
 
   // Generalized Poisson data
   if (!rtdbjson["nwpw"]["generalized_poisson"].is_null()) 
   {
      auto gpoissonjson = rtdbjson["nwpw"]["generalized_poisson"];
      if (gpoissonjson["on"].is_boolean())            pgpoisson_on = gpoissonjson["on"];
      if (gpoissonjson["relax_dielec"].is_boolean())  pgpoisson_relax_dielec = gpoissonjson["relax_dielec"];
      if (gpoissonjson["cube_dielec"].is_boolean())   pgpoisson_cube_dielec = gpoissonjson["cube_dielec"];
      if (gpoissonjson["dielec"].is_number_float())   pgpoisson_dielec = gpoissonjson["dielec"];
      if (gpoissonjson["filter"].is_number_float())   pgpoisson_filter = gpoissonjson["filter"];
      if (gpoissonjson["rho0"].is_number_float())     pgpoisson_rho0 = gpoissonjson["rho0"];
      if (gpoissonjson["beta"].is_number_float())     pgpoisson_beta = gpoissonjson["beta"];
      if (gpoissonjson["rhomin"].is_number_float())   pgpoisson_rhomin = gpoissonjson["rhomin"];
      if (gpoissonjson["rhomax"].is_number_float())   pgpoisson_rhomax = gpoissonjson["rhomax"];
      if (gpoissonjson["rmin"].is_number_float())     pgpoisson_rmin = gpoissonjson["rmin"];
      if (gpoissonjson["rmax"].is_number_float())     pgpoisson_rmax = gpoissonjson["rmax"];
      if (gpoissonjson["rcut_ion"].is_number_float()) pgpoisson_rcut_ion = gpoissonjson["rcut_ion"];
      if (gpoissonjson["alpha"].is_number_float())    pgpoisson_alpha = gpoissonjson["alpha"];
      if (gpoissonjson["model"].is_number_integer())  pgpoisson_model = gpoissonjson["model"];
      if (gpoissonjson["maxit"].is_number_integer())  pgpoisson_maxit = gpoissonjson["maxit"];
   }

   // staged_gpu_fft
   if (!rtdbjson["nwpw"]["staged_gpu_fft"].is_null()) 
   {
      auto gstaged = rtdbjson["nwpw"]["staged_gpu_fft"];
      if (gstaged["on"].is_boolean()) pstaged_gpu_fft = gstaged["on"];
   }

   // fft_container_size
   if (rtdbjson["nwpw"]["fft_container_size"].is_number_integer()) 
      pfft_container_size = rtdbjson["nwpw"]["fft_container_size"];


 
   // Nose data
   pnose_on = false;
   pnose_restart = false;
   pnose_mchain = 0;
   pnose_mchain = 0;
   pnose_te = 0.0;
   pnose_tr = 0.0;
   pnose_pe = 0.0;
   pnose_pr = 0.0;
   pnose_ne_chain = 0.0;
   pnose_eke0 = 0.0;
   if (!rtdbjson["nwpw"]["car-parrinello"]["Nose-Hoover"].is_null()) 
   {
 
      auto nosejson = rtdbjson["nwpw"]["car-parrinello"]["Nose-Hoover"];
     
      if (nosejson["on"].is_boolean())      pnose_on = nosejson["on"];
      if (nosejson["restart"].is_boolean()) pnose_restart = nosejson["restart"];
     
      if (nosejson["Mchain"].is_number_integer()) pnose_mchain = nosejson["Mchain"];
      if (nosejson["Nchain"].is_number_integer()) pnose_nchain = nosejson["Nchain"];
     
      if (nosejson["Te"].is_number_float()) pnose_te = nosejson["Te"];
      if (nosejson["Tr"].is_number_float()) pnose_tr = nosejson["Tr"];
      if (nosejson["Pe"].is_number_float()) pnose_pe = nosejson["Pe"];
      if (nosejson["Pr"].is_number_float()) pnose_pr = nosejson["Pr"];
     
      if (nosejson["Ne_chain"].is_number_float()) pnose_ne_chain = nosejson["Ne_chain"];
      if (nosejson["eke0"].is_number_float())     pnose_eke0 = nosejson["eke0"];
     
      if (!nosejson["Xem"].is_null()) {
         size_t nu = nosejson["Xem"].size();
         for (size_t i = 0; i < nu; ++i)
            pnose_xem.push_back(nosejson["Xem"][i]);
      }
      if (!nosejson["Xe0"].is_null()) {
         size_t nu = nosejson["Xe0"].size();
         for (size_t i=0; i<nu; ++i)
            pnose_xe0.push_back(nosejson["Xe0"][i]);
      }
      if (!nosejson["Xe1"].is_null()) {
         size_t nu = nosejson["Xe1"].size();
         for (size_t i=0; i<nu; ++i)
            pnose_xe1.push_back(nosejson["Xe1"][i]);
      }
      if (!nosejson["Qe"].is_null()) {
         size_t nu = nosejson["Qe"].size();
         for (size_t i=0; i<nu; ++i)
            pnose_qe.push_back(nosejson["Qe"][i]);
      }
     
      if (!nosejson["Xrm"].is_null()) {
         size_t nu = nosejson["Xrm"].size();
         for (size_t i=0; i<nu; ++i)
            pnose_xrm.push_back(nosejson["Xrm"][i]);
      }
      if (!nosejson["Xr0"].is_null()) {
         size_t nu = nosejson["Xr0"].size();
         for (size_t i=0; i<nu; ++i)
            pnose_xr0.push_back(nosejson["Xr0"][i]);
      }
      if (!nosejson["Xr1"].is_null()) {
         size_t nu = nosejson["Xr1"].size();
         for (size_t i=0; i<nu; ++i)
            pnose_xr1.push_back(nosejson["Xr1"][i]);
      }
      if (!nosejson["Qr"].is_null()) {
         size_t nu = nosejson["Qr"].size();
         for (size_t i=0; i<nu; ++i)
            pnose_qr.push_back(nosejson["Qr"][i]);
      }
   }
 
   psa_on = false;
   psa_decay[0] = 1.0;
   psa_decay[1] = 1.0;
   if (rtdbjson["nwpw"]["car-parrinello"]["sa_decay"][0].is_number_float()) {
     psa_on = true;
     psa_decay[0] = rtdbjson["nwpw"]["car-parrinello"]["SA_decay"][0];
   }
   if (rtdbjson["nwpw"]["car-parrinello"]["sa_decay"][1].is_number_float()) {
     psa_on = true;
     psa_decay[1] = rtdbjson["nwpw"]["car-parrinello"]["SA_decay"][1];
   }
 
   if (!rtdbjson["driver"]["maxiter"].is_null()) {
     pdriver_maxiter = rtdbjson["driver"]["maxiter"];
   }
   if (!rtdbjson["driver"]["lmbfgs_size"].is_null()) {
     pdriver_lmbfgs_size = rtdbjson["driver"]["lmbfgs_size"];
   }
   if (!rtdbjson["driver"]["gmax"].is_null()) {
     pdriver_gmax = rtdbjson["driver"]["gmax"];
   }
   if (!rtdbjson["driver"]["grms"].is_null()) {
     pdriver_grms = rtdbjson["driver"]["grms"];
   }
   if (!rtdbjson["driver"]["xmax"].is_null()) {
     pdriver_xmax = rtdbjson["driver"]["xmax"];
   }
   if (!rtdbjson["driver"]["xrms"].is_null()) {
     pdriver_xrms = rtdbjson["driver"]["xrms"];
   }
   if (!rtdbjson["driver"]["trust"].is_null()) {
     pdriver_trust = rtdbjson["driver"]["trust"];
   }
}

void Control2::add_permanent_dir(char fname[]) 
{
   char stmp[256];
 
   if (strlen(ppermanent_dir) > 0) {
      strcpy(stmp, fname);
      strcpy(fname, ppermanent_dir);
      strcat(fname, "/");
      strcat(fname, stmp);
   }
}

void Control2::add_scratch_dir(char fname[]) 
{
   char stmp[256];
 
   if (strlen(pscratch_dir) > 0) {
      strcpy(stmp, fname);
      strcpy(fname, pscratch_dir);
      strcat(fname, "/");
      strcat(fname, stmp);
   }
}

int Control2::number_cubefiles() 
{
   json rtdbjson = json::parse(myrtdbstring);
   int num_cubefiles = 0;
   if (!rtdbjson["nwpw"]["dplot"].is_null())
      num_cubefiles = rtdbjson["nwpw"]["dplot"].size();
   return num_cubefiles;
}

int Control2::cubetype_cubefiles(const int i) 
{
   int cubetype = 0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["nwpw"]["dplot"].is_null()) 
   {
      int count = 0;
      std::string key = "";
      for (auto &el : rtdbjson["nwpw"]["dplot"].items()) {
         if (count == i)
            key = mystring_lowercase(el.key());
        ++count;
      }
      
      if (mystring_contains(key, "orbital-"))
        cubetype = std::stoi(mystring_split(key, "-")[1]);
      if (mystring_contains(key, "orbital2-"))             cubetype = -99;
      if (mystring_contains(key, "density-total"))         cubetype = -1;
      if (mystring_contains(key, "density-diff"))          cubetype = -2;
      if (mystring_contains(key, "density-alpha"))         cubetype = -3;
      if (mystring_contains(key, "density-up"))            cubetype = -3;
      if (mystring_contains(key, "density-beta"))          cubetype = -4;
      if (mystring_contains(key, "density-down"))          cubetype = -4;
      if (mystring_contains(key, "density-laplacian"))     cubetype = -5;
      if (mystring_contains(key, "density-potential"))     cubetype = -6;
      if (mystring_contains(key, "density-fullpotential")) cubetype = -7;
      if (mystring_contains(key, "elf"))                   cubetype = -8;
   }
   return cubetype;
}

std::string Control2::cubename_cubefiles(const int i) 
{
   std::string cubename = "";
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["nwpw"]["dplot"].is_null()) {
 
     int count = 0;
     for (auto &el : rtdbjson["nwpw"]["dplot"].items()) {
       if (count == i)
         cubename = el.value();
       ++count;
     }
   }
   return cubename;
}

std::string Control2::cubekey_cubefiles(const int i) 
{
   std::string cubekey = "";
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["nwpw"]["dplot"].is_null()) 
   {
      int count = 0;
      for (auto &el : rtdbjson["nwpw"]["dplot"].items()) 
      {
         if (count == i) cubekey = el.key();
         ++count;
      }
   }
   return cubekey;
}

double Control2::position_tolerance_cubefiles() 
{
   double tol = 0.0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["nwpw"]["dplot"]["position_tolerance"].is_null())
     tol = rtdbjson["nwpw"]["dplot"]["position_tolerance"];
   return tol;
}

int Control2::ncell_cubefiles(const int i) 
{
   int n = 0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["nwpw"]["dplot"]["ncell"].is_null())
     n = rtdbjson["nwpw"]["dplot"]["ncell"][i];
   
   return n;
}

double Control2::origin_cubefiles(const int i) 
{
   double origin = 0.0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["nwpw"]["dplot"]["ncell"].is_null())
     origin = rtdbjson["nwpw"]["dplot"]["origin"][i];
 
   return origin;
}


// bond constraints
int Control2::nhb_bond()
{
   int nhb = 0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["bond"].is_null())
      nhb = rtdbjson["constraints"]["bond"].size();
   return nhb;
}
int Control2::i0_bond(const int i)
{
   int i0=0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["bond"][i]["i0"].is_null())
      i0 = rtdbjson["constraints"]["bond"][i]["i0"];
   return i0;
}
int Control2::j0_bond(const int i)
{
   int j0=0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["bond"][i]["j0"].is_null())
      j0 = rtdbjson["constraints"]["bond"][i]["j0"];
   return j0;
}
double Control2::Kspring0_bond(const int i)
{
   double Kspring0=0.0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["bond"][i]["Kspring0"].is_null())
      Kspring0 = rtdbjson["constraints"]["bond"][i]["Kspring0"];
   return Kspring0;
}
double Control2::R0_bond(const int i)
{
   double R0=1.0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["bond"][i]["R0"].is_null())
      R0 = rtdbjson["constraints"]["bond"][i]["R0"];
   return R0;
}

// cbond constraints
int Control2::nhcb_cbond()
{
   int nhcb = 0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["cbond"].is_null())
      nhcb = rtdbjson["constraints"]["cbond"].size();
   return nhcb;
}
int Control2::i0_cbond(const int i)
{
   int i0=0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["cbond"][i]["i0"].is_null())
      i0 = rtdbjson["constraints"]["cbond"][i]["i0"];
   return i0;
}
int Control2::j0_cbond(const int i)
{
   int j0=0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["cbond"][i]["j0"].is_null())
      j0 = rtdbjson["constraints"]["cbond"][i]["j0"];
   return j0;
}
int Control2::k0_cbond(const int i)
{
   int k0=0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["cbond"][i]["k0"].is_null())
      k0 = rtdbjson["constraints"]["cbond"][i]["k0"];
   return k0;
}
int Control2::l0_cbond(const int i)
{
   int l0=0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["cbond"][i]["l0"].is_null())
      l0 = rtdbjson["constraints"]["cbond"][i]["l0"];
   return l0;
}
double Control2::Kspring0_cbond(const int i)
{
   double Kspring0=0.0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["cbond"][i]["Kspring0"].is_null())
      Kspring0 = rtdbjson["constraints"]["cbond"][i]["Kspring0"];
   return Kspring0;
}
double Control2::Rij0_cbond(const int i)
{
   double Rij0=1.0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["cbond"][i]["Rij0"].is_null())
      Rij0 = rtdbjson["constraints"]["cbond"][i]["Rij0"];
   return Rij0;
}
double Control2::Rkl0_cbond(const int i)
{
   double Rkl0=1.0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["cbond"][i]["Rkl0"].is_null())
      Rkl0 = rtdbjson["constraints"]["cbond"][i]["Rkl0"];
   return Rkl0;
}





// angle constraints
int Control2::nha_angle()
{
   int nha = 0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["angle"].is_null())
      nha = rtdbjson["constraints"]["angle"].size();
   return nha;
}
int Control2::i0_angle(const int i)
{
   int i0=0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["angle"][i]["i0"].is_null())
      i0 = rtdbjson["constraints"]["angle"][i]["i0"];
   return i0;
}
int Control2::j0_angle(const int i)
{
   int j0=0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["angle"][i]["j0"].is_null())
      j0 = rtdbjson["constraints"]["angle"][i]["j0"];
   return j0;
}
int Control2::k0_angle(const int i)
{
   int k0=0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["angle"][i]["k0"].is_null())
      k0 = rtdbjson["constraints"]["angle"][i]["k0"];
   return k0;
}
double Control2::Kspring0_angle(const int i)
{
   double Kspring0=0.0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["angle"][i]["Kspring0"].is_null())
      Kspring0 = rtdbjson["constraints"]["angle"][i]["Kspring0"];
   return Kspring0;
}
double Control2::Theta0_angle(const int i)
{
   double Theta0=1.0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["angle"][i]["Theta0"].is_null())
      Theta0 = rtdbjson["constraints"]["angle"][i]["Theta0"];
   return Theta0;
}



// bondings constraints
//### Adding reaction ###
//# reaction_type    = AB + C --> AC + B
//# reaction_indexes = 1 2 6
//# reaction_gamma   = -2.400
//##### current gamma=-2.400 ######
//constraints
//   clear
//   spring bondings 1.000000 -2.400 1.0 1 2 -1.0 1 6
//end
int Control2::nhc_bondings()
{
   int nhc = 0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["bondings"].is_null())
      nhc = rtdbjson["constraints"]["bondings"].size();
   return nhc;
}

std::vector<double> Control2::coef_bondings(const int i)
{
   json rtdbjson = json::parse(myrtdbstring);
   return rtdbjson["constraints"]["bondings"][i]["coef"];
}

std::vector<int> Control2::indx_bondings(const int i)
{
   json rtdbjson = json::parse(myrtdbstring);
   return  rtdbjson["constraints"]["bondings"][i]["indexes"];
}
double Control2::Kspring0_bondings(const int i)
{
   double Kspring0=1.0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["bondings"][i]["Kspring0"].is_null())
      Kspring0 = rtdbjson["constraints"]["bondings"][i]["Kspring0"];
   return Kspring0;
}
double Control2::gamma0_bondings(const int i)
{
   double gamma0=-1.0;
   json rtdbjson = json::parse(myrtdbstring);
   if (!rtdbjson["constraints"]["bondings"][i]["gamma0"].is_null())
      gamma0 = rtdbjson["constraints"]["bondings"][i]["gamma0"];
   return gamma0;
}


} // namespace pwdft

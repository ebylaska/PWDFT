#include        <iostream>
#include        <cmath>
#include        <cstdlib>
#include        <string>

using namespace std;


#include	"json.hpp"
#include	"rtdb.hpp"
#include	"Control2.hpp"
#include	"Parallel.hpp"

#include	"parsestring.hpp"

using json = nlohmann::json;




/***********************************
 *                                 *
 *         factor_count2           *
 *                                 *
 ***********************************/
static int factor_count2(int n, int m)
{  
   int f = 0;
   int nn = n;
   while ((nn%m)==0)
   {  
      nn /= m;
      ++f;
   }
   return f;
}
static int dum_ipow(int a, int n)
{  
   int npa = 1;
   for (auto i=0; i<n; ++i) npa *= a;
   return npa;
}

/***********************************
 *                                 *
 *         control_set_ngrid       *
 *                                 *
 ***********************************/
/* return n so that it is a multiple of 2,3,5,7 
*/
static int control_set_ngrid(double x, bool mult2)
{
   int nx = (int) floor(x + 0.5);
   if (mult2 && ((nx%2)!=0)) ++nx;

   int nf2 = factor_count2(nx,2);
   int nf3 = factor_count2(nx,3);
   int nf5 = factor_count2(nx,5);
   int nf7 = factor_count2(nx,7);
   int ntest = dum_ipow(2,nf2) * dum_ipow(3,nf3) * dum_ipow(5,nf5) * dum_ipow(7,nf7);
   while (nx != ntest)
   {
      ++nx;
      if (mult2) ++nx;
      nf2 = factor_count2(nx,2);
      nf3 = factor_count2(nx,3);
      nf5 = factor_count2(nx,5);
      nf7 = factor_count2(nx,7);
      ntest = dum_ipow(2,nf2) * dum_ipow(3,nf3) * dum_ipow(5,nf5) * dum_ipow(7,nf7);
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
static void control_ngrid_default(double *unita, double ecut, int mapping, int *ngrid)
{
   double gx,gy,gz,xh,yh,zh,unitg[9];
   double twopi = 8.0*atan(1.0);
   unitg[0] = unita[4]*unita[8] - unita[5]*unita[7];
   unitg[1] = unita[5]*unita[6] - unita[3]*unita[8];
   unitg[2] = unita[3]*unita[7] - unita[4]*unita[6];
   unitg[3] = unita[7]*unita[2] - unita[8]*unita[1];
   unitg[4] = unita[8]*unita[0] - unita[6]*unita[2];
   unitg[5] = unita[6]*unita[1] - unita[7]*unita[0];
   unitg[6] = unita[1]*unita[5] - unita[2]*unita[4];
   unitg[7] = unita[2]*unita[3] - unita[0]*unita[5];
   unitg[8] = unita[0]*unita[4] - unita[1]*unita[3];
   double volume = unita[0]*unitg[0] + unita[1]*unitg[1] + unita[2]*unitg[2];
   for (int i=0; i<9; ++i)
      unitg[i] *= twopi/volume;

   gx = unitg[0];
   gy = unitg[1];
   gz = unitg[2];
   xh = sqrt(2.00*ecut/(gx*gx + gy*gy + gz*gz))+0.5;

   gx = unitg[3];
   gy = unitg[4];
   gz = unitg[5];
   yh = sqrt(2.00*ecut/(gx*gx + gy*gy + gz*gz))+0.5;

   gx = unitg[6];
   gy = unitg[7];
   gz = unitg[8];
   zh = sqrt(2.00*ecut/(gx*gx + gy*gy + gz*gz))+0.5;

   ngrid[0] = control_set_ngrid(2.0*xh,true);
   ngrid[1] = control_set_ngrid(2.0*yh,true);
   ngrid[2] = control_set_ngrid(2.0*zh,true);
   if (mapping==1)
   {
      if (ngrid[1]>ngrid[2])
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

Control2::Control2(const int np0, const string rtdbstring)
{
   myrtdbstring  = rtdbstring;
   json rtdbjson = json::parse(rtdbstring);

   int  np = np0;
   bool is_cpmd;

   ptotal_ion_charge = 0;
   pne[0] = 0;
   pne[1] = 0;

   ptask = 0;
   if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),"energy"))           ptask = 1;
   if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),"gradient"))         ptask = 2;
   if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),"optimize"))         ptask = 3;
   if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),"freq"))             ptask = 4;
   if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),"steepest_descent")) ptask = 5;
   if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),"car-parrinello"))   ptask = 6;
   if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),"noit_"))            ptask *= -1;

   ptotal_ion_charge = -1.0;
   ptotal_charge = 0.0;
   if (rtdbjson["charge"].is_number_float()) ptotal_charge = rtdbjson["charge"];

   /* get parallel mappings */
   pmapping = 1;
   if (rtdbjson["nwpw"]["mapping"].is_number_integer()) pmapping = rtdbjson["nwpw"]["mapping"];

   /* set mapping1d */
   pmapping1d = 1;
   if (rtdbjson["nwpw"]["mapping1d"].is_number_integer()) pmapping1d = rtdbjson["nwpw"]["mapping1d"];

   /* qsize */
   pqsize = 4;
   if (rtdbjson["nwpw"]["pfft3_qsize"].is_number_integer()) pqsize = rtdbjson["nwpw"]["pfft3_qsize"];

   /* get np_dimensions */
   //if (!myrtdb.get("nwpw:np_dimensions",rtdb_int,3,np_dimensions))
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

      pnp_dimensions[2] = 1;
      if (pnp_dimensions[1]<1) pnp_dimensions[1] = 1;
      if (pnp_dimensions[0]<1) pnp_dimensions[0] = 1;

      /* reset np_dimensions[1] if it is not a  multiple of np2 */
      np = np/pnp_dimensions[2];
      while ( ((np%pnp_dimensions[1])!=0) && (pnp_dimensions[1]>1) )
            pnp_dimensions[1] = pnp_dimensions[1] - 1;
      pnp_dimensions[0] = np/pnp_dimensions[1];
    }

   pbalance = 1;
   if (rtdbjson["nwpw"]["nobalance"].is_boolean()) 
      if (rtdbjson["nwpw"]["nobalance"]) pbalance = 0;

   pgeometry_optimize = false;
   if (rtdbjson["nwpw"]["steepest_descent"]["geometry_optimize"].is_boolean()) 
      pgeometry_optimize = rtdbjson["nwpw"]["steepest_descent"]["geometry_optimize"];

   psp_library_dir = "";
   if (rtdbjson["psp_library_dir"].is_string()) 
      psp_library_dir = rtdbjson["psp_library_dir"];

   if (!rtdbjson["nwpw"]["pseudopotentials"].is_null()) 
      for (auto& el : rtdbjson["nwpw"]["pseudopotentials"].items()) {
         psp_libraries[el.key()] = el.value();
      }

   pprint_level = 2;
   if (rtdbjson["print"].is_string()) 
   {
      if (mystring_contains(mystring_lowercase(rtdbjson["print"]),"debug"))  pprint_level = 4;
      if (mystring_contains(mystring_lowercase(rtdbjson["print"]),"high"))   pprint_level = 3;
      if (mystring_contains(mystring_lowercase(rtdbjson["print"]),"medium")) pprint_level = 2;
      if (mystring_contains(mystring_lowercase(rtdbjson["print"]),"low"))    pprint_level = 1;
      if (mystring_contains(mystring_lowercase(rtdbjson["print"]),"off"))    pprint_level = 0;
      if (mystring_contains(mystring_lowercase(rtdbjson["print"]),"none"))   pprint_level = 0;
   }

   string permanent_dir_str = "";
   string scratch_dir_str   = "";
   if (rtdbjson["permanent_dir"].is_string()) 
      permanent_dir_str = rtdbjson["permanent_dir"];

   if (rtdbjson["scratch_dir"].is_string()) 
      scratch_dir_str = rtdbjson["scratch_dir"];


   string input_movecs  = "eric.movecs";  
   string output_movecs = "eric.movecs";  
   string input_v_movecs  = "eric.vmovecs";  
   string output_v_movecs = "eric.vmovecs";  
   if (rtdbjson["dbname"].is_string()) 
   {
       string dbname = rtdbjson["dbname"];
       input_movecs  = dbname + ".movecs";  
       output_movecs = dbname + ".movecs";  
       input_v_movecs  = dbname + ".vmovecs";  
       output_v_movecs = dbname + ".vmovecs";  
   }
   // read from nwpw block
   if (rtdbjson["nwpw"]["input_wavefunction_filename"].is_string())  input_movecs = rtdbjson["nwpw"]["input_wavefunction_filename"];
   if (rtdbjson["nwpw"]["output_wavefunction_filename"].is_string()) output_movecs = rtdbjson["nwpw"]["output_wavefunction_filename"];


   // from car-parrinello block
   if (ptask==6) if (rtdbjson["nwpw"]["car-parrinello"]["input_wavefunction_filename"].is_string()) input_movecs = rtdbjson["nwpw"]["car-parrinello"]["input_wavefunction_filename"];
   if (ptask==6) if (rtdbjson["nwpw"]["car-parrinello"]["output_wavefunction_filename"].is_string()) output_movecs = rtdbjson["nwpw"]["car-parrinello"]["output_wavefunction_filename"];
   if (ptask==6) if (rtdbjson["nwpw"]["car-parrinello"]["input_v_wavefunction_filename"].is_string()) input_v_movecs = rtdbjson["nwpw"]["car-parrinello"]["input_v_wavefunction_filename"];
   if (ptask==6) if (rtdbjson["nwpw"]["car-parrinello"]["output_v_wavefunction_filename"].is_string()) output_v_movecs = rtdbjson["nwpw"]["car-parrinello"]["output_v_wavefunction_filename"];

   // from steepest_descent block
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["input_wavefunction_filename"].is_string()) input_movecs = rtdbjson["nwpw"]["steepest_descent"]["input_wavefunction_filename"];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["output_wavefunction_filename"].is_string()) output_movecs = rtdbjson["nwpw"]["steepest_descent"]["output_wavefunction_filename"];
 
   if (permanent_dir_str.size()>0)
   { 
      if (input_movecs[0]!='/')  input_movecs  = permanent_dir_str + "/" + input_movecs;
      if (output_movecs[0]!='/') output_movecs = permanent_dir_str + "/" + output_movecs;
      if (input_v_movecs[0]!='/')  input_v_movecs  = permanent_dir_str + "/" + input_v_movecs;
      if (output_v_movecs[0]!='/') output_v_movecs = permanent_dir_str + "/" + output_v_movecs;
   }

   strcpy(pinput_movecs_filename,const_cast<char*>(input_movecs.data()));
   strcpy(poutput_movecs_filename,const_cast<char*>(output_movecs.data()));
   strcpy(pinput_v_movecs_filename,const_cast<char*>(input_v_movecs.data()));
   strcpy(poutput_v_movecs_filename,const_cast<char*>(output_v_movecs.data()));
   strcpy(ppermanent_dir,const_cast<char*>(permanent_dir_str.data()));
   strcpy(pscratch_dir,const_cast<char*>(scratch_dir_str.data()));
   //strcpy(ppsp_library_dir,const_cast<char*>(psp_library_dir_str.data()));


   pfake_mass = 400000.0;
   if (rtdbjson["nwpw"]["fake_mass"].is_number_float()) pfake_mass = rtdbjson["nwpw"]["fake_mass"];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["fake_mass"].is_number_float()) pfake_mass = rtdbjson["nwpw"]["steepest_descent"]["fake_mass"];
   if (ptask==6) if (rtdbjson["nwpw"]["car-parrinello"]["fake_mass"].is_number_float())   pfake_mass = rtdbjson["nwpw"]["car-parrinello"]["fake_mass"];


   ptime_step = 5.8;
   if (rtdbjson["nwpw"]["time_step"].is_number_float()) ptime_step = rtdbjson["nwpw"]["time_step"];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["time_step"].is_number_float()) ptime_step = rtdbjson["nwpw"]["steepest_descent"]["time_step"];
   if (ptask==6) if (rtdbjson["nwpw"]["car-parrinello"]["time_step"].is_number_float())   ptime_step = rtdbjson["nwpw"]["car-parrinello"]["time_step"];

   ploop[0]=10 ; ploop[1]=100; 
   if (rtdbjson["nwpw"]["loop"][0].is_number_integer()) ploop[0] = rtdbjson["nwpw"]["loop"][0];
   if (rtdbjson["nwpw"]["loop"][1].is_number_integer()) ploop[1] = rtdbjson["nwpw"]["loop"][1];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["loop"][0].is_number_integer()) ploop[0] = rtdbjson["nwpw"]["steepest_descent"]["loop"][0];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["loop"][1].is_number_integer()) ploop[1] = rtdbjson["nwpw"]["steepest_descent"]["loop"][1];
   if (ptask==6) if (rtdbjson["nwpw"]["car-parrinello"]["loop"][0].is_number_integer()) ploop[0] = rtdbjson["nwpw"]["car-parrinello"]["loop"][0];
   if (ptask==6) if (rtdbjson["nwpw"]["car-parrinello"]["loop"][1].is_number_integer()) ploop[1] = rtdbjson["nwpw"]["car-parrinello"]["loop"][1];


   ptolerances[0] = 1.0e-9; ptolerances[1] = 1.0e-9; ptolerances[2] = 1.0e-4; 
   if (rtdbjson["nwpw"]["tolerances"][0].is_number_float()) ptolerances[0] = rtdbjson["nwpw"]["tolerances"][0];
   if (rtdbjson["nwpw"]["tolerances"][1].is_number_float()) ptolerances[1] = rtdbjson["nwpw"]["tolerances"][1];
   if (rtdbjson["nwpw"]["tolerances"][2].is_number_float()) ptolerances[2] = rtdbjson["nwpw"]["tolerances"][2];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["tolerances"][0].is_number_float()) ptolerances[0] = rtdbjson["nwpw"]["steepest_descent"]["tolerances"][0];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["tolerances"][1].is_number_float()) ptolerances[1] = rtdbjson["nwpw"]["steepest_descent"]["tolerances"][1];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["tolerances"][2].is_number_float()) ptolerances[2] = rtdbjson["nwpw"]["steepest_descent"]["tolerances"][2];

   pecut=9000.0;
   pwcut=9000.0;
   if (rtdbjson["nwpw"]["cutoff"][0].is_number_float()) pwcut = rtdbjson["nwpw"]["cutoff"][0];
   if (rtdbjson["nwpw"]["cutoff"][1].is_number_float()) pecut = rtdbjson["nwpw"]["cutoff"][1];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["cutoff"][0].is_number_float()) pwcut = rtdbjson["nwpw"]["steepest_descent"]["cutoff"][0];
   if (ptask==5) if (rtdbjson["nwpw"]["steepest_descent"]["cutoff"][1].is_number_float()) pecut = rtdbjson["nwpw"]["steepest_descent"]["cutoff"][1];

   prcut=0.0;
   if (rtdbjson["nwpw"]["ewald_rcut"].is_number_float()) prcut = rtdbjson["nwpw"]["ewald_rcut"];
   pncut=1;
   if (rtdbjson["nwpw"]["ewald_ncut"].is_number_integer()) pncut = rtdbjson["nwpw"]["ewald_ncut"];
   pmultiplicity=1;
   if (rtdbjson["nwpw"]["mult"].is_number_integer()) pmultiplicity = rtdbjson["nwpw"]["mult"];
   pispin=1;
   if (rtdbjson["nwpw"]["ispin"].is_number_integer()) pispin = rtdbjson["nwpw"]["ispin"];

   xcstring = "";
   if (rtdbjson["nwpw"]["xc"].is_string()) xcstring = rtdbjson["nwpw"]["xc"];


   punita[0] = 20.0; punita[1] =  0.0; punita[2] =  0.0;
   punita[3] =  0.0; punita[4] = 20.0; punita[5] =  0.0;
   punita[6] =  0.0; punita[7] =  0.0; punita[8] = 20.0;
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][0].is_number_float()) punita[0] = rtdbjson["nwpw"]["simulation_cell"]["unita"][0];
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][1].is_number_float()) punita[1] = rtdbjson["nwpw"]["simulation_cell"]["unita"][1];
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][2].is_number_float()) punita[2] = rtdbjson["nwpw"]["simulation_cell"]["unita"][2];

   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][3].is_number_float()) punita[3] = rtdbjson["nwpw"]["simulation_cell"]["unita"][3];
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][4].is_number_float()) punita[4] = rtdbjson["nwpw"]["simulation_cell"]["unita"][4];
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][5].is_number_float()) punita[5] = rtdbjson["nwpw"]["simulation_cell"]["unita"][5];

   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][6].is_number_float()) punita[6] = rtdbjson["nwpw"]["simulation_cell"]["unita"][6];
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][7].is_number_float()) punita[7] = rtdbjson["nwpw"]["simulation_cell"]["unita"][7];
   if (rtdbjson["nwpw"]["simulation_cell"]["unita"][8].is_number_float()) punita[8] = rtdbjson["nwpw"]["simulation_cell"]["unita"][8];


   pngrid[0] = -1; pngrid[1] = -1; pngrid[2] = -1;
   if (rtdbjson["nwpw"]["simulation_cell"]["ngrid"][0].is_number_integer()) pngrid[0] = rtdbjson["nwpw"]["simulation_cell"]["ngrid"][0];
   if (rtdbjson["nwpw"]["simulation_cell"]["ngrid"][1].is_number_integer()) pngrid[1] = rtdbjson["nwpw"]["simulation_cell"]["ngrid"][1];
   if (rtdbjson["nwpw"]["simulation_cell"]["ngrid"][2].is_number_integer()) pngrid[2] = rtdbjson["nwpw"]["simulation_cell"]["ngrid"][2];
   if ((pngrid[0]<0) || (pngrid[1]<0) || (pngrid[2]<0))
   {
      if (pecut >5000.0)
      {
         pngrid[0] = 32; pngrid[1] = 32; pngrid[2] = 32;
      }
      else
         control_ngrid_default(punita, pecut, pmapping, pngrid);
   }


   pewald_grid[0] = pngrid[0]; pewald_grid[1] = pngrid[1]; pewald_grid[2] = pngrid[2];
   if (rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][0].is_number_integer()) pewald_grid[0] = rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][0];
   if (rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][1].is_number_integer()) pewald_grid[1] = rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][1];
   if (rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][2].is_number_integer()) pewald_grid[2] = rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][2];


   pminimizer = 1;
   if (rtdbjson["nwpw"]["minimizer"].is_number_integer()) pminimizer = rtdbjson["nwpw"]["minimizer"];

}

void Control2::add_permanent_dir(char fname[])
{
   char stmp[256]; 

   if (strlen(ppermanent_dir)>0)
   {
      strcpy(stmp,fname);
      strcpy(fname,ppermanent_dir);
      strcat(fname,"/");
      strcat(fname,stmp);
   }
}

void Control2::add_scratch_dir(char fname[])
{
   char stmp[256]; 

   if (strlen(pscratch_dir)>0)
   {
      strcpy(stmp,fname);
      strcpy(fname,pscratch_dir);
      strcat(fname,"/");
      strcat(fname,stmp);
   }
}



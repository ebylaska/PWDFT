#include        <iostream>
#include        <cmath>
#include        <cstdlib>
#include        <string>

using namespace std;


#include	"json.hpp"
#include	"rtdb.hpp"
#include	"Control2.hpp"

#include	"parsestring.hpp"

using json = nlohmann::json;


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

   ptask = 0;
   if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),"energy"))           ptask = 1;
   if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),"gradient"))         ptask = 2;
   if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),"optimize"))         ptask = 3;
   if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),"freq"))             ptask = 4;
   if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),"steepest_descent")) ptask = 5;
   if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]),"car-parrinello"))   ptask = 6;


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


   pngrid[0] = 32; pngrid[1] = 32; pngrid[2] = 32;
   if (rtdbjson["nwpw"]["simulation_cell"]["ngrid"][0].is_number_integer()) pngrid[0] = rtdbjson["nwpw"]["simulation_cell"]["ngrid"][0];
   if (rtdbjson["nwpw"]["simulation_cell"]["ngrid"][1].is_number_integer()) pngrid[1] = rtdbjson["nwpw"]["simulation_cell"]["ngrid"][1];
   if (rtdbjson["nwpw"]["simulation_cell"]["ngrid"][2].is_number_integer()) pngrid[2] = rtdbjson["nwpw"]["simulation_cell"]["ngrid"][2];


   pewald_grid[0] = pngrid[0]; pewald_grid[1] = pngrid[1]; pewald_grid[2] = pngrid[2];
   if (rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][0].is_number_integer()) pewald_grid[0] = rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][0];
   if (rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][1].is_number_integer()) pewald_grid[1] = rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][1];
   if (rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][2].is_number_integer()) pewald_grid[2] = rtdbjson["nwpw"]["simulation_cell"]["ewald_ngrid"][2];

}

void Control2::add_permanent_dir(char fname[])
{
   char stmp[80]; 

   if (strlen(ppermanent_dir)>0)
   {
      strcpy(stmp,fname);
      strcpy(fname,ppermanent_dir);
      strcat(fname,"/");
      strcat(fname,stmp);
   }
}




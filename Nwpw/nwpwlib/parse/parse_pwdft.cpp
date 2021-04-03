#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

#include "json.hpp"
#include "parsestring.hpp"


using json = nlohmann::json;

using namespace std;

json periodic_table_mass = json::parse("{ \"H\"  : 1.008, \"He\" : 4.0026, \"Li\" : 7.016, \"Be\" : 9.01218, \"B\"  : 11.00931, \"C\"  : 12.0, \"N\"  : 14.00307, \"O\"  : 15.99491, \"F\"  : 18.9984, \"Ne\" : 19.99244, \"Na\" : 22.9898, \"Mg\" : 23.98504, \"Al\" : 26.98154, \"Si\" : 27.97693, \"P\"  : 30.97376, \"S\"  : 31.97207, \"Cl\" : 34.96885, \"Ar\" : 39.9624, \"K\"  : 38.96371, \"Ca\" : 39.96259, \"Sc\" : 44.95592, \"Ti\" : 45.948, \"V\"  : 50.9440, \"Cr\" : 51.9405, \"Mn\" : 54.9381, \"Fe\" : 55.9349, \"Co\" : 58.9332, \"Ni\" : 57.9353, \"Cu\" : 62.9298, \"Zn\" : 63.9291, \"Ga\" : 68.9257, \"Ge\" : 73.9219, \"As\" : 74.9216, \"Se\" : 78.9183, \"Br\" : 79.9165, \"Kr\" : 83.912, \"Rb\" : 84.9117, \"Sr\" : 87.9056, \"Y\"  : 88.9054, \"Zr\" : 89.9043, \"Nb\" : 92.9060, \"Mo\" : 97.9055, \"Tc\" : 97.9072, \"Ru\" : 101.9037, \"Rh\" : 102.9048, \"Pd\" : 105.9032, \"Ag\" : 106.90509, \"Cd\" : 113.9036, \"In\" : 114.9041, \"Sn\" : 117.9018, \"Sb\" : 120.9038, \"Te\" : 129.9067, \"I\"  : 126.9004, \"Xe\" : 131.9042, \"Cs\" : 132.9051, \"Ba\" : 137.9050, \"La\" : 138.9061, \"Ce\" : 139.9053, \"Pr\" : 140.9074, \"Nd\" : 143.9099, \"Pm\" : 144.9128, \"Sm\" : 151.9195, \"Eu\" : 152.920, \"Gd\" : 157.9241, \"Tb\" : 159.9250, \"Dy\" : 163.9288, \"Ho\" : 164.9303, \"Er\" : 165.930, \"Tm\" : 168.9344, \"Yb\" : 173.9390, \"Lu\" : 174.9409, \"Hf\" : 179.9468, \"Ta\" : 180.948, \"W\"  : 183.9510, \"Re\" : 186.9560, \"Os\" : 189.9586, \"Ir\" : 192.9633, \"Pt\" : 194.9648, \"Au\" : 196.9666, \"Hg\" : 201.9706, \"Tl\" : 204.9745, \"Pb\" : 207.9766, \"Bi\" : 208.9804, \"Po\" : 209.9829, \"At\" : 210.9875, \"Rn\" : 222.0175, \"Fr\" : 223.0198, \"Ra\" : 226.0254, \"Ac\" : 227.0278, \"Th\" : 232.0382, \"Pa\" : 231.0359, \"U\"  : 238.0508, \"Np\" : 237.0482, \"Pu\" : 244.0642, \"Am\" : 243.0614, \"Cm\" : 247.0704, \"Bk\" : 247.0703, \"Cf\" : 251.0796, \"Es\" : 252.0829, \"Fm\" : 257.0950, \"Md\" : 258.0986, \"No\" : 259.1009, \"Lr\" : 262.1100, \"Rf\" : 261.1087, \"Ha\" : 262.1138, \"Sg\" : 266.1219, \"Bh\" : 262.1229, \"Hs\" : 267.1318, \"Mt\" : 268.1388 }");


json periodic_table_Z = json::parse("{ \"H\"  : 1, \"He\" : 2, \"Li\" : 3, \"Be\" : 4, \"B\"  : 5, \"C\"  : 6, \"N\"  : 7, \"O\"  : 8, \"F\"  : 9, \"Ne\" : 10, \"Na\" : 11, \"Mg\" : 12, \"Al\" : 13, \"Si\" : 14, \"P\"  : 15, \"S\"  : 16, \"Cl\" : 17, \"Ar\" : 18, \"K\"  : 19, \"Ca\" : 20, \"Sc\" : 21, \"Ti\" : 22, \"V\"  : 23, \"Cr\" : 24, \"Mn\" : 25, \"Fe\" : 26, \"Co\" : 27, \"Ni\" : 28, \"Cu\" : 29, \"Zn\" : 30, \"Ga\" : 31, \"Ge\" : 32, \"As\" : 33, \"Se\" : 34, \"Br\" : 35, \"Kr\" : 36, \"Rb\" : 37, \"Sr\" : 38, \"Y\"  : 39, \"Zr\" : 40, \"Nb\" : 41, \"Mo\" : 42, \"Tc\" : 43, \"Ru\" : 44, \"Rh\" : 45, \"Pd\" : 46, \"Ag\" : 47, \"Cd\" : 48, \"In\" : 49, \"Sn\" : 50, \"Sb\" : 51, \"Te\" : 52, \"I\"  : 53, \"Xe\" : 54, \"Cs\" : 55, \"Ba\" : 56, \"La\" : 57, \"Ce\" : 58, \"Pr\" : 59, \"Nd\" : 60, \"Pm\" : 61, \"Sm\" : 62, \"Eu\" : 63, \"Gd\" : 64, \"Tb\" : 65, \"Dy\" : 66, \"Ho\" : 67, \"Er\" : 68, \"Tm\" : 69, \"Yb\" : 70, \"Lu\" : 71, \"Hf\" : 72, \"Ta\" : 73, \"W\"  : 74, \"Re\" : 75, \"Os\" : 76, \"Ir\" : 77, \"Pt\" : 78, \"Au\" : 79, \"Hg\" : 80, \"Tl\" : 81, \"Pb\" : 82, \"Bi\" : 83, \"Po\" : 84, \"At\" : 85, \"Rn\" : 86, \"Fr\" : 87, \"Ra\" : 88, \"Ac\" : 89, \"Th\" : 90, \"Pa\" : 91, \"U\"  : 92, \"Np\" : 93, \"Pu\" : 94, \"Am\" : 95, \"Cm\" : 96, \"Bk\" : 97, \"Cf\" : 98, \"Es\" : 99, \"Fm\" : 100, \"Md\" : 101, \"No\" : 102, \"Lr\" : 103, \"Rf\" : 104, \"Ha\" : 105, \"Sg\" : 106, \"Bh\" : 107, \"Hs\" : 108, \"Mt\" : 109 }");

/**************************************************
 *                                                *
 *                parse_geometry                  *
 *                                                *
 **************************************************/

static json parse_geometry(json geom, int *curptr, vector<string> lines)
{
   json geomjson;

   int cur = *curptr;
   int center = 1;
   int autoz = 0;
   int autosym = 0;
   //double angs_to_au = 1.0/0.52917715;
   double angs_to_au = 1.88972598858;
   double conv = angs_to_au;
   vector<string> ss;
   string geometry = "geometry";

   ss = mystring_split0(lines[cur]);
   if (ss.size()>1) 
   {
      int nogeomname =  mystring_contains(mystring_lowercase(ss[1]),"au")
                     || mystring_contains(mystring_lowercase(ss[1]),"a.u.")
                     || mystring_contains(mystring_lowercase(ss[1]),"an")
                     || mystring_contains(mystring_lowercase(ss[1]),"nm")
                     || mystring_contains(mystring_lowercase(ss[1]),"na")
                     || mystring_contains(mystring_lowercase(ss[1]),"pm")
                     || mystring_contains(mystring_lowercase(ss[1]),"pi")
                     || mystring_contains(mystring_lowercase(ss[1]),"center")
                     || mystring_contains(mystring_lowercase(ss[1]),"autoz")
                     || mystring_contains(mystring_lowercase(ss[1]),"autosym");
      if (!nogeomname) geometry = ss[1];
   }


   if (mystring_contains(mystring_lowercase(lines[cur])," au"))   conv = 1.0;
   if (mystring_contains(mystring_lowercase(lines[cur])," a.u.")) conv = 1.0;
   if (mystring_contains(mystring_lowercase(lines[cur])," bo"))   conv = 1.0;
   if (mystring_contains(mystring_lowercase(lines[cur])," an"))   conv = angs_to_au;
   if (mystring_contains(mystring_lowercase(lines[cur])," nm"))   conv = 10.0*angs_to_au;
   if (mystring_contains(mystring_lowercase(lines[cur])," na"))   conv = 10.0*angs_to_au;
   if (mystring_contains(mystring_lowercase(lines[cur])," pm"))   conv = 0.01*angs_to_au;
   if (mystring_contains(mystring_lowercase(lines[cur])," pi"))   conv = 0.01*angs_to_au;

   if (mystring_contains(mystring_lowercase(lines[cur]),"nocenter"))    center = 0;
   else if (mystring_contains(mystring_lowercase(lines[cur]),"center")) center = 1;
   if (mystring_contains(mystring_lowercase(lines[cur]),"noautoz"))     autoz = 0;
   else if (mystring_contains(mystring_lowercase(lines[cur]),"autoz"))  autoz = 1;
   if (mystring_contains(mystring_lowercase(lines[cur]),"noautosym"))    autosym = 0;
   else if (mystring_contains(mystring_lowercase(lines[cur]),"autosym")) autosym = 1;

   geomjson["conv"]    = conv;
   geomjson["center"]  = center;
   geomjson["autoz"]   = autoz;
   geomjson["autosym"] = autosym;

   int endcount = 1;
   vector<string> symbols;
   vector<double> coords;
   vector<double> velocities;
   vector<double> masses;
   vector<double> charges;
   ++cur;
   int nion = 0;
   double mm;
   string line;
   while (endcount>0)
   {
      line = lines[cur];
      
      ss = mystring_split0(lines[cur]);
      symbols.push_back(mystring_capitalize(ss[0]));
      coords.push_back(std::stod(ss[1])*conv);
      coords.push_back(std::stod(ss[2])*conv);
      coords.push_back(std::stod(ss[3])*conv);

      double vx = 0.0;
      double vy = 0.0;
      double vz = 0.0;
      if (ss.size()>6)
      {
         bool ffail = false;
         try {vx = std::stod(ss[4])*conv;} catch(std::exception& ia) {vx = 0.0; ffail=true;}
         try {vy = std::stod(ss[5])*conv;} catch(std::exception& ia) {vy = 0.0; ffail=true;}
         try {vz = std::stod(ss[6])*conv;} catch(std::exception& ia) {vz = 0.0; ffail=true;}
         if (ffail) { vx=0.0; vy=0.0; vz=0.0;}
      }
      velocities.push_back(vx);
      velocities.push_back(vy);
      velocities.push_back(vz);

      mm = periodic_table_mass[mystring_capitalize(ss[0])];
      if (mystring_contains(mystring_lowercase(line),"mass"))
         mm = std::stod(mystring_split0(mystring_split(line,"mass")[1])[0]);
      masses.push_back(mm);

      mm = (double) periodic_table_Z[mystring_capitalize(ss[0])];
      if (mystring_contains(mystring_lowercase(line),"charge"))
         mm = std::stod(mystring_split0(mystring_split(line,"charge")[1])[0]);
      charges.push_back(mm);

      ++nion;
      ++cur;
      if (mystring_contains(mystring_lowercase(lines[cur]),"end"))
         --endcount;
   }

   if (center)
   {
      int ii;
      double xcm=0.0;
      double ycm=0.0;
      double zcm=0.0;
      for (ii=0; ii<nion; ++ii)
      {
         xcm += coords[3*ii];
         ycm += coords[3*ii+1];
         zcm += coords[3*ii+2];
      }
      xcm /= ((double) nion);
      ycm /= ((double) nion);
      zcm /= ((double) nion);
      for (int ii=0; ii<nion; ++ii)
      {
         coords[3*ii]   -= xcm;
         coords[3*ii+1] -= ycm;
         coords[3*ii+2] -= zcm;
      }
   }

   geomjson["symbols"]    = symbols;
   geomjson["coords"]     = coords;
   geomjson["velocities"] = velocities;
   geomjson["nion"]       = nion;
   geomjson["masses"]     = masses;
   geomjson["charges"]    = charges;

   *curptr = cur;

   geom[geometry] = geomjson;

   return geom;
}
/**************************************************
 *                                                *
 *              parse_pseudopotentials             *
 *                                                *
 **************************************************/
static json parse_pseudopotentials(json pseudopotentials, int *curptr, vector<string> lines)
{
   int cur = *curptr;
   int endcount = 1;
   ++cur;
   string line;
   vector<string> ss;

   while (endcount>0)
   {
      line = lines[cur];
      if (mystring_contains(lines[cur],"library"))
      {
         ss = mystring_split0(lines[cur]);
         if (ss.size()>2)
            if (mystring_contains(ss[1],"library"))
                pseudopotentials[ss[0]] = ss[2];
      }
      ++cur;
      if (mystring_contains(lines[cur],"end"))
         --endcount;
   }
   *curptr = cur;
   return pseudopotentials;
}

/**************************************************
 *                                                *
 *                parse_simulation_cell           *
 *                                                *
 **************************************************/

static json parse_simulation_cell(json celljson, int *curptr, vector<string> lines)
{
   int cur = *curptr;
   int endcount = 1;
   ++cur;
   string line;
   vector<string> ss;

   while (endcount>0)
   {
      line = mystring_lowercase(lines[cur]);
      if (mystring_contains(line,"lattice_vectors"))
      {
         vector<double> unita = {0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0};
         ++cur;
         line = mystring_lowercase(lines[cur]);
         ss = mystring_split0(line);
         if (ss.size()>2) 
         {
             unita[0] = std::stod(ss[0]);
             unita[1] = std::stod(ss[1]);
             unita[2] = std::stod(ss[2]);
         }
         ++cur;
         line = mystring_lowercase(lines[cur]);
         ss = mystring_split0(line);
         if (ss.size()>2) 
         {
             unita[3] = std::stod(ss[0]);
             unita[4] = std::stod(ss[1]);
             unita[5] = std::stod(ss[2]);
         }
         ++cur;
         line = mystring_lowercase(lines[cur]);
         ss = mystring_split0(line);
         if (ss.size()>2) 
         {
             unita[6] = std::stod(ss[0]);
             unita[7] = std::stod(ss[1]);
             unita[8] = std::stod(ss[2]);
         }
         celljson["unita"] = unita;
      }
      else if (mystring_contains(line,"lattice"))
      {
         vector<double> unita = {0.0,0.0,0.0, 0.0,0.0,0.0, 0.0,0.0,0.0};
         celljson["unita"] = unita;
      }
      else if (mystring_contains(line,"ngrid"))
      {
         ss = mystring_split0(line);
         if (ss.size()>3) 
         {
            vector<int> ngrid = {std::stoi(ss[1]),std::stoi(ss[2]),std::stoi(ss[3])};
            celljson["ngrid"] = ngrid;
         }
      }
      else if (mystring_contains(line,"sc"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) 
         {
            double a = std::stod(ss[1]);
            vector<double> unita = {a,0.0,0.0, 0.0,a,0.0, 0.0,0.0,a};
            celljson["unita"] = unita;
         }
      }
      else if (mystring_contains(line,"fcc"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) 
         {
            double a = 0.5*std::stod(ss[1]);
            vector<double> unita = {a,a,0.0, a,0.0,a, 0.0,a,a};
            celljson["unita"] = unita;
         }
      }
      else if (mystring_contains(line,"bcc"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) 
         {
            double a = 0.5*std::stod(ss[1]);
            vector<double> unita = {-a,a,a, a,-a,a, a,a,-a};
            celljson["unita"] = unita;
         }
      }
      else if (mystring_contains(line,"box_delta"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) 
         {
            double a = std::stod(ss[1]);
            celljson["box_delta"] = a;
         }
      }
      else if (mystring_contains(line,"box_orient"))
      {
         celljson["box_orient"] = true;
      }
      else if (mystring_contains(line,"box_different_lengths"))
      {
         celljson["box_different_lengths"] = true;
      }


      ++cur;
      if (mystring_contains(lines[cur],"end"))
         --endcount;
   }

   *curptr = cur;

   return celljson;
}

/**************************************************
 *                                                *
 *                parse_steepest_descent          *
 *                                                *
 **************************************************/

static json parse_steepest_descent(json sdjson, int *curptr, vector<string> lines)
{
   int cur = *curptr;
   int endcount = 1;
   ++cur;
   string line;
   vector<string> ss;

   while (endcount>0)
   {
      line = mystring_lowercase(lines[cur]);

      if (mystring_contains(line,"loop"))
      {
         vector<int> loop = {1, 1};
         ss = mystring_split0(line);
         if (ss.size()>1) loop[0] = std::stoi(ss[1]);
         if (ss.size()>2) loop[1] = std::stoi(ss[2]);
         sdjson["loop"] = loop;
      }
      else if (mystring_contains(line,"xc"))
      {
         sdjson["xc"] = mystring_trim(mystring_split(line,"xc")[1]);
      }
      else if (mystring_contains(line,"geometry_optimize"))
      {
         sdjson["geometry_optimize"] = true;
      }
      else if (mystring_contains(line,"input_wavefunction_filename"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) sdjson["input_wavefunction_filename"] = ss[1];
      }
      else if (mystring_contains(line,"output_wavefunction_filename"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) sdjson["output_wavefunction_filename"] = ss[1];
      }
      else if (mystring_contains(line,"time_step"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) sdjson["time_step"] = std::stod(ss[1]);
      }
      else if (mystring_contains(line,"fake_mass"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) sdjson["fake_mass"] = std::stod(ss[1]);
      }
      else if (mystring_contains(line,"cutoff"))
      {
         ss = mystring_split0(line);
         if (ss.size()==2) sdjson["cutoff"] = {std::stod(ss[1]),2*std::stod(ss[1])};
         if (ss.size()>2)  sdjson["cutoff"] = {std::stod(ss[1]),std::stod(ss[2])};
      }
      else if (mystring_contains(line,"tolerances"))
      {
         ss = mystring_split0(line);
         if (ss.size()==2) sdjson["tolerances"] = {std::stod(ss[1]),std::stod(ss[1]),1.0e-4};
         if (ss.size()==3) sdjson["tolerances"] = {std::stod(ss[1]),std::stod(ss[2]),1.0e-4};
         if (ss.size()>3)  sdjson["tolerances"] = {std::stod(ss[1]),std::stod(ss[2]),std::stod(ss[3])};
      }

      ++cur;
      if (mystring_contains(lines[cur],"end"))
         --endcount;
   }

   *curptr = cur;

   return sdjson;
}

/**************************************************
 *                                                *
 *                parse_car_parrinello            *
 *                                                *
 **************************************************/

static json parse_car_parrinello(json cpmdjson, int *curptr, vector<string> lines)
{
   int cur = *curptr;
   int endcount = 1;
   ++cur;
   string line;
   vector<string> ss;

   while (endcount>0)
   {
      line = mystring_lowercase(lines[cur]);
      if (mystring_contains(line,"loop"))
      {
         vector<int> loop = {1, 1};
         ss = mystring_split0(line);
         if (ss.size()>1) loop[0] = std::stoi(ss[1]);
         if (ss.size()>2) loop[1] = std::stoi(ss[2]);
         cpmdjson["loop"] = loop;
      }
      else if (mystring_contains(line,"xc"))
      {
         cpmdjson["xc"] = mystring_trim(mystring_split(line,"xc")[1]);
      }
      else if (mystring_contains(line,"input_wavefunction_filename"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) cpmdjson["input_wavefunction_filename"] = ss[1];
      }
      else if (mystring_contains(line,"output_wavefunction_filename"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) cpmdjson["output_wavefunction_filename"] = ss[1];
      }
      else if (mystring_contains(line,"input_v_wavefunction_filename"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) cpmdjson["input_v_wavefunction_filename"] = ss[1];
      }
      else if (mystring_contains(line,"output_v_wavefunction_filename"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) cpmdjson["output_v_wavefunction_filename"] = ss[1];
      }
      else if (mystring_contains(line,"time_step"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) cpmdjson["time_step"] = std::stod(ss[1]);
      }
      else if (mystring_contains(line,"fake_mass"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) cpmdjson["fake_mass"] = std::stod(ss[1]);
      }
      else if (mystring_contains(line,"cutoff"))
      {
         ss = mystring_split0(line);
         if (ss.size()==2) cpmdjson["cutoff"] = {std::stod(ss[1]),2*std::stod(ss[1])};
         if (ss.size()>2)  cpmdjson["cutoff"] = {std::stod(ss[1]),std::stod(ss[2])};
      }
      else if (mystring_contains(line,"scaling"))
      {
         ss = mystring_split0(line);
         if (ss.size()==2) cpmdjson["scaling"] = {std::stod(ss[1]),std::stod(ss[1])};
         if (ss.size()>2)  cpmdjson["scaling"] = {std::stod(ss[1]),std::stod(ss[2])};
      }
      else if (mystring_contains(line,"sa_decay"))
      {
         ss = mystring_split0(line);
         if (ss.size()==2) cpmdjson["sa_decay"] = {std::stod(ss[1]),std::stod(ss[1])};
         if (ss.size()>2)  cpmdjson["sa_decay"] = {std::stod(ss[1]),std::stod(ss[2])};
      }
      else if (mystring_contains(line,"xyz_filename"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) cpmdjson["xyz_filename"] = ss[1];
      }
      else if (mystring_contains(line,"emotion_filename"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) cpmdjson["emotion_filename"] = ss[1];
      }
      else if (mystring_contains(line,"ion_motion_filename"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) cpmdjson["ion_motion_filename"] = ss[1];
      }
      else if (mystring_contains(line,"eigmotion_filename"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) cpmdjson["eigmotion_filename"] = ss[1];
      }
      else if (mystring_contains(line,"omotion_filename"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) cpmdjson["omotion_filename"] = ss[1];
      }
      else if (mystring_contains(line,"hmotion_filename"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) cpmdjson["hmotion_filename"] = ss[1];
      }
      else if (mystring_contains(line,"fei"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) cpmdjson["fei"] = ss[1];
      }
      else if (mystring_contains(line,"dipole_motion"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) cpmdjson["dipole_motion"] = ss[1];
      }
      else if (mystring_contains(line,"intitial_velocities"))
      {
         ss = mystring_split0(line);
         if (ss.size()==2) 
            cpmdjson["initial_velocities"] = {std::stod(ss[1]),12345};
         else if (ss.size()>2) 
            cpmdjson["initial_velocities"] = {std::stod(ss[1]),std::stoi(ss[2])};
          else
            cpmdjson["initial_velocities"] = {298.15,12345};
      }




/*
      else if (mystring_contains(line,"xc"))
      else if (mystring_contains(line,"tolerances"))
      else if (mystring_contains(line,"geometry_optimize"))
      else if (mystring_contains(line,"energy"))
      else if (mystring_contains(line,"temperature"))
      else if (mystring_contains(line,"initial_velocities"))
      else if (mystring_contains(line,"fei_quench"))
      else if (mystring_contains(line,"sa_decay"))
*/


      ++cur;
      if (mystring_contains(lines[cur],"end"))
         --endcount;
   }

   *curptr = cur;

   return cpmdjson;
}



/**************************************************
 *                                                *
 *                parse_nwpw                      *
 *                                                *
 **************************************************/

static json parse_nwpw(json nwpwjson, int *curptr, vector<string> lines)
{
   //json nwpwjson;

   int cur = *curptr;
   int endcount = 1;
   ++cur;
   string line;
   vector<string> ss;

   while (endcount>0)
   {
      line = mystring_lowercase(lines[cur]);

/*
      if (mystring_contains(line,"steepest_descent"))
      if (mystring_contains(line,"car-parrinello"))
*/
      if (mystring_contains(line,"simulation_cell"))
      {
         if  (nwpwjson["simulation_cell"].is_null())
         {
            json simulation_cell; 
            nwpwjson["simulation_cell"] = simulation_cell;
         }
         *curptr = cur;
         nwpwjson["simulation_cell"] = parse_simulation_cell(nwpwjson["simulation_cell"],curptr,lines);
         cur = *curptr;
      }
      else if (mystring_contains(line,"pseudopotentials"))
      {
         if  (nwpwjson["pseudopotentials"].is_null())
         {
            json pseudopotentials; 
            nwpwjson["pseudopotentials"] = pseudopotentials;
         }
         *curptr = cur;
         nwpwjson["pseudopotentials"] = parse_pseudopotentials(nwpwjson["pseudopotentials"],curptr,lines);
         cur = *curptr;
      }
      else if (mystring_contains(line,"steepest_descent"))
      {
         if  (nwpwjson["steepest_descent"].is_null())
         {
            json steepest_descent;
            nwpwjson["steepest_descent"] = steepest_descent;
         }
         *curptr = cur;
         nwpwjson["steepest_descent"] = parse_steepest_descent(nwpwjson["steepest_descent"],curptr,lines);
         cur = *curptr;
      }
      else if (mystring_contains(line,"car-parrinello"))
      {
         if  (nwpwjson["car-parrinello"].is_null())
         {
            json car_parrinello;
            nwpwjson["car-parrinello"] = car_parrinello;
         }
         *curptr = cur;
         nwpwjson["car-parrinello"] = parse_car_parrinello(nwpwjson["car-parrinello"],curptr,lines);
         cur = *curptr;
      }
      else if (mystring_contains(line,"nobalance"))
      {
         nwpwjson["nobalance"] = true;
      }
      else if (mystring_contains(line,"mapping"))
      {
         ss = mystring_split0(line);
         if (ss.size()>1) nwpwjson["mapping"] = std::stoi(ss[1]);
      }
      else if (mystring_contains(line,"np_dimensions"))
      {
         ss = mystring_split0(line);
         if (ss.size()>2) nwpwjson["np_dimensions"] = {std::stoi(ss[1]),std::stoi(ss[2])};
         if (ss.size()>3) nwpwjson["np_dimensions"] = {std::stoi(ss[1]),std::stoi(ss[2]),std::stoi(ss[3])};
      }
      else if (mystring_contains(line,"loop"))
      {
         vector<int> loop;
         loop.push_back(1);
         loop.push_back(1);
         ss = mystring_split0(line);
         if (ss.size()>1) loop[0] = std::stoi(ss[1]);
         if (ss.size()>2) loop[1] = std::stoi(ss[2]);
         nwpwjson["loop"] = loop;
      }
      else if (mystring_contains(line,"xc"))
      {
         nwpwjson["xc"] = mystring_trim(mystring_split(line,"xc")[1]);
      }
      else if (mystring_contains(line,"cutoff"))
      {
         ss = mystring_split0(line);
         if (ss.size()==2) nwpwjson["cutoff"] = {std::stod(ss[1]),2*std::stod(ss[1])};
         if (ss.size()>2)  nwpwjson["cutoff"] = {std::stod(ss[1]),std::stod(ss[2])};
      }
      else if (mystring_contains(line,"intitial_velocities"))
      {
         ss = mystring_split0(line);
         if (ss.size()==2) 
            nwpwjson["initial_velocities"] = {std::stod(ss[1]),12345};
         else if (ss.size()>2) 
            nwpwjson["initial_velocities"] = {std::stod(ss[1]),std::stoi(ss[2])};
          else
            nwpwjson["initial_velocities"] = {298.15,12345};
      }
      else if (mystring_contains(line,"cg"))
      {
         if (mystring_contains(line,"stiefel"))
             nwpwjson["minimizer"] = 4;
         else
            nwpwjson["minimizer"] = 1;
      }
      else if (mystring_contains(line,"lmbfgs"))
      {
         if (mystring_contains(line,"stiefel"))
            nwpwjson["minimizer"] = 7;
         else
            nwpwjson["minimizer"] = 2;
      }
      else if (mystring_contains(line,"scf"))
      {
         if (mystring_contains(line,"potential"))
            nwpwjson["minimizer"] = 5;
         else
            nwpwjson["minimizer"] = 8;
      }
      else if (mystring_contains(line,"vectors"))
      {
         if (mystring_contains(line," input"))
            nwpwjson["input_wavefunction_filename"] = mystring_split0(mystring_trim(mystring_split(line," input")[1]))[0];
         
         if (mystring_contains(line," output"))
            nwpwjson["output_wavefunction_filename"] = mystring_split0(mystring_trim(mystring_split(line," output")[1]))[0];

         if (mystring_contains(line,"vinput"))
            nwpwjson["input_v_wavefunction_filename"] = mystring_split0(mystring_trim(mystring_split(line,"vinput")[1]))[0];

         if (mystring_contains(line,"voutput"))
            nwpwjson["output_v_wavefunction_filename"] = mystring_split0(mystring_trim(mystring_split(line,"voutput")[1]))[0];
      }

      ++cur;
      if (mystring_contains(lines[cur],"end"))
         --endcount;
   }
   //std::cout << "NWPWJSON= " << nwpwjson.dump() << std::endl;

   *curptr = cur;

   return nwpwjson;
}




/**************************************************
 *                                                *
 *                parse_rtdbjson                  *
 *                                                *
 **************************************************/

json parse_rtdbjson(json rtdb)
{
   vector<string> lines = rtdb["nwinput_lines"];
   int n   = rtdb["nwinput_nlines"];
   int cur = rtdb["nwinput_cur"];

   bool foundtask = false;
   while ((cur<n) && (!foundtask))
   {
      if (mystring_contains(mystring_lowercase(lines[cur]),"unset "))
      {
         vector<string> ss = mystring_split0(lines[cur]);
         if (ss.size()>1) auto count_erase = rtdb.erase(ss[1]);
      }
      else if (mystring_contains(mystring_lowercase(lines[cur]),"set "))
      {
         vector<string> ss = mystring_split0(lines[cur]);
         if (ss.size()>2) rtdb[ss[1]] = ss[2];
      }
      else if (mystring_contains(mystring_lowercase(lines[cur]),"start"))
      {
         rtdb["dbname"] = mystring_trim(mystring_split(mystring_split(lines[cur],"start")[1],"\n")[0]);
      }
      else if (mystring_contains(mystring_lowercase(lines[cur]),"geometry"))
      {
         rtdb["geometries"] = parse_geometry(rtdb["geometries"], &cur,lines);
      }
      else if (mystring_contains(mystring_lowercase(lines[cur]),"title"))
      {
         rtdb["title"] = mystring_trim(mystring_ireplace(mystring_split(mystring_ireplace(lines[cur],"TITLE","title"),"title")[1],"\"",""));
      }
      else if (mystring_contains(mystring_lowercase(lines[cur]),"charge"))
      {
         rtdb["charge"] = std::stoi(mystring_trim(mystring_split(mystring_split(lines[cur],"charge")[1],"\n")[0]));
      }
      else if (mystring_contains(mystring_lowercase(lines[cur]),"nwpw"))
      {
         rtdb["nwpw"] = parse_nwpw(rtdb["nwpw"],&cur,lines);
      }
      else if (mystring_contains(mystring_lowercase(lines[cur]),"driver"))
      {
         std::cout << "driver: " << lines[cur] << std::endl;
      }
      else if (mystring_contains(mystring_lowercase(lines[cur]),"task"))
      {
         rtdb["current_task"] = lines[cur];
         foundtask = true;
      }
      else if (mystring_contains(mystring_lowercase(lines[cur]),"print"))
      {
         rtdb["print"] = mystring_trim(mystring_split(mystring_split(lines[cur],"print")[1],"\n")[0]);
      }

      ++cur;
   }
   rtdb["nwinput_cur"] = cur;
   rtdb["foundtask"]   = foundtask;

   return rtdb;
}


/**************************************************
 *                                                *
 *               parse_nwinput                    *
 *                                                *
 **************************************************/

std::string parse_nwinput(std::string nwinput)
{
   // fetch the permanent_dir and scratch_dir
   string permanent_dir = ".";
   string scratch_dir   = ".";
   string psp_library_dir = "";
   if (mystring_contains(mystring_lowercase(nwinput),"permanent_dir"))
   {
      if (!mystring_contains(mystring_trim(mystring_split(mystring_split(nwinput,"permanent_dir")[0],"\n").back()),"#"))
         permanent_dir  = mystring_rtrim_slash(mystring_trim(mystring_split(mystring_split(nwinput,"permanent_dir")[1],"\n")[0]));
   }
   if (mystring_contains(mystring_lowercase(nwinput),"scratch_dir"))
   {
      if (!mystring_contains(mystring_trim(mystring_split(mystring_split(nwinput,"scratch_dir")[0],"\n").back()),"#"))
         scratch_dir  = mystring_rtrim_slash(mystring_trim(mystring_split(mystring_split(nwinput,"scratch_dir")[1],"\n")[0]));
   }
   if (mystring_contains(mystring_lowercase(nwinput),"psp_library_dir"))
   {
      if (!mystring_contains(mystring_trim(mystring_split(mystring_split(nwinput,"psp_library_dir")[0],"\n").back()),"#"))
         psp_library_dir  = mystring_rtrim_slash(mystring_trim(mystring_split(mystring_split(nwinput,"psp_library_dir")[1],"\n")[0]));
   }

   // fetch the dbname
   string dbname = "nwchemex";
   if (mystring_contains(mystring_lowercase(nwinput),"start"))
      dbname  = mystring_trim(mystring_split(mystring_split(nwinput,"start")[1],"\n")[0]);

   json rtdb;
   if (mystring_contains(mystring_lowercase(nwinput),"restart"))
   {
      // read a JSON file
      string dbname0 = permanent_dir + "/" + dbname + ".json";
      std::ifstream ifile(dbname0);
      ifile >> rtdb;
   }
   else
   {
      // intialize the rtdb structure
      json nwpw,geometries;
      rtdb["nwpw"]     = nwpw;
      rtdb["geometries"] = geometries;
   }

   // set the dbname
   rtdb["dbname"] = dbname;

   // set the permanent_dir and scratch_dir
   rtdb["permanent_dir"] = permanent_dir;
   rtdb["scratch_dir"]   = scratch_dir;
   rtdb["psp_library_dir"] = psp_library_dir;


   // split nwinput into lines
   vector<string> lines = mystring_split(nwinput,"\n");

   // Remove comments
   for (auto i = lines.begin(); i != lines.end(); ++i)
      *i = mystring_split(*i,"#")[0];


   rtdb["nwinput_lines"]  = lines;
   rtdb["nwinput_nlines"] = lines.size();
   rtdb["nwinput_cur"]    = 0;

   rtdb = parse_rtdbjson(rtdb);
  

   return rtdb.dump();
}


/**************************************************
 *                                                *
 *               parse_rtdbstring                 *
 *                                                *
 **************************************************/

std::string parse_rtdbstring(std::string rtdbstring)
{
   auto rtdb =  json::parse(rtdbstring);

   rtdb = parse_rtdbjson(rtdb);
   
   return rtdb.dump();
}


/**************************************************
 *                                                *
 *                 parse_task                     *
 *                                                *
 **************************************************/

int parse_task(std::string rtdbstring)
{
   auto rtdb =  json::parse(rtdbstring);
   int task = 0;
   if (rtdb["foundtask"])
   {
      // Look for pspw jobs
      if (mystring_contains(mystring_lowercase(rtdb["current_task"]),"pspw"))
      {
         if (mystring_contains(mystring_lowercase(rtdb["current_task"]),"energy"))           task = 1;
         if (mystring_contains(mystring_lowercase(rtdb["current_task"]),"gradient"))         task = 2;
         if (mystring_contains(mystring_lowercase(rtdb["current_task"]),"optimize"))         task = 3;
         if (mystring_contains(mystring_lowercase(rtdb["current_task"]),"freq"))             task = 4;
         if (mystring_contains(mystring_lowercase(rtdb["current_task"]),"steepest_descent")) task = 5;
         if (mystring_contains(mystring_lowercase(rtdb["current_task"]),"car-parrinello"))   task = 6;
      }
   }

   return task;
}

/**************************************************
 *                                                *
 *       parse_input_wavefunction_filename        *
 *                                                *
 **************************************************/
std::string parse_input_wavefunction_filename(std::string rtdbstring)
{
   auto rtdb =  json::parse(rtdbstring);

   std::string filename  = "eric.movecs";
   if (rtdb["dbname"].is_string())
   {
       string dbname = rtdb["dbname"];
       filename = dbname + ".movecs";
   }
   // read from nwpw block
   if (rtdb["nwpw"]["input_wavefunction_filename"].is_string()) filename = rtdb["nwpw"]["input_wavefunction_filename"];

   return filename;
}

/**************************************************
 *                                                *
 *           parse_gen_lowlevel_rtdbstrs          *
 *                                                *
 **************************************************/
std::vector<std::string> parse_gen_lowlevel_rtdbstrs(std::string rtdbstring)
{
   auto rtdbjson =  json::parse(rtdbstring);

   double ecut0=100.0;
   double wcut0=50.0;
   if (rtdbjson["nwpw"]["cutoff"][0].is_number_float()) wcut0 = rtdbjson["nwpw"]["cutoff"][0];
   if (rtdbjson["nwpw"]["cutoff"][1].is_number_float()) ecut0 = rtdbjson["nwpw"]["cutoff"][1];

   int steps;
   double dx;
   if (wcut0<5.0) {
      steps = 0; dx = 1.0;
   } else if (wcut0<=20.0) {
      steps = 2; dx = 1.0/2.0;
   } else if (wcut0<=30.0) {
      steps = 3; dx = 1.0/3.0;
   } else {
      steps = 4; dx = 1.0/4.0;
   }

   std::vector<std::string> gen_rtdbs;

   for (int i=1; i<steps; ++i)
   {
      rtdbjson["nwpw"]["cutoff"][0] = i*dx*wcut0;
      rtdbjson["nwpw"]["cutoff"][1] = i*dx*ecut0;
      rtdbjson["current_task"]      = "energy";

      gen_rtdbs.push_back(rtdbjson.dump());
   }

   return gen_rtdbs;
}


/**************************************************
 *                                                *
 *                 parse_write                    *
 *                                                *
 **************************************************/

void parse_write(std::string rtdbstring)
{
   auto rtdbjson =  json::parse(rtdbstring);
   std::string pdir = rtdbjson["permanent_dir"];
   std::string dbname0 = rtdbjson["dbname"];
   std::cout << "writing rtdbjson = " <<  pdir + "/" + dbname0 + ".json" << std::endl;
   std::ofstream ofile(pdir + "/" + dbname0 + ".json");
   ofile << std::setw(4) << rtdbjson << std::endl;
}


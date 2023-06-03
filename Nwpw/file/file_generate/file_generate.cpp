#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

#include "iofmt.hpp"
#include "util_date.hpp"

//#include	"control.hpp"
//#include "mpi.h"


#include "json.hpp"
using json = nlohmann::json;

#include "parsestring.hpp"

namespace pwdft {

/******************************************
 *                                        *
 *               file_generate            *
 *                                        *
 ******************************************/
int file_generate(std::string &rtdbstring)
{
   auto rtdbjson = json::parse(rtdbstring);

   std::cout << "current_task=" << rtdbjson["current_task"] << std::endl;

   // XYZ file generation
   if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]), "xyz")) 
   {
      //std::vector<std::string> ss;
      auto ss = mystring_split0(rtdbjson["current_task"]);
      std::string xyz_filename = "eric.xyz";
      if (ss.size() > 1) 
          xyz_filename = ss[ss.size()-1];

      std::string geomname = "geometry";
      if (rtdbjson["geometry"].is_string())
         geomname = rtdbjson["geometry"];

      json geomjson = rtdbjson["geometries"][geomname];
      json simujson = rtdbjson["nwpw"]["simulation_cell"];
      auto nion    = geomjson["nion"];
      auto symbols = geomjson["symbols"];
      double AACONV = 0.529177;


      double unita[9];
      unita[0] = 20.0;
      unita[1] = 0.0;
      unita[2] = 0.0;
      unita[3] = 0.0;
      unita[4] = 20.0;
      unita[5] = 0.0;
      unita[6] = 0.0;
      unita[7] = 0.0;
      unita[8] = 20.0;
 
      if (geomjson["unita"][0].is_number_float()) unita[0] = geomjson["unita"][0];
      if (geomjson["unita"][1].is_number_float()) unita[1] = geomjson["unita"][1];
      if (geomjson["unita"][2].is_number_float()) unita[2] = geomjson["unita"][2];
 
      if (geomjson["unita"][3].is_number_float()) unita[3] = geomjson["unita"][3];
      if (geomjson["unita"][4].is_number_float()) unita[4] = geomjson["unita"][4];
      if (geomjson["unita"][5].is_number_float()) unita[5] = geomjson["unita"][5];
 
      if (geomjson["unita"][6].is_number_float()) unita[6] = geomjson["unita"][6];
      if (geomjson["unita"][7].is_number_float()) unita[7] = geomjson["unita"][7];
      if (geomjson["unita"][8].is_number_float()) unita[8] = geomjson["unita"][8];

      if (simujson["unita"][0].is_number_float()) unita[0] = simujson["unita"][0];
      if (simujson["unita"][1].is_number_float()) unita[1] = simujson["unita"][1];
      if (simujson["unita"][2].is_number_float()) unita[2] = simujson["unita"][2];
 
      if (simujson["unita"][3].is_number_float()) unita[3] = simujson["unita"][3];
      if (simujson["unita"][4].is_number_float()) unita[4] = simujson["unita"][4];
      if (simujson["unita"][5].is_number_float()) unita[5] = simujson["unita"][5];
 
      if (simujson["unita"][6].is_number_float()) unita[6] = simujson["unita"][6];
      if (simujson["unita"][7].is_number_float()) unita[7] = simujson["unita"][7];
      if (simujson["unita"][8].is_number_float()) unita[8] = simujson["unita"][8];


      //std::ofstream xyz_stream(xyz_filename,std::ios::app);
      std::ofstream xyz_stream;
      if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]), " new")) 
         xyz_stream.open(xyz_filename);
      else
         xyz_stream.open(xyz_filename, std::ios_base::app);

      xyz_stream << nion << std::endl;
      xyz_stream << std::fixed << std::setprecision(6) << std::setw(12)
                 << AACONV*unita[0] << std::setw(12) << AACONV*unita[1] << std::setw(12) << AACONV*unita[2] << std::setw(12)
                 << AACONV*unita[3] << std::setw(12) << AACONV*unita[4] << std::setw(12) << AACONV*unita[5] << std::setw(12)
                 << AACONV*unita[6] << std::setw(12) << AACONV*unita[7] << std::setw(12) << AACONV*unita[8] << std::endl;

      for (auto ii=0; ii<nion; ++ii) 
      {
         std::string sym1 = symbols[ii];
         if (sym1.length() == 1) sym1 += " ";
         sym1 += "    ";
         double xx = geomjson["coords"][3*ii];
         double yy = geomjson["coords"][3*ii+1];
         double zz = geomjson["coords"][3*ii+2];
         xyz_stream << sym1 << std::fixed << std::setprecision(6) << std::setw(12)
                    << AACONV * xx << " " << std::setw(12)
                    << AACONV * yy << " " << std::setw(12)
                    << AACONV * zz << std::endl;
      }
   }

   // CIF file generation
   else if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]), "cif")) 
   {
      double shift = 0.5;
      if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]), "shift")) 
          shift = std::stod(mystring_split0(mystring_split(rtdbjson["current_task"],"shift")[1])[0]);

      std::cout << " - cell shifted (shift " << shift << ")" << std::endl;

      //std::vector<std::string> ss;
      auto ss = mystring_split0(rtdbjson["current_task"]);
      std::string cif_filename = "eric.cif";
      if (ss.size() > 1)
          cif_filename = ss[ss.size()-1];

      std::string geomname = "geometry";
      if (rtdbjson["geometry"].is_string())
         geomname = rtdbjson["geometry"];

      json geomjson = rtdbjson["geometries"][geomname];
      json simujson = rtdbjson["nwpw"]["simulation_cell"];
      auto nion    = geomjson["nion"];
      auto symbols = geomjson["symbols"];
      double AACONV = 0.529177;


      double a[9];
      a[0] = 20.0;
      a[1] = 0.0;
      a[2] = 0.0;
      a[3] = 0.0;
      a[4] = 20.0;
      a[5] = 0.0;
      a[6] = 0.0;
      a[7] = 0.0;
      a[8] = 20.0;

      if (geomjson["unita"][0].is_number_float()) a[0] = geomjson["unita"][0];
      if (geomjson["unita"][1].is_number_float()) a[1] = geomjson["unita"][1];
      if (geomjson["unita"][2].is_number_float()) a[2] = geomjson["unita"][2];

      if (geomjson["unita"][3].is_number_float()) a[3] = geomjson["unita"][3];
      if (geomjson["unita"][4].is_number_float()) a[4] = geomjson["unita"][4];
      if (geomjson["unita"][5].is_number_float()) a[5] = geomjson["unita"][5];

      if (geomjson["unita"][6].is_number_float()) a[6] = geomjson["unita"][6];
      if (geomjson["unita"][7].is_number_float()) a[7] = geomjson["unita"][7];
      if (geomjson["unita"][8].is_number_float()) a[8] = geomjson["unita"][8];

      if (simujson["unita"][0].is_number_float()) a[0] = simujson["unita"][0];
      if (simujson["unita"][1].is_number_float()) a[1] = simujson["unita"][1];
      if (simujson["unita"][2].is_number_float()) a[2] = simujson["unita"][2];

      if (simujson["unita"][3].is_number_float()) a[3] = simujson["unita"][3];
      if (simujson["unita"][4].is_number_float()) a[4] = simujson["unita"][4];
      if (simujson["unita"][5].is_number_float()) a[5] = simujson["unita"][5];

      if (simujson["unita"][6].is_number_float()) a[6] = simujson["unita"][6];
      if (simujson["unita"][7].is_number_float()) a[7] = simujson["unita"][7];
      if (simujson["unita"][8].is_number_float()) a[8] = simujson["unita"][8];

      double b[9];
      b[0] = a[4]*a[8] - a[5]*a[7];
      b[1] = a[5]*a[6] - a[3]*a[8];
      b[2] = a[3]*a[7] - a[4]*a[6];
 
      b[3] = a[7]*a[2] - a[8]*a[1];
      b[4] = a[8]*a[0] - a[6]*a[2];
      b[5] = a[6]*a[1] - a[7]*a[0];
 
      b[6] = a[1]*a[5] - a[2]*a[4];
      b[7] = a[2]*a[3] - a[0]*a[5];
      b[8] = a[0]*a[4] - a[1]*a[3];
      double volume = a[0]*b[0] + a[1]*b[1]+ a[2]*b[2];
 
      for (auto i=0; i<9; ++i)
         b[i] /= volume;
 
      // determine a,b,c,alpha,beta,gmma
      double pi = 4.0*std::atan(1.0);
      double aa = std::sqrt(a[0]*a[0] + a[1]*a[1] +a[2]*a[2]);
      double bb = std::sqrt(a[3]*a[3] + a[4]*a[4] +a[5]*a[5]);
      double cc = std::sqrt(a[6]*a[6] + a[7]*a[7] +a[8]*a[8]);
 
      double da = (a[3]-a[6])*(a[3]-a[6])
                + (a[4]-a[7])*(a[4]-a[7])
                + (a[5]-a[8])*(a[5]-a[8]);
      double alpha = (bb*bb + cc*cc - da)/(2.0*bb*cc);
      alpha = std::acos(alpha)*180.0/pi;
 
      double db = (a[6]-a[0])*(a[6]-a[0])
                + (a[7]-a[1])*(a[7]-a[1])
                + (a[8]-a[2])*(a[8]-a[2]);
      double beta = (cc*cc + aa*aa - db)/(2.0*cc*aa);
      beta = std::acos(beta)*180.0/pi;
 
      double dg = (a[0]-a[3])*(a[0]-a[3])
                + (a[1]-a[4])*(a[1]-a[4])
                + (a[2]-a[5])*(a[2]-a[5]);
      double gmma = (aa*aa + bb*bb - dg)/(2.0*aa*bb);
      gmma = std::acos(gmma)*180.0/pi;


      std::ofstream cif_stream;
      if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]), " new"))
         cif_stream.open(cif_filename);
      else
         cif_stream.open(cif_filename, std::ios_base::app);

      cif_stream << "data_nwchemex_pwdft" << std::endl << std::endl;
      cif_stream << "_audit_creation_date   " << util_date() << std::endl;
      cif_stream << "_audit_creation_method    generated by pwdft module of NWChemEx" << std::endl;
      cif_stream<< std::endl << std::endl;
      cif_stream << "_cell_length_a   " << Ffmt(16,4) << aa * 0.529177 << std::endl;
      cif_stream << "_cell_length_b   " << Ffmt(16,4) << bb * 0.529177 << std::endl;
      cif_stream << "_cell_length_c   " << Ffmt(16,4) << cc * 0.529177 << std::endl;
      cif_stream << "_cell_angle_alpha" << Ffmt(16,4) << alpha << std::endl;
      cif_stream << "_cell_angle_beta " << Ffmt(16,4) << beta  << std::endl;
      cif_stream << "_cell_angle_gamma" << Ffmt(16,4) << gmma  << std::endl;
      cif_stream << std::endl;
      cif_stream << "_symmetry_space_group_name_H-M     P1  " << std::endl;
      cif_stream << std::endl;
      cif_stream << "loop_" << std::endl;
      cif_stream << "_atom_site_type_symbol" << std::endl;
      //stream << 1242 FORMAT('_atom_site_label')
      cif_stream << "_atom_site_fract_x" << std::endl;
      cif_stream << "_atom_site_fract_y" << std::endl;
      cif_stream << "_atom_site_fract_z" << std::endl;
      for (auto ii=0; ii<nion; ++ii)
      {
         std::string sym1 = symbols[ii];
         if (sym1.length() == 1) sym1 += " ";
         sym1 += "      ";
         double xx = geomjson["coords"][3*ii];
         double yy = geomjson["coords"][3*ii+1];
         double zz = geomjson["coords"][3*ii+2];

         double f0 = b[0]*xx + b[1]*yy + b[2]*zz + shift;
         double f1 = b[3]*xx + b[4]*yy + b[5]*zz + shift;
         double f2 = b[6]*xx + b[7]*yy + b[8]*zz + shift;
  
         cif_stream << sym1 << std::fixed << std::setprecision(6) 
                    << std::setw(14) << f0 
                    << std::setw(14) << f1 
                    << std::setw(14) << f2 
                    << std::endl;
      }
   }
   // ion_motion file generation
   else if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]), "ion_motion")) 
   {
      double time = 0.0;
      if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]), "time")) 
          time = std::stod(mystring_split0(mystring_split(rtdbjson["current_task"],"time")[1])[0]);

      std::cout << " - ion_motion time = " << time << std::endl;

      //std::vector<std::string> ss;
      auto ss = mystring_split0(rtdbjson["current_task"]);
      std::string ion_motion_filename = "eric.ion_motion";
      if (ss.size() > 1)
          ion_motion_filename = ss[ss.size()-1];

      std::string geomname = "geometry";
      if (rtdbjson["geometry"].is_string())
         geomname = rtdbjson["geometry"];

      json geomjson = rtdbjson["geometries"][geomname];
      json simujson = rtdbjson["nwpw"]["simulation_cell"];
      int  nion    = geomjson["nion"];
      auto symbols = geomjson["symbols"];

      double a1x = 20.0;
      double a1y = 0.0;
      double a1z = 0.0;

      double a2x = 0.0;
      double a2y = 20.0;
      double a2z = 0.0;

      double a3x = 0.0;
      double a3y = 0.0;
      double a3z = 20.0;

      if (geomjson["unita"][0].is_number_float()) a1x = geomjson["unita"][0];
      if (geomjson["unita"][1].is_number_float()) a1y = geomjson["unita"][1];
      if (geomjson["unita"][2].is_number_float()) a1z = geomjson["unita"][2];

      if (geomjson["unita"][3].is_number_float()) a2x = geomjson["unita"][3];
      if (geomjson["unita"][4].is_number_float()) a2y = geomjson["unita"][4];
      if (geomjson["unita"][5].is_number_float()) a2z = geomjson["unita"][5];

      if (geomjson["unita"][6].is_number_float()) a3x = geomjson["unita"][6];
      if (geomjson["unita"][7].is_number_float()) a3y = geomjson["unita"][7];
      if (geomjson["unita"][8].is_number_float()) a3z = geomjson["unita"][8];

      if (simujson["unita"][0].is_number_float()) a1x = simujson["unita"][0];
      if (simujson["unita"][1].is_number_float()) a1y = simujson["unita"][1];
      if (simujson["unita"][2].is_number_float()) a1z = simujson["unita"][2];

      if (simujson["unita"][3].is_number_float()) a2x = simujson["unita"][3];
      if (simujson["unita"][4].is_number_float()) a2y = simujson["unita"][4];
      if (simujson["unita"][5].is_number_float()) a2z = simujson["unita"][5];

      if (simujson["unita"][6].is_number_float()) a3x = simujson["unita"][6];
      if (simujson["unita"][7].is_number_float()) a3y = simujson["unita"][7];
      if (simujson["unita"][8].is_number_float()) a3z = simujson["unita"][8];

      double bx1 = (a2y*a3z - a2z*a3y);
      double by1 = (a2z*a3x - a2x*a3z);
      double bz1 = (a2x*a3y - a2y*a3x);
      double omega = std::fabs(a1x*bx1 + a1y*by1 + a1z*bz1);


      std::ofstream ion_motion_stream;
      if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]), " new"))
         ion_motion_stream.open(ion_motion_filename);
      else
         ion_motion_stream.open(ion_motion_filename, std::ios_base::app);

      ion_motion_stream << Efmt(15,6) << time << Ifmt(6) << nion << Efmt(15,6) << omega
                        << Efmt(15,6) << a1x << Efmt(15,6) << a1y << Efmt(15,6) << a1z
                        << Efmt(15,6) << a2x << Efmt(15,6) << a2y << Efmt(15,6) << a2z
                        << Efmt(15,6) << a3x << Efmt(15,6) << a3y << Efmt(15,6) << a3z
                        << std::endl;;

      for (auto ii=0; ii<nion; ++ii) 
      {
         std::string sym1 = symbols[ii];
         double xx = geomjson["coords"][3*ii];
         double yy = geomjson["coords"][3*ii+1];
         double zz = geomjson["coords"][3*ii+2];
         double vx = geomjson["velocities"][3*ii];
         double vy = geomjson["velocities"][3*ii+1];
         double vz = geomjson["velocities"][3*ii+2];
         ion_motion_stream << Ifmt(6) << ii+1 << " " 
                           << Lfmt(3) << sym1 << Lfmt(5) << sym1
                           << Efmt(15,6) << xx << Efmt(15,6) << yy << Efmt(15,6) << zz
                           << Efmt(15,6) << vx << Efmt(15,6) << vy << Efmt(15,6) << vz
                           << std::endl;
      }
   }

   return 0;
}

} // namespace pwdft

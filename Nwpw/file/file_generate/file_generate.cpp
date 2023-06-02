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
      auto nion    = geomjson["nion"];
      auto symbols = geomjson["symbols"];
      double AACONV = 0.529177;


      //std::ofstream xyz_stream(xyz_filename,std::ios::app);
      std::ofstream xyz_stream;
      if (mystring_contains(mystring_lowercase(rtdbjson["current_task"]), " new")) 
         xyz_stream.open(xyz_filename);
      else
         xyz_stream.open(xyz_filename, std::ios_base::app);

      xyz_stream << nion << std::endl << std::endl;
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

   return 0;
}

} // namespace pwdft

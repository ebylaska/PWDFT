/* Brillouin.cpp -
   Author - Eric Bylaska
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <sstream>


#include "json.hpp"
using json = nlohmann::json;

#include "iofmt.hpp"
#include "Control2.hpp"
#include "Lattice.hpp"
#include "Brillouin.hpp"


namespace pwdft {


Brillouin::Brillouin(std::string rtdbstring,  Lattice *mylattice, Control2& control) 
{
 
   auto rtdbjson = json::parse(rtdbstring);
 
 
   json brillouinjson = rtdbjson["nwpw"]["brillouin_zone"];

   // Optional Monkhorst-Pack provenance (nice for printing / reproducibility)
   if (brillouinjson.contains("monkhorst-pack") && brillouinjson["monkhorst-pack"].is_array())
   {
      auto mp = brillouinjson["monkhorst-pack"];
      if (mp.size() > 0) nkx = mp[0].get<int>();
      if (mp.size() > 1) nky = mp[1].get<int>();
      if (mp.size() > 2) nkz = mp[2].get<int>();
   }
   if (brillouinjson.contains("monkhorst-pack-shift") && brillouinjson["monkhorst-pack-shift"].is_array())
   {
      auto sh = brillouinjson["monkhorst-pack-shift"];
      if (sh.size() > 0) skx = sh[0].get<double>();
      if (sh.size() > 1) sky = sh[1].get<double>();
      if (sh.size() > 2) skz = sh[2].get<double>();
   }
 
   bool nobrillread = false;
   if (brillouinjson["kvectors"].is_null())
   {
      nbrillouin  = 1;
      nbrillouin0 = 1;
      nobrillread = true;
   }
   else
   {
      // Effective list we actually use
      nbrillouin  = static_cast<int>(brillouinjson["kvectors"].size());
 
      // Default: no reduction info known
      nbrillouin0 = nbrillouin;
 
      // Optional metadata (preferred)
      if (brillouinjson.contains("num_kpoints_original") &&
          !brillouinjson["num_kpoints_original"].is_null())
          nbrillouin0 = brillouinjson["num_kpoints_original"].get<int>();
 
      if (brillouinjson.contains("num_kpoints_effective") &&
          !brillouinjson["num_kpoints_effective"].is_null())
      {
          int neff = brillouinjson["num_kpoints_effective"].get<int>();
          // sanity: if it disagrees with kvectors.size(), kvectors wins
          // (but keep the check for debugging)
          // if (neff != nbrillouin) ... optional warning
      }
   }

   weight   = new double[nbrillouin];
   kvector  = new double[3*nbrillouin];
   ksvector = new double[3*nbrillouin];

   if (nobrillread)
   {
      kvector[0] = ksvector[0] = 0.0;
      kvector[1] = ksvector[1] = 0.0;
      kvector[2] = ksvector[2] = 0.0;
      weight[0]  = 1.0;
   }
   else
      for (auto nb=0; nb<nbrillouin; ++nb) 
      {
         ksvector[3*nb]   = brillouinjson["kvectors"][nb][0];
         ksvector[3*nb+1] = brillouinjson["kvectors"][nb][1];
         ksvector[3*nb+2] = brillouinjson["kvectors"][nb][2];
         weight[nb]       = brillouinjson["kvectors"][nb][3];
      
         kvector[3*nb]   = ksvector[3*nb]  *mylattice->unitg(0,0)
                         + ksvector[3*nb+1]*mylattice->unitg(0,1)
                         + ksvector[3*nb+2]*mylattice->unitg(0,2);
      
         kvector[3*nb+1] = ksvector[3*nb]  *mylattice->unitg(1,0)
                         + ksvector[3*nb+1]*mylattice->unitg(1,1)
                         + ksvector[3*nb+2]*mylattice->unitg(1,2);
      
         kvector[3*nb+2] = ksvector[3*nb]  *mylattice->unitg(2,0)
                         + ksvector[3*nb+1]*mylattice->unitg(2,1)
                         + ksvector[3*nb+2]*mylattice->unitg(2,2);
      }
}

/*******************************************
 *                                         *
 *        Brillouin::print_zone            *
 *                                         *
 *******************************************/
std::string Brillouin::print_zone()
{
   std::stringstream stream;

   std::ios init(NULL);
   init.copyfmt(stream);

   //stream << "      number of zone points = " << Ifmt(3) << nbrillouin << std::endl;
   //stream << "      number of zone points = " << Ifmt(3) << nbrillouin << " (reduced from " << nbrillouin0 << " zone points without symmetry)" << std::endl;
   if (nkx > 0 && nky > 0 && nkz > 0)
   {
      stream << "      monkhorst-pack grid = "
             << nkx << " " << nky << " " << nkz;
 
      // Only print shift if it’s non-trivial (or always, your choice)
      if (std::fabs(skx) > 0.0 || std::fabs(sky) > 0.0 || std::fabs(skz) > 0.0)
         stream << "   shift = " << skx << " " << sky << " " << skz;
 
      stream << "\n";
   }

   if (nbrillouin0 != nbrillouin)
      stream << "      number of zone points = " << Ifmt(3) << nbrillouin
             << " (reduced from " << nbrillouin0 << " zone points without symmetry)\n";
   else
      stream << "      number of zone points = " << Ifmt(3) << nbrillouin
             << " (no symmetry reduction)\n";

   for (auto nb=0; nb<nbrillouin; ++nb)
   {
      stream << "      weight = " << Ffmt(8,3) << weight[nb]
             << " ks = <" << Ffmt(8,3) << ksvector[3*nb]   << " "
                          << Ffmt(8,3) << ksvector[3*nb+1] << " "
                          << Ffmt(8,3) << ksvector[3*nb+2] << "> "
             << " k = <" << Ffmt(8,3) << kvector[3*nb]   << " "
                         << Ffmt(8,3) << kvector[3*nb+1] << " "
                         << Ffmt(8,3) << kvector[3*nb+2] << "> " 
             << std::endl;
   }
   return stream.str();
}

/*******************************************
 *                                         *
 *        Brillouin::print_zone_point      *
 *                                         *
 *******************************************/
std::string Brillouin::print_zone_point(const int nb)
{
   std::stringstream stream;

   std::ios init(NULL);
   init.copyfmt(stream);

   stream << " Brillouin zone point: " << Ifmt(3) << nb+1 << std::endl;
   stream << "   weight = " << Ffmt(8,3) << weight[nb] << std::endl
          << "   k = <" << Ffmt(8,3) << ksvector[3*nb]   << " "
                       << Ffmt(8,3) << ksvector[3*nb+1] << " "
                       << Ffmt(8,3) << ksvector[3*nb+2] << "> . <b1,b2,b3>" << std::endl 
          << "     = <" << Ffmt(8,3) << kvector[3*nb]   << " "
                      << Ffmt(8,3) << kvector[3*nb+1] << " "
                      << Ffmt(8,3) << kvector[3*nb+2] << "> " << std::endl;
   return stream.str();
}





} // namespace pwdft

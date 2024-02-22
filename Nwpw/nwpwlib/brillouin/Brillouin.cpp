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
 
   bool nobrillread = false;
   if (brillouinjson["kvectors"].is_null())
   {
       nbrillouin = 1;
       nobrillread = true;
   }
   else
      nbrillouin = brillouinjson["kvectors"].size();

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

   stream << "      number of zone points = " << Ifmt(3) << nbrillouin << std::endl;
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

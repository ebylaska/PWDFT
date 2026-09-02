#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
//
//#include "iofmt.hpp"
#include "Parallel.hpp"
#include "util_date.hpp"
#include "mpi.h"

//#include "gdevice.hpp"


#include "json.hpp"
using json = nlohmann::json;

using minimizer_function = int (*)(MPI_Comm, std::string&, std::ostream&);

namespace pwdft {

static bool read_unita( const json& value, std::array<double, 9>& unita)
{
   if (!value.is_array() || value.size() != 9)
       return false;

   for (int i = 0; i < 9; ++i)
   {
      if (!value[i].is_number())
          return false;
      unita[i] = value[i].get<double>();
   }

   return true;
}

static void write_unita(json& value, const std::array<double, 9>& unita)
{
   value = json::array();

   for (double x : unita)
      value.push_back(x);
}

static double unita_relative_difference(const std::array<double, 9>& current, const std::array<double, 9>& frozen)
{
   double difference_squared = 0.0;
   double frozen_squared = 0.0;

   for (int i = 0; i < 9; ++i)
   {
      const double difference = current[i] - frozen[i];
      difference_squared += difference*difference;
      frozen_squared     += frozen[i]*frozen[i];
   }

   const double denominator = std::max(1.0, std::sqrt(frozen_squared));
   return std::sqrt(difference_squared)/denominator;
}


static bool update_unita_frozen(std::string& rtdbstring, std::ostream& coutput, const double lattice_tolerance)
{
   if (lattice_tolerance < 0.0)
   {
       coutput << "driver_optimizer: negative lattice tolerance\n";
       return false;
   }

   json rtdbjson;

   try
   {
       rtdbjson = json::parse(rtdbstring);
   }
   catch (const json::exception& ex)
   {
       coutput
           << "driver_optimizer: invalid RTDB JSON: "
           << ex.what()
           << '\n';

       return false;
   }

   const std::string geomname =
       (rtdbjson.contains("geometry") &&
        rtdbjson["geometry"].is_string())
           ? rtdbjson["geometry"].get<std::string>()
           : "geometry";

   if (!rtdbjson.contains("geometries") ||
       !rtdbjson["geometries"].is_object() ||
       !rtdbjson["geometries"].contains(geomname))
   {
       coutput
           << "driver_optimizer: geometry '"
           << geomname
           << "' not found\n";

       return false;
   }

   const json& geometry =
       rtdbjson["geometries"][geomname];

   std::array<double, 9> current_unita{};

   if (!read_unita(
           geometry.value("unita", json{}),
           current_unita))
   {
       coutput
           << "driver_optimizer: invalid current geometry lattice\n";

       return false;
   }

   /*
    * If the NWPW simulation-cell lattice exists, use it as the
    * current physical lattice. This preserves the intended precedence.
    */
   if (rtdbjson.contains("nwpw") &&
       rtdbjson["nwpw"].is_object() &&
       rtdbjson["nwpw"].contains("simulation_cell") &&
       rtdbjson["nwpw"]["simulation_cell"].is_object())
   {
       const json& simulation_cell =
           rtdbjson["nwpw"]["simulation_cell"];

       std::array<double, 9> simulation_unita{};

       if (read_unita(
               simulation_cell.value("unita", json{}),
               simulation_unita))
       {
           current_unita = simulation_unita;
       }
   }

   json& simulation_cell =
       rtdbjson["nwpw"]["simulation_cell"];

   if (!simulation_cell.is_object())
       simulation_cell = json::object();

   std::array<double, 9> frozen_unita{};

   const bool has_frozen_unita =
       read_unita(
           simulation_cell.value(
               "unita_frozen",
               json{}),
           frozen_unita);

   if (!has_frozen_unita)
   {
       write_unita(
           simulation_cell["unita_frozen"],
           current_unita);

       coutput
           << "driver_optimizer: initialized unita_frozen\n";
   }
   else
   {
       const double difference =
           unita_relative_difference(
               current_unita,
               frozen_unita);

       if (difference > lattice_tolerance)
       {
           write_unita(
               simulation_cell["unita_frozen"],
               current_unita);

           coutput
               << "driver_optimizer: resetting unita_frozen\n"
               << "  relative lattice change = "
               << difference
               << "\n"
               << "  tolerance               = "
               << lattice_tolerance
               << '\n';
       }
   }

   /*
    * Critical: pass the modified RTDB back to the caller/minimizer.
    */
   rtdbstring = rtdbjson.dump();

   return true;
}

/******************************************
 *                                        *
 *            driver_optimizer            *
 *                                        *
 ******************************************/
/**
 * @brief Common top-level optimization driver.
 *
 * The driver performs shared workflow-level setup and indirectly invokes
 * the selected electronic-structure backend through a function callback.
 *
 * The callback may be:
 *
 *     pspw_minimizer
 *     band_minimizer
 *
 * For unit-cell optimization, the RTDB variable `unita_frozen` stores
 * the reference lattice used to define the fixed numerical grid and
 * plane-wave basis support. The current lattice may change during the
 * optimization, while `unita_frozen` remains fixed for the optimization
 * stage. This prevents repeated grid/basis changes and follows the
 * NWChem-style driver architecture.
 */
int driver_optimizer(MPI_Comm comm_world0, std::string &rtdbstring, std::ostream &coutput, minimizer_function minimizer)
{
   if (minimizer == nullptr)
   {
      coutput << "driver_optimizer: null minimizer\n";
      return 1;
   }

   Parallel myparallel(comm_world0);
   bool oprint = myparallel.is_master();

   if (oprint) 
   {
      std::ios_base::sync_with_stdio();
      constexpr int width = 78;

      coutput << '\n'
              << std::string(width, '=')
              << '\n'
              << "                         PWDFT OPTIMIZATION DRIVER\n"
              << std::string(width, '=')
              << '\n'
              << '\n'
              << "  Architecture           : NWChem-style driver dispatch\n"
              << "  Role                   : top-level orchestration layer\n"
              << "  Backend                : PSPW or band minimizer callback\n"
              << "  Implementation         : NorthwestEx C++ driver\n"
              << "  Method                 : Grassmann/Stiefel manifold\n"
              << "  Cell optimization      : RTDB unita_frozen reference lattice\n"
              << "  Numerical grid         : fixed during optimization stage\n"
              << "  Date                   : "
              << util_date()
              << '\n'
              << '\n'
              << "  The current lattice may change during unit-cell optimization.\n"
              << "  The RTDB variable unita_frozen preserves the reference lattice\n"
              << "  used to define the numerical grid and basis support.\n"
              << "  For large unitcell changes the unita_frozen is changed.\n"
              << "  The reference lattice used to establish the numerical grid and\n"
              << "  basis-support policy. It remains unchanged during the current \n"
              << "  optimization stage but is not the physical lattice used for the\n"
              << "  current energy, force, or stress evaluation.\n"
              << '\n'
              << std::string(width, '-')
              << '\n';
   }
  
   /*
    * Common driver-level work goes here.
    *
    * For now, if the driver is only dispatching, call the
    * selected minimizer directly.
    */

   // add unita_frozen in rtdb
   double lattice_tolerance = 0.01;
   if (!update_unita_frozen(rtdbstring,coutput,lattice_tolerance))
      return 1;

    
   return minimizer(comm_world0, rtdbstring, coutput);

}


} // namespace pwdft 


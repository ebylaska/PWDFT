
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <iomanip>
//
//#include "iofmt.hpp"
#include "Parallel.hpp"
#include "Control2.hpp"
#include "Lattice.hpp"
#include "Ion.hpp"
#include "util_date.hpp"
#include "mpi.h"

//#include "gdevice.hpp"


struct SymmetryInfo {
    std::string space_group_name = "unknown";
    std::string type = "unknown";
    int group_order = -1;
    bool is_primitive = false;
    bool is_cubic = false;

    std::string system = "unknown"; // now 'system' defines your symmetry constraints!

    // Add more fields if needed

    // Helper for test
    bool has_symmetry() const {
        return (space_group_name != "unknown" && group_order > 0);
    }
};



#include "json.hpp"
using json = nlohmann::json;

using minimizer_function = int (*)(MPI_Comm, std::string&, std::ostream&);

namespace pwdft {

static void print_lattice_state( const json& rtdbjson, std::ostream& coutput, const std::string& tag)
{
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
            << tag
            << "Geometry '"
            << geomname
            << "' is missing\n";

        return;
    }

    const json& geometry =
        rtdbjson["geometries"].at(geomname);

    const json* current_unita =
        geometry.contains("unita") &&
        geometry["unita"].is_array()
            ? &geometry["unita"]
            : nullptr;

    const json* frozen_unita =
        nullptr;

    if (rtdbjson.contains("nwpw") &&
        rtdbjson["nwpw"].is_object() &&
        rtdbjson["nwpw"].contains("simulation_cell") &&
        rtdbjson["nwpw"]["simulation_cell"].is_object())
    {
        const json& simulation_cell =
            rtdbjson["nwpw"]["simulation_cell"];

        if (simulation_cell.contains("unita_frozen") &&
            simulation_cell["unita_frozen"].is_array())
        {
            frozen_unita =
                &simulation_cell["unita_frozen"];
        }
    }

    if (current_unita != nullptr &&
        current_unita->size() == 9)
    {
        coutput
            << tag
            << "Current unita:\n";

        for (int j = 0; j < 3; ++j)
        {
            coutput  << std::setprecision(10)
                << tag
                << "  "
                << (*current_unita)[3*j + 0].get<double>()
                << " "
                << (*current_unita)[3*j + 1].get<double>()
                << " "
                << (*current_unita)[3*j + 2].get<double>()
                << '\n';
        }
    }

    if (frozen_unita != nullptr &&
        frozen_unita->size() == 9)
    {
        coutput
            << tag
            << "Frozen unita_frozen:\n";

        for (int j = 0; j < 3; ++j)
        {
            coutput  << std::setprecision(10)
                << tag
                << "  "
                << (*frozen_unita)[3*j + 0].get<double>()
                << " "
                << (*frozen_unita)[3*j + 1].get<double>()
                << " "
                << (*frozen_unita)[3*j + 2].get<double>()
                << '\n';
        }
    }

    if (rtdbjson.contains("driver") &&
        rtdbjson["driver"].is_object() &&
        rtdbjson["driver"].contains("numerical_grid") &&
        rtdbjson["driver"]["numerical_grid"].is_object())
    {
        const json& grid = rtdbjson["driver"]["numerical_grid"];

        coutput
            << tag
            << "FFT grid = "
            << grid.value("nx", -1)
            << " x "
            << grid.value("ny", -1)
            << " x "
            << grid.value("nz", -1)
            << " waves0 = "
            << grid.value("nwave0", -1)
            << " waves1 = "
            << grid.value("nwave1", -1)
            << " npack0 = "
            << grid.value("npack0", -1)
            << " npack1 = "
            << grid.value("npack1", -1)
            << '\n';
    }
    else
    {
        coutput << tag << "Numerical grid: unavailable\n";
    }

}

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


static bool update_unita_frozen(std::string& rtdbstring, std::ostream& coutput, const double lattice_tolerance, const bool oprint)
{
   if (lattice_tolerance < 0.0)
   {
        if (oprint) coutput << "driver_optimizer: negative lattice tolerance\n";
       return false;
   }

   json rtdbjson;

   try
   {
       rtdbjson = json::parse(rtdbstring);
   }
   catch (const json::exception& ex)
   {
       if (oprint) coutput << "driver_optimizer: invalid RTDB JSON: " << ex.what() << '\n';
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
       if (oprint) coutput << "driver_optimizer: geometry '" << geomname << "' not found\n";
       return false;
   }

   const json& geometry =
       rtdbjson["geometries"][geomname];

   std::array<double, 9> current_unita{};

   if (!read_unita(
           geometry.value("unita", json{}),
           current_unita))
   {
       if (oprint) coutput << "driver_optimizer: invalid current geometry lattice\n";

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

       if (oprint) coutput << "driver_optimizer: initialized unita_frozen\n";
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

           if (oprint) coutput << "driver_optimizer: resetting unita_frozen\n"
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
 *            compute_egs_values          *
 *                                        *
 ******************************************/
/**
 * @brief Compute and collect energy, gradient, and stress data.
 *
 * Prepares an RTDB request for the selected calculation type, invokes the
 * supplied PSPW or band minimizer callback, and extracts the resulting
 * quantities from the backend-specific RTDB section.
 *
 * The calculation selected by @p option is:
 *
 *   - @c 1: energy
 *   - @c 2: gradient
 *   - @c 3: stress
 *
 * For stress calculations, the NWPW `includestress` option is enabled before
 * invoking the minimizer. The callback updates its RTDB string argument with
 * the calculation results. The updated RTDB string is copied back into
 * @p rtdbstring before this function returns.
 *
 * The returned JSON object contains the following fields when provided by
 * the selected backend:
 *
 *   - @c energy
 *   - @c gradient
 *   - @c stress
 *   - @c stress_sym
 *   - @c lstress
 *
 * The backend is selected from the returned RTDB. A `pspw` section is checked
 * first, followed by a `band` section.
 *
 * @param[in] option
 *     Calculation type: 1 for energy, 2 for gradient, or 3 for stress.
 *
 * @param[in] comm_world0
 *     MPI communicator used by the minimizer.
 *
 * @param[in] minimizer
 *     Callback function for the PSPW or band minimizer. The callback must
 *     have the signature:
 *
 *         int(MPI_Comm, std::string&, std::ostream&)
 *
 * @param[in,out] rtdbstring
 *     RTDB JSON string containing the input state. It is replaced with the
 *     updated RTDB returned by the minimizer.
 *
 * @param[in] coutput
 *     Output stream used by the minimizer and diagnostic messages.
 *
 * @return
 *     A JSON object containing the extracted energy, gradient, and stress
 *     results. Missing backend fields are returned as null JSON values.
 *
 * @throws std::invalid_argument
 *     If @p minimizer is null or @p option is not in the range 1--3.
 *
 * @throws nlohmann::json::exception
 *     If the input or returned RTDB string is not valid JSON.
 *
 * @throws std::runtime_error
 *     If the minimizer fails or neither PSPW nor band results are present
 *     in the returned RTDB.
 */
static json compute_egs_values(const int option, 
                               MPI_Comm comm_world0, 
                               minimizer_function minimizer, 
                               std::string& rtdbstring, 
                               std::ostream& coutput)
{
   if (minimizer == nullptr)
      throw std::invalid_argument("compute_egs_values: null minimizer");
   if (option < 1 || option > 3)
      throw std::invalid_argument("compute_egs_values: invalid option");

   // Prepare the input RTDB for this evaluation.
   json request = json::parse(rtdbstring);

   if (option == 1)
   {
      request["current_task"] = "task pspw energy";
      request["nwpw"]["includestress"] = false;
   }
   else if (option == 2)
   {
      request["current_task"] = "task pspw gradient";
      request["nwpw"]["includestress"] = false;
   }
   else if (option == 3)
   {
      request["current_task"] = "task pspw stress";
      request["nwpw"]["includestress"] = true;
   }

   request["driver"]["cell_optimization"] = true;
   request["driver"]["use_frozen_lattice"] = true;

   const std::string tag = "@";

   print_lattice_state(request, coutput, tag);



   // The callback modifies this string by reference.
   std::string rtdbstring1 = request.dump();

   const int ierr = minimizer(comm_world0, rtdbstring1, coutput);

   if (ierr != 0)
      throw std::runtime_error("compute_egs_values: minimizer failed with error " + std::to_string(ierr));


   // Parse the UPDATED RTDB returned by the minimizer.
   const json output = json::parse(rtdbstring1);

   const json* backend = nullptr;

   if (output.contains("pspw") && output["pspw"].is_object())
      backend = &output["pspw"];
   else if (output.contains("band") && output["band"].is_object())
      backend = &output["band"];
   else
      throw std::runtime_error("compute_egs_values: neither pspw nor band results found");

   const json& backend_result = *backend;
   json result = json::object();

   result["energy"]     = backend_result.value("energy", json{});
   result["gradient"]   = backend_result.value("gradient", json{});
   result["stress"]     = backend_result.value("stress", json{});
   result["stress_sym"] = backend_result.value("stress_sym", json{});
   result["lstress"] = backend_result.value("lstress", json{});

   // Preserve the most recent backend RTDB for the caller.
   rtdbstring = std::move(rtdbstring1);

   return result;
}


/******************************************
 *                                        *
 *          scale_cubic_cell              *
 *                                        *
 ******************************************/
static void scale_cubic_cell(std::string& rtdbstring, const double scale)
{
    json rtdbjson = json::parse(rtdbstring);

    const std::string geomname =
        rtdbjson.contains("geometry") &&
        rtdbjson["geometry"].is_string()
            ? rtdbjson["geometry"].get<std::string>()
            : "geometry";

    json& geometry =
        rtdbjson["geometries"][geomname];

    for (int i = 0; i < 9; ++i)
    {
        geometry["unita"][i] =
            geometry["unita"][i].get<double>() *
            scale;
    }

    if (geometry.contains("coords") &&
        geometry["coords"].is_array())
    {
        for (std::size_t i = 0;
             i < geometry["coords"].size();
             ++i)
        {
            geometry["coords"][i] =
                geometry["coords"][i].get<double>() *
                scale;
        }
    }

    /*
     * If a separate current simulation-cell lattice exists,
     * keep it consistent with geometry.unita.
     */
    if (rtdbjson.contains("nwpw") &&
        rtdbjson["nwpw"].is_object() &&
        rtdbjson["nwpw"].contains("simulation_cell") &&
        rtdbjson["nwpw"]["simulation_cell"].is_object())
    {
        json& cell =
            rtdbjson["nwpw"]["simulation_cell"];

        if (cell.contains("unita") &&
            cell["unita"].is_array() &&
            cell["unita"].size() == 9)
        {
            for (int i = 0; i < 9; ++i)
            {
                cell["unita"][i] =
                    cell["unita"][i].get<double>() *
                    scale;
            }
        }
    }

    rtdbstring = rtdbjson.dump();
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

   const bool master = myparallel.is_master();

   // Add unita_frozen in rtdb
   // This driver is the unit-cell/geometry optimization path,
   // so establish or validate unita_frozen before Control2 reads the RTDB.
   //constexpr double lattice_tolerance = 1.0e-2;
   constexpr double lattice_tolerance = 0.05;

   if (!update_unita_frozen(rtdbstring, coutput, lattice_tolerance, master))
      return 1;


   Control2 control(myparallel.np(),rtdbstring);

   bool hprint = master && control.print_level("high");
   bool oprint = master && control.print_level("medium");
   bool lprint = master && control.print_level("low");

   /* reset Parallel base_stdio_print = lprint */
   myparallel.base_stdio_print = lprint;
   std::string tag  = "@@";

   // Fetch the space groupd symmetry if it exists
   SymmetryInfo symmetry_info;
   json parse_json = json::parse(rtdbstring);
   if (parse_json.contains("effective_symmetry"))
   {    
      const json& effective_symmetry = parse_json.at("effective_symmetry");
      symmetry_info.space_group_name = effective_symmetry.value("name", "unknown");
      symmetry_info.type = effective_symmetry.value("type", "unknown");
      symmetry_info.group_order = effective_symmetry.value("order", -1);
      symmetry_info.is_primitive = effective_symmetry.value("primitive", false);
     
      // Simple cubic detection, expand as needed
      symmetry_info.is_cubic = (symmetry_info.space_group_name.find("Fd-3m") != std::string::npos) ||
                               (symmetry_info.group_order == 192);

      int sgnum = symmetry_info.group_order;
      if      (sgnum >= 1   && sgnum <= 2)   symmetry_info.system = "triclinic";
      else if (sgnum >= 3   && sgnum <= 15)  symmetry_info.system = "monoclinic";
      else if (sgnum >= 16  && sgnum <= 74)  symmetry_info.system = "orthorhombic";
      else if (sgnum >= 75  && sgnum <= 142) symmetry_info.system = "tetragonal";
      else if (sgnum >= 143 && sgnum <= 167) symmetry_info.system = "trigonal";
      else if (sgnum >= 168 && sgnum <= 194) symmetry_info.system = "hexagonal";
      else if (sgnum >= 195 && sgnum <= 230) symmetry_info.system = "cubic";
      if (symmetry_info.space_group_name.find("Fd-3m") != std::string::npos)
         symmetry_info.system = "cubic";
   }
   

   if (oprint) 
   {
      std::ios_base::sync_with_stdio();
      constexpr int width = 78;

      coutput << '\n'
       << tag << std::string(width, '=') << '\n'
       << tag << "                         PWDFT OPTIMIZATION DRIVER\n"
       << tag << std::string(width, '=') << '\n'
       << tag << '\n'
       << tag << "  Architecture           : NWChem-style driver dispatch\n"
       << tag << "  Role                   : top-level orchestration layer\n"
       << tag << "  Backend                : PSPW or band minimizer callback\n"
       << tag << "  Implementation         : NorthwestEx C++ driver\n"
       << tag << "  Method                 : Grassmann/Stiefel manifold\n";
      if (symmetry_info.has_symmetry())
      {
         coutput << tag << "  Symmetry information:" << std::endl;
         coutput << tag << "    Space group name:   " << symmetry_info.space_group_name << std::endl;
         coutput << tag << "    Symmetry type:      " << symmetry_info.type << std::endl;
         coutput << tag << "    Group order:        " << symmetry_info.group_order << std::endl;
         coutput << tag << "    Primitive cell:     " << (symmetry_info.is_primitive ? "true" : "false") << std::endl;
         coutput << tag << "    Crystal system:     " << symmetry_info.system <<std::endl;
         // Add more fields if needed
      } 
      else 
      {
         coutput << tag << "  No symmetry information detected." << std::endl;
      }
       
      coutput 
       << tag << "  Cell optimization      : RTDB unita_frozen reference lattice\n"
       << tag << "  Numerical grid         : fixed during optimization stage\n"
       << tag << "  Lattice tolerance      : " << lattice_tolerance << '\n' 
       << tag << "  Date                   : " << util_date() << '\n'
       << tag << '\n'
       << tag << "  The current lattice may change during unit-cell optimization.\n"
       << tag << "  RTDB variable unita_frozen stores the reference lattice used\n"
       << tag << "  to establish the numerical grid and basis-support policy.\n"
       << tag << "  It is reset when the relative lattice change exceeds the\n"
       << tag << "  configured tolerance or when a new optimization stage begins.\n"
       << tag << "  Current energy, force, and stress evaluations use the current\n"
       << tag << "  physical lattice, not unita_frozen.\n"
       << tag << '\n'
       << tag << std::string(width, '-') << '\n';
   }

  
   // Common driver-level work goes here.
   //  - For now, if the driver is only dispatching, call the
   //  - selected minimizer directly.

   // Add unita_frozen in rtdb

   // Relative Frobenius-norm tolerance for resetting unita_frozen.
   //  - For an isotropic lattice scaling, 1.0e-2 corresponds approximately
   //  - to a 1% change in the lattice constant.
    
   if (oprint) coutput << tag <<  "Initial Stress Calcultions" << std::endl;
   json result = compute_egs_values(3,comm_world0,minimizer,rtdbstring, coutput);


   /*
   const json& effective_symmetry = json::parse(rtdbstring).at("effective_symmetry");
   if (oprint) {
      coutput << tag << "Symmetry information:" << std::endl;
      if (effective_symmetry.contains("name"))
          coutput << tag << "  Space group name:  " << effective_symmetry.at("name").get<std::string>() << std::endl;
      if (effective_symmetry.contains("type"))
          coutput << tag << "  Symmetry type:     " << effective_symmetry.at("type").get<std::string>() << std::endl;
      if (effective_symmetry.contains("order"))
          coutput << tag << "  Group order:       " << effective_symmetry.at("order").get<int>() << std::endl;
      if (effective_symmetry.contains("num_centering"))
          coutput << tag << "  Centerings:        " << effective_symmetry.at("num_centering").get<int>() << std::endl;
      if (effective_symmetry.contains("tolerance"))
          coutput << tag << "  Tolerance:         " << effective_symmetry.at("tolerance").get<double>() << std::endl;
      if (effective_symmetry.contains("primitive"))
          coutput << tag << "  Primitive cell:    " << (effective_symmetry.at("primitive").get<bool>() ? "true" : "false") << std::endl;
      if (effective_symmetry.contains("coords_type"))
          coutput << tag << "  Coordinates:       " << effective_symmetry.at("coords_type").get<std::string>() << std::endl;
      if (effective_symmetry.contains("translation_type"))
          coutput << tag << "  Translation type:  " << effective_symmetry.at("translation_type").get<std::string>() << std::endl;
     
      // Print how many symmetry operations (ops)
      if (effective_symmetry.contains("ops"))
          coutput << tag << "  Symmetry operations: " << effective_symmetry.at("ops").size() << std::endl;
     
      // If you want, print first few symmetry operations (rotation/translation)
      if (effective_symmetry.contains("ops")) {
          int nprint = std::min(3, static_cast<int>(effective_symmetry.at("ops").size()));
          coutput << tag << "  First " << nprint << " symmetry operations:" << std::endl;
          for (int i = 0; i < nprint; ++i) {
              const auto& op = effective_symmetry.at("ops").at(i);
              coutput << tag << "    R = [";
              for (int r = 0; r < 3; ++r) {
                  for (int c = 0; c < 3; ++c) {
                      coutput << op.at("R").at(r).at(c).get<double>();
                      if (c < 2) coutput << ", ";
                  }
                  if (r < 2) coutput << " | ";
              }
              coutput << "]  ";
              coutput << "t = [";
              for (int t = 0; t < 3; ++t) {
                  coutput << op.at("t").at(t).get<double>();
                  if (t < 2) coutput << ", ";
              }
              coutput << "]" << std::endl;
          }
      }
     
      // Print fingerprint if present
      if (effective_symmetry.contains("sym_fingerprint"))
          coutput << tag << "  Symmetry fingerprint: " << effective_symmetry.at("sym_fingerprint").get<std::string>() << std::endl;
     
      if (effective_symmetry.contains("source"))
          coutput << tag << "  Symmetry source:      " << effective_symmetry.at("source").get<std::string>() << std::endl;
   }
   */
  

   const double energy = result.at("energy").get<double>();
   const json& lstress = result.at("lstress");

   // can you implement a simple lattice optimization first using lstress

   if (oprint)
   {
      coutput << std::fixed << std::setprecision(10)
              << tag << "Current energy = " << energy << '\n'
              << tag << "dE/da     = " << lstress.at(0).get<double>() << '\n'
              << tag << "dE/db     = " << lstress.at(1).get<double>() << '\n'
              << tag << "dE/dc     = " << lstress.at(2).get<double>() << '\n'
              << tag << "dE/dalpha = " << lstress.at(3).get<double>() << '\n'
              << tag << "dE/dbeta  = " << lstress.at(4).get<double>() << '\n'
              << tag << "dE/dgamma = " << lstress.at(5).get<double>() << '\n';
   }



   // Finite-difference check of the isotropic lattice derivative.
   //  - This evaluates: - dE/ds ≈ [E((1+delta)A) - E((1-delta)A)] / (2*delta)
   //  - where A is the current direct lattice and s is an isotropic scale factor.

   constexpr double fd_delta = 1.0e-2;

   std::string plus_rtdb = rtdbstring;
   std::string minus_rtdb = rtdbstring;

   scale_cubic_cell(plus_rtdb, 1.0 + fd_delta);
   scale_cubic_cell(minus_rtdb, 1.0 - fd_delta);

   json plus_result  = compute_egs_values(1, comm_world0, minimizer, plus_rtdb, coutput);
   json minus_result = compute_egs_values(1, comm_world0, minimizer, minus_rtdb, coutput);

   const double eplus = plus_result.at("energy").get<double>();
   const double eminus = minus_result.at("energy").get<double>();

   const double finite_difference = (eplus - eminus) / (2.0 * fd_delta);

   json current_json = json::parse(rtdbstring);

   const std::string geomname =
       current_json.contains("geometry") &&
       current_json["geometry"].is_string()
           ? current_json["geometry"].get<std::string>()
           : "geometry";


   const json& current_unita = current_json["geometries"][geomname]["unita"];

   const double a = current_unita.at(0).get<double>();
   const double dE_da = lstress.at(0).get<double>();
   const double dE_db = lstress.at(1).get<double>();
   const double dE_dc = lstress.at(2).get<double>();

   const double analytic_dE_dscale = a * (dE_da + dE_db + dE_dc);

   if (oprint)
   {
      coutput << std::setprecision(12)
              << tag
              << "Finite-difference check:\n"
              << tag
              << "  E(+delta)          = "
              << eplus
              << '\n'
              << tag
              << "  E(-delta)          = "
              << eminus
              << '\n'
              << tag
              << "  dE/dscale FD       = "
              << finite_difference
              << '\n'
              << tag
              << "  dE/dscale lstress  = "
              << analytic_dE_dscale
              << '\n'
              << tag 
              << "  dE_da = "
              << dE_da
              << '\n';
      coutput << tag << "\n";
      coutput << tag << "\n";
   }


   // Start optimization Here!!!!
   int lstep = 0;
   double lenergy = 0.0;

   double step    = 0.0025;
  
   constexpr double minimum_step = 1.0e-5;
   constexpr double minimum_gradient = 1.0e-4;
   bool converged = false;

   //constexpr double minimum_step = 1.0e-5; // control stuff
   constexpr int max_steps = 25;

   for (int istep=0; istep<max_steps; ++istep)
   {
      json current_result = compute_egs_values(3, comm_world0, minimizer, rtdbstring, coutput);

      const double current_energy = current_result.at("energy").get<double>();
      const json& lstress = current_result.at("lstress");

    const double dE_dcell =
        (
            lstress.at(0).get<double>() +
            lstress.at(1).get<double>() +
            lstress.at(2).get<double>()
        ) / 3.0;

    if (oprint)
    {
        coutput << std::defaultfloat << std::setprecision(10)
            << tag
            << "Cell step "
            << istep
            << " current energy = "
            << current_energy
            << " dE/dcell = "
            << dE_dcell
            << " step = "
            << step
            << '\n';
    }

   std::string expanded_rtdb = rtdbstring;
   std::string contracted_rtdb = rtdbstring;

   scale_cubic_cell(expanded_rtdb,   1.0 + step);
   scale_cubic_cell(contracted_rtdb, 1.0 - step);

   json expanded_result = compute_egs_values(1,comm_world0, minimizer, expanded_rtdb, coutput);
   json contracted_result = compute_egs_values(1, comm_world0, minimizer, contracted_rtdb, coutput);

   const double expanded_energy = expanded_result.at("energy").get<double>();
   const double contracted_energy = contracted_result.at("energy").get<double>();

   if (oprint) {
      coutput << "\n"
              << tag << "----------------------------------------------\n"
              << tag << " PWDFT Lattice Optimization                   \n"
              << tag << "----------------------------------------------\n"
              << tag << " Step        : " << istep << '\n'
              << tag << " Energy      : " << std::fixed << std::setprecision(10) << current_energy << " Hartree\n"
              << tag << " Lattice a   : " << std::fixed << std::setprecision(6) << a << " Bohr"
              << " (" << std::fixed << std::setprecision(3) << a * 0.529177 << " Å)\n"
              << tag << " dE/da       : " << lstress.at(0).get<double>() << '\n'
              << tag << " dE/db       : " << lstress.at(1).get<double>() << '\n'
              << tag << " dE/dc       : " << lstress.at(2).get<double>() << '\n'
              << tag << " Step size   : " << step << '\n';
     
      // Show action taken
      if (expanded_energy < current_energy && expanded_energy <= contracted_energy) {
          coutput << tag << " Action      : Expanded lattice, accepted.\n";
      } else if (contracted_energy < current_energy) {
          coutput << tag << " Action      : Contracted lattice, accepted.\n";
      } else {
          coutput << tag << " Action      : No improvement, step rejected (minimal lattice change).\n";
      }
      coutput << tag << "----------------------------------------------\n";
   }


   double grad_norm = std::sqrt(std::pow(lstress.at(0).get<double>(), 2) +
                                std::pow(lstress.at(1).get<double>(), 2) +
                                std::pow(lstress.at(2).get<double>(), 2));

   // Convergence check
   if ((step < minimum_step) && (grad_norm < minimum_gradient))
   {
      converged = true;
      break;
   }

if (expanded_energy < current_energy && expanded_energy <= contracted_energy)
{
    rtdbstring = std::move(expanded_rtdb);

    if (oprint)
    {
        coutput  << std::defaultfloat << std::setprecision(10)
            << tag
            << "Accepted expansion, step = " 
            << istep
            << ", energy = "
            << expanded_energy
            << '\n';
    }
}
else if (contracted_energy < current_energy)
{
    rtdbstring =
        std::move(contracted_rtdb);

    if (oprint)
    {
        coutput  << std::defaultfloat << std::setprecision(10)
            << tag
            << "Accepted contraction, step = "  
            << istep
            << ", energy = "
            << contracted_energy
            << '\n';
       coutput << std::defaultfloat << std::setprecision(10)
            << tag
            << "curent_unita = "  
            << current_unita << '\n';
    }
}
else
{
    step *= 0.5;

    if (oprint)
    {
        coutput  << std::defaultfloat << std::setprecision(10)
            << tag
            << "Rejected both directions, step = "
            << istep
            << '\n';
    }

    if (step < minimum_step)
    {
       break;
    }
}

   lstep = istep;
   lenergy = current_energy;
}


// After the optimization loop (use the final geometry and its energy/stress)
json final_result = compute_egs_values(3, comm_world0, minimizer, rtdbstring, coutput);
// Use 'rtdbstring' as possibly updated in the last step

const double final_energy = final_result.at("energy").get<double>();
const json& final_lstress = final_result.at("lstress");

// Extract the geometry from rtdbstring itself (not initial current_unita!)
json final_json = json::parse(rtdbstring);
const json& final_unita = final_json["geometries"][geomname]["unita"];
const double final_a = final_unita.at(0).get<double>();



// After the optimization loop:
if (oprint)
{
    coutput << "\n";
    coutput << tag << "==============================================\n";
    coutput << tag << " PWDFT Lattice Optimization COMPLETE\n";
    coutput << tag << "==============================================\n";

    coutput << tag << " Final lattice parameter (a): "
            << std::fixed << std::setprecision(6) << final_a
            << " Bohr = "
            << std::fixed << std::setprecision(3) << final_a * 0.529177
            << " Å\n";
    coutput << tag << " Minimum energy (total): "
            << std::fixed << std::setprecision(8) << final_energy
            << " Hartree\n";
    coutput << tag << " Gradients at minimum: "
            << "dE/da = " << std::fixed << std::setprecision(5) << final_lstress.at(0).get<double>()
            << ", dE/db = " << final_lstress.at(1).get<double>()
            << ", dE/dc = " << final_lstress.at(2).get<double>() << '\n';
    coutput << tag << " Optimization steps taken: " << lstep+1 << '\n';

    std::string status = (step < minimum_step) ? "Converged" : "Stopped (max steps reached)";
    coutput << tag << " Status: " << status << '\n';
    coutput << tag << "==============================================\n";
}


   return 0;

}


} // namespace pwdft 


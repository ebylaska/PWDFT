
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

//#include "lattice_minimizer.hpp"
#include "lattice_common.hpp"


//#include "gdevice.hpp"




#include "json.hpp"
using json = nlohmann::json;


namespace pwdft {


using minimizer_function = electronic_minimizer;



enum class RelaxCombinedTask {
    // --- 1. Classical Ground State Structural Modes ---
    GeometryOnly  = 0, // Frozen cell box; adjust inner coordinates
    LatticeOnly   = 1, // Frozen atom coordinates; scale bounding box
    Both          = 2, // Co-optimize atoms and bounding box in tandem

};



/*******************************************
 *                                         *
 *          update_unita_frozen            *
 *                                         *
 *******************************************/
/**
 * Initialize or refresh the frozen reference lattice used by NWPW.
 *
 * The routine parses the serialized RTDB, locates the active geometry, and
 * determines the current physical lattice. The geometry lattice is used by
 * default; a valid `nwpw.simulation_cell.unita` takes precedence when present.
 *
 * The current lattice is stored in
 *
 *     nwpw.simulation_cell.unita_frozen
 *
 * if no valid frozen lattice exists. Also, if a frozen lattice is already present,
 * it is replaced only when its relative difference from the current lattice
 * is strictly greater than `lattice_tolerance`. The modified RTDB is then
 * serialized back into `rtdbstring`, including cases where no reset was
 * required.
 *
 * Diagnostic messages are written to `coutput` only when `oprint` is true.
 *
 * @param[in,out] rtdbstring
 *     Serialized RTDB JSON. On successful return, it contains the initialized
 *     or checked `unita_frozen` reference lattice.
 *
 * @param[out] coutput
 *     Stream used for validation errors and lattice-update diagnostics.
 *
 * @param[in] lattice_tolerance
 *     Maximum permitted relative change between the current and frozen
 *     lattices. A negative value is invalid.
 *
 * @param[in] oprint
 *     Enables diagnostic output when true.
 *
 * @return
 *     `true` if the RTDB was parsed and processed successfully; `false` if
 *     the tolerance is invalid, the JSON cannot be parsed, the active
 *     geometry is missing, or the current geometry lattice is invalid.
 *
 * @note
 *     A return value of `true` does not necessarily mean that
 *     `unita_frozen` changed; it also indicates a successful check for which
 *     the existing frozen lattice remained within tolerance.
 */
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

   const json& geometry = rtdbjson["geometries"][geomname];

   std::array<double, 9> current_unita{};

   if (!read_unita( geometry.value("unita", json{}), current_unita))
   {
      if (oprint) coutput << "driver_optimizer: invalid current geometry lattice\n";
      return false;
   }

   // If the NWPW simulation-cell lattice exists, use it as the
   //    current physical lattice. This preserves the intended precedence.
   if (rtdbjson.contains("nwpw") &&
       rtdbjson["nwpw"].is_object() &&
       rtdbjson["nwpw"].contains("simulation_cell") &&
       rtdbjson["nwpw"]["simulation_cell"].is_object())
   {
      const json& simulation_cell = rtdbjson["nwpw"]["simulation_cell"];

      std::array<double, 9> simulation_unita{};

      if (read_unita( simulation_cell.value("unita", json{}), simulation_unita))
      {
         current_unita = simulation_unita;
      }
   }

   json& simulation_cell = rtdbjson["nwpw"]["simulation_cell"];

   if (!simulation_cell.is_object())
      simulation_cell = json::object();

   std::array<double, 9> frozen_unita{};

   const bool has_frozen_unita = read_unita(simulation_cell.value( "unita_frozen", json{}), frozen_unita);

   if (!has_frozen_unita)
   {
      write_unita(simulation_cell["unita_frozen"], current_unita);
      if (oprint) coutput << "driver_optimizer: initialized unita_frozen\n";
   }
   else
   {
      const double difference = unita_relative_difference(current_unita, frozen_unita);

      if (difference > lattice_tolerance)
      {
         write_unita(simulation_cell["unita_frozen"], current_unita);

         if (oprint) 
            coutput << "driver_optimizer: resetting unita_frozen\n"
                    << "  relative lattice change = "
                    << difference
                    << "\n"
                    << "  tolerance               = "
                    << lattice_tolerance
                    << '\n';
      }
   }

   // Critical: pass the modified RTDB back to the caller/minimizer.
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

   const bool master = myparallel.is_master();

   // Add unita_frozen in rtdb
   //  - This driver is the unit-cell/geometry optimization path,
   //  - so establish or validate unita_frozen before Control2 reads the RTDB.
   constexpr double lattice_tolerance = 0.05;

   if (!update_unita_frozen(rtdbstring, coutput, lattice_tolerance, master))
      return 1;


   Control2 control(myparallel.np(),rtdbstring);

   bool hprint = master && control.print_level("high");
   bool oprint = master && control.print_level("medium");
   bool lprint = master && control.print_level("low");

   // Initialize relax_type to a safe default fallback
   RelaxCombinedTask driver_relax_type = static_cast<RelaxCombinedTask>(control.driver_relax_type());



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
      symmetry_info.ita_number = effective_symmetry.value("ita_number", -1);

      // Read symmetry operations if present.
      if (effective_symmetry.contains("ops") && effective_symmetry["ops"].is_array())
      {
         symmetry_info.ops.reserve(effective_symmetry["ops"].size());
         for (const auto& op : effective_symmetry["ops"])
         {
            FracSymOp s{};
        
            const auto& R = op.at("R");
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    s.R[3*i + j] = R.at(i).at(j).get<double>();
        
            const auto& t = op.at("t");
            for (int i = 0; i < 3; ++i)
                s.t[i] = t.at(i).get<double>();
        
            symmetry_info.ops.push_back(s);
         }
      }
     
      int sgnum = symmetry_info.ita_number;
      if      (sgnum >= 1   && sgnum <= 2)   symmetry_info.system = "triclinic";
      else if (sgnum >= 3   && sgnum <= 15)  symmetry_info.system = "monoclinic";
      else if (sgnum >= 16  && sgnum <= 74)  symmetry_info.system = "orthorhombic";
      else if (sgnum >= 75  && sgnum <= 142) symmetry_info.system = "tetragonal";
      else if (sgnum >= 143 && sgnum <= 167) symmetry_info.system = "trigonal";
      else if (sgnum >= 168 && sgnum <= 194) symmetry_info.system = "hexagonal";
      else if (sgnum >= 195 && sgnum <= 230) symmetry_info.system = "cubic";
      else  symmetry_info.system = "unknown";
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
              << tag << "  Architecture                   : NWChem-style driver dispatch\n"
              << tag << "  Role                           : top-level orchestration layer\n"
              << tag << "  Backend                        : PSPW or band minimizer callback\n"
              << tag << "  Implementation                 : NorthwestEx C++ driver\n"
              << tag << "  Method                         : Grassmann/Stiefel manifold\n" 
              << tag << std::endl;
      if (symmetry_info.has_symmetry())
      {
         coutput << tag << "  Symmetry information:            " << std::endl;
         coutput << tag << "      Space group name           : " << symmetry_info.space_group_name << std::endl;
         coutput << tag << "      Space group number         : " << symmetry_info.ita_number << " (ITC)" <<  std::endl;
         coutput << tag << "      Symmetry type              : " << symmetry_info.type << std::endl;
         coutput << tag << "      Group order                : " << symmetry_info.group_order << std::endl;
         coutput << tag << "      Primitive cell             : " << (symmetry_info.is_primitive ? "true" : "false") << std::endl;
         coutput << tag << "      Crystal system             : " << symmetry_info.system <<std::endl;
         coutput << tag << std::endl;
         // Add more fields if needed
      } 
      else 
      {
         coutput << tag << "  No symmetry information detected." << std::endl;
      }
       
      //coutput << tag << "  Cell optimization      : RTDB unita_frozen reference lattice\n";
      if (driver_relax_type == RelaxCombinedTask::GeometryOnly) coutput << tag << "  Geometry only optimization     : RTDB unita_frozen reference lattice\n";
      if (driver_relax_type == RelaxCombinedTask::LatticeOnly)  coutput << tag << "  Lattice only optimization      : RTDB unita_frozen reference lattice\n";
      if (driver_relax_type == RelaxCombinedTask::Both)         coutput << tag << "  Geometry and Cell optimization : RTDB unita_frozen reference lattice\n";

      coutput << tag << "  Numerical grid                 : fixed during optimization stage\n"
              << tag << "  Lattice tolerance              : " << lattice_tolerance << '\n' 
              << tag << "  Date                           : " << util_date() << '\n'
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

  
   /*
   // Common driver-level work goes here.
   //  - For now, if the driver is only dispatching, call the
   //  - selected minimizer directly.
   // Add unita_frozen in rtdb
   // Relative Frobenius-norm tolerance for resetting unita_frozen.
   //  - For an isotropic lattice scaling, 1.0e-2 corresponds approximately
   //  - to a 1% change in the lattice constant.
   if (oprint) coutput << tag <<  "Initial Stress Calculations" << std::endl;
   json result         = compute_egs_values(3,comm_world0,minimizer,rtdbstring, coutput);
   const double energy = result.at("energy").get<double>();
   const json& lstress = result.at("lstress");

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
   */


   // ---------------------------------------------------------------
   // Dispatch to the crystal-system-specific lattice minimizer.
   //
   // pick_lattice_minimizer never returns nullptr: unknown systems fall
   // back to general_lattice_minimizer, which currently evaluates once
   // and returns without modifying the lattice.
   // ---------------------------------------------------------------
   lattice_minimizer lm = pick_lattice_minimizer(symmetry_info.system);

   // define ctx lattice optimization controls
   LatticeContext ctx;
   ctx.oprint           = oprint;
   ctx.max_steps        = control.driver_lattice_maxiter();
   ctx.minimum_gradient = control.driver_lattice_gmax();
   ctx.initial_step     = control.driver_lattice_step();
   ctx.minimum_step     = control.driver_lattice_xmin();

   const int ierr = lm(comm_world0, rtdbstring, coutput, minimizer, ctx);
   if (ierr != 0)
   {
      coutput << tag << "lattice minimizer returned " << ierr << '\n';
      return ierr;
   }

   return 0;

}


} // namespace pwdft 


// cubic_lattice_minimizer.cpp
//
// Safeguarded one-dimensional quasi-Newton optimizer for cubic cells
// (a = b = c, all angles 90 degrees).
//
// The optimization coordinate is
//
//     x = log(a)
//
// so a trial displacement dx scales the complete cubic cell by exp(dx).
// This guarantees that the lattice parameter remains positive.
//
// Assuming lstress contains the lattice derivatives
//
//     dE/da, dE/db, dE/dc,
//
// the derivative along the constrained cubic path is
//
//     dE/da|cubic = dE/da + dE/db + dE/dc
//
// and the derivative with respect to x = log(a) is
//
//     dE/dx = a * dE/da|cubic.
//
// After two accepted points, the code uses a secant approximation to the
// one-dimensional Newton step. The proposed step is restricted by a trust
// radius and safeguarded by energy backtracking.
//
// IMPORTANT: If lstress contains physical stress rather than direct lattice
// derivatives, cubic_gradient() must be changed to include the appropriate
// sign, volume, and coordinate-transformation factors.

#include "lattice_minimizer.hpp"
#include "lattice_common.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <string>
#include <utility>

#include <mpi.h>

#include "json.hpp"

namespace pwdft {

using json = nlohmann::json;

namespace {

constexpr double bohr_to_angstrom = 0.529177;

/*
 * Return the active geometry name.
 */
static std::string active_geometry_name(const json& rtdb)
{
    if (rtdb.contains("geometry") && rtdb["geometry"].is_string())
        return rtdb["geometry"].get<std::string>();

    return "geometry";
}

/*
 * Extract the cubic lattice parameter from geometry.unita.
 *
 * This follows the convention already used by the original minimizer:
 * the first element of unita is the cubic lattice parameter a.
 */
static double cubic_lattice_parameter(const std::string& rtdbstring)
{
    const json rtdb = json::parse(rtdbstring);
    const std::string geomname = active_geometry_name(rtdb);

    return rtdb.at("geometries")
        .at(geomname)
        .at("unita")
        .at(0)
        .get<double>();
}

/*
 * Derivative along the constrained path a = b = c.
 *
 * This is correct when lstress[0], lstress[1], and lstress[2] are,
 * respectively, dE/da, dE/db, and dE/dc.
 */
static double cubic_gradient(const json& lstress)
{
    return lstress.at(0).get<double>()
         + lstress.at(1).get<double>()
         + lstress.at(2).get<double>();
}

} // namespace



/*****************************************
 *                                       *
 *        cubic_lattice_minimizer        *
 *                                       *
 *****************************************/
/**
 * Optimize an isotropic cubic unit cell using a safeguarded 1D secant method.
 *
 * The cubic constraint keeps a = b = c with fixed cell angles, reducing the
 * lattice optimization to one degree of freedom.  The optimization coordinate
 * is x = log(a), so a displacement dx scales the cell by exp(dx) and preserves
 * a positive lattice parameter.
 *
 * Energy and lattice derivatives are obtained through `compute_egs_values()`.
 * Assuming `lstress[0..2]` contain dE/da, dE/db, and dE/dc, respectively, the
 * derivative along the constrained cubic path is
 *
 *     dE/da|cubic = dE/da + dE/db + dE/dc,
 *
 * and the derivative in the logarithmic coordinate is
 *
 *     dE/dx = a * dE/da|cubic.
 *
 * After two accepted points, a secant approximation to the 1D Newton step is
 * used.  The proposed step is limited by a trust radius and must be a descent
 * direction.  Trial cells are accepted only when they lower the energy;
 * otherwise, the step is backtracked.  If the secant estimate is unavailable
 * or unreliable, the routine falls back to a bounded step opposite the
 * gradient.
 *
 * The accepted lattice is written back to `rtdbstring`.  Cell scaling is
 * performed by `scale_cubic_cell()`, which updates the physical lattice and
 * associated coordinates.  Electronic energy and stress evaluations are
 * delegated to the supplied `minimizer` callback.
 *
 * Convergence is determined from the magnitude of the constrained cubic
 * derivative and the tolerances in `ctx`.  Diagnostic output is controlled by
 * `ctx.oprint` and prefixed with `ctx.tag`.
 *
 * @param comm
 *     MPI communicator used by the electronic-structure calculations.
 *
 * @param rtdbstring
 *     Serialized RTDB JSON state.  On entry, it contains the initial cell; on
 *     return, it contains the final accepted cubic cell.
 *
 * @param coutput
 *     Output stream used for optimization diagnostics and electronic
 *     minimizer output.
 *
 * @param minimizer
 *     Electronic minimizer callback used to evaluate energy and lattice
 *     derivatives for each cell.
 *
 * @param ctx
 *     Lattice-optimization settings, including the initial trust radius,
 *     minimum step, gradient tolerance, maximum number of steps, and output
 *     controls.
 *
 * @return
 *     Zero on completion.  Convergence status is currently reported through
 *     the diagnostic output rather than through the return value.
 *
 * @note
 *     The derivative transformation above is valid only if `lstress` contains
 *     direct lattice derivatives.  If it instead contains physical stress
 *     tensor components, the appropriate sign, volume, and chain-rule factors
 *     must be applied.
 */
int cubic_lattice_minimizer(MPI_Comm comm,
                            std::string& rtdbstring,
                            std::ostream& coutput,
                            electronic_minimizer minimizer,
                            const LatticeContext& ctx)
{
   // Retain the original minimizer's local diagnostic behavior.
   const bool oprint      = ctx.oprint;
   const std::string& tag = ctx.tag;

   const int max_steps           = ctx.max_steps;
   const double minimum_step     = ctx.minimum_step;
   const double minimum_gradient = ctx.minimum_gradient;

   // Here "step" is a trust radius in x = log(a). For small values it is
   //    approximately the maximum relative change in the lattice parameter.
   double step = std::abs(ctx.initial_step);

   if (!std::isfinite(step) || step <= 0.0)
       step = 1.0e-2;

   bool converged = false;
   int steps_taken = 0;

   // Previous accepted point used by the secant curvature estimate.
   bool have_previous_point = false;
   double previous_x = 0.0;
   double previous_gradient_x = 0.0;

   if (oprint)
   {
       coutput << '\n'
               << tag << "==============================================\n"
               << tag << " PWDFT cubic lattice optimization\n"
               << tag << " Safeguarded secant method in log(a)\n"
               << tag << "==============================================\n";
   }

   for (int istep=0; istep<max_steps; ++istep)
   {
      const double a = cubic_lattice_parameter(rtdbstring);

      if (!std::isfinite(a) || a <= 0.0)
      {
         coutput << tag << "ERROR: invalid cubic lattice parameter a = " << a << '\n';
         break;
      }

      const double x = std::log(a);

      // Evaluate both energy and lattice derivatives at the current accepted cell.
      const json current_result   = compute_egs_values(3, comm, minimizer, rtdbstring, coutput);
      const double current_energy = current_result.at("energy").get<double>();
      const json& lstress = current_result.at("lstress");
      const double dE_da  = lstress.at(0).get<double>();
      const double dE_db  = lstress.at(1).get<double>();
      const double dE_dc  = lstress.at(2).get<double>();

      // Gradient along a = b = c and then with respect to x = log(a).
      const double gradient_a = cubic_gradient(lstress);
      const double gradient_x = a * gradient_a;

      if (oprint)
      {
         coutput << '\n'
                 << tag << "----------------------------------------------\n"
                 << tag << " Step        : " << istep << '\n'
                 << tag << " Energy      : "
                        << std::fixed << std::setprecision(10)
                        << current_energy << " Hartree\n"
                 << tag << " Lattice a   : "
                        << std::fixed << std::setprecision(6)
                        << a << " Bohr ("
                        << std::fixed << std::setprecision(3)
                        << a * bohr_to_angstrom << " A)\n"
                 << tag << " dE/da       : "
                        << std::defaultfloat << std::setprecision(10)
                        << dE_da << '\n'
                 << tag << " dE/db       : " << dE_db << '\n'
                 << tag << " dE/dc       : " << dE_dc << '\n'
                 << tag << " Cubic dE/da : " << gradient_a << '\n'
                 << tag << " dE/dlog(a)  : " << gradient_x << '\n'
                 << tag << " Trust radius: " << step << '\n';
      }

      if (!std::isfinite(current_energy) || !std::isfinite(gradient_a) || !std::isfinite(gradient_x))
      {
         if (oprint)
         {
            coutput << tag << " Action      : non-finite energy or gradient; stop.\n"
                    << tag << "----------------------------------------------\n";
         }
         break;
      }

      // Use the derivative with respect to the physical cubic lattice
      //    parameter for convergence, so minimum_gradient retains the same
      //    general units as the reported lattice derivatives.
      if (std::abs(gradient_a) < minimum_gradient)
      {
         converged = true;

         if (oprint)
         {
             coutput << tag << " Action      : gradient converged.\n"
                     << tag << "----------------------------------------------\n";
         }
         break;
      }

      double dx = 0.0;
      bool used_secant = false;

      // In one dimension, this is the quasi-Newton/secant step:
      //    dx = -g(x) (x - x_previous) / (g(x) - g_previous).
      if (have_previous_point)
      {
         const double delta_x               = x - previous_x;
         const double delta_gradient        = gradient_x - previous_gradient_x;
         const double gradient_scale        = std::max({1.0, std::abs(gradient_x), std::abs(previous_gradient_x)});
         const double denominator_tolerance = 100.0 * std::numeric_limits<double>::epsilon() * gradient_scale;

         if (std::isfinite(delta_x) &&
             std::isfinite(delta_gradient) &&
             std::abs(delta_x) > 0.0 &&
             std::abs(delta_gradient) > denominator_tolerance)
         {
            dx          = -gradient_x * delta_x / delta_gradient;
            used_secant = std::isfinite(dx);
         }
      }

      // On the first iteration, or if the secant denominator is too small,
      //    take a trust-radius step opposite to the gradient.
      if (!used_secant)
         dx = -std::copysign(step, gradient_x);

      // Restrict the proposed displacement to the current trust radius.
      dx = std::clamp(dx, -step, step);

      // Reject a non-finite or non-descent secant direction before doing
      //    any energy evaluations.
      if (!std::isfinite(dx) || gradient_x * dx >= 0.0)
      {
         dx = -std::copysign(step, gradient_x);
         used_secant = false;
      }

      // Avoid a zero displacement caused by numerical cancellation.
      if (std::abs(dx) < std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(x)))
      {
         dx = -std::copysign(step, gradient_x);
         used_secant = false;
      }

      std::string accepted_rtdb;
      bool accepted          = false;
      double accepted_dx     = 0.0;
      double accepted_energy = current_energy;

      double trial_dx        = dx;

      // Safeguarded backtracking. Only energy is requested at trial cells;
      //    stress is evaluated at the next accepted outer iteration.
      constexpr int maximum_backtracks = 12;

      for (int iback = 0; iback < maximum_backtracks; ++iback)
      {
         if (std::abs(trial_dx) < minimum_step)
            break;

         const double scale_factor = std::exp(trial_dx);

         if (!std::isfinite(scale_factor) || scale_factor <= 0.0)
         {
             trial_dx *= 0.5;
             continue;
         }

         std::string trial_rtdb = rtdbstring;
         scale_cubic_cell(trial_rtdb, scale_factor);

         const json trial_result   = compute_egs_values(1, comm, minimizer, trial_rtdb, coutput);
         const double trial_energy = trial_result.at("energy").get<double>();

         // Monotonic acceptance is deliberately used instead of a strict
         // Armijo condition because electronic minimization noise can make
         // an Armijo test unnecessarily restrictive.
         if (std::isfinite(trial_energy) && trial_energy < current_energy)
         {
            accepted        = true;
            accepted_dx     = trial_dx;
            accepted_energy = trial_energy;
            accepted_rtdb   = std::move(trial_rtdb);
            break;
         }

         trial_dx *= 0.5;
      }

      if (accepted)
      {
         // Preserve the old accepted point. At the next outer iteration,
         //    rtdbstring will represent the new point, giving the two points
         //    needed for a secant curvature estimate.
         previous_x          = x;
         previous_gradient_x = gradient_x;
         have_previous_point = true;

         rtdbstring = std::move(accepted_rtdb);
         ++steps_taken;

         // If backtracking substantially reduced the proposal, tighten
         //    the trust radius. Otherwise allow it to grow moderately.
         const double accepted_magnitude = std::abs(accepted_dx);

         if (accepted_magnitude < 0.5*step)
         {
            step = std::max(minimum_step, 2.0*accepted_magnitude);
         }
         else
         {
            const double maximum_step = std::max(std::abs(ctx.initial_step), minimum_step);
            step = std::min(1.5 * step, maximum_step);
         }

         if (oprint)
         {
            coutput << tag << " Method      : "
                           << (used_secant
                                   ? "secant/quasi-Newton"
                                   : "gradient fallback")
                           << '\n'
                    << tag << " Action      : accepted\n"
                    << tag << " Delta log(a): "
                           << std::defaultfloat
                           << std::setprecision(10)
                           << accepted_dx << '\n'
                    << tag << " Scale factor: "
                           << std::exp(accepted_dx) << '\n'
                    << tag << " New energy  : "
                           << std::fixed
                           << std::setprecision(10)
                           << accepted_energy << " Hartree\n"
                    << tag
                           << "----------------------------------------------\n";
         }
      }
      else
      {
         // The local secant curvature estimate was not useful. Reset its
         //    history and reduce the trust radius.
         have_previous_point = false;
         step *= 0.5;

         if (oprint)
         {
            coutput << tag << " Method      : "
                           << (used_secant
                                   ? "secant/quasi-Newton"
                                   : "gradient fallback")
                           << '\n'
                    << tag
                           << " Action      : rejected; halve trust radius.\n"
                    << tag << " New radius  : "
                           << std::defaultfloat
                           << std::setprecision(10)
                           << step << '\n'
                    << tag
                           << "----------------------------------------------\n";
         }

         if (step < minimum_step)
            break;
      }
   }

   // Recompute energy and lattice derivatives at the final accepted cell.
   const json final_result   = compute_egs_values(3, comm, minimizer, rtdbstring, coutput);
   const double final_energy = final_result.at("energy").get<double>();
   const json& final_lstress = final_result.at("lstress");
   const double final_dE_da = final_lstress.at(0).get<double>();
   const double final_dE_db = final_lstress.at(1).get<double>();
   const double final_dE_dc = final_lstress.at(2).get<double>();

   const double final_gradient_a = final_dE_da + final_dE_db + final_dE_dc;
   const double final_a          = cubic_lattice_parameter(rtdbstring);

   // A final evaluation may satisfy the convergence criterion even if that
   //   point was accepted on the last permitted iteration.
   if (std::isfinite(final_gradient_a) && std::abs(final_gradient_a) < minimum_gradient)
   {
      converged = true;
   }

   if (oprint)
   {
      coutput << '\n'
              << tag << "==============================================\n"
              << tag << " PWDFT cubic lattice optimization COMPLETE\n"
              << tag << "==============================================\n"
              << tag << " Final lattice parameter (a): "
                     << std::fixed << std::setprecision(6)
                     << final_a << " Bohr = "
                     << std::fixed << std::setprecision(3)
                     << final_a * bohr_to_angstrom << " A\n"
              << tag << " Final accepted energy (total): "
                     << std::fixed << std::setprecision(8)
                     << final_energy << " Hartree\n"
              << tag << " Gradients at minimum: "
                     << "dE/da = "
                     << std::defaultfloat << std::setprecision(10)
                     << final_dE_da
                     << ", dE/db = " << final_dE_db
                     << ", dE/dc = " << final_dE_dc << '\n'
              << tag << " Constrained cubic gradient: "
                     << final_gradient_a << '\n'
              << tag << " Accepted lattice steps: "
                     << steps_taken << '\n'
              << tag << " Status: "
                     << (converged
                             ? "Converged"
                             : "Stopped before gradient convergence")
                     << '\n'
              << tag << "==============================================\n";
   }

   return 0;
}

} // namespace pwdft

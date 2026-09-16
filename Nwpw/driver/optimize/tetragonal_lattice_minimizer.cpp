// tetragonal_lattice_minimizer.cpp
//
// Two-dimensional quasi-Newton (BFGS) optimizer for tetragonal cells
// (a = b, all angles 90 degrees; c independent).
//
// Free parameters: a and c. Optimization coordinates:
//
//     x = (log(a), log(c))
//
// so trial displacements scale the cell anisotropically as
//
//     a_new = a * exp(dx_a),   c_new = c * exp(dx_c)
//
// which keeps both lattice parameters positive.
//
// The energy derivatives along the tetragonal constraint are
//
//     dE/da|tet = lstress[0] + lstress[1]     (a and b move together)
//     dE/dc|tet = lstress[2]
//
// and the log-coordinate gradient is
//
//     g = (a * dE/da|tet, c * dE/dc|tet)
//
// A 2x2 inverse-Hessian approximation B is maintained by BFGS and updated
// from each accepted secant pair (dx, dg). The proposed Newton step is
//
//     dx = -B * g
//
// clipped to a trust radius and safeguarded by energy backtracking.

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

// ---------------------------------------------------------------------------
// 2x2 symmetric inverse-Hessian approximation
// ---------------------------------------------------------------------------

struct Hessian2 {
    double B00 = 1.0;
    double B01 = 0.0;
    double B11 = 1.0;

    void reset_to_identity()
    {
        B00 = 1.0;
        B01 = 0.0;
        B11 = 1.0;
    }

    std::pair<double, double> apply(double g0, double g1) const
    {
        return { B00 * g0 + B01 * g1,
                 B01 * g0 + B11 * g1 };
    }

    // BFGS update: B <- (I - rho*s*y^T) B (I - rho*y*s^T) + rho*s*s^T
    // Returns false if the curvature condition y.s > 0 fails.
    bool bfgs_update(double s0, double s1, double y0, double y1)
    {
        const double ys = y0 * s0 + y1 * s1;

        if (!(ys > 1.0e-12 * (std::abs(y0) + std::abs(y1)) *
                          (std::abs(s0) + std::abs(s1)) + 1.0e-16))
            return false;

        const double rho = 1.0 / ys;

        // Compute V = I - rho * y * s^T  (2x2)
        const double V00 = 1.0 - rho * y0 * s0;
        const double V01 =       -rho * y0 * s1;
        const double V10 =       -rho * y1 * s0;
        const double V11 = 1.0 - rho * y1 * s1;

        // t = V * B  (V is not symmetric in general)
        const double t00 = V00 * B00 + V01 * B01;
        const double t01 = V00 * B01 + V01 * B11;
        const double t10 = V10 * B00 + V11 * B01;
        const double t11 = V10 * B01 + V11 * B11;

        // B_new = t * V^T + rho * s * s^T
        const double u00 = t00 * V00 + t01 * V01 + rho * s0 * s0;
        const double u01 = t00 * V10 + t01 * V11 + rho * s0 * s1;
        const double u11 = t10 * V10 + t11 * V11 + rho * s1 * s1;

        B00 = u00;
        B01 = 0.5 * (u01 + u01);   // symmetrize (t*V^T is already symmetric
        B11 = u11;                 //  for symmetric B and our construction)
        B01 = u01;

        return true;
    }
};

// ---------------------------------------------------------------------------
// Small helpers
// ---------------------------------------------------------------------------

static double active_cubic_or_tetragonal_gradient_a(const json& lstress)
{
    // For both cubic and tetragonal, the a/b direction moves both a and b.
    return lstress.at(0).get<double>() + lstress.at(1).get<double>();
}

} // namespace



int tetragonal_lattice_minimizer(MPI_Comm comm,
                                 std::string& rtdbstring,
                                 std::ostream& coutput,
                                 electronic_minimizer minimizer,
                                 const LatticeContext& ctx)
{
    const bool oprint      = ctx.oprint;
    const std::string& tag = ctx.tag;

    const int    max_steps        = ctx.max_steps;
    const double minimum_step     = ctx.minimum_step;
    const double minimum_gradient = ctx.minimum_gradient;

    double trust_radius = std::abs(ctx.initial_step);
    if (!std::isfinite(trust_radius) || trust_radius <= 0.0)
        trust_radius = 1.0e-2;

    // Log-coordinate state used by the BFGS update.
    bool have_previous_point = false;
    double previous_xa = 0.0, previous_xc = 0.0;
    double previous_ga = 0.0, previous_gc = 0.0;

    Hessian2 B;

    bool converged = false;
    int  steps_taken = 0;

    if (oprint)
    {
        coutput << '\n'
                << tag << "==============================================\n"
                << tag << " PWDFT tetragonal lattice optimization\n"
                << tag << " 2D quasi-Newton (BFGS) in log(a), log(c)\n"
                << tag << "==============================================\n";
    }

    for (int istep = 0; istep < max_steps; ++istep)
    {
        // -- Read current (a, c) ------------------------------------------
        const auto [a, c] = read_tetragonal_lattice(rtdbstring);

        if (!(a > 0.0) || !(c > 0.0))
        {
            coutput << tag << "ERROR: invalid tetragonal cell a=" << a
                    << " c=" << c << '\n';
            break;
        }

        const double xa = std::log(a);
        const double xc = std::log(c);

        // -- Evaluate energy and gradient at current cell -----------------
        const json   current_result   = compute_egs_values(3, comm, minimizer,
                                                           rtdbstring, coutput);
        const double current_energy   = current_result.at("energy").get<double>();
        const json&  lstress          = current_result.at("lstress");

        const double dE_da_a  = lstress.at(0).get<double>();
        const double dE_da_b  = lstress.at(1).get<double>();
        const double dE_dc_c  = lstress.at(2).get<double>();

        const double gradient_a = dE_da_a + dE_da_b;   // constrained tetragonal
        const double gradient_c = dE_dc_c;

        // Log-coordinate gradient
        const double ga = a * gradient_a;
        const double gc = c * gradient_c;

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
                    << tag << " Lattice c   : "
                           << std::fixed << std::setprecision(6)
                           << c << " Bohr ("
                           << std::fixed << std::setprecision(3)
                           << c * bohr_to_angstrom << " A)\n"
                    << tag << " dE/da|tet   : "
                           << std::defaultfloat << std::setprecision(10)
                           << gradient_a << '\n'
                    << tag << " dE/dc|tet   : "
                           << gradient_c << '\n'
                    << tag << " dE/dlog(a)  : " << ga << '\n'
                    << tag << " dE/dlog(c)  : " << gc << '\n'
                    << tag << " Trust radius: " << trust_radius << '\n';
        }

        if (!std::isfinite(current_energy) ||
            !std::isfinite(ga) || !std::isfinite(gc))
        {
            if (oprint)
                coutput << tag << " Action      : non-finite energy or gradient; stop.\n"
                        << tag << "----------------------------------------------\n";
            break;
        }

        // -- Convergence --------------------------------------------------
        if (std::abs(gradient_a) < minimum_gradient &&
            std::abs(gradient_c) < minimum_gradient)
        {
            converged = true;
            if (oprint)
                coutput << tag << " Action      : gradient converged.\n"
                        << tag << "----------------------------------------------\n";
            break;
        }

        // -- BFGS update from the previous accepted step ------------------
        bool updated_hessian = false;
        if (have_previous_point)
        {
            const double s_a = xa - previous_xa;
            const double s_c = xc - previous_xc;
            const double y_a = ga - previous_ga;
            const double y_c = gc - previous_gc;

            updated_hessian = B.bfgs_update(s_a, s_c, y_a, y_c);
        }

        // -- Proposed step ------------------------------------------------
        std::pair<double, double> dx;
        bool used_bfgs = false;

        if (have_previous_point && updated_hessian)
        {
            dx = B.apply(ga, gc);
            dx.first  = -dx.first;
            dx.second = -dx.second;
            used_bfgs = true;
        }
        else
        {
            // Fallback: steepest descent
            dx.first  = -ga;
            dx.second = -gc;
        }

        // Reject non-descent directions; fall back to steepest descent
        if (ga * dx.first + gc * dx.second >= 0.0 ||
            !std::isfinite(dx.first) || !std::isfinite(dx.second))
        {
            dx.first  = -ga;
            dx.second = -gc;
            used_bfgs = false;
        }

        // Trust region: clip ||dx||_2 <= trust_radius
        const double dx_norm = std::hypot(dx.first, dx.second);
        if (dx_norm > trust_radius && dx_norm > 0.0)
        {
            const double scale = trust_radius / dx_norm;
            dx.first  *= scale;
            dx.second *= scale;
        }

        if (!std::isfinite(dx.first) || !std::isfinite(dx.second) ||
            std::hypot(dx.first, dx.second) <
                std::numeric_limits<double>::epsilon())
        {
            if (oprint)
                coutput << tag << " Action      : zero step; stop.\n"
                        << tag << "----------------------------------------------\n";
            break;
        }

        // -- Backtracking line search -------------------------------------
        bool   accepted          = false;
        double accepted_dx_a     = 0.0;
        double accepted_dx_c     = 0.0;
        double accepted_energy   = current_energy;
        std::string accepted_rtdb;

        double trial_dx_a = dx.first;
        double trial_dx_c = dx.second;

        constexpr int maximum_backtracks = 12;

        for (int iback = 0; iback < maximum_backtracks; ++iback)
        {
            const double trial_norm = std::hypot(trial_dx_a, trial_dx_c);
            if (trial_norm < minimum_step)
                break;

            const double a_trial = std::exp(xa + trial_dx_a);
            const double c_trial = std::exp(xc + trial_dx_c);

            if (!std::isfinite(a_trial) || !std::isfinite(c_trial) ||
                !(a_trial > 0.0) || !(c_trial > 0.0))
            {
                trial_dx_a *= 0.5;
                trial_dx_c *= 0.5;
                continue;
            }

            std::string trial_rtdb = rtdbstring;
            try
            {
                set_tetragonal_cell(trial_rtdb, a_trial, c_trial);
            }
            catch (const std::exception& e)
            {
                if (oprint)
                    coutput << tag << " set_tetragonal_cell failed: "
                            << e.what() << "; halving step.\n";
                trial_dx_a *= 0.5;
                trial_dx_c *= 0.5;
                continue;
            }

            const json trial_result =
                compute_egs_values(1, comm, minimizer, trial_rtdb, coutput);
            const double trial_energy = trial_result.at("energy").get<double>();

            if (std::isfinite(trial_energy) && trial_energy < current_energy)
            {
                accepted         = true;
                accepted_dx_a    = trial_dx_a;
                accepted_dx_c    = trial_dx_c;
                accepted_energy  = trial_energy;
                accepted_rtdb    = std::move(trial_rtdb);
                break;
            }

            trial_dx_a *= 0.5;
            trial_dx_c *= 0.5;
        }

        // -- Accept or reject ---------------------------------------------
        if (accepted)
        {
            previous_xa = xa;
            previous_xc = xc;
            previous_ga = ga;
            previous_gc = gc;
            have_previous_point = true;

            rtdbstring = std::move(accepted_rtdb);
            ++steps_taken;

            const double accepted_norm = std::hypot(accepted_dx_a, accepted_dx_c);

            if (accepted_norm < 0.5 * trust_radius)
                trust_radius = std::max(minimum_step, 2.0 * accepted_norm);
            else
                trust_radius = std::min(1.5 * trust_radius,
                                        std::max(std::abs(ctx.initial_step),
                                                 minimum_step));

            if (oprint)
            {
                coutput << tag << " Method      : "
                        << (used_bfgs ? "BFGS/quasi-Newton" : "gradient fallback")
                        << '\n'
                        << tag << " Action      : accepted\n"
                        << tag << " Delta log(a): "
                        << std::defaultfloat << std::setprecision(10)
                        << accepted_dx_a << '\n'
                        << tag << " Delta log(c): " << accepted_dx_c << '\n'
                        << tag << " Scale a     : " << std::exp(accepted_dx_a) << '\n'
                        << tag << " Scale c     : " << std::exp(accepted_dx_c) << '\n'
                        << tag << " New energy  : "
                        << std::fixed << std::setprecision(10)
                        << accepted_energy << " Hartree\n"
                        << tag << "----------------------------------------------\n";
            }
        }
        else
        {
            have_previous_point = false;
            B.reset_to_identity();
            trust_radius *= 0.5;

            if (oprint)
            {
                coutput << tag << " Method      : "
                        << (used_bfgs ? "BFGS/quasi-Newton" : "gradient fallback")
                        << '\n'
                        << tag << " Action      : rejected; halve trust radius.\n"
                        << tag << " New radius  : "
                        << std::defaultfloat << std::setprecision(10)
                        << trust_radius << '\n'
                        << tag << "----------------------------------------------\n";
            }

            if (trust_radius < minimum_step)
                break;
        }
    }

    // -- Final report -----------------------------------------------------
    const json   final_result   = compute_egs_values(3, comm, minimizer,
                                                     rtdbstring, coutput);
    const double final_energy   = final_result.at("energy").get<double>();
    const json&  final_lstress  = final_result.at("lstress");

    const double final_dE_da_a = final_lstress.at(0).get<double>();
    const double final_dE_da_b = final_lstress.at(1).get<double>();
    const double final_dE_dc_c = final_lstress.at(2).get<double>();

    const double final_gradient_a = final_dE_da_a + final_dE_da_b;
    const double final_gradient_c = final_dE_dc_c;

    const auto [final_a, final_c] = read_tetragonal_lattice(rtdbstring);

    if (std::isfinite(final_gradient_a) &&
        std::isfinite(final_gradient_c) &&
        std::abs(final_gradient_a) < minimum_gradient &&
        std::abs(final_gradient_c) < minimum_gradient)
    {
        converged = true;
    }

    if (oprint)
    {
        coutput << '\n'
                << tag << "==============================================\n"
                << tag << " PWDFT tetragonal lattice optimization COMPLETE\n"
                << tag << "==============================================\n"
                << tag << " Final lattice a: "
                       << std::fixed << std::setprecision(6) << final_a
                       << " Bohr = "
                       << std::fixed << std::setprecision(3)
                       << final_a * bohr_to_angstrom << " A\n"
                << tag << " Final lattice c: "
                       << std::fixed << std::setprecision(6) << final_c
                       << " Bohr = "
                       << std::fixed << std::setprecision(3)
                       << final_c * bohr_to_angstrom << " A\n"
                << tag << " Minimum energy (total): "
                       << std::fixed << std::setprecision(8)
                       << final_energy << " Hartree\n"
                << tag << " Gradients at minimum:\n"
                << tag << "   dE/da|tet = " << final_gradient_a << '\n'
                << tag << "   dE/dc|tet = " << final_gradient_c << '\n'
                << tag << " Accepted lattice steps: " << steps_taken << '\n'
                << tag << " Status: "
                       << (converged ? "Converged"
                                     : "Stopped before gradient convergence")
                       << '\n'
                << tag << "==============================================\n";
    }

    return 0;
}

} // namespace pwdft

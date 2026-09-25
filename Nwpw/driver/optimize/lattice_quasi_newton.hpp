// lattice_quasi_newton.hpp
//
// Templated N-dimensional quasi-Newton engine for lattice optimization.
//
// All per-system lattice minimizers share the same skeleton:
//
//   1. read cell parameters (in whatever coordinates they chose)
//   2. evaluate energy and gradient via compute_egs_values
//   3. test convergence, propose a step (BFGS or fallback)
//   4. clip to a trust radius, backtrack until energy improves
//   5. update the inverse-Hessian B and the trust radius
//
// Only the coordinate system, the cell read/write, and the gradient
// mapping differ between systems. Those are supplied by the caller as
// lambdas; everything else lives here.

#pragma once

#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <string>

#include <mpi.h>

#include "lattice_minimizer.hpp"
#include "lattice_common.hpp"
#include "json.hpp"

namespace pwdft {

using json = nlohmann::json;

inline constexpr double bohr_to_angstrom = 0.529177;

// ---------------------------------------------------------------------------
// N-dimensional symmetric inverse-Hessian approximation
// ---------------------------------------------------------------------------

template <int N>
struct HessianN
{
    std::array<std::array<double, N>, N> B{};

    HessianN(double initial_diag = 0.1)
    {
        reset_to_scaled_identity(initial_diag);
    }

    void reset_to_scaled_identity(double d)
    {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                B[i][j] = (i == j) ? d : 0.0;
    }

    std::array<double, N> apply(const std::array<double, N>& g) const
    {
        std::array<double, N> out{};
        for (int i = 0; i < N; ++i)
        {
            double s = 0.0;
            for (int j = 0; j < N; ++j)
                s += B[i][j] * g[j];
            out[i] = s;
        }
        return out;
    }

    // BFGS update. Returns false if the curvature condition y.s > 0 fails.
    bool bfgs_update(const std::array<double, N>& s,
                     const std::array<double, N>& y)
    {
        double ys = 0.0;
        double y_abs = 0.0;
        double s_abs = 0.0;
        for (int i = 0; i < N; ++i)
        {
            ys    += y[i] * s[i];
            y_abs += std::abs(y[i]);
            s_abs += std::abs(s[i]);
        }

        if (!(ys > 1.0e-12 * y_abs * s_abs + 1.0e-16))
            return false;

        const double rho = 1.0 / ys;

        // V = I - rho * y * s^T
        std::array<std::array<double, N>, N> V{};
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                V[i][j] = ((i == j) ? 1.0 : 0.0) - rho * y[i] * s[j];

        // t = V * B
        std::array<std::array<double, N>, N> t{};
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
            {
                double sum = 0.0;
                for (int k = 0; k < N; ++k)
                    sum += V[i][k] * B[k][j];
                t[i][j] = sum;
            }

        // B <- t * V^T + rho * s * s^T
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
            {
                double sum = 0.0;
                for (int k = 0; k < N; ++k)
                    sum += t[i][k] * V[j][k];
                sum += rho * s[i] * s[j];
                B[i][j] = sum;
            }

        return true;
    }

    // Keep B positive-definite after an update. Diagonal floor plus a
    // conservative off-diagonal clamp (stronger than the strict PD condition
    // for N > 2, but safe).
    void enforce_positive_definite(double B_min = 1.0e-2)
    {
        for (int i = 0; i < N; ++i)
            if (B[i][i] < B_min)
                B[i][i] = B_min;

        const double denom = std::max(1, N - 1);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
            {
                if (i == j) continue;
                const double bound =
                    std::sqrt(B[i][i] * B[j][j]) * 0.95 / denom;
                if (B[i][j] >  bound) B[i][j] =  bound;
                if (B[i][j] < -bound) B[i][j] = -bound;
            }
    }

    // Solve B z = s (B positive-definite), return s . z.
    // Returns 0 if B is numerically singular.
    double inverse_quadratic_form(const std::array<double, N>& s) const
    {
        std::array<std::array<double, N>, N> A = B;
        std::array<double, N> z = s;

        for (int k = 0; k < N; ++k)
        {
            const double pivot = A[k][k];
            if (std::abs(pivot) < 1.0e-14)
                return 0.0;
            for (int i = k + 1; i < N; ++i)
            {
                const double factor = A[i][k] / pivot;
                for (int j = k; j < N; ++j)
                    A[i][j] -= factor * A[k][j];
                z[i] -= factor * z[k];
            }
        }

        for (int i = N - 1; i >= 0; --i)
        {
            double sum = z[i];
            for (int j = i + 1; j < N; ++j)
                sum -= A[i][j] * z[j];
            z[i] = sum / A[i][i];
        }

        double result = 0.0;
        for (int i = 0; i < N; ++i)
            result += s[i] * z[i];
        return result;
    }
};

// ---------------------------------------------------------------------------
// Generic quasi-Newton driver
//
// Caller supplies:
//   read_x(rtdb)          -> std::array<double,N>   (current coords)
//   write_x(rtdb, x)      -> void                   (modify rtdb to coords x)
//   compute_g(x, lstress) -> std::array<double,N>   (gradient in coords)
//   print_step(out, tag, istep, x, g, E, trust)
//
// Returns 0 on completion. Convergence status is reported via output.
// ---------------------------------------------------------------------------

template <int N, typename ReadX, typename WriteX, typename ComputeG, typename PrintStep, typename PrintFinal>
int run_quasi_newton(std::string& rtdbstring,
                     MPI_Comm comm,
                     electronic_minimizer minimizer,
                     std::ostream& coutput,
                     const LatticeContext& ctx,
                     const std::string& system_name,
                     ReadX&& read_x,
                     WriteX&& write_x,
                     ComputeG&& compute_g,
                     PrintStep&& print_step,
                     PrintFinal&& print_final)
{
    const bool oprint = ctx.oprint;
    const std::string& tag = ctx.tag;

    const int    max_steps        = ctx.max_steps;
    const double minimum_step     = ctx.minimum_step;
    const double minimum_gradient = ctx.minimum_gradient;

    double trust_radius = std::abs(ctx.initial_step);
    if (!std::isfinite(trust_radius) || trust_radius <= 0.0)
        trust_radius = 1.0e-2;
    const double trust_upper_bound =
        std::max(10.0 * std::abs(ctx.initial_step), minimum_step);

    HessianN<N> B(0.1);

    bool have_previous_point = false;
    std::array<double, N> previous_x{};
    std::array<double, N> previous_g{};

    bool converged   = false;
    int  steps_taken = 0;

    if (oprint)
    {
        coutput << tag << "==============================================================================\n"
                << tag << " PWDFT " << system_name << " lattice optimization\n"
                << tag << " " << N << "D quasi-Newton (BFGS) in log coordinates\n"
                << tag << "==============================================================================\n";
    }

    for (int istep = 0; istep < max_steps; ++istep)
    {
        const std::array<double, N> x = read_x(rtdbstring);

        const json   current_result = compute_egs_values(3, comm, minimizer, rtdbstring, coutput);
        const double current_energy = current_result.at("energy").get<double>();
        const json&  lstress        = current_result.at("lstress");

        const std::array<double, N> g = compute_g(x, lstress);

        if (oprint)
            print_step(coutput, tag, istep, x, g, current_energy, trust_radius);

        bool finite = std::isfinite(current_energy);
        for (int i = 0; i < N; ++i)
            finite = finite && std::isfinite(x[i]) && std::isfinite(g[i]);
        if (!finite)
        {
            if (oprint)
                coutput << tag << " Action        : non-finite; stop.\n"
                        << tag << "----------------------------------------------\n";
            break;
        }

        double max_g = 0.0;
        for (int i = 0; i < N; ++i)
            max_g = std::max(max_g, std::abs(g[i]));
        if (max_g < minimum_gradient)
        {
            converged = true;
            if (oprint)
                coutput << tag << " Action        : gradient converged.\n"
                        << tag << "----------------------------------------------\n";
            break;
        }

        // BFGS update from previous accepted step
        bool updated_hessian = false;
        if (have_previous_point)
        {
            std::array<double, N> s{}, y{};
            double s_norm2 = 0.0, y_norm2 = 0.0, g_norm2 = 0.0;
            for (int i = 0; i < N; ++i)
            {
                s[i] = x[i] - previous_x[i];
                y[i] = g[i] - previous_g[i];
                s_norm2 += s[i] * s[i];
                y_norm2 += y[i] * y[i];
                g_norm2 += g[i] * g[i];
            }
            const double s_norm = std::sqrt(s_norm2);
            const double y_norm = std::sqrt(y_norm2);
            const double g_norm = std::sqrt(g_norm2);

            const bool usable_pair = (s_norm > 1.0e-6) && (y_norm > 1.0e-3 * g_norm);

            if (usable_pair)
            {
                updated_hessian = B.bfgs_update(s, y);
                if (updated_hessian)
                    B.enforce_positive_definite();
            }
        }

        // Propose step
        std::array<double, N> dx{};
        bool used_bfgs = false;
        if (have_previous_point && updated_hessian)
        {
            dx = B.apply(g);
            for (int i = 0; i < N; ++i) dx[i] = -dx[i];
            used_bfgs = true;
        }
        else
        {
            for (int i = 0; i < N; ++i) dx[i] = -g[i];
        }

        // Reject non-descent
        double g_dot_dx = 0.0;
        for (int i = 0; i < N; ++i) g_dot_dx += g[i] * dx[i];
        if (g_dot_dx >= 0.0)
        {
            for (int i = 0; i < N; ++i) dx[i] = -g[i];
            used_bfgs = false;
        }

        // Trust region clip
        double dx_norm2 = 0.0;
        for (int i = 0; i < N; ++i) dx_norm2 += dx[i] * dx[i];
        double dx_norm = std::sqrt(dx_norm2);
        if (dx_norm > trust_radius && dx_norm > 0.0)
        {
            const double scale = trust_radius / dx_norm;
            for (int i = 0; i < N; ++i) dx[i] *= scale;
            dx_norm = trust_radius;
        }

        if (dx_norm < std::numeric_limits<double>::epsilon())
        {
            if (oprint)
                coutput << tag << " Action        : zero step; stop.\n"
                        << tag << "----------------------------------------------\n";
            break;
        }

        // Backtracking line search
        bool accepted = false;
        std::array<double, N> accepted_dx{};
        double accepted_energy = current_energy;
        std::string accepted_rtdb;

        std::array<double, N> trial_dx = dx;
        constexpr int maximum_backtracks = 12;

        for (int iback = 0; iback < maximum_backtracks; ++iback)
        {
            double trial_norm2 = 0.0;
            for (int i = 0; i < N; ++i) trial_norm2 += trial_dx[i] * trial_dx[i];
            if (std::sqrt(trial_norm2) < minimum_step)
                break;

            std::array<double, N> x_trial{};
            for (int i = 0; i < N; ++i) x_trial[i] = x[i] + trial_dx[i];

            std::string trial_rtdb = rtdbstring;
            bool wrote_ok = false;
            try
            {
                write_x(trial_rtdb, x_trial);
                wrote_ok = true;
            }
            catch (const std::exception& e)
            {
                if (oprint)
                    coutput << tag << " write_x failed: " << e.what()
                            << "; halving step.\n";
            }
            if (!wrote_ok)
            {
                for (int i = 0; i < N; ++i) trial_dx[i] *= 0.5;
                continue;
            }

            const json trial_result =
                compute_egs_values(1, comm, minimizer, trial_rtdb, coutput);
            const double trial_energy = trial_result.at("energy").get<double>();

            if (std::isfinite(trial_energy) && trial_energy < current_energy)
            {
                accepted        = true;
                accepted_dx     = trial_dx;
                accepted_energy = trial_energy;
                accepted_rtdb   = std::move(trial_rtdb);
                break;
            }

            for (int i = 0; i < N; ++i) trial_dx[i] *= 0.5;
        }

        if (accepted)
        {
            previous_x = x;
            previous_g = g;
            have_previous_point = true;

            rtdbstring = std::move(accepted_rtdb);
            ++steps_taken;

            // Trust-region update via rho
            const double sHs = B.inverse_quadratic_form(accepted_dx);

            double g_dot_s = 0.0;
            for (int i = 0; i < N; ++i) g_dot_s += g[i] * accepted_dx[i];
            const double predicted_reduction = -g_dot_s - 0.5 * sHs;
            const double actual_reduction    = current_energy - accepted_energy;

            double rho = 0.0;
            if (predicted_reduction > 1.0e-14)
                rho = actual_reduction / predicted_reduction;

            double accepted_norm2 = 0.0;
            for (int i = 0; i < N; ++i) accepted_norm2 += accepted_dx[i] * accepted_dx[i];
            const double accepted_norm = std::sqrt(accepted_norm2);

            if (rho < 0.25)
                trust_radius *= 0.5;
            else if (rho > 0.75 && accepted_norm > 0.9 * trust_radius)
                trust_radius = std::min(2.0 * trust_radius, trust_upper_bound);

            if (oprint)
                coutput << tag << " Method        : "
                        << (used_bfgs ? "BFGS/quasi-Newton" : "gradient fallback")
                        << '\n'
                        << tag << " Action        : accepted\n"
                        << tag << " rho           : " << rho << '\n'
                        << tag << " Trust radius  : " << trust_radius << '\n'
                        << tag << " New energy    : "
                        << std::fixed << std::setprecision(10)
                        << accepted_energy << " Hartree\n"
                        << tag << "----------------------------------------------\n";
        }
        else
        {
            trust_radius *= 0.5;

            if (oprint)
                coutput << tag << " Method        : "
                        << (used_bfgs ? "BFGS/quasi-Newton" : "gradient fallback")
                        << '\n'
                        << tag << " Action        : rejected; halve trust radius.\n"
                        << tag << " New radius    : "
                        << std::defaultfloat << std::setprecision(10)
                        << trust_radius << '\n'
                        << tag << "----------------------------------------------\n";

            if (trust_radius < minimum_step)
                break;
        }
    }

    // Final report
    const json   final_result  = compute_egs_values(3, comm, minimizer, rtdbstring, coutput);
    const double final_energy  = final_result.at("energy").get<double>();
    const json&  final_lstress = final_result.at("lstress");
    const std::array<double, N> final_x = read_x(rtdbstring);
    const std::array<double, N> final_g = compute_g(final_x, final_lstress);

    double final_max_g = 0.0;
    for (int i = 0; i < N; ++i)
        final_max_g = std::max(final_max_g, std::abs(final_g[i]));
    if (final_max_g < minimum_gradient)
        converged = true;

    if (oprint)
    {
        coutput << '\n'
                << tag << "==============================================\n"
                << tag << " PWDFT " << system_name
                       << " lattice optimization COMPLETE\n"
                << tag << "==============================================\n";

        print_final(coutput, tag, final_x, final_energy);

        coutput << tag << " Final energy   : "
                       << std::fixed << std::setprecision(8)
                       << final_energy << " Hartree\n"
                << tag << " Max |gradient| : "
                       << std::defaultfloat << std::setprecision(10)
                       << final_max_g << '\n'
                << tag << " Accepted steps : " << steps_taken << '\n'
                << tag << " Status         : "
                       << (converged ? "Converged"
                                     : "Stopped before gradient convergence")
                       << '\n'
                << tag << "==============================================\n";
    }

    return 0;
}

} // namespace pwdft

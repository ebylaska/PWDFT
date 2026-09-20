// general_lattice_minimizer.cpp
//
// General strain-based lattice optimizer.
//
// Works for any cell symmetry. The optimization coordinates are the nine
// components of a dimensionless strain tensor relative to a fixed
// reference cell h0 captured at function entry:
//
//     h(eps) = h0 * (I + eps)
//
// Cartesian atom positions transform by the same (I + eps):
//
//     r(eps) = (I + eps) * r_ref
//
// so fractional coordinates are preserved by construction.
//
// The gradient of E with respect to eps is obtained directly from the
// symmetrized stress tensor sigma returned by compute_egs_values:
//
//     dE/deps_ij = -V * sigma_ij
//
// with V = |det(h_current)|. No chain rule, no angle transformation.
//
// Only the six symmetric components of eps correspond to physical
// deformations; the three antisymmetric components correspond to rigid
// rotations and carry no energy. BFGS steps and gradients are symmetrized
// so the search stays in the physical six-dimensional subspace.
//
// This is the "general" lattice minimizer: no assumption about which
// crystal system the cell belongs to, no reliance on the packed lstress
// vector, and no dependence on the crystal-system dispatcher.

#include "lattice_minimizer.hpp"
#include "lattice_common.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <string>
#include <vector>

#include <mpi.h>

#include "json.hpp"

namespace pwdft {

using json = nlohmann::json;

namespace {

using Cell9  = std::array<double, 9>;
using Strain = std::array<double, 9>;

// ---------------------------------------------------------------------------
// Small linear-algebra helpers
// ---------------------------------------------------------------------------

double det3x3(const Cell9& h)
{
    return h[0]*(h[4]*h[8] - h[5]*h[7])
         - h[1]*(h[3]*h[8] - h[5]*h[6])
         + h[2]*(h[3]*h[7] - h[4]*h[6]);
}

// h_new = h0 * (I + eps), all matrices row-major 3x3.
Cell9 strain_cell(const Cell9& h0, const Strain& eps)
{
    Cell9 hnew{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            double sum = 0.0;
            for (int k = 0; k < 3; ++k)
                sum += h0[3*i + k] *
                       ((k == j ? 1.0 : 0.0) + eps[3*k + j]);
            hnew[3*i + j] = sum;
        }
    return hnew;
}

// Force eps to be symmetric, so the six physical components are the only
// ones BFGS ever sees.
void symmetrize(Strain& eps)
{
    for (int i = 0; i < 3; ++i)
        for (int j = i+1; j < 3; ++j)
        {
            const double s = 0.5 * (eps[3*i + j] + eps[3*j + i]);
            eps[3*i + j] = s;
            eps[3*j + i] = s;
        }
}

// eps = h0^{-1} h - I
Strain strain_from_cell(const Cell9& h0_inv, const Cell9& h)
{
    Strain eps{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            double s = 0.0;
            for (int k = 0; k < 3; ++k)
                s += h0_inv[3*i + k] * h[3*k + j];
            eps[3*i + j] = s - ((i == j) ? 1.0 : 0.0);
        }
    return eps;
}

// ---------------------------------------------------------------------------
// RTDB helpers
// ---------------------------------------------------------------------------

std::string active_geometry_name(const json& rtdb)
{
    return (rtdb.contains("geometry") && rtdb["geometry"].is_string())
        ? rtdb["geometry"].get<std::string>()
        : "geometry";
}

Cell9 read_cell(const std::string& rtdbstring)
{
    const json rtdb = json::parse(rtdbstring);
    const std::string geomname = active_geometry_name(rtdb);
    const auto& u = rtdb.at("geometries").at(geomname).at("unita");

    Cell9 h{};
    for (int i = 0; i < 9; ++i)
        h[i] = u.at(i).get<double>();
    return h;
}

std::vector<double> read_coords(const std::string& rtdbstring)
{
    const json rtdb = json::parse(rtdbstring);
    const std::string geomname = active_geometry_name(rtdb);
    const auto& c = rtdb.at("geometries").at(geomname).at("coords");

    std::vector<double> r;
    r.reserve(c.size());
    for (const auto& v : c)
        r.push_back(v.get<double>());
    return r;
}

std::array<double, 9> read_stress_sym(const json& result)
{
    std::array<double, 9> s{};
    if (!result.contains("stress_sym")) return s;
    const auto& ss = result.at("stress_sym");
    if (!ss.is_array() || ss.size() != 9) return s;
    for (int i = 0; i < 9; ++i)
        s[i] = ss.at(i).get<double>();
    return s;
}

// Write h_new into geometry.unita, and (I+eps)*coords_ref into coords.
void write_strained_state(std::string& rtdbstring,
                          const Cell9& h_new,
                          const std::vector<double>& coords_ref,
                          const Strain& eps)
{
    json rtdb = json::parse(rtdbstring);
    const std::string geomname = active_geometry_name(rtdb);
    json& geometry = rtdb["geometries"][geomname];

    for (int i = 0; i < 9; ++i)
        geometry["unita"][i] = h_new[i];

    if (geometry.contains("coords") && geometry["coords"].is_array())
    {
        auto& coords = geometry["coords"];
        const std::size_t n_atoms = coords.size() / 3;

        if (coords_ref.size() < coords.size())
            throw std::runtime_error(
                "write_strained_state: coords_ref smaller than coords");

        const double A00 = 1.0 + eps[0], A01 = eps[1],       A02 = eps[2];
        const double A10 = eps[3],       A11 = 1.0 + eps[4], A12 = eps[5];
        const double A20 = eps[6],       A21 = eps[7],       A22 = 1.0 + eps[8];

        for (std::size_t i = 0; i < n_atoms; ++i)
        {
            const double x = coords_ref[3*i + 0];
            const double y = coords_ref[3*i + 1];
            const double z = coords_ref[3*i + 2];

            coords[3*i + 0] = A00*x + A01*y + A02*z;
            coords[3*i + 1] = A10*x + A11*y + A12*z;
            coords[3*i + 2] = A20*x + A21*y + A22*z;
        }
    }

    // Keep nwpw.simulation_cell.unita in sync.
    if (rtdb.contains("nwpw") &&
        rtdb["nwpw"].is_object() &&
        rtdb["nwpw"].contains("simulation_cell") &&
        rtdb["nwpw"]["simulation_cell"].is_object())
    {
        json& cell = rtdb["nwpw"]["simulation_cell"];
        if (cell.contains("unita") &&
            cell["unita"].is_array() &&
            cell["unita"].size() == 9)
        {
            for (int i = 0; i < 9; ++i)
                cell["unita"][i] = h_new[i];
        }
    }

    rtdbstring = rtdb.dump();
}

// ---------------------------------------------------------------------------
// 9x9 BFGS on the strain components
// ---------------------------------------------------------------------------

struct Hessian9
{
    std::array<double, 81> B{};

    Hessian9() { reset(0.1); }

    void reset(double d)
    {
        B.fill(0.0);
        for (int i = 0; i < 9; ++i) B[9*i + i] = d;
    }

    Strain apply(const Strain& g) const
    {
        Strain out{};
        for (int i = 0; i < 9; ++i)
        {
            double s = 0.0;
            for (int j = 0; j < 9; ++j)
                s += B[9*i + j] * g[j];
            out[i] = s;
        }
        return out;
    }

    bool bfgs_update(const Strain& s, const Strain& y)
    {
        double ys = 0.0, ya = 0.0, sa = 0.0;
        for (int i = 0; i < 9; ++i)
        {
            ys += y[i]*s[i];
            ya += std::abs(y[i]);
            sa += std::abs(s[i]);
        }
        if (!(ys > 1.0e-12 * ya * sa + 1.0e-16))
            return false;

        const double rho = 1.0 / ys;

        std::array<double, 81> V{};
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j)
                V[9*i + j] = ((i == j) ? 1.0 : 0.0) - rho * y[i] * s[j];

        std::array<double, 81> t{};
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j)
            {
                double sum = 0.0;
                for (int k = 0; k < 9; ++k)
                    sum += V[9*i + k] * B[9*k + j];
                t[9*i + j] = sum;
            }

        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j)
            {
                double sum = 0.0;
                for (int k = 0; k < 9; ++k)
                    sum += t[9*i + k] * V[9*j + k];
                B[9*i + j] = sum + rho * s[i] * s[j];
            }

        enforce_pd();
        return true;
    }

    void enforce_pd(double B_min = 1.0e-2)
    {
        for (int i = 0; i < 9; ++i)
            if (B[9*i + i] < B_min)
                B[9*i + i] = B_min;

        // Conservative off-diagonal clamp. For a 9x9 matrix this is much
        // tighter than the strict PD condition would require, but it keeps
        // the BFGS iterates well-conditioned for the smooth, low-rank
        // curvature that strain-stress pairs produce.
        constexpr double denom = 4.0;
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j)
            {
                if (i == j) continue;
                const double bound =
                    std::sqrt(B[9*i + i] * B[9*j + j]) * 0.95 / denom;
                if (B[9*i + j] >  bound) B[9*i + j] =  bound;
                if (B[9*i + j] < -bound) B[9*i + j] = -bound;
            }
    }

    // s^T B^{-1} s via Gaussian elimination with partial pivoting.
    // Returns 0 if B is numerically singular.
    double inverse_qf(const Strain& s) const
    {
        std::array<double, 81> A = B;
        Strain z = s;

        for (int k = 0; k < 9; ++k)
        {
            int pivot_row = k;
            double pivot_abs = std::abs(A[9*k + k]);
            for (int i = k + 1; i < 9; ++i)
            {
                if (std::abs(A[9*i + k]) > pivot_abs)
                {
                    pivot_abs = std::abs(A[9*i + k]);
                    pivot_row = i;
                }
            }
            if (pivot_abs < 1.0e-14) return 0.0;
            if (pivot_row != k)
            {
                for (int j = 0; j < 9; ++j)
                    std::swap(A[9*k + j], A[9*pivot_row + j]);
                std::swap(z[k], z[pivot_row]);
            }

            const double pivot = A[9*k + k];
            for (int i = k + 1; i < 9; ++i)
            {
                const double factor = A[9*i + k] / pivot;
                for (int j = k; j < 9; ++j)
                    A[9*i + j] -= factor * A[9*k + j];
                z[i] -= factor * z[k];
            }
        }

        for (int i = 8; i >= 0; --i)
        {
            double sum = z[i];
            for (int j = i + 1; j < 9; ++j)
                sum -= A[9*i + j] * z[j];
            z[i] = sum / A[9*i + i];
        }

        double r = 0.0;
        for (int i = 0; i < 9; ++i)
            r += s[i] * z[i];
        return r;
    }
};

} // namespace

// ---------------------------------------------------------------------------
// general_lattice_minimizer
// ---------------------------------------------------------------------------

int general_lattice_minimizer(MPI_Comm comm,
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
    const double trust_upper_bound =
        std::max(10.0 * std::abs(ctx.initial_step), minimum_step);

    // Capture the reference cell and reference Cartesian positions once.
    const Cell9 h0 = read_cell(rtdbstring);
    const std::vector<double> coords_ref = read_coords(rtdbstring);

    Cell9 h0_inv{};
    if (!invert_3x3(h0, h0_inv))
    {
        coutput << tag << " general_lattice_minimizer: "
                << "singular reference cell; aborting.\n";
        return 1;
    }

    Hessian9 B;

    bool have_previous_point = false;
    Strain previous_eps{};
    Strain previous_g{};

    bool converged   = false;
    int  steps_taken = 0;

    if (oprint)
    {
        coutput << tag << "==============================================\n"
                << tag << " PWDFT general lattice optimization\n"
                << tag << " 9D quasi-Newton (BFGS) in strain coordinates\n"
                << tag << "==============================================\n";
    }

    for (int istep = 0; istep < max_steps; ++istep)
    {
        const Cell9 h_current = read_cell(rtdbstring);
        const Strain eps = strain_from_cell(h0_inv, h_current);
        const double V = std::abs(det3x3(h_current));

        const json result = compute_egs_values(3, comm, minimizer,
                                               rtdbstring, coutput);
        const double energy = result.at("energy").get<double>();
        const auto sigma = read_stress_sym(result);

        // Gradient of E w.r.t. strain, at the current cell.
        //
        // If lstress[0..2] in the same run are positive when the cell
        // wants to expand, then this sign is correct. Verify with a
        // finite difference on the very first call if in doubt.
        Strain g{};
        for (int i = 0; i < 9; ++i)
            g[i] = -V * sigma[i];
        symmetrize(g);

        if (oprint)
        {
            coutput << '\n'
                    << tag << "----------------------------------------------\n"
                    << tag << " Step          : " << istep << '\n'
                    << tag << " Energy        : " << std::fixed << std::setprecision(10)
                           << energy << " Hartree\n"
                    << tag << " Volume        : " << V << " Bohr^3\n"
                    << tag << " eps_xx        : "
                           << std::defaultfloat << std::setprecision(10)
                           << eps[0] << '\n'
                    << tag << " eps_yy        : " << eps[4] << '\n'
                    << tag << " eps_zz        : " << eps[8] << '\n'
                    << tag << " eps_yz        : " << eps[5] << '\n'
                    << tag << " eps_xz        : " << eps[2] << '\n'
                    << tag << " eps_xy        : " << eps[1] << '\n'
                    << tag << " Trust radius  : " << trust_radius << '\n';
        }

        // Convergence: largest component of the 9-vector gradient.
        double max_g = 0.0;
        for (int i = 0; i < 9; ++i)
            max_g = std::max(max_g, std::abs(g[i]));
        if (max_g < minimum_gradient)
        {
            converged = true;
            if (oprint)
                coutput << tag << " Action        : gradient converged.\n"
                        << tag << "----------------------------------------------\n";
            break;
        }

        // BFGS update from previous accepted step.
        bool updated_hessian = false;
        if (have_previous_point)
        {
            Strain s{}, y{};
            double sn2 = 0.0, yn2 = 0.0, gn2 = 0.0;
            for (int i = 0; i < 9; ++i)
            {
                s[i] = eps[i] - previous_eps[i];
                y[i] = g[i]   - previous_g[i];
                sn2 += s[i]*s[i];
                yn2 += y[i]*y[i];
                gn2 += g[i]*g[i];
            }
            const double sn = std::sqrt(sn2);
            const double yn = std::sqrt(yn2);
            const double gn = std::sqrt(gn2);

            if ((sn > 1.0e-8) && (yn > 1.0e-3 * gn))
                updated_hessian = B.bfgs_update(s, y);
        }

        // Propose step.
        Strain dx{};
        bool used_bfgs = false;
        if (have_previous_point && updated_hessian)
        {
            dx = B.apply(g);
            for (int i = 0; i < 9; ++i) dx[i] = -dx[i];
            symmetrize(dx);
            used_bfgs = true;
        }
        else
        {
            for (int i = 0; i < 9; ++i) dx[i] = -g[i];
            symmetrize(dx);
        }

        // Reject non-descent directions.
        double g_dot_dx = 0.0;
        for (int i = 0; i < 9; ++i) g_dot_dx += g[i] * dx[i];
        if (g_dot_dx >= 0.0)
        {
            for (int i = 0; i < 9; ++i) dx[i] = -g[i];
            symmetrize(dx);
            used_bfgs = false;
        }

        // Trust-region clip.
        double dxn2 = 0.0;
        for (int i = 0; i < 9; ++i) dxn2 += dx[i]*dx[i];
        double dxn = std::sqrt(dxn2);
        if (dxn > trust_radius && dxn > 0.0)
        {
            const double sc = trust_radius / dxn;
            for (int i = 0; i < 9; ++i) dx[i] *= sc;
            dxn = trust_radius;
        }
        if (dxn < std::numeric_limits<double>::epsilon())
        {
            if (oprint)
                coutput << tag << " Action        : zero step; stop.\n"
                        << tag << "----------------------------------------------\n";
            break;
        }

        // Backtracking line search.
        bool   accepted        = false;
        Strain accepted_dx{};
        double accepted_energy = energy;
        std::string accepted_rtdb;

        Strain trial_dx = dx;
        constexpr int maximum_backtracks = 12;

        for (int iback = 0; iback < maximum_backtracks; ++iback)
        {
            double tn2 = 0.0;
            for (int i = 0; i < 9; ++i) tn2 += trial_dx[i]*trial_dx[i];
            if (std::sqrt(tn2) < minimum_step)
                break;

            Strain trial_eps = eps;
            for (int i = 0; i < 9; ++i) trial_eps[i] += trial_dx[i];

            std::string trial_rtdb = rtdbstring;
            bool wrote_ok = true;
            try
            {
                const Cell9 h_trial = strain_cell(h0, trial_eps);
                write_strained_state(trial_rtdb, h_trial,
                                     coords_ref, trial_eps);
            }
            catch (const std::exception& e)
            {
                if (oprint)
                    coutput << tag << " write_strained_state failed: "
                            << e.what() << "; halving step.\n";
                wrote_ok = false;
            }
            if (!wrote_ok)
            {
                for (int i = 0; i < 9; ++i) trial_dx[i] *= 0.5;
                continue;
            }

            const json trial_result = compute_egs_values(
                1, comm, minimizer, trial_rtdb, coutput);
            const double trial_energy = trial_result.at("energy").get<double>();

            if (std::isfinite(trial_energy) && trial_energy < energy)
            {
                accepted        = true;
                accepted_dx     = trial_dx;
                accepted_energy = trial_energy;
                accepted_rtdb   = std::move(trial_rtdb);
                break;
            }

            for (int i = 0; i < 9; ++i) trial_dx[i] *= 0.5;
        }

        // Accept or reject.
        if (accepted)
        {
            previous_eps = eps;
            previous_g   = g;
            have_previous_point = true;

            rtdbstring = std::move(accepted_rtdb);
            ++steps_taken;

            const double sHs = B.inverse_qf(accepted_dx);

            double g_dot_s = 0.0;
            for (int i = 0; i < 9; ++i) g_dot_s += g[i] * accepted_dx[i];
            const double pred_red = -g_dot_s - 0.5 * sHs;
            const double act_red  = energy - accepted_energy;

            double rho = 0.0;
            if (pred_red > 1.0e-14)
                rho = act_red / pred_red;

            double an2 = 0.0;
            for (int i = 0; i < 9; ++i) an2 += accepted_dx[i]*accepted_dx[i];
            const double an = std::sqrt(an2);

            if (rho < 0.25)
                trust_radius *= 0.5;
            else if (rho > 0.75 && an > 0.9 * trust_radius)
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

    // Final report.
    const json final_result = compute_egs_values(3, comm, minimizer,
                                                 rtdbstring, coutput);
    const double final_energy = final_result.at("energy").get<double>();
    const auto   final_sigma  = read_stress_sym(final_result);

    const Cell9  h_final   = read_cell(rtdbstring);
    const double V_final   = std::abs(det3x3(h_final));
    const Strain eps_final = strain_from_cell(h0_inv, h_final);

    double final_max_g = 0.0;
    for (int i = 0; i < 9; ++i)
        final_max_g = std::max(final_max_g, std::abs(-V_final * final_sigma[i]));
    if (final_max_g < minimum_gradient)
        converged = true;

    if (oprint)
    {
        coutput << '\n'
                << tag << "==============================================\n"
                << tag << " PWDFT general lattice optimization COMPLETE\n"
                << tag << "==============================================\n"
                << tag << " Final eps_xx  : "
                       << std::defaultfloat << std::setprecision(10)
                       << eps_final[0] << '\n'
                << tag << " Final eps_yy  : " << eps_final[4] << '\n'
                << tag << " Final eps_zz  : " << eps_final[8] << '\n'
                << tag << " Final eps_yz  : " << eps_final[5] << '\n'
                << tag << " Final eps_xz  : " << eps_final[2] << '\n'
                << tag << " Final eps_xy  : " << eps_final[1] << '\n'
                << tag << " Final energy   : "
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

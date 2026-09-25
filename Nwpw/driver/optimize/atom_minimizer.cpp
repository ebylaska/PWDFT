// atom_minimizer.cpp
//
// Atom-coordinate optimizer for PWDFT.
//
// Optimization variable: fractional coordinates of all atoms.
// Engine: L-BFGS (limited-memory BFGS).
// Symmetry: gradient projection and position symmetrization using the
// ops passed in AtomContext.

#include "atom_minimizer.hpp"
#include "lattice_common.hpp"    // for invert_3x3

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

#include "json.hpp"

namespace pwdft {

using json = nlohmann::json;

namespace {

// ---------------------------------------------------------------------------
// RTDB helpers
// ---------------------------------------------------------------------------

std::string active_geometry_name(const json& rtdb)
{
   return (rtdb.contains("geometry") && rtdb["geometry"].is_string())
        ? rtdb["geometry"].get<std::string>()
        : "geometry";
}

std::array<double, 9> read_unita(const std::string& rtdbstring)
{
   const json rtdb = json::parse(rtdbstring);
   const std::string geomname = active_geometry_name(rtdb);
   const auto& u = rtdb.at("geometries").at(geomname).at("unita");
   std::array<double, 9> A{};

   for (int i=0; i<9; ++i) 
      A[i] = u.at(i).get<double>();

   return A;
}

std::vector<double> read_coords_cart(const std::string& rtdbstring)
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

std::vector<double> read_gradient_cart(const json& result)
{
   if (!result.contains("gradient")) 
      return {};

   const auto& g = result.at("gradient");

   if (!g.is_array()) 
      return {};

   std::vector<double> out;
   out.reserve(g.size());
   for (const auto& v : g) 

   {
      if (!v.is_number()) 
         return {};
      out.push_back(v.get<double>());
   }

   return out;
}

void write_coords_cart(std::string& rtdbstring,
                       const std::vector<double>& coords_cart)
{
   json rtdb = json::parse(rtdbstring);
   const std::string geomname = active_geometry_name(rtdb);
   json& geometry = rtdb["geometries"][geomname];

   geometry["coords"] = json::array();
   for (double x : coords_cart) 
      geometry["coords"].push_back(x);

   rtdbstring = rtdb.dump();
}

// ---------------------------------------------------------------------------
// Coordinate transforms
//
// Convention: unita is row-major 3x3, rows are lattice vectors:
//     a_row = (A[3*row+0], A[3*row+1], A[3*row+2])
// Cartesian r and fractional f relate by r_col = sum_row f_row * A[row][col].
// So r = f * A  (row-vector times matrix).
// ---------------------------------------------------------------------------

FracCoord cart_to_frac(const std::array<double, 3>& r,
                       const std::array<double, 9>& A_inv)
{
   return {r[0]*A_inv[0] + r[1]*A_inv[3] + r[2]*A_inv[6],
           r[0]*A_inv[1] + r[1]*A_inv[4] + r[2]*A_inv[7],
           r[0]*A_inv[2] + r[1]*A_inv[5] + r[2]*A_inv[8]
   };
}

std::array<double, 3> frac_to_cart(const FracCoord& f,
                                   const std::array<double, 9>& A)
{
   return {f[0]*A[0] + f[1]*A[3] + f[2]*A[6],
           f[0]*A[1] + f[1]*A[4] + f[2]*A[7],
           f[0]*A[2] + f[1]*A[5] + f[2]*A[8]
   };
}

FracCoord grad_cart_to_frac(const std::array<double, 3>& gc,
                            const std::array<double, 9>& A)
{
   return {A[0]*gc[0] + A[1]*gc[1] + A[2]*gc[2],
           A[3]*gc[0] + A[4]*gc[1] + A[5]*gc[2],
           A[6]*gc[0] + A[7]*gc[1] + A[8]*gc[2]
   };
}

// ---------------------------------------------------------------------------
// Symmetry projection
// ---------------------------------------------------------------------------

/***************************************
 *                                     *
 *      project_gradient_symmetry      *
 *                                     *
 ***************************************/
void project_gradient_symmetry(std::vector<FracCoord>& g,
                               const AtomPermutation& perm,
                               const std::vector<FracSymOp>& ops)
{
   if (!perm.valid() || ops.empty()) return;

   const int n_atoms = perm.n_atoms();
   const int n_ops   = perm.n_ops();
   if (n_atoms == 0 || n_ops == 0) return;

   std::vector<FracCoord> g_out(n_atoms, FracCoord{0,0,0});

   for (int op = 0; op < n_ops; ++op)
   {
      const FracSymOp& s = ops[op];
      for (int i=0; i<n_atoms; ++i)
      {
         const int j = perm.perm[op][i];
         const FracCoord rotated = s.apply_transpose(g[j]);
         for (int k=0; k<3; ++k) 
            g_out[i][k] += rotated[k];
      }
   }

   const double inv_nops = 1.0/static_cast<double>(n_ops);
   for (int i = 0; i < n_atoms; ++i)
      for (int k = 0; k < 3; ++k)
         g[i][k] = g_out[i][k] * inv_nops;
}


/***************************************
 *                                     *
 *         symmetrize_positions        *
 *                                     *
 ***************************************/
void symmetrize_positions(std::vector<FracCoord>& f, const AtomPermutation& perm, const std::vector<FracSymOp>& ops)
{
   if (!perm.valid() || ops.empty()) return;

   const int n_atoms = perm.n_atoms();
   const int n_ops   = perm.n_ops();
   if (n_atoms == 0 || n_ops == 0) return;

   std::vector<FracCoord> f_out(n_atoms, FracCoord{0,0,0});

   for (int op=0; op<n_ops; ++op)
   {
      const FracSymOp& s = ops[op];
      for (int i=0; i<n_atoms; ++i)
      {
         const int j = perm.perm[op][i];
         const FracCoord image = wrap_frac(s.apply(f[j]));
         for (int k=0; k<3; ++k) 
         {
            // Round to the nearest integer cell translation step
            double diff = image[k] - f[i][k];

            if (std::abs(diff) < 1.0e-2)
               f_out[i][k] += image[k];
            else
               f_out[i][k] += f[i][k];
         }
      }
   }

   const double inv_nops = 1.0/static_cast<double>(n_ops);
   for (int i=0; i<n_atoms; ++i)
      for (int k=0; k<3; ++k)
         f[i][k] = f_out[i][k] * inv_nops;
}


// ---------------------------------------------------------------------------
// L-BFGS engine
// ---------------------------------------------------------------------------

struct LBFGS {
    int N;
    int m;
    int n_stored  = 0;
    int next_slot = 0;

    std::vector<std::vector<double>> s_hist;
    std::vector<std::vector<double>> y_hist;
    std::vector<double> rho_hist;

    LBFGS(int N_, int m_) : N(N_), m(m_)
    {
       s_hist.resize(m, std::vector<double>(N, 0.0));
       y_hist.resize(m, std::vector<double>(N, 0.0));
       rho_hist.resize(m, 0.0);
    }

    void update(const std::vector<double>& s, const std::vector<double>& y)
    {
       double ys = 0.0;
       for (int i=0; i<N; ++i) 
          ys += y[i] * s[i];
       if (!(ys > 1.0e-16)) return;

       s_hist[next_slot] = s;
       y_hist[next_slot] = y;
       rho_hist[next_slot] = 1.0/ys;

       next_slot = (next_slot + 1) % m;
       if (n_stored < m) ++n_stored;
    }

    // Two-loop recursion: returns H * g where H is the current inverse
    // Hessian approximation. Step direction is -H*g.
    std::vector<double> apply(const std::vector<double>& g) const
    {
        std::vector<double> q = g;
        std::vector<double> alpha(n_stored, 0.0);

        // First loop: newest to oldest
        for (int k = 0; k < n_stored; ++k)
        {
            const int idx = (next_slot - 1 - k + m) % m;
            double sq = 0.0;
            for (int i = 0; i < N; ++i) sq += s_hist[idx][i] * q[i];
            alpha[k] = rho_hist[idx] * sq;
            for (int i = 0; i < N; ++i) q[i] -= alpha[k] * y_hist[idx][i];
        }

        // Initial Hessian scaling H0 = gamma * I
        double gamma = 1.0;
        if (n_stored > 0)
        {
           const int last = (next_slot - 1 + m) % m;
           double yy = 0.0, sy = 0.0;
           for (int i=0; i<N; ++i) 
           {
              yy += y_hist[last][i] * y_hist[last][i];
              sy += s_hist[last][i] * y_hist[last][i];
           }
           if (yy > 1.0e-16) gamma = sy / yy;
        }
        for (int i=0; i<N; ++i) 
           q[i] *= gamma;

        // Second loop: oldest to newest
        for (int k=n_stored-1; k>=0; --k)
        {
           const int idx = (next_slot - 1 - k + m) % m;
           double yq = 0.0;

           for (int i=0; i<N; ++i) 
              yq += y_hist[idx][i] * q[i];

           const double beta = rho_hist[idx] * yq;
           const double coef = alpha[k] - beta;

           for (int i=0; i<N; ++i) 
              q[i] += coef * s_hist[idx][i];
        }

        return q;
    }
};

} // namespace

// ---------------------------------------------------------------------------
// compute_atom_permutation
// ---------------------------------------------------------------------------

AtomPermutation compute_atom_permutation(const std::vector<FracCoord>& coords_frac, const std::vector<FracSymOp>& ops, double tol)
{
   const int n_atoms = static_cast<int>(coords_frac.size());
   const int n_ops   = static_cast<int>(ops.size());

   if (n_atoms == 0)
      throw std::runtime_error("compute_atom_permutation: empty coordinate list");

   if (n_ops == 0)
      throw std::runtime_error("compute_atom_permutation: empty symmetry op list");

   std::vector<FracCoord> f_ref(n_atoms);
   for (int i=0; i<n_atoms; ++i)
      f_ref[i] = wrap_frac(coords_frac[i]);

   const double tol_sq = tol * tol;

   AtomPermutation result;
   result.perm.assign(n_ops, std::vector<int>(n_atoms, -1));

   for (int g=0; g<n_ops; ++g)
   {
      const FracSymOp& op = ops[g];

      for (int i = 0; i < n_atoms; ++i)
      {
         const FracCoord f_image = wrap_frac(op.apply(f_ref[i]));

         int    best_j  = -1;
         double best_d2 = std::numeric_limits<double>::max();

         for (int j = 0; j < n_atoms; ++j)
         {
            const double d2 = frac_distance_sq(f_image, f_ref[j]);
            if (d2 < best_d2) {
                best_d2 = d2;
                best_j  = j;
            }
         }

         if (best_d2 > tol_sq)
         {
            std::ostringstream oss;
            oss << "compute_atom_permutation: op " << g
                << " maps atom " << i
                << " to no atom within tolerance " << tol
                << " (nearest distance = " << std::sqrt(best_d2) << ")";
            throw std::runtime_error(oss.str());
         }

         result.perm[g][i] = best_j;
      }
   }

   return result;
}

// ---------------------------------------------------------------------------
// atom_minimizer
// ---------------------------------------------------------------------------

int atom_minimizer(MPI_Comm comm,
                   std::string& rtdbstring,
                   std::ostream& coutput,
                   electronic_minimizer minimizer,
                   const AtomContext& ctx)
{
   coutput << "Im in atom_minimizer A" << std::endl;

   const bool oprint      = ctx.oprint;
   const std::string& tag = ctx.tag;

   const int    max_steps        = ctx.max_steps;
   const int    lbfgs_memory     = std::max(3, ctx.lbfgs_memory);
   const double minimum_gradient = ctx.minimum_gradient;

   double trust_radius = std::abs(ctx.initial_step);
   if (!std::isfinite(trust_radius) || trust_radius <= 0.0)
      trust_radius = 1.0e-1;

   const double trust_upper_bound = 10.0 * trust_radius;

   // --- read cell and coords ---
   const std::array<double, 9> A = read_unita(rtdbstring);
   std::array<double, 9> A_inv{};
   if (!invert_3x3(A, A_inv))
   {
      coutput << tag << " atom_minimizer: singular cell; aborting.\n";
      return 1;
   }

   const std::vector<double> coords_cart_0 = read_coords_cart(rtdbstring);
   if (coords_cart_0.empty() || coords_cart_0.size() % 3 != 0)
   {
      coutput << tag << " atom_minimizer: malformed coords; aborting.\n";
      return 1;
   }

   const int n_atoms = static_cast<int>(coords_cart_0.size() / 3);
   const int N       = 3*n_atoms;

   // --- initial fractional positions ---
   std::vector<FracCoord> f(n_atoms);
   for (int i=0; i<n_atoms; ++i)
   {
      const std::array<double, 3> r{coords_cart_0[3*i+0], coords_cart_0[3*i+1], coords_cart_0[3*i+2] };
      f[i] = cart_to_frac(r, A_inv);
   }

   // --- permutation if symmetry is active ---
   AtomPermutation perm;
   const bool want_symmetry = ctx.use_symmetry && !ctx.ops.empty();

   if (want_symmetry)
   {
      try 
      {
         perm = compute_atom_permutation(f, ctx.ops);
      }
      catch (const std::exception& e) 
      {
         coutput << tag << " atom_minimizer: " << e.what() << '\n'
                 << tag << " atom_minimizer: continuing without symmetry.\n";
      }
   }

   const bool symmetry_active = want_symmetry && perm.valid();

   if (oprint)
   {
      coutput << '\n'
              << tag << "==============================================\n"
              << tag << " PWDFT atom optimization\n"
              << tag << " L-BFGS in fractional coordinates\n"
              << tag << " n_atoms = " << n_atoms
              << "  dimension = " << N
              << "  memory = "     << lbfgs_memory << '\n'
              << tag << " symmetry = "
              << (symmetry_active ? "on" : "off")
              << " (ops read = " << ctx.ops.size() << ")\n"
              << tag << "==============================================\n";

   }

   LBFGS lbfgs(N, lbfgs_memory);

   std::vector<double> x_prev(N, 0.0);
   std::vector<double> g_prev(N, 0.0);
   bool have_previous = false;

   bool converged   = false;
   int  steps_taken = 0;

   double energy0 = 0.0;
   for (int istep = 0; istep < max_steps; ++istep)
   {
      // --- fractional -> cartesian, flatten into vectors ---
      std::vector<double> x(N);
      std::vector<double> coords_cart(3 * n_atoms);
      for (int i = 0; i < n_atoms; ++i)
      {
         for (int k=0; k<3; ++k) 
            x[3*i + k] = f[i][k];

         const auto r = frac_to_cart(f[i], A);
         coords_cart[3*i + 0] = r[0];
         coords_cart[3*i + 1] = r[1];
         coords_cart[3*i + 2] = r[2];
      }

      // --- write and evaluate ---
      std::string rtdb_eval = rtdbstring;
      write_coords_cart(rtdb_eval, coords_cart);

      const json result = compute_egs_values(2, comm, minimizer, rtdb_eval, coutput);
      const double energy = result.at("energy").get<double>();

      const std::vector<double> g_cart = read_gradient_cart(result);
      if (static_cast<int>(g_cart.size()) != N)
      {
         coutput << tag << " atom_minimizer: gradient size mismatch (" << g_cart.size() << " vs " << N << "); aborting.\n";
         return 1;
      }

      // --- cartesian gradient -> fractional ---
      std::vector<FracCoord> g_frac(n_atoms);
      for (int i = 0; i < n_atoms; ++i)
      {
         const std::array<double, 3> gc{ g_cart[3*i + 0], g_cart[3*i + 1], g_cart[3*i + 2] };

         // NOTE: if RTDB "gradient" is -dE/dr (a force), negate here.
         g_frac[i] = grad_cart_to_frac(gc, A);
      }

      if (symmetry_active)
         project_gradient_symmetry(g_frac, perm, ctx.ops);

      std::vector<double> g(N);
      for (int i = 0; i < n_atoms; ++i)
         for (int k = 0; k < 3; ++k)
            g[3*i + k] = g_frac[i][k];

      // --- diagnostics ---
      double gmax = 0.0, grms2 = 0.0;
      for (int i = 0; i < N; ++i)
      {
         gmax = std::max(gmax, std::abs(g[i]));
         grms2 += g[i] * g[i];
      }
      const double grms = std::sqrt(grms2) / static_cast<double>(n_atoms);

      if (oprint)
      {
         coutput << '\n'
                 << tag << "----------------------------------------------\n"
                 << tag << " Step          : " << istep << '\n'
                 << tag << " Energy        : " << std::fixed << std::setprecision(10)
                        << energy << " Hartree\n"
                 << tag << " Gmax          : "
                        << std::defaultfloat << std::setprecision(6)
                        << gmax << '\n'
                 << tag << " Grms          : " << grms << '\n'
                 << tag << " Trust radius  : " << trust_radius << '\n';
         
         // Refactored Compact Version
         //coutput << tag << " Step: " << istep 
         //        << " | E: " << std::fixed << std::setprecision(10) << energy 
         //        << " | Gmax: " << std::defaultfloat << std::setprecision(6) << gmax 
         //        << " | Grms: " << grms 
         //        << " | TR: " << trust_radius << '\n';


      }

      // --- convergence ---
      if (gmax < minimum_gradient)
      {
         converged = true;
         if (oprint)
            coutput << tag << " Action        : gradient converged.\n"
                    << tag << "----------------------------------------------\n";
         break;
      }

      // --- L-BFGS update from previous accepted step ---
      if (have_previous)
      {
         std::vector<double> s(N), y(N);
         for (int i = 0; i < N; ++i) 
         {
            s[i] = x[i] - x_prev[i];
            y[i] = g[i] - g_prev[i];
         }

         double sn2 = 0.0, yn2 = 0.0, gn2 = 0.0;
         for (int i = 0; i < N; ++i) 
         {
            sn2 += s[i]*s[i];
            yn2 += y[i]*y[i];
            gn2 += g[i]*g[i];
         }
         const double sn = std::sqrt(sn2);
         const double yn = std::sqrt(yn2);
         const double gn = std::sqrt(gn2);

         if (sn > 1.0e-8 && yn > 1.0e-3 * gn)
            lbfgs.update(s, y);
      }

      // --- propose step ---
      std::vector<double> dx(N);
      bool used_lbfgs = false;
      if (have_previous && lbfgs.n_stored > 0)
      {
         const std::vector<double> Hg = lbfgs.apply(g);
         for (int i=0; i<N; ++i) 
            dx[i] = -Hg[i];
         used_lbfgs = true;
      }
      else
      {
         for (int i=0; i<N; ++i) 
            dx[i] = -g[i];
      }

      // Descent check
      double g_dot_dx = 0.0;
      for (int i = 0; i < N; ++i) g_dot_dx += g[i] * dx[i];
      if (g_dot_dx >= 0.0)
      {
         for (int i=0; i<N; ++i) 
            dx[i] = -g[i];

         used_lbfgs = false;
      }

      // Trust-region clip
      double dxn2 = 0.0;
      for (int i=0; i<N; ++i) 
         dxn2 += dx[i] * dx[i];

      double dxn = std::sqrt(dxn2);
      if (dxn > trust_radius && dxn > 0.0)
      {
         const double sc = trust_radius / dxn;
         for (int i=0; i<N; ++i) 
            dx[i] *= sc;
         dxn = trust_radius;
      }
      if (dxn < std::numeric_limits<double>::epsilon())
      {
         if (oprint)
            coutput << tag << " Action        : zero step; stop.\n"
                    << tag << "----------------------------------------------\n";
         break;
      }

      // --- backtracking line search ---
      bool accepted = false;
      std::vector<double> accepted_dx;
      double accepted_energy = energy;
      std::vector<FracCoord> accepted_f;

      std::vector<double> trial_dx = dx;
      constexpr int max_backtracks = 12;


      for (int iback=0; iback<max_backtracks; ++iback)
      {
         double tn2 = 0.0;
         for (int i=0; i<N; ++i) 
            tn2 += trial_dx[i]*trial_dx[i];

         if (std::sqrt(tn2) < 1.0e-8) break;

         std::vector<FracCoord> f_trial(n_atoms);
         for (int i=0; i<n_atoms; ++i)
            for (int k=0; k<3; ++k)
               f_trial[i][k] = f[i][k] + trial_dx[3*i + k];

         if (symmetry_active)
            symmetrize_positions(f_trial, perm, ctx.ops);

         std::vector<double> coords_trial(3 * n_atoms);
         for (int i=0; i<n_atoms; ++i)
         {
            const auto r = frac_to_cart(f_trial[i], A);
            coords_trial[3*i + 0] = r[0];
            coords_trial[3*i + 1] = r[1];
            coords_trial[3*i + 2] = r[2];
         }

         std::string trial_rtdb = rtdbstring;
         write_coords_cart(trial_rtdb, coords_trial);

         const json trial_result = compute_egs_values(1, comm, minimizer, trial_rtdb, coutput);
         const double trial_energy = trial_result.at("energy").get<double>();

         if (std::isfinite(trial_energy) && trial_energy < energy)
         {
            accepted        = true;
            accepted_dx     = trial_dx;
            accepted_energy = trial_energy;
            accepted_f      = std::move(f_trial);
            rtdbstring      = std::move(trial_rtdb);
            break;
         }

         for (int i=0; i<N; ++i) 
            trial_dx[i] *= 0.5;
      }

      // --- accept or reject ---
      if (accepted)
      {
         x_prev = x;
         g_prev = g;
         have_previous = true;

         f = std::move(accepted_f);
         ++steps_taken;

         const std::vector<double> Hdx = lbfgs.apply(accepted_dx);
         double sHs = 0.0;
         for (int i = 0; i < N; ++i) sHs += accepted_dx[i] * Hdx[i];

         double g_dot_s = 0.0;
         for (int i=0; i<N; ++i) 
            g_dot_s += g[i] * accepted_dx[i];

         const double pred_red = -g_dot_s - 0.5 * sHs;
         const double act_red  = energy - accepted_energy;

         double rho = 0.0;
         if (pred_red > 1.0e-14) 
            rho = act_red / pred_red;

         double an2 = 0.0;
         for (int i=0; i<N; ++i) 
            an2 += accepted_dx[i]*accepted_dx[i];

         const double an = std::sqrt(an2);

         if (rho < 0.25)
            trust_radius *= 0.5;
         else if (rho > 0.75 && an > 0.9 * trust_radius)
            trust_radius = std::min(2.0 * trust_radius, trust_upper_bound);

         if (oprint)
            coutput << tag << " Method        : " << (used_lbfgs ? "L-BFGS" : "gradient fallback") << '\n'
                    << tag << " Action        : accepted\n"
                    << tag << " rho           : " << rho << '\n'
                    << tag << " Trust radius  : " << trust_radius << '\n'
                    << tag << " New energy    : " << std::fixed << std::setprecision(10) << accepted_energy << " Hartree\n"
                    << tag << " Delta energy  : " << std::fixed << std::setprecision(10) << accepted_energy - energy << " Hartree\n"
                    << tag << "----------------------------------------------\n";
         /*
         if (oprint) {
            coutput << "@@ Step " << istep 
                    << " | Energy: " << std::scientific << energy 
                    << " | Method: " << (used_lbfgs ? "L-BFGS" : "gradient fallback")
                    << " | Action: " << "accepted"
                    << " | rho: " << rho 
                    << " | Trust radius: " << trust_radius 
                    << " | New energy: " << std::fixed << std::setprecision(10) << accepted_energy << " Hartree\n";
        }
        */

      }
      else
      {
         trust_radius *= 0.5;
         if (oprint)
            coutput << tag << " Method        : " << (used_lbfgs ? "L-BFGS" : "gradient fallback") << '\n'
                    << tag << " Action        : rejected; halve trust radius.\n"
                    << tag << " New radius    : " << std::defaultfloat << std::setprecision(10) << trust_radius << '\n'
                    << tag << "----------------------------------------------\n";

         /*
         if (oprint) {
            coutput << "@@ Step " << istep 
                    << " | Energy: " << std::scientific << energy 
                    << " | Method: " << (used_lbfgs ? "L-BFGS" : "gradient fallback")
                    << " | Action: " << "rejected; halve trust radius"
                    << " | Trust radius: " << trust_radius 
                    << " | New energy: " << std::fixed << std::setprecision(10) << accepted_energy << " Hartree\n";
         }
         */

         if (trust_radius < 1.0e-8) break;
      }
   }

   // --- final write-back ---
   std::vector<double> coords_final(3 * n_atoms);
   for (int i=0; i<n_atoms; ++i)
   {
      const auto r = frac_to_cart(f[i], A);
      coords_final[3*i + 0] = r[0];
      coords_final[3*i + 1] = r[1];
      coords_final[3*i + 2] = r[2];
   }
   write_coords_cart(rtdbstring, coords_final);

   if (oprint)
   {
      coutput << '\n'
              << tag << "==============================================\n"
              << tag << " PWDFT atom optimization COMPLETE\n"
              << tag << "==============================================\n"
              << tag << " Accepted steps : " << steps_taken << '\n'
              << tag << " Status         : " << (converged ? "Converged" : "Stopped before gradient convergence") << '\n'
              << tag << "==============================================\n";
   }

   return 0;
}

} // namespace pwdft

// alternating_optimizer.cpp
//
// Alternating atom + lattice optimizer.
//
// Each outer iteration:
//   1. atom_minimizer   (L-BFGS, fractional coords, symmetry projected)
//   2. lattice minimizer (per-system BFGS on log/angle coordinates)
//   3. energy check
//
// The two phases use the same electronic minimizer callback and the same
// RTDB, so no state is lost between them. Symmetry ops are re-derived by
// the atom phase from the current RTDB each call, so a cell change in
// phase 2 is picked up correctly by the next atom phase.

#include "alternating_optimizer.hpp"
#include "lattice_common.hpp"    // for compute_egs_values

#include <cmath>
#include <iomanip>
#include <ostream>
#include <string>

#include <mpi.h>

#include "json.hpp"

namespace pwdft {

using json = nlohmann::json;

namespace {

double evaluate_energy(MPI_Comm comm,
                       electronic_minimizer minimizer,
                       std::string& rtdbstring,
                       std::ostream& coutput)
{
    const json r = compute_egs_values(1, comm, minimizer, rtdbstring, coutput);
    return r.at("energy").get<double>();
}

} // namespace

int alternating_optimizer(MPI_Comm comm,
                          std::string& rtdbstring,
                          std::ostream& coutput,
                          electronic_minimizer minimizer,
                          const AtomContext& atom_ctx,
                          const LatticeContext& lattice_ctx,
                          lattice_minimizer lm,
                          int    max_outer,
                          double energy_tol)
{
    if (lm == nullptr)
    {
        coutput << "alternating_optimizer: null lattice minimizer\n";
        return 1;
    }

    const bool oprint = atom_ctx.oprint && lattice_ctx.oprint;
    const std::string& tag = atom_ctx.tag;

    if (oprint)
    {
        coutput << '\n'
                << tag << "==============================================\n"
                << tag << " PWDFT alternating geometry + lattice optimization\n"
                << tag << " Alternating atom and cell phases\n"
                << tag << " outer max = " << max_outer
                << "  E tol = " << energy_tol << '\n'
                << tag << "==============================================\n";
    }

    // Baseline energy before any phase runs.
    double E_prev = evaluate_energy(comm, minimizer, rtdbstring, coutput);

    bool converged = false;
    int  outer     = 0;

    for (outer = 0; outer < max_outer; ++outer)
    {
        if (oprint)
        {
            coutput << '\n'
                    << tag << "----------------------------------------------\n"
                    << tag << " Outer iteration " << outer << '\n'
                    << tag << " Starting energy : " << std::fixed
                           << std::setprecision(10) << E_prev
                           << " Hartree\n"
                    << tag << "----------------------------------------------\n";
        }

        // --- atom phase ---
        if (oprint)
            coutput << '\n'
                    << tag << ">>>>> Atom phase (outer " << outer << ")\n";

        const int ierr_atom = atom_minimizer(comm, rtdbstring, coutput,
                                             minimizer, atom_ctx);
        if (ierr_atom != 0)
        {
            coutput << tag << " combined_optimizer: atom phase failed ierr="
                    << ierr_atom << '\n';
            return ierr_atom;
        }

        // --- cell phase ---
        if (oprint)
            coutput << '\n'
                    << tag << ">>>>> Cell phase (outer " << outer << ")\n";

        const int ierr_cell = lm(comm, rtdbstring, coutput,
                                 minimizer, lattice_ctx);
        if (ierr_cell != 0)
        {
            coutput << tag << " combined_optimizer: cell phase failed ierr="
                    << ierr_cell << '\n';
            return ierr_cell;
        }

        // --- convergence check ---
        const double E_now = evaluate_energy(comm, minimizer, rtdbstring, coutput);
        const double dE    = std::abs(E_now - E_prev);

        if (oprint)
        {
            coutput << '\n'
                    << tag << " Outer iteration " << outer
                           << " complete.\n"
                    << tag << " Energy after phase : "
                           << std::fixed << std::setprecision(10)
                           << E_now << " Hartree\n"
                    << tag << " |Delta E|          : "
                           << std::defaultfloat << std::setprecision(3)
                           << dE << " Hartree\n";
        }

        if (dE < energy_tol)
        {
            converged = true;
            E_prev = E_now;
            if (oprint)
                coutput << tag << " Action             : converged\n";
            break;
        }

        E_prev = E_now;
    }

    if (oprint)
    {
        coutput << '\n'
                << tag << "==============================================\n"
                << tag << " PWDFT alternating optimization COMPLETE\n"
                << tag << "==============================================\n"
                << tag << " Outer iterations : "
                       << (converged ? outer : outer + 1) << '\n'
                << tag << " Final energy     : "
                       << std::fixed << std::setprecision(10)
                       << E_prev << " Hartree\n"
                << tag << " Status           : "
                       << (converged ? "Converged"
                                     : "Stopped (max outer reached)")
                       << '\n'
                << tag << "==============================================\n";
    }

    return 0;
}

} // namespace pwdft

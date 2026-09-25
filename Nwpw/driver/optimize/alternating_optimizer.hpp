// alternating_optimizer.hpp
#pragma once

#include <mpi.h>
#include <ostream>
#include <string>

#include "atom_minimizer.hpp"
#include "lattice_minimizer.hpp"

namespace pwdft {

// Alternating atom/lattice relaxation.
//
// Runs atom_minimizer to convergence at fixed cell, then the dispatched
// lattice minimizer to convergence at fixed atoms, and repeats until the
// energy change between consecutive outer iterations falls below
// energy_tol, or max_outer is reached.
//
// The two phases share nothing but the RTDB string; each drives the same
// electronic minimizer callback. This mirrors the NWChem driver pattern
// of alternating geometry and cell optimization.
//
// Returns 0 on successful completion (including non-converged runs that
// hit max_outer), non-zero if either phase reports a fatal error.

int alternating_optimizer(MPI_Comm comm,
                          std::string& rtdbstring,
                          std::ostream& coutput,
                          electronic_minimizer minimizer,
                          const AtomContext& atom_ctx,
                          const LatticeContext& lattice_ctx,
                          lattice_minimizer lm,
                          int    max_outer  = 10,
                          double energy_tol = 1.0e-5);

} // namespace pwdft

// atom_minimizer.hpp
#pragma once

#include <array>
#include <mpi.h>
#include <ostream>
#include <string>
#include <vector>
#include <cmath>

#include "lattice_minimizer.hpp"   // for FracSymOp, SymmetryInfo

namespace pwdft {

using FracCoord = std::array<double, 3>;

struct AtomContext {
    bool        oprint = false;
    std::string tag    = "@@";

    int         max_steps        = 50;
    int         lbfgs_memory     = 10;
    double      minimum_gradient = 1.0e-3;
    double      initial_step     = 0.1;

    bool        use_symmetry     = true;

    // Symmetry ops for projection. Populated by the driver from
    // symmetry_info.ops. Empty list (or use_symmetry=false) means no
    // projection is applied.
    std::vector<FracSymOp> ops;
};

struct AtomContext {
    bool        oprint = false;
    std::string tag    = "@@";

    int         max_steps        = 50;
    int         lbfgs_memory     = 10;
    double      minimum_gradient = 1.0e-3;
    double      initial_step     = 0.1;

    bool        use_symmetry     = true;
};

// For each (op, atom), the index of the atom that the op maps it to.
// perm[op][atom] = permuted atom index.
struct AtomPermutation {
    std::vector<std::vector<int>> perm;

    int n_ops()  const { return static_cast<int>(perm.size()); }
    int n_atoms() const {
        return perm.empty() ? 0 : static_cast<int>(perm[0].size());
    }

    bool valid() const { return !perm.empty() && !perm[0].empty(); }
};

// Compute the atom permutation for a given fractional coordinate set
// and list of symmetry operations. Throws std::runtime_error if any
// op/atom pair cannot be matched within tolerance.
AtomPermutation compute_atom_permutation(
    const std::vector<FracCoord>& coords_frac,
    const std::vector<FracSymOp>& ops,
    double tol = 1.0e-4);

// Wrap a fractional coordinate into [0, 1).
inline FracCoord wrap_frac(const FracCoord& f) {
    FracCoord r;
    for (int i = 0; i < 3; ++i) {
        r[i] = f[i] - std::floor(f[i]);
    }
    return r;
}

// Minimum-image squared distance between two fractional coordinates.
// Both inputs should be wrapped into [0, 1).
inline double frac_distance_sq(const FracCoord& a, const FracCoord& b) {
    double d2 = 0.0;
    for (int i = 0; i < 3; ++i) {
        double d = a[i] - b[i];
        d -= std::round(d);            // minimum image
        d2 += d * d;
    }
    return d2;
}

// Atom optimizer entry point. Same signature shape as the lattice
// minimizers: takes the RTDB by reference, returns 0 on success.
int atom_minimizer(MPI_Comm comm,
                   std::string& rtdbstring,
                   std::ostream& coutput,
                   electronic_minimizer minimizer,
                   const AtomContext& ctx);

} // namespace pwdft

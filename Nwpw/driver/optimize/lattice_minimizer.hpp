// lattice_minimizer.hpp
#pragma once

#include <mpi.h>
#include <ostream>
#include <string>
#include <array>
#include <vector>

namespace pwdft {

// ---------------------------------------------------------------------------
// Minimizer callback types
//
// Two distinct signatures, deliberately not unified:
//
//   electronic_minimizer — the PSPW/band minimizer. Takes (comm, rtdb, out),
//                          returns 0 on success. This is what the driver
//                          receives and what the lattice minimizers call
//                          (indirectly, via compute_egs_values).
//
//   lattice_minimizer    — a cell-shape optimizer. Takes the same three
//                          arguments plus the electronic minimizer it should
//                          drive for energy/stress evaluations.
// ---------------------------------------------------------------------------

using electronic_minimizer = int (*)(MPI_Comm, std::string&, std::ostream&);


// ---------------------------------------------------------------------------
// LatticeContext
//
// Per-call context passed from the driver to whichever lattice minimizer
// pick_lattice_minimizer selects. Must be defined before the lattice_minimizer
// typedef below, since the typedef names it.
// ---------------------------------------------------------------------------

struct LatticeContext {
    bool        oprint = false;
    std::string tag    = "@@";

    int         max_steps        = 25;
    double      initial_step     = 0.0025;
    double      minimum_step     = 1.0e-5;
    double      minimum_gradient = 1.0e-4;
};

using lattice_minimizer    = int (*)(MPI_Comm,
                                     std::string&,
                                     std::ostream&,
                                     electronic_minimizer,
                                     const LatticeContext&);



// ---------------------------------------------------------------------------
// FracSymOp
//
// One symmetry operation of the effective space group, in the same form the
// RTDB "effective_symmetry.ops" block uses:
//
//     f' = R * f + t
//
// where f is a fractional coordinate 3-vector, R is a 3x3 integer-valued
// rotation (stored row-major), and t is a fractional translation.
//
// This is the primitive action on fractional coordinates. It's all the atom
// optimizer needs for projecting gradients and symmetrizing positions.
// ---------------------------------------------------------------------------

struct FracSymOp {
    std::array<double, 9> R{};   // row-major 3x3
    std::array<double, 3> t{};   // fractional translation

    // Apply to a fractional coordinate: f' = R * f + t
    std::array<double, 3> apply(const std::array<double, 3>& f) const
    {
        return {
            R[0]*f[0] + R[1]*f[1] + R[2]*f[2] + t[0],
            R[3]*f[0] + R[4]*f[1] + R[5]*f[2] + t[1],
            R[6]*f[0] + R[7]*f[1] + R[8]*f[2] + t[2]
        };
    }

    // Apply the transpose: g' = R^T * g   (used for projecting gradients,
    // since the gradient transforms with the inverse rotation).
    std::array<double, 3> apply_transpose(const std::array<double, 3>& g) const
    {
        return {
            R[0]*g[0] + R[3]*g[1] + R[6]*g[2],
            R[1]*g[0] + R[4]*g[1] + R[7]*g[2],
            R[2]*g[0] + R[5]*g[1] + R[8]*g[2]
        };
    }

    bool is_identity(double tol = 1.0e-8) const
    {
        for (int i = 0; i < 9; ++i)
            if (std::abs(R[i] - ((i % 4 == 0) ? 1.0 : 0.0)) > tol)
                return false;
        for (int i = 0; i < 3; ++i)
            if (std::abs(t[i]) > tol)
                return false;
        return true;
    }
};



// ---------------------------------------------------------------------------
// SymmetryInfo
//
// Plain data struct describing the effective symmetry of the current cell.
// Populated once by the driver from the RTDB "effective_symmetry" block, then
// passed to pick_lattice_minimizer to select a strategy, and to the atom
// optimizer for symmetry projection.
//
// Members mirror the JSON keys:
//   name       -> space_group_name
//   type       -> type
//   order      -> group_order       (number of symmetry operations)
//   primitive  -> is_primitive
//   ita_number -> ita_number        (1..230; -1 if absent)
//   ops        -> ops               (array of {R, t})
//
// "system" is derived here from ita_number and is the field the dispatch
// actually switches on.
//
// "ops" is empty when no symmetry block is present in the RTDB. When
// non-empty, it contains one entry per operation of the effective space
// group, in the same order the SCF used them. The atom optimizer uses this
// list for gradient projection and position symmetrization.
// ---------------------------------------------------------------------------

struct SymmetryInfo {
    std::string space_group_name = "unknown";
    std::string type             = "unknown";
    int         ita_number       = -1;
    int         group_order      = -1;
    bool        is_primitive     = false;
    bool        is_cubic         = false;
    std::string system           = "unknown"; // triclinic..cubic, or "unknown"

    std::vector<FracSymOp> ops;

    bool has_symmetry() const {
        return (space_group_name != "unknown" && group_order > 0);
    }

    // Whether we have enough information to project with symmetry.
    bool has_ops() const {
        return !ops.empty();
    }
};



// ---------------------------------------------------------------------------
// Per-system lattice minimizer entry points
//
// Each has the lattice_minimizer signature. Declare one line here when you
// add a new .cpp implementation, and add a branch in pick_lattice_minimizer.
// Do not declare a symbol here until its definition exists, or you will trade
// a compile error for a link error.
// ---------------------------------------------------------------------------

int general_lattice_minimizer(MPI_Comm,
                              std::string&,
                              std::ostream&,
                              electronic_minimizer,
                              const LatticeContext&);

int cubic_lattice_minimizer(MPI_Comm,
                            std::string&,
                            std::ostream&,
                            electronic_minimizer,
                            const LatticeContext&);

int tetragonal_lattice_minimizer(MPI_Comm,        // <-- add this
                                 std::string&,
                                 std::ostream&,
                                 electronic_minimizer,
                                 const LatticeContext&);

int hexagonal_lattice_minimizer(MPI_Comm,        // <-- add this
                                std::string&,
                                std::ostream&,
                                electronic_minimizer,
                                const LatticeContext&);

int orthorhombic_lattice_minimizer(MPI_Comm,
                                   std::string&,
                                   std::ostream&,
                                   electronic_minimizer,
                                   const LatticeContext&);

int monoclinic_lattice_minimizer(MPI_Comm,
                                 std::string&,
                                 std::ostream&,
                                 electronic_minimizer,
                                 const LatticeContext&);

int trigonal_lattice_minimizer(MPI_Comm, std::string&, 
                               std::ostream&,
                               electronic_minimizer, 
                               const LatticeContext&);

int triclinic_lattice_minimizer(MPI_Comm,
                                std::string&,
                                std::ostream&,
                                electronic_minimizer,
                                const LatticeContext&);


// ---------------------------------------------------------------------------
// Dispatch
//
// Returns a lattice_minimizer based on the crystal-system string produced by
// the driver from SymmetryInfo::ita_number. Falls back to
// general_lattice_minimizer when the system is unknown or unrecognized.
// Never returns nullptr.
// ---------------------------------------------------------------------------

lattice_minimizer pick_lattice_minimizer(const std::string& system);

} // namespace pwdft

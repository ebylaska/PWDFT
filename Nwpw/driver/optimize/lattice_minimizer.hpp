// lattice_minimizer.hpp
#pragma once

#include <mpi.h>
#include <ostream>
#include <string>

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
// SymmetryInfo
//
// Plain data struct describing the effective symmetry of the current cell.
// Populated once by the driver from the RTDB "effective_symmetry" block, then
// passed to pick_lattice_minimizer to select a strategy.
//
// Members mirror the JSON keys:
//   name       -> space_group_name
//   type       -> type
//   order      -> group_order       (number of symmetry operations)
//   primitive  -> is_primitive
//   ita_number -> ita_number        (1..230; -1 if absent)
//
// "system" is derived here from ita_number and is the field the dispatch
// actually switches on.
// ---------------------------------------------------------------------------

struct SymmetryInfo {
    std::string space_group_name = "unknown";
    std::string type             = "unknown";
    int         ita_number       = -1;
    int         group_order      = -1;
    bool        is_primitive     = false;
    bool        is_cubic         = false;
    std::string system           = "unknown"; // triclinic..cubic, or "unknown"

    bool has_symmetry() const {
        return (space_group_name != "unknown" && group_order > 0);
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

// Declare these as you implement them:
// int tetragonal_lattice_minimizer(...);
// int orthorhombic_lattice_minimizer(...);
// int hexagonal_lattice_minimizer(...);
// int trigonal_lattice_minimizer(...);
// int monoclinic_lattice_minimizer(...);
// int triclinic_lattice_minimizer(...);

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

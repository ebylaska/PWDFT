// general_lattice_minimizer.cpp
//
// Fallback lattice minimizer for cases where the crystal system is unknown
// or unrecognized by pick_lattice_minimizer. It performs a single energy
// evaluation at the current cell and returns, leaving the lattice unchanged.
//
// The intended long-term role of this function is to host a generic
// (triclinic) line search over all six cell parameters, so that any cell can
// be optimized even when the space group is not identified. Until that is
// implemented, returning after one evaluation is the safest behavior: the
// caller's post-dispatch code (energy extraction, final reporting) still sees
// a valid, consistent RTDB, and no incorrect lattice deformation is applied.
//
// Signature matches lattice_minimizer from lattice_minimizer.hpp:
//
//     int (MPI_Comm, std::string&, std::ostream&, electronic_minimizer)

#include "lattice_minimizer.hpp"

#include <ostream>
#include <string>

#include <mpi.h>

namespace pwdft {

int general_lattice_minimizer(MPI_Comm comm,
                              std::string& rtdbstring,
                              std::ostream& coutput,
                              electronic_minimizer minimizer,
                              const LatticeContext& ctx)
{
    coutput << "@@ general_lattice_minimizer: no crystal-system-specific "
               "optimizer available; skipping lattice optimization.\n";

    // Intentionally do not modify rtdbstring.
    //
    // When this becomes a real triclinic optimizer, the structure will be:
    //
    //   1. read a, b, c, alpha, beta, gamma from rtdbstring
    //   2. build a LatticeDoF (see lattice_common.hpp) for the 6 parameters
    //   3. run_lattice_line_search(dof, geomname, rtdbstring, comm, coutput,
    //                              tag, minimizer)
    //
    // The argument is named in the header; it is commented out here to
    // silence the unused-parameter warning while the body is a stub. Restore
    // the name when the body starts using it.

    return 0;
}

} // namespace pwdft

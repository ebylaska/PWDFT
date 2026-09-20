#include "lattice_minimizer.hpp"

namespace pwdft {

lattice_minimizer pick_lattice_minimizer(const std::string& system)
{
    if (system == "cubic") return cubic_lattice_minimizer;
    if (system == "tetragonal")  return tetragonal_lattice_minimizer;
    if (system == "hexagonal") return hexagonal_lattice_minimizer;
    if (system == "trigonal") return trigonal_lattice_minimizer;
    if (system == "orthorhombic") return orthorhombic_lattice_minimizer;
    if (system == "monoclinic") return monoclinic_lattice_minimizer;
    if (system == "triclinic") return triclinic_lattice_minimizer;
    // ... add as implemented
    return general_lattice_minimizer;
}

} // namespace pwdft

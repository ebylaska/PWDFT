#include "lattice_minimizer.hpp"

namespace pwdft {

lattice_minimizer pick_lattice_minimizer(const std::string& system)
{
    if (system == "cubic") return cubic_lattice_minimizer;
    if (system == "tetragonal")  return tetragonal_lattice_minimizer;
    if (system == "hexagonal") return hexagonal_lattice_minimizer;
    // ... add as implemented
    return general_lattice_minimizer;
}

} // namespace pwdft

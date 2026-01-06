#ifndef _APPLY_KPOINT_SYMMETRY_PRUNING_HPP_
#define _APPLY_KPOINT_SYMMETRY_PRUNING_HPP_


#include "json.hpp"
using json = nlohmann::json;

namespace pwdft {

/**
 * @brief Apply symmetry-based pruning to Monkhorst–Pack k-points.
 *
 * Reads effective symmetry information from the RTDB JSON and reduces
 * the full Monkhorst–Pack k-point set to the irreducible Brillouin zone
 * using point-group rotations and time-reversal symmetry.
 *
 * The routine:
 *  - Requires `effective_symmetry` to be present and authoritative
 *  - Operates purely in fractional k-space
 *  - Preserves total k-point weight exactly
 *  - Records both original and effective k-point counts
 *
 * This function does **not** recompute symmetry or lattice information;
 * it consumes symmetry that has already been resolved and frozen by
 * `resolve_symmetry_and_cell()`.
 *
 * If no symmetry is present, or if k-points are absent, the input RTDB
 * string is returned unchanged.
 *
 * @param rtdbstr Serialized RTDB JSON string.
 * @return Updated RTDB JSON string with pruned k-points and bookkeeping.
 *
 * @note
 *  - This implementation uses strict point-group + time-reversal symmetry.
 *  - Legacy NWChem behavior may produce fewer k-points; PWDFT preserves
 *    full k-star structure required for modern methods.
 */

std::string apply_kpoint_symmetry_pruning(std::string rtdbstr);

} // namespace pwdft


#endif // _APPLY_KPOINT_SYMMETRY_PRUNING_HPP_


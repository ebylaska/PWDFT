#ifndef _RESOLVE_SYMMETRY_CELL_HPP_
#define _RESOLVE_SYMMETRY_CELL_HPP_

#include "json.hpp"
using json = nlohmann::json;

namespace pwdft {

/**
 * Resolve simulation cell and symmetry into an effective, self-consistent form.
 *
 * - Applies precedence rules for simulation_cell vs geometry
 * - Handles autosym and explicit symmetry
 * - Performs conventional → primitive conversion if requested
 * - Writes results into rtdbjson["effective_symmetry"]
 * - Updates rtdbjson["nwpw"]["simulation_cell"]
 *
 * This function does NOT modify the user's geometry intent.
 */
std::string resolve_symmetry_and_cell(std::string rtdbstring);

} // namespace pwdft

#endif


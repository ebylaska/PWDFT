#ifndef _PWD_FT_UTIL_PATHS_HPP_
#define _PWD_FT_UTIL_PATHS_HPP_

#include <string>

namespace pwdft {

/**
 * Resolve the base NWPW data directory.
 *
 * Resolution order:
 *   1. $NWPW_DATADIR environment variable
 *   2. Installed data directory (CMAKE_INSTALL_PREFIX/share/nwpw)
 *   3. Source-tree fallback (developer builds)
 */
std::string resolve_datadir();

/**
 * Resolve the pseudopotential library directory.
 *
 * Resolution order:
 *   1. $NWPW_LIBRARY or $NWCHEM_NWPW_LIBRARY
 *   2. Installed data directory + "/libraryps"
 *   3. Source-tree fallback
 */
std::string resolve_libraryps();

/**
 * Resolve the van der Waals kernel directory.
 */
std::string resolve_vdw();

/**
 * Resolve the symmetry database directory.
 */
std::string resolve_symmetry();


/**
 * Resolve the spacegroups.dat dattabase.
 */
std::string resolve_spacegroups_dat();

} // namespace pwdft

#endif


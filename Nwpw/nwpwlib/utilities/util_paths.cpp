/** @defgroup nwpw_utilities Utility Infrastructure
 *  Core utility services used across NWPW/PWDFT.
 */
/**
 * @file util_paths.cpp
 * @brief Runtime data path resolution utilities for NWPW / PWDFT.
 *
 * This file provides centralized logic for resolving filesystem paths
 * to shared runtime data used by the NWPW/PWDFT codebase, including:
 *
 *  - Pseudopotential libraries
 *  - van der Waals kernel data
 *  - Symmetry databases (e.g., spacegroups.dat)
 *
 * Path resolution follows a strict and documented precedence:
 *
 *  1. Environment variable overrides (e.g., NWPW_DATADIR, NWPW_LIBRARY)
 *  2. Installed shared data directories (via CMake install layout)
 *  3. Source-tree fallbacks for developer builds
 *
 * This design allows:
 *
 *  - Relocatable binary installations
 *  - Clean separation of build-time vs runtime paths
 *  - Consistent behavior across developer, CI, and HPC deployments
 *
 * All filesystem queries are intentionally centralized here to avoid
 * hard-coded paths and duplicated logic elsewhere in the codebase.
 *
 * @note This file is intended to be dependency-light and safe to use
 *       from low-level components such as symmetry, pseudopotentials,
 *       and initialization routines.
 *
 * @author
 *   Eric J. Bylaska
 *
 * @ingroup nwpw_utilities
 */

#include "util_paths.hpp"

//#include "NwpwUtilConfig.hpp"
#include "NwpwConfig.h"

#include <cstdlib>
#include <filesystem>


/** @ingroup nwpw_utilities */
namespace pwdft {

/******************************************
 *                                        *
 *               exists                   *
 *                                        *
 ******************************************/
/**
 * @brief Check whether a filesystem path exists.
 *
 * Lightweight helper used throughout path-resolution logic to probe
 * for installed runtime data directories (e.g., pseudopotentials,
 * symmetry tables) before falling back to source-tree defaults.
 *
 * @param p Filesystem path to test.
 * @return true if the path exists, false otherwise.
 */
static bool exists(const std::string& p)
{
    return std::filesystem::exists(p);
}


/******************************************
 *                                        *
 *            resolve_datadir             *
 *                                        *
 ******************************************/
/**
 * @brief Resolve the root shared data directory for NWPW/PWDFT.
 *
 * Resolution order:
 * 1. Environment variable `NWPW_DATADIR`
 * 2. Installed shared data directory (from build configuration)
 * 3. Source-tree directory (developer fallback)
 *
 * This function defines the canonical runtime data root used by all
 * other path-resolution helpers.
 *
 * @return Absolute path to the resolved data directory.
 */
std::string resolve_datadir()
{
    if (const char* env = std::getenv("NWPW_DATADIR"))
        return env;

    if (exists(Nwpw_DATADIR))
        return Nwpw_DATADIR;

    return Nwpw_SOURCE_DIR;
}


/******************************************
 *                                        *
 *            resolve_libraryps           *
 *                                        *
 ******************************************/
/**
 * @brief Resolve the pseudopotential library directory.
 *
 * Resolution order:
 * 1. Environment variable `NWPW_LIBRARY`
 * 2. Environment variable `NWCHEM_NWPW_LIBRARY`
 * 3. Installed data directory (`$NWPW_DATADIR/libraryps`)
 * 4. Source-tree fallback
 *
 * @return Absolute path to the pseudopotential library.
 */
std::string resolve_libraryps()
{
    if (const char* env = std::getenv("NWPW_LIBRARY"))
        return env;

    if (const char* env = std::getenv("NWCHEM_NWPW_LIBRARY"))
        return env;

    auto installed = resolve_datadir() + "/libraryps";
    if (exists(installed))
        return installed;

    return Nwpw_LIBRARYPS_Default;
}

//std::string resolve_vdw()
//{
//    return resolve_libraryps() + "/VDW";
//}

/******************************************
 *                                        *
 *            resolve_vdw                 *
 *                                        *
 ******************************************/
/**
 * @brief Resolve the van der Waals (vdW) kernel data directory.
 *
 * Resolution order:
 * 1. Installed data directory (`$NWPW_DATADIR/vdw`)
 * 2. Source-tree fallback
 *
 * @return Absolute path to the vdW data directory.
 */
std::string resolve_vdw()
{
    auto installed = resolve_datadir() + "/vdw";
    if (exists(installed))
        return installed;

    return Nwpw_LIBRARYVDW_Default;
}


/******************************************
 *                                        *
 *            resolve_symmetry            *
 *                                        *
 ******************************************/
/**
 * @brief Resolve the symmetry data directory.
 *
 * Resolution order:
 * 1. Installed data directory (`$NWPW_DATADIR/symmetry`)
 * 2. Source-tree fallback
 *
 * @return Absolute path to the symmetry data directory.
 */
std::string resolve_symmetry()
{
    auto installed = resolve_datadir() + "/symmetry";
    if (exists(installed))
        return installed;

    return Nwpw_SYMMETRY_Default;
}


/******************************************
 *                                        *
 *        resolve_spacegroups_dat         *
 *                                        *
 ******************************************/
/**
 * @brief Resolve the path to the space group database file.
 *
 * Returns the fully qualified path to `spacegroups.dat`, using
 * the resolved symmetry data directory. This is the canonical
 * access point for symmetry database loading.
 *
 * @return Absolute path to `spacegroups.dat`.
 */
std::string resolve_spacegroups_dat()
{
    return resolve_symmetry() + "/spacegroups.dat";
}


} // namespace pwdft


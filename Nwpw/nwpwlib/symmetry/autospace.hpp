#ifndef _AUTOSPACE_HPP_
#define _AUTOSPACE_HPP_

#pragma once

#include <string>
#include <vector>


namespace pwdft {

/**
 * Result of space-group auto-detection.
 * This is a *proposal*, not yet authoritative.
 */
struct AutoSpaceResult
{
    bool success = false;

    std::string group_name;     // e.g. "Fm-3m"
    int group_number = 0;       // e.g. 225

    double tolerance_used = 0.0;

    bool primitive_suggested = false; // True if detector believes primitive cell is more appropriate; caller decides whether to actually reduce the cell.
    bool lattice_consistent = true;

    std::string method;         // "spglib", "internal", "none"
    std::string message;        // diagnostics / failure reason
};

/**
 * Attempt to detect a space group from lattice + atoms.
 *
 * Inputs are cartesian coordinates in the *current* cell.
 */
AutoSpaceResult detect_space_group(
    const double unita[9],
    const std::vector<std::string>& symbols,
    const std::vector<double>& coords_xyz,
    double tolerance
);

} // namespace pwdft


#endif

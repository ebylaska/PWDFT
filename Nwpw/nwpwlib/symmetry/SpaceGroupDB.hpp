#ifndef _SPACEGROUPDB_HPP_
#define _SPACEGROUPDB_HPP_

#pragma once
#include <string>
#include <vector>
#include <unordered_map>
#include <cstddef>
#include <mutex>

namespace pwdft {

/**
 * @brief Single symmetry operation (R | t)
 *
 * Represents an affine symmetry transformation:
 *   r' = R r + t
 */
struct SymOp {
    double R[3][3];   ///< Rotation matrix
    double t[3];      ///< Translation vector
};

/**
 * @brief Crystallographic space group definition
 */
struct SpaceGroup {
    std::string name;          ///< Hermann–Mauguin symbol
    int number;                ///< ITA space-group number
    int num_ops;               ///< Number of symmetry operations
    std::vector<SymOp> ops;    ///< Symmetry operators
};

/**
 * @brief Thread-safe database of crystallographic space groups.
 *
 * SpaceGroupDB supports exactly-once loading from a data file and
 * provides lock-free, read-only access after initialization.
 *
 * The database is immutable once loaded.
 */
class SpaceGroupDB {
public:
    //SpaceGroupDB() = default;
    SpaceGroupDB() = delete;
    explicit SpaceGroupDB(const std::string& filename);


    /** Access space group by index. */
    const SpaceGroup& at(size_t idx) const;

    /** Retrieve space groups by ITA number. */
    const std::vector<const SpaceGroup*>& by_number(int number) const;

    /** Number of loaded space groups. */
    size_t size() const;

    /* destructor */

private:
    /** Load database from file (thread-safe, exactly-once). */
    void load(const std::string& filename);

    void load_impl(const std::string& filename);

    std::once_flag load_flag_;
    bool loaded_ = false;

    std::vector<SpaceGroup> groups_;
    std::unordered_map<int, std::vector<const SpaceGroup*>> number_index_;
};

/**
 * @brief Access the global SpaceGroup database.
 *
 * Lazily loads spacegroups.dat on first use.
 * Thread-safe initialization (C++11 static initialization).
 *
 * @return const reference to global SpaceGroupDB
 */
const SpaceGroupDB& spacegroup_db();

} // namespace pwdft

#endif


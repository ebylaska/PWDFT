
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "SpaceGroupDB.hpp"
#include "util_paths.hpp"

namespace pwdft {


/*******************************************
 *                                         *
 *           spacegroup_db                 *
 *                                         *
 *******************************************/
 /**
 * @brief Global accessor for the crystallographic space-group database.
 *
 * Returns a reference to a lazily initialized, immutable `SpaceGroupDB`
 * containing all crystallographic space groups defined in
 * `spacegroups.dat`.
 *
 * The database is initialized on first use via a function-local static
 * object. This guarantees:
 *
 *  - Thread-safe initialization (per C++11 static initialization rules)
 *  - Exactly-once loading of the space-group database
 *  - No static initialization order dependencies across translation units
 *
 * The location of `spacegroups.dat` is resolved at runtime using
 * `resolve_spacegroups_dat()`, allowing seamless operation in both
 * installed and developer (source-tree) configurations.
 *
 * This accessor should be used wherever symmetry information is required
 * (e.g., lattice construction, symmetry reduction, k-point generation),
 * rather than instantiating `SpaceGroupDB` directly.
 *
 * @return Const reference to the global `SpaceGroupDB` instance.
 */
const SpaceGroupDB& spacegroup_db()
{
   static SpaceGroupDB db(resolve_spacegroups_dat());
   return db;
}



/*******************************************
 *                                         *
 *       SpaceGroupDB::SpaceGroupDB        *
 *                                         *
 *******************************************/
/**
 * @brief Construct and immediately load a space-group database.
 *
 * @param filename Path to the spacegroups.dat file.
 *
 * Equivalent to default construction followed by load(filename).
 * Throws on file I/O failure or parse error.
 */
SpaceGroupDB::SpaceGroupDB(const std::string& filename)
{
    load(filename);
}


/*******************************************
 *                                         *
 *           SpaceGroupDB::load            *
 *                                         *
 *******************************************/
/**
 * @brief Load space-group definitions from a data file.
 *
 * Reads crystallographic space-group data from the specified file
 * and populates the internal database. The file is expected to contain
 * space-group entries in the standard NWPW format:
 *
 *   name  ITA_number  num_ops
 *   (3×4 symmetry matrices for each operation)
 *
 * This function may only be called once per object instance.
 *
 * @param filename Path to the spacegroups.dat file.
 *
 * @throws std::logic_error if called more than once.
 * @throws std::runtime_error if the file cannot be opened.
 */
void SpaceGroupDB::load(const std::string& filename)
{
    std::call_once(load_flag_, &SpaceGroupDB::load_impl, this, filename);
}

void SpaceGroupDB::load_impl(const std::string& filename)
{
    std::ifstream in(filename);
    if (!in)
        throw std::runtime_error("Failed to open spacegroups.dat");

    groups_.clear();
    number_index_.clear();

    while (true) {
        SpaceGroup sg;

        if (!(in >> sg.name))
            break; // EOF

        in >> sg.number;
        in >> sg.num_ops;

        sg.ops.resize(sg.num_ops);

        for (int op = 0; op < sg.num_ops; ++op) {
            for (int row = 0; row < 3; ++row) {
                in >> sg.ops[op].R[row][0]
                   >> sg.ops[op].R[row][1]
                   >> sg.ops[op].R[row][2]
                   >> sg.ops[op].t[row];
            }
        }

        groups_.push_back(std::move(sg));
    }

    for (const auto& sg : groups_) {
        number_index_[sg.number].push_back(&sg);
    }

    loaded_ = true;
}



/*******************************************
 *                                         *
 *           SpaceGroupDB::at              *
 *                                         *
 *******************************************/
/**
 * @brief Access a space group by database index.
 *
 * Provides read-only access to the space group stored at the given
 * index in the internal database.
 *
 * @param idx Zero-based index into the space-group list.
 *
 * @return Const reference to the requested SpaceGroup.
 *
 * @throws std::logic_error if the database has not been loaded.
 * @throws std::out_of_range if idx is invalid.
 */
const SpaceGroup& SpaceGroupDB::at(size_t idx) const
{
    if (!loaded_)
        throw std::logic_error("SpaceGroupDB used before load()");
    return groups_.at(idx);
}


/*******************************************
 *                                         *
 *         SpaceGroupDB::by_number         *
 *                                         *
 *******************************************/
/**
 * @brief Retrieve space groups by ITA number.
 *
 * Returns all space-group entries matching the given International
 * Tables for Crystallography (ITA) number. Multiple entries may exist
 * due to different settings or origin choices.
 *
 * @param number ITA space-group number (1–230).
 *
 * @return Vector of pointers to matching SpaceGroup objects.
 *         Returns an empty vector if no matches are found.
 *
 * @note Returned pointers remain valid for the lifetime of the database.
 */
const std::vector<const SpaceGroup*>& SpaceGroupDB::by_number(int number) const
{
    static const std::vector<const SpaceGroup*> empty;
    auto it = number_index_.find(number);
    return (it != number_index_.end()) ? it->second : empty;
}


/*******************************************
 *                                         *
 *           SpaceGroupDB::size            *
 *                                         *
 *******************************************/
/**
 * @brief Return the number of space groups in the database.
 *
 * @return Total number of loaded space-group entries.
 */
size_t SpaceGroupDB::size() const
{
    return groups_.size();
}



} // namespace pwdft


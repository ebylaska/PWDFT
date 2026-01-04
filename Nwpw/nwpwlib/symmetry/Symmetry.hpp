#ifndef _SYMMETRY_HPP_
#define _SYMMETRY_HPP_

#pragma once
#include <string>
#include <vector>
#include <memory>

#include "json.hpp"
using json = nlohmann::json;
#include "SpaceGroupDB.hpp"   // SymOp, SpaceGroup
#include "PointGroupGenerators.hpp"

namespace pwdft {

//class Lattice;   // forward declaration

/**
 * @brief Crystallographic or molecular symmetry container.
 *
 * Symmetry encapsulates either:
 *  - a crystallographic space group (from SpaceGroupDB), or
 *  - a molecular point group (generated on demand).
 *
 * The class is immutable after construction and may be safely
 * shared across lattice, ion, and Brillouin-zone objects.
 */
class Symmetry {
public:
    /** Construct trivial symmetry (C1 / P1). */
    Symmetry();

    /** Construct from space-group name (e.g. "P4bm"). */
    explicit Symmetry(const std::string& spacegroup_name);

    /** Construct from ITA space-group number. */
    explicit Symmetry(int ita_number);

    /** Construct from point-group symbol (e.g. "D6h", "Ih"). */
    static Symmetry from_point_group(const std::string& pg_symbol);

    /** Construct from json. */
    static Symmetry from_json(const nlohmann::json& j);


    /** Number of symmetry operators. */
    size_t order() const;

    /** Access symmetry operators. */
    const std::vector<SymOp>& operators() const;

    /** Return symmetry name (HM or Schoenflies). */
    const std::string& name() const;

    /** True if this is a space group. */
    bool is_space_group() const;

    /** True if this is a point group. */
    bool is_point_group() const;

    /** Number of centering lattice points (space groups only). */
    int num_centering() const;

    /** Apply conventional → primitive transformation. */
    //void convert_to_primitive(Lattice& lattice);

    bool is_trivial() const { return ops_.size() == 1; }


private:
    enum class Type {
        Trivial,
        PointGroup,
        SpaceGroup
    };

    Type type_ = Type::Trivial;

    std::string name_;
    int ita_number_ = 0;

    int num_centering_ = 1;

    std::vector<SymOp> ops_;
};


} // namespace pwdft

#endif


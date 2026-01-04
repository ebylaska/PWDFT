#ifndef _POINTGROUPGENERATORS_HPP_
#define _POINTGROUPGENERATORS_HPP_


#pragma once
#include <vector>
#include <string>
#include "SpaceGroupDB.hpp"  // for SymOp

namespace pwdft {

/**
 * @class PointGroupGenerators
 * @brief Programmatic generators for three-dimensional crystallographic
 *        and molecular point groups.
 *
 * The @c PointGroupGenerators class provides a complete, table-free
 * implementation of all standard three-dimensional point groups using
 * analytic symmetry generators rather than lookup tables.
 *
 * Supported point-group families include:
 *   - Cyclic groups: @f$ C_n @f$, @f$ C_{nv} @f$, @f$ C_{nh} @f$
 *   - Dihedral groups: @f$ D_n @f$, @f$ D_{nh} @f$, @f$ D_{nd} @f$
 *   - Improper rotation groups: @f$ S_{2n} @f$
 *   - Trivial and inversion groups: @f$ C_1 @f$, @f$ C_s @f$, @f$ C_i @f$
 *   - Polyhedral groups: @f$ T @f$, @f$ T_d @f$, @f$ T_h @f$,
 *     @f$ O @f$, @f$ O_h @f$, @f$ I @f$, @f$ I_h @f$
 *
 * All symmetry operations are returned as @c SymOp objects representing
 * affine transformations @f$ (\mathbf{R} | \mathbf{t}) @f$, with
 * zero translation vectors appropriate for point-group symmetry acting
 * about the origin.
 *
 * The implementation is based on a minimal set of primitive symmetry
 * constructors:
 *   - Identity
 *   - Proper rotations about arbitrary axes
 *   - Mirror reflections through planes
 *   - Inversion through the origin
 *
 * Higher-order point groups are constructed through systematic
 * group-theoretic composition of these primitives, ensuring:
 *   - Mathematical correctness
 *   - Group closure
 *   - Absence of redundant operators
 *   - Exact operator counts for each group
 *
 * This class is intentionally stateless and provides only static
 * member functions. It is designed to serve as the point-group
 * symmetry backend for higher-level objects such as @c Symmetry,
 * @c Lattice, @c BrillouinZone, and k-point reduction logic.
 *
 * @note
 *     Unlike space groups, which are loaded from tabulated data,
 *     point groups are generated analytically to avoid duplication,
 *     reduce data maintenance, and enable precise control over
 *     symmetry construction and validation.
 *
 * @see SymOp
 * @see SpaceGroupDB
 */

class PointGroupGenerators {
public:
    static std::vector<SymOp> generate(const std::string& symbol);

private:
    static std::vector<SymOp> Cn(int n);
    static std::vector<SymOp> Cnv(int n);
    static std::vector<SymOp> Cnh(int n);

    static std::vector<SymOp> Dn(int n);
    static std::vector<SymOp> Dnh(int n);
    static std::vector<SymOp> Dnd(int n);

    static std::vector<SymOp> S2n(int n);

    static std::vector<SymOp> Cs();
    static std::vector<SymOp> Ci();

    static std::vector<SymOp> T();
    static std::vector<SymOp> Td();
    static std::vector<SymOp> Th();

    static std::vector<SymOp> O();
    static std::vector<SymOp> Oh();

    static std::vector<SymOp> I();
    static std::vector<SymOp> Ih();

    static SymOp rotation_about_axis(const double axis[3], double angle_rad);

    static SymOp identity();
    static SymOp inversion();
    static SymOp mirror(const double normal[3]);
};

}

#endif

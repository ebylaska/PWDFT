/* PointGroupGenerators.cpp  -
*/


#include <fstream>
#include <sstream>
#include <stdexcept>
#include <array>
#include <cmath>


#include "PointGroupGenerators.hpp"

namespace pwdft {


/*******************************************
 *                                         *
 *      PointGroupGenerators::generate     *
 *                                         *
 *******************************************/
 /**
 * @brief Generate symmetry operators for a finite 3D point group.
 *
 * Constructs the complete set of symmetry operations corresponding to a
 * three-dimensional point group specified by its Schoenflies symbol
 * (e.g., C2v, D3h, Oh, Ih).
 *
 * This routine serves as a single dispatch point for point-group symmetry
 * generation.  The Schoenflies symbol is parsed and mapped to an analytic
 * generator function.  No lookup tables are used; all point groups are
 * constructed algorithmically from fundamental symmetry elements such as
 * rotations, mirrors, and inversion.
 *
 * Supported point groups include:
 *   - Trivial groups:        C1, Cs, Ci
 *   - Cyclic groups:         Cn, Cnv, Cnh
 *   - Dihedral groups:       Dn, Dnh, Dnd
 *   - Improper rotations:    S2n
 *   - Polyhedral groups:     T, Td, Th
 *                            O, Oh
 *                            I, Ih
 *
 * The returned symmetry operators are represented as @c SymOp objects
 * with zero translation vectors, appropriate for point-group symmetry
 * acting about the origin.
 *
 * This function is intended to be called by higher-level symmetry
 * management code (e.g., the @c Symmetry class) and cleanly separates
 * symbolic symmetry selection from group-theoretic construction.
 *
 * @param symbol
 *     Schoenflies symbol identifying the desired point group.
 *
 * @return
 *     A vector containing the complete set of symmetry operators for the
 *     specified point group.
 *
 * @throws std::runtime_error
 *     If the symbol is malformed, invalid, or corresponds to an
 *     unsupported point group.
 */
std::vector<SymOp> PointGroupGenerators::generate(const std::string& symbol)
{
    // Trivial groups
    if (symbol == "C1")
        return { identity() };

    if (symbol == "Cs")
        return Cs();

    if (symbol == "Ci")
        return Ci();

    // Icosahedral / cubic / tetrahedral
    if (symbol == "T")
        return T();
    if (symbol == "Td")
        return Td();
    if (symbol == "Th")
        return Th();

    if (symbol == "O")
        return O();
    if (symbol == "Oh")
        return Oh();

    if (symbol == "I")
        return I();
    if (symbol == "Ih")
        return Ih();

    // Axial groups: parse prefix + integer
    auto parse_n = [&](size_t pos) -> int {
        int n = 0;
        for (; pos < symbol.size(); ++pos) {
            if (!std::isdigit(symbol[pos])) break;
            n = 10*n + (symbol[pos] - '0');
        }
        return n;
    };

    // Cn, Cnv, Cnh
    if (symbol[0] == 'C') {
        int n = parse_n(1);
        if (n <= 0)
            throw std::runtime_error("Invalid Cn group: " + symbol);

        if (symbol.size() == 1 + std::to_string(n).size())
            return Cn(n);

        if (symbol.back() == 'v')
            return Cnv(n);

        if (symbol.back() == 'h')
            return Cnh(n);
    }

    // Dn, Dnh, Dnd
    if (symbol[0] == 'D') {
        int n = parse_n(1);
        if (n <= 0)
            throw std::runtime_error("Invalid Dn group: " + symbol);

        if (symbol.size() == 1 + std::to_string(n).size())
            return Dn(n);

        if (symbol.substr(symbol.size()-2) == "nh")
            return Dnh(n);

        if (symbol.substr(symbol.size()-2) == "nd")
            return Dnd(n);
    }

    // Improper rotation groups S2n
    if (symbol[0] == 'S') {
        int n = parse_n(1);
        if (n <= 0 || n % 2 != 0)
            throw std::runtime_error("Invalid S2n group: " + symbol);

        return S2n(n/2);
    }

    throw std::runtime_error("Unsupported point group: " + symbol);
}


/*******************************************
 *                                         *
 *      PointGroupGenerators:Cn            *
 *                                         *
 *******************************************/
/**
 * @brief Generate the cyclic point group Cn.
 *
 * Constructs the cyclic point group @f$ C_n @f$, consisting of all proper
 * rotations about a single principal axis by integer multiples of
 * @f$ 2\pi/n @f$.
 *
 * The symmetry operations are generated analytically as rotations about
 * the Cartesian z-axis, which is taken as the principal symmetry axis.
 * All operations are proper rotations (determinant +1) and have zero
 * translation vectors, appropriate for point-group symmetry.
 *
 * The resulting group contains exactly @f$ n @f$ symmetry operations:
 * the identity and @f$ n-1 @f$ rotations.
 *
 * @param n
 *     Order of the cyclic group. Must satisfy @f$ n \ge 1 @f$.
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ C_n @f$.
 *
 * @throws std::invalid_argument
 *     If @p n is less than 1.
 */
std::vector<SymOp> PointGroupGenerators::Cn(int n)
{
    if (n < 1)
        throw std::invalid_argument("Cn: n must be >= 1");

    std::vector<SymOp> ops;
    ops.reserve(n);

    // principal rotation axis (z)
    const double axis[3] = {0.0, 0.0, 1.0};

    for (int k = 0; k < n; ++k) {
        double angle = 2.0 * M_PI * k / n;
        SymOp op = rotation_about_axis(axis, angle);

        // ensure zero translation (point group)
        op.t[0] = 0.0;
        op.t[1] = 0.0;
        op.t[2] = 0.0;

        ops.push_back(op);
    }

    return ops;
}

/*******************************************
 *                                         *
 *      PointGroupGenerators::Cnv          *
 *                                         *
 *******************************************/
/**
 * @brief Generate the cyclic point group Cnv.
 *
 * Constructs the point group @f$ C_{nv} @f$, consisting of the cyclic
 * rotation group @f$ C_n @f$ augmented by @f$ n @f$ vertical mirror planes.
 *
 * The construction proceeds in two steps:
 *   1. Generate all proper rotations of @f$ C_n @f$ about the principal
 *      symmetry axis (taken to be the Cartesian z-axis).
 *   2. Add @f$ n @f$ vertical mirror planes containing the z-axis, with
 *      mirror normals lying in the xy-plane and evenly spaced by
 *      @f$ \pi/n @f$.
 *
 * All symmetry operations are represented as @c SymOp objects with zero
 * translation vectors, appropriate for point-group symmetry acting about
 * the origin.
 *
 * The resulting group contains exactly @f$ 2n @f$ symmetry operations.
 *
 * @param n
 *     Order of the principal rotation axis. Must satisfy @f$ n \ge 1 @f$.
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ C_{nv} @f$.
 *
 * @throws std::invalid_argument
 *     If @p n is less than 1.
 */
std::vector<SymOp> PointGroupGenerators::Cnv(int n)
{
    if (n < 1)
        throw std::invalid_argument("Cnv: n must be >= 1");

    std::vector<SymOp> ops;

    // 1) Start with Cn rotations
    auto cn = Cn(n);
    ops.insert(ops.end(), cn.begin(), cn.end());

    // 2) Add n vertical mirror planes
    //    Planes contain z-axis, normals lie in xy-plane
    for (int k = 0; k < n; ++k) {
        double angle = M_PI * k / n;

        double normal[3] = {
            std::cos(angle),
            std::sin(angle),
            0.0
        };

        ops.push_back(mirror(normal));
    }

    return ops;
}


/*******************************************
 *                                         *
 *      PointGroupGenerators::Cnh          *
 *                                         *
 *******************************************/
/**
 * @brief Generate the cyclic point group Cnh.
 *
 * Constructs the point group @f$ C_{nh} @f$, consisting of the cyclic
 * rotation group @f$ C_n @f$ augmented by a horizontal mirror plane
 * perpendicular to the principal rotation axis.
 *
 * The construction proceeds as follows:
 *   1. Generate all proper rotations of @f$ C_n @f$ about the principal
 *      symmetry axis (taken to be the Cartesian z-axis).
 *   2. Construct the horizontal mirror plane @f$ \sigma_h @f$ coincident
 *      with the xy-plane.
 *   3. Form the improper symmetry operations by combining the horizontal
 *      mirror with each rotation: @f$ \sigma_h C_n^k @f$.
 *
 * All symmetry operations are represented as @c SymOp objects with zero
 * translation vectors, appropriate for point-group symmetry acting about
 * the origin.
 *
 * The resulting group contains exactly @f$ 2n @f$ symmetry operations.
 *
 * @param n
 *     Order of the principal rotation axis. Must satisfy @f$ n \ge 1 @f$.
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ C_{nh} @f$.
 *
 * @throws std::invalid_argument
 *     If @p n is less than 1.
 */
std::vector<SymOp> PointGroupGenerators::Cnh(int n)
{
    if (n < 1)
        throw std::invalid_argument("Cnh: n must be >= 1");

    std::vector<SymOp> ops;

    // 1) Start with Cn rotations
    auto cn = Cn(n);
    ops.insert(ops.end(), cn.begin(), cn.end());

    // 2) Horizontal mirror plane (xy-plane)
    //    normal = z
    //SymOp sigma_h = mirror({0.0, 0.0, 1.0});
    const double z_axis[3] = {0.0, 0.0, 1.0};
    SymOp sigma_h = mirror(z_axis);

    // 3) Combine σ_h with all rotations: σ_h * Cn^k
    for (const auto& r : cn) {
        SymOp op{};

        // R = σ_h * R_k
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                op.R[i][j] = 0.0;
                for (int k = 0; k < 3; ++k)
                    op.R[i][j] += sigma_h.R[i][k] * r.R[k][j];
            }

        // point group → zero translation
        op.t[0] = op.t[1] = op.t[2] = 0.0;

        ops.push_back(op);
    }

    return ops;
}



/*******************************************
 *                                         *
 *      PointGroupGenerators::Dn           *
 *                                         *
 *******************************************/
/**
 * @brief Generate the dihedral point group Dn.
 *
 * Constructs the dihedral point group @f$ D_n @f$, consisting of the
 * cyclic rotation group @f$ C_n @f$ together with @f$ n @f$ twofold
 * rotations about axes perpendicular to the principal rotation axis.
 *
 * The construction proceeds as follows:
 *   1. Generate all proper rotations of @f$ C_n @f$ about the principal
 *      symmetry axis (taken to be the Cartesian z-axis).
 *   2. Add @f$ n @f$ additional twofold (C2) rotations about axes lying
 *      in the xy-plane and passing through the origin.  These axes are
 *      evenly spaced by angles of @f$ \pi/n @f$.
 *
 * All symmetry operations are proper rotations (determinant +1) and are
 * represented as @c SymOp objects with zero translation vectors,
 * appropriate for point-group symmetry acting about the origin.
 *
 * The resulting group contains exactly @f$ 2n @f$ symmetry operations.
 *
 * @param n
 *     Order of the principal rotation axis. Must satisfy @f$ n \ge 2 @f$.
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ D_n @f$.
 *
 * @throws std::invalid_argument
 *     If @p n is less than 2.
 */
std::vector<SymOp> PointGroupGenerators::Dn(int n)
{
    if (n < 2)
        throw std::invalid_argument("Dn: n must be >= 2");

    std::vector<SymOp> ops;

    // 1) Add Cn rotations about z
    auto cn = Cn(n);
    ops.insert(ops.end(), cn.begin(), cn.end());

    // 2) Add n perpendicular C2 rotations
    //    Axes lie in xy-plane, angles k*pi/n
    for (int k = 0; k < n; ++k) {
        double angle = M_PI * k / n;

        double axis[3] = {
            std::cos(angle),
            std::sin(angle),
            0.0
        };

        // C2 = rotation by pi about axis
        SymOp c2 = rotation_about_axis(axis, M_PI);

        ops.push_back(c2);
    }

    return ops;
}


/*******************************************
 *                                         *
 *      PointGroupGenerators::Dnh          *
 *                                         *
 *******************************************/
/**
 * @brief Generate the dihedral point group Dnh.
 *
 * Constructs the dihedral point group @f$ D_{nh} @f$, consisting of the
 * proper dihedral rotations of @f$ D_n @f$ augmented by a horizontal
 * mirror plane perpendicular to the principal rotation axis.
 *
 * The construction proceeds as follows:
 *   1. Generate all proper rotations of the dihedral group @f$ D_n @f$,
 *      including the principal @f$ C_n @f$ rotations about the z-axis
 *      and the @f$ n @f$ perpendicular twofold rotations in the xy-plane.
 *   2. Construct the horizontal mirror plane @f$ \sigma_h @f$ coincident
 *      with the xy-plane.
 *   3. Form the improper symmetry operations by combining the horizontal
 *      mirror with each element of @f$ D_n @f$: @f$ \sigma_h D_n @f$.
 *
 * All symmetry operations are represented as @c SymOp objects with zero
 * translation vectors, appropriate for point-group symmetry acting about
 * the origin.
 *
 * The resulting group contains exactly @f$ 4n @f$ symmetry operations.
 *
 * @param n
 *     Order of the principal rotation axis. Must satisfy @f$ n \ge 2 @f$.
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ D_{nh} @f$.
 *
 * @throws std::invalid_argument
 *     If @p n is less than 2.
 */
std::vector<SymOp> PointGroupGenerators::Dnh(int n)
{
    if (n < 2)
        throw std::invalid_argument("Dnh: n must be >= 2");

    std::vector<SymOp> ops;

    // 1) Start with Dn
    auto dn = Dn(n);
    ops.insert(ops.end(), dn.begin(), dn.end());

    // 2) Horizontal mirror plane (xy-plane)
    //SymOp sigma_h = mirror({0.0, 0.0, 1.0});
    const double z_axis[3] = {0.0, 0.0, 1.0};
    SymOp sigma_h = mirror(z_axis);

    // 3) Form the coset σ_h * Dn
    for (const auto& r : dn) {
        SymOp op{};

        // R = σ_h * R_dn
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                op.R[i][j] = 0.0;
                for (int k = 0; k < 3; ++k)
                    op.R[i][j] += sigma_h.R[i][k] * r.R[k][j];
            }

        // point group → zero translation
        op.t[0] = 0.0;
        op.t[1] = 0.0;
        op.t[2] = 0.0;

        ops.push_back(op);
    }

    return ops;
}


/*******************************************
 *                                         *
 *      PointGroupGenerators::Dnd          *
 *                                         *
 *******************************************/
/**
 * @brief Generate the dihedral point group Dnd.
 *
 * Constructs the dihedral point group @f$ D_{nd} @f$, consisting of the
 * proper dihedral rotations of @f$ D_n @f$ augmented by @f$ n @f$ diagonal
 * mirror planes.
 *
 * The construction proceeds as follows:
 *   1. Generate all proper rotations of the dihedral group @f$ D_n @f$,
 *      including the principal @f$ C_n @f$ rotations about the z-axis
 *      and the @f$ n @f$ perpendicular twofold rotations in the xy-plane.
 *   2. Construct @f$ n @f$ diagonal mirror planes @f$ \sigma_d @f$ that
 *      contain the principal z-axis and bisect the angles between
 *      adjacent perpendicular C2 axes.
 *   3. Form the improper symmetry operations by combining each diagonal
 *      mirror plane with all elements of @f$ D_n @f$: @f$ \sigma_d D_n @f$.
 *
 * All symmetry operations are represented as @c SymOp objects with zero
 * translation vectors, appropriate for point-group symmetry acting about
 * the origin.
 *
 * The resulting group contains exactly @f$ 4n @f$ symmetry operations.
 *
 * @param n
 *     Order of the principal rotation axis. Must satisfy @f$ n \ge 2 @f$.
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ D_{nd} @f$.
 *
 * @throws std::invalid_argument
 *     If @p n is less than 2.
 */
std::vector<SymOp> PointGroupGenerators::Dnd(int n)
{
    if (n < 2)
        throw std::invalid_argument("Dnd: n must be >= 2");

    std::vector<SymOp> ops;

    // 1) Start with Dn
    auto dn = Dn(n);
    ops.insert(ops.end(), dn.begin(), dn.end());

    // 2) Diagonal mirror planes σ_d
    //    Planes contain z-axis, bisect angles between C2 axes
    for (int k = 0; k < n; ++k) {
        double angle = M_PI / (2.0 * n) + M_PI * k / n;

        double normal[3] = {
            std::cos(angle),
            std::sin(angle),
            0.0
        };

        SymOp sigma_d = mirror(normal);

        // 3) Form coset σ_d * Dn
        for (const auto& r : dn) {
            SymOp op{};

            // R = σ_d * R_dn
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) {
                    op.R[i][j] = 0.0;
                    for (int m = 0; m < 3; ++m)
                        op.R[i][j] += sigma_d.R[i][m] * r.R[m][j];
                }

            // point group → zero translation
            op.t[0] = 0.0;
            op.t[1] = 0.0;
            op.t[2] = 0.0;

            ops.push_back(op);
        }
    }

    return ops;
}


/*******************************************
 *                                         *
 *      PointGroupGenerators::S2n          *
 *                                         *
 *******************************************/
/**
 * @brief Generate the improper rotation point group S2n.
 *
 * Constructs the improper rotation point group @f$ S_{2n} @f$, consisting
 * of a cyclic subgroup of proper rotations augmented by improper
 * rotations formed from a rotation followed by a horizontal reflection.
 *
 * The construction proceeds as follows:
 *   1. Generate the proper cyclic rotations @f$ C_n @f$ about the
 *      principal symmetry axis (taken to be the Cartesian z-axis).
 *   2. Construct the improper rotation generator
 *      @f$ S_{2n} = \sigma_h C_{2n} @f$, where @f$ C_{2n} @f$ is a rotation
 *      by @f$ \pi/n @f$ about the z-axis and @f$ \sigma_h @f$ is reflection
 *      through the horizontal (xy) plane.
 *   3. Form the remaining improper elements by combining the generator
 *      with each proper rotation: @f$ S_{2n} C_n^k @f$.
 *
 * All symmetry operations are represented as @c SymOp objects with zero
 * translation vectors, appropriate for point-group symmetry acting about
 * the origin.
 *
 * The resulting group contains exactly @f$ 2n @f$ symmetry operations:
 * @f$ n @f$ proper rotations and @f$ n @f$ improper rotations.
 *
 * @param n
 *     Order parameter defining the group @f$ S_{2n} @f$. Must satisfy
 *     @f$ n \ge 1 @f$.
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ S_{2n} @f$.
 *
 * @throws std::invalid_argument
 *     If @p n is less than 1.
 */
std::vector<SymOp> PointGroupGenerators::S2n(int n)
{
    if (n < 1)
        throw std::invalid_argument("S2n: n must be >= 1");

    std::vector<SymOp> ops;

    // 1) Proper rotations: Cn
    auto cn = Cn(n);
    ops.insert(ops.end(), cn.begin(), cn.end());

    // 2) Generator: S2n = σ_h * C_(2n)
    //    C_(2n) is rotation by pi/n about z
    //SymOp c2n = rotation_about_axis({0.0, 0.0, 1.0}, M_PI / n);
    //SymOp sigma_h = mirror({0.0, 0.0, 1.0});
    const double z_axis[3] = {0.0, 0.0, 1.0};
    SymOp c2n = rotation_about_axis(z_axis, M_PI / n);
    SymOp sigma_h = mirror(z_axis);

    SymOp s2n{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            s2n.R[i][j] = 0.0;
            for (int k = 0; k < 3; ++k)
                s2n.R[i][j] += sigma_h.R[i][k] * c2n.R[k][j];
        }
    s2n.t[0] = s2n.t[1] = s2n.t[2] = 0.0;

    // 3) Improper elements: S2n * Cn^k
    for (const auto& r : cn) {
        SymOp op{};

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) {
                op.R[i][j] = 0.0;
                for (int k = 0; k < 3; ++k)
                    op.R[i][j] += s2n.R[i][k] * r.R[k][j];
            }

        op.t[0] = op.t[1] = op.t[2] = 0.0;
        ops.push_back(op);
    }

    return ops;
}


/*******************************************
 *                                         *
 *      PointGroupGenerators::Cs           *
 *                                         *
 *******************************************/
/**
 * @brief Generate the mirror point group Cs.
 *
 * Constructs the point group @f$ C_s @f$, consisting of the identity
 * operation and a single mirror reflection.
 *
 * By convention, the mirror plane is taken to be the horizontal
 * (xy) plane, with its normal aligned along the Cartesian z-axis.
 * This choice fixes a consistent orientation for downstream symmetry
 * usage while remaining fully general for point-group symmetry.
 *
 * All symmetry operations are represented as @c SymOp objects with zero
 * translation vectors, appropriate for point-group symmetry acting about
 * the origin.
 *
 * The resulting group contains exactly two symmetry operations:
 * the identity and one mirror reflection.
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ C_s @f$.
 */
std::vector<SymOp> PointGroupGenerators::Cs()
{
    std::vector<SymOp> ops;
    ops.reserve(2);

    // Identity
    ops.push_back(identity());

    // Single mirror plane (by convention: xy-plane, normal = z)
    //ops.push_back(mirror({0.0, 0.0, 1.0}));
    const double z_axis[3] = {0.0, 0.0, 1.0};
    ops.push_back(mirror(z_axis));

    return ops;
}


/*******************************************
 *                                         *
 *      PointGroupGenerators::Ci           *
 *                                         *
 *******************************************/
/**
 * @brief Generate the inversion point group Ci.
 *
 * Constructs the point group @f$ C_i @f$, consisting of the identity
 * operation and inversion through the origin.
 *
 * The inversion operation maps a point @f$ \mathbf{r} @f$ to
 * @f$ -\mathbf{r} @f$.  No mirror planes or rotational symmetries are
 * present beyond the identity.
 *
 * All symmetry operations are represented as @c SymOp objects with zero
 * translation vectors, appropriate for point-group symmetry acting about
 * the origin.
 *
 * The resulting group contains exactly two symmetry operations:
 * the identity and inversion.
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ C_i @f$.
 */
std::vector<SymOp> PointGroupGenerators::Ci()
{
    std::vector<SymOp> ops;
    ops.reserve(2);

    // Identity
    ops.push_back(identity());

    // Inversion through the origin
    ops.push_back(inversion());

    return ops;
}



/*******************************************
 *                                         *
 *      PointGroupGenerators::T            *
 *                                         *
 *******************************************/
/**
 * @brief Generate the proper tetrahedral point group T.
 *
 * Constructs the proper tetrahedral point group @f$ T @f$, consisting
 * exclusively of proper rotational symmetries of a regular tetrahedron.
 * No mirror planes or inversion operations are included.
 *
 * The construction uses the standard Cartesian embedding of the
 * tetrahedral group:
 *   - The identity operation.
 *   - Three twofold (C2) rotations about the Cartesian x, y, and z axes.
 *   - Eight threefold (C3) rotations about the four body diagonals of a
 *     cube, corresponding to the axes passing through opposite vertices
 *     of a regular tetrahedron.
 *
 * All symmetry operations are proper rotations (determinant +1) and are
 * represented as @c SymOp objects with zero translation vectors,
 * appropriate for point-group symmetry acting about the origin.
 *
 * The resulting group contains exactly 12 symmetry operations.
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ T @f$.
 */
std::vector<SymOp> PointGroupGenerators::T()
{
    std::vector<SymOp> ops;
    ops.reserve(12);

    // Identity
    ops.push_back(identity());

    // C2 rotations about x, y, z
    {
        double x_axis[3] = {1.0, 0.0, 0.0};
        double y_axis[3] = {0.0, 1.0, 0.0};
        double z_axis[3] = {0.0, 0.0, 1.0};

        ops.push_back(rotation_about_axis(x_axis, M_PI));
        ops.push_back(rotation_about_axis(y_axis, M_PI));
        ops.push_back(rotation_about_axis(z_axis, M_PI));
    }

    // C3 rotations about body diagonals
    const double inv = 1.0 / std::sqrt(3.0);
    const double axes[4][3] = {
        {  inv,  inv,  inv},
        { -inv, -inv,  inv},
        { -inv,  inv, -inv},
        {  inv, -inv, -inv}
    };

    for (const auto& a : axes) {
        ops.push_back(rotation_about_axis(a,  2.0 * M_PI / 3.0));
        ops.push_back(rotation_about_axis(a, -2.0 * M_PI / 3.0));
    }

    return ops;
}


/*******************************************
 *                                         *
 *      PointGroupGenerators::Td           *
 *                                         *
 *******************************************/
/**
 * @brief Generate the tetrahedral point group Td.
 *
 * Constructs the tetrahedral point group @f$ T_d @f$, consisting of the
 * proper rotational symmetries of a regular tetrahedron augmented by
 * diagonal mirror planes.
 *
 * The construction proceeds as follows:
 *   1. Generate all proper rotational symmetry operations of the
 *      tetrahedral group @f$ T @f$.
 *   2. Add six diagonal mirror planes @f$ \sigma_d @f$ corresponding to
 *      planes that bisect the edges of the tetrahedron.  These mirror
 *      planes have normals lying along the face diagonals of a cube.
 *
 * All symmetry operations are represented as @c SymOp objects with zero
 * translation vectors, appropriate for point-group symmetry acting about
 * the origin.
 *
 * The resulting group contains exactly 24 symmetry operations:
 *   - 12 proper rotations
 *   - 12 improper operations (mirror reflections)
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ T_d @f$.
 */
std::vector<SymOp> PointGroupGenerators::Td()
{
    std::vector<SymOp> ops;

    // Start with all proper tetrahedral rotations
    std::vector<SymOp> T_ops = T();
    ops.insert(ops.end(), T_ops.begin(), T_ops.end());

    // Add 6 diagonal mirror planes (σ_d)
    const double inv = 1.0 / std::sqrt(2.0);
    const double normals[6][3] = {
        {  inv,  inv,  0.0},
        {  inv, -inv,  0.0},
        {  inv,  0.0,  inv},
        {  inv,  0.0, -inv},
        {  0.0,  inv,  inv},
        {  0.0,  inv, -inv}
    };

    for (const auto& n : normals) {
        ops.push_back(mirror(n));
    }

    return ops;
}


/*******************************************
 *                                         *
 *      PointGroupGenerators::Th           *
 *                                         *
 *******************************************/
/**
 * @brief Generate the tetrahedral point group Th.
 *
 * Constructs the tetrahedral point group @f$ T_h @f$, consisting of the
 * proper rotational symmetries of a regular tetrahedron augmented by
 * inversion symmetry.
 *
 * The construction proceeds as follows:
 *   1. Generate all proper rotational symmetry operations of the
 *      tetrahedral group @f$ T @f$.
 *   2. Form the improper symmetry operations by composing each proper
 *      rotation with inversion: @f$ iT @f$.
 *
 * No mirror planes are present in this group; all improper symmetry
 * operations arise solely from inversion combined with proper rotations.
 *
 * All symmetry operations are represented as @c SymOp objects with zero
 * translation vectors, appropriate for point-group symmetry acting about
 * the origin.
 *
 * The resulting group contains exactly 24 symmetry operations:
 *   - 12 proper rotations
 *   - 12 improper operations (inversion × rotation)
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ T_h @f$.
 */
std::vector<SymOp> PointGroupGenerators::Th()
{
    std::vector<SymOp> ops;

    // Start with proper tetrahedral rotations
    std::vector<SymOp> T_ops = T();
    ops.insert(ops.end(), T_ops.begin(), T_ops.end());

    // Add inversion * each rotation
    SymOp inv = inversion();

    for (const auto& op : T_ops) {
        SymOp composed{};

        // R = (-I) * R
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                composed.R[i][j] = -op.R[i][j];

        // t = 0 (point group)
        for (int i = 0; i < 3; ++i)
            composed.t[i] = 0.0;

        ops.push_back(composed);
    }

    return ops;
}


/*******************************************
 *                                         *
 *      PointGroupGenerators::O            *
 *                                         *
 *******************************************/
/**
 * @brief Generate the proper octahedral point group O.
 *
 * Constructs the proper octahedral point group @f$ O @f$, consisting
 * exclusively of proper rotational symmetries of a cube or regular
 * octahedron.  No mirror planes or inversion symmetry are included.
 *
 * The construction uses a standard Cartesian embedding and includes:
 *   - The identity operation.
 *   - Three fourfold (C4) rotation axes coincident with the Cartesian
 *     x, y, and z axes, including rotations by ±90° and 180°.
 *   - Four threefold (C3) rotation axes along the body diagonals of
 *     the cube, corresponding to rotations by ±120°.
 *   - Six twofold (C2) rotation axes passing through the centers of
 *     opposite edges of the cube.
 *
 * All symmetry operations are proper rotations (determinant +1) and are
 * represented as @c SymOp objects with zero translation vectors,
 * appropriate for point-group symmetry acting about the origin.
 *
 * The resulting group contains exactly 24 symmetry operations.
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ O @f$.
 */
std::vector<SymOp> PointGroupGenerators::O()
{
    std::vector<SymOp> ops;
    ops.reserve(24);

    // Identity
    ops.push_back(identity());

    // Cartesian axes
    const double x_axis[3] = {1.0, 0.0, 0.0};
    const double y_axis[3] = {0.0, 1.0, 0.0};
    const double z_axis[3] = {0.0, 0.0, 1.0};

    // C4 rotations about x, y, z
    const double angles[3] = { M_PI/2.0, M_PI, -M_PI/2.0 };
    for (double a : angles) {
        ops.push_back(rotation_about_axis(x_axis, a));
        ops.push_back(rotation_about_axis(y_axis, a));
        ops.push_back(rotation_about_axis(z_axis, a));
    }

    // Body-diagonal C3 rotations
    const double inv3 = 1.0 / std::sqrt(3.0);
    const double diagonals[4][3] = {
        {  inv3,  inv3,  inv3},
        { -inv3, -inv3,  inv3},
        { -inv3,  inv3, -inv3},
        {  inv3, -inv3, -inv3}
    };

    for (const auto& a : diagonals) {
        ops.push_back(rotation_about_axis(a,  2.0 * M_PI / 3.0));
        ops.push_back(rotation_about_axis(a, -2.0 * M_PI / 3.0));
    }

    // C2 rotations about edge-center axes
    const double inv2 = 1.0 / std::sqrt(2.0);
    const double edge_axes[6][3] = {
        { inv2,  inv2,  0.0},
        { inv2, -inv2,  0.0},
        { inv2,  0.0,  inv2},
        { inv2,  0.0, -inv2},
        { 0.0,  inv2,  inv2},
        { 0.0,  inv2, -inv2}
    };

    for (const auto& a : edge_axes) {
        ops.push_back(rotation_about_axis(a, M_PI));
    }

    return ops;
}


/*******************************************
 *                                         *
 *      PointGroupGenerators::Oh           *
 *                                         *
 *******************************************/
/**
 * @brief Generate the full octahedral point group Oh.
 *
 * Constructs the full octahedral point group @f$ O_h @f$, consisting of
 * the proper rotational symmetries of the octahedral group @f$ O @f$
 * augmented by inversion symmetry.
 *
 * The construction proceeds as follows:
 *   1. Generate all proper rotational symmetry operations of the
 *      octahedral group @f$ O @f$.
 *   2. Form the improper symmetry operations by composing each proper
 *      rotation with inversion: @f$ iO @f$.
 *
 * Mirror planes and improper rotations (e.g., S4) arise naturally from
 * the combination of inversion with proper rotations and are not
 * constructed explicitly.
 *
 * All symmetry operations are represented as @c SymOp objects with zero
 * translation vectors, appropriate for point-group symmetry acting about
 * the origin.
 *
 * The resulting group contains exactly 48 symmetry operations:
 *   - 24 proper rotations
 *   - 24 improper operations (inversion × rotation)
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ O_h @f$.
 */
std::vector<SymOp> PointGroupGenerators::Oh()
{
    std::vector<SymOp> ops;

    // Proper octahedral rotations
    std::vector<SymOp> O_ops = O();
    ops.insert(ops.end(), O_ops.begin(), O_ops.end());

    // Inversion * each rotation
    for (const auto& op : O_ops) {
        SymOp composed{};

        // R = (-I) * R
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                composed.R[i][j] = -op.R[i][j];

        // Point group: no translation
        for (int i = 0; i < 3; ++i)
            composed.t[i] = 0.0;

        ops.push_back(composed);
    }

    return ops;
}


/*******************************************
 *                                         *
 *      PointGroupGenerators::I            *
 *                                         *
 *******************************************/
/**
 * @brief Generate the proper icosahedral point group I.
 *
 * Constructs the proper icosahedral point group @f$ I @f$, consisting
 * exclusively of proper rotational symmetries of a regular icosahedron
 * (or equivalently, a regular dodecahedron).  No mirror planes or
 * inversion symmetry are included.
 *
 * The construction uses a standard Cartesian embedding of the
 * icosahedral group and includes:
 *   - The identity operation.
 *   - Fifteen twofold (C2) rotations about axes passing through the
 *     midpoints of opposite edges.
 *   - Twenty threefold (C3) rotations about axes passing through the
 *     centers of opposite triangular faces.
 *   - Twenty-four fivefold (C5) rotations about axes passing through
 *     opposite vertices, expressed using the golden ratio
 *     @f$ \varphi = (1+\sqrt{5})/2 @f$.
 *
 * All symmetry operations are proper rotations (determinant +1) and are
 * represented as @c SymOp objects with zero translation vectors,
 * appropriate for point-group symmetry acting about the origin.
 *
 * The resulting group contains exactly 60 symmetry operations.
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ I @f$.
 */
std::vector<SymOp> PointGroupGenerators::I()
{
    std::vector<SymOp> ops;
    ops.reserve(60);

    const double phi = 0.5 * (1.0 + std::sqrt(5.0));
    const double pi  = M_PI;

    // Identity
    ops.push_back(identity());

    /* -------------------------
     *  2-fold rotations (C2)
     * ------------------------- */
    {
        const double axes[][3] = {
            { 1, 0, 0 }, { -1, 0, 0 },
            { 0, 1, 0 }, { 0,-1, 0 },
            { 0, 0, 1 }, { 0, 0,-1 }
        };

        for (const auto& a : axes)
            ops.push_back(rotation_about_axis(a, pi));
    }

    /* -------------------------
     *  3-fold rotations (C3)
     * ------------------------- */
    {
        const double axes[][3] = {
            {  1,  1,  1 }, { -1,  1,  1 },
            {  1, -1,  1 }, {  1,  1, -1 },
            { -1, -1,  1 }, { -1,  1, -1 },
            {  1, -1, -1 }, { -1, -1, -1 }
        };

        for (const auto& a : axes) {
            ops.push_back(rotation_about_axis(a,  2.0*pi/3.0));
            ops.push_back(rotation_about_axis(a, -2.0*pi/3.0));
        }
    }

    /* -------------------------
     *  5-fold rotations (C5)
     * ------------------------- */
    {
        const double axes[][3] = {
            { 0,  1,  phi }, { 0, -1,  phi },
            { 0,  1, -phi }, { 0, -1, -phi },

            { 1,  phi, 0 }, { -1,  phi, 0 },
            { 1, -phi, 0 }, { -1, -phi, 0 },

            { phi, 0,  1 }, { -phi, 0,  1 },
            { phi, 0, -1 }, { -phi, 0, -1 }
        };

        for (const auto& a : axes) {
            ops.push_back(rotation_about_axis(a,  2.0*pi/5.0));
            ops.push_back(rotation_about_axis(a,  4.0*pi/5.0));
            ops.push_back(rotation_about_axis(a, -2.0*pi/5.0));
            ops.push_back(rotation_about_axis(a, -4.0*pi/5.0));
        }
    }

    return ops;
}


/*******************************************
 *                                         *
 *      PointGroupGenerators::Ih           *
 *                                         *
 *******************************************/
/**
 * @brief Generate the full icosahedral point group Ih.
 *
 * Constructs the full icosahedral point group @f$ I_h @f$, which includes
 * both proper and improper symmetry operations of a regular icosahedron
 * (or dodecahedron).  This group is obtained by augmenting the proper
 * icosahedral rotation group @f$ I @f$ with inversion symmetry.
 *
 * The construction proceeds as follows:
 *   - Generate all 60 proper rotational symmetries of the icosahedron
 *     via @c PointGroupGenerators::I().
 *   - For each proper rotation @f$ R @f$, generate an improper operation
 *     by composing inversion with the rotation, yielding @f$ -R @f$.
 *
 * The resulting group consists of:
 *   - 60 proper rotations (determinant +1)
 *   - 60 improper operations formed by inversion × rotation
 *     (determinant −1)
 *
 * All symmetry operations are represented as @c SymOp objects with
 * zero translation vectors, appropriate for point-group symmetry
 * acting about the origin.
 *
 * The total number of symmetry operations is exactly 120.
 *
 * @return
 *     A vector of @c SymOp objects representing the complete set of
 *     symmetry operations for the point group @f$ I_h @f$.
 */
std::vector<SymOp> PointGroupGenerators::Ih()
{
    std::vector<SymOp> ops;

    // Proper icosahedral rotations (60)
    std::vector<SymOp> I_ops = I();
    ops.insert(ops.end(), I_ops.begin(), I_ops.end());

    // Inversion * each rotation → improper elements
    for (const auto& op : I_ops) {
        SymOp composed{};

        // R = (-I) * R
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                composed.R[i][j] = -op.R[i][j];

        // Point group → zero translation
        for (int i = 0; i < 3; ++i)
            composed.t[i] = 0.0;

        ops.push_back(composed);
    }

    return ops;  // total = 120
}



/*******************************************
 *                                         *
 * PointGroupGenerators:rotation_about_axis*
 *                                         *
 *******************************************/
/**
 * @brief Construct a proper rotation about an arbitrary axis.
 *
 * Generates a symmetry operation corresponding to a right-handed
 * rotation by a specified angle about an arbitrary axis passing
 * through the origin.
 *
 * The rotation matrix is constructed using Rodrigues' rotation
 * formula.  The input axis is normalized internally; the direction
 * of rotation follows the right-hand rule.
 *
 * This routine produces a proper rotation with determinant +1
 * and zero translation, appropriate for point-group symmetry
 * operations.
 *
 * @param axis
 *     A 3-element array specifying the rotation axis direction.
 *     The axis need not be normalized but must be nonzero.
 *
 * @param angle_rad
 *     Rotation angle in radians. Positive angles correspond to
 *     right-handed rotations about the given axis.
 *
 * @return
 *     A @c SymOp representing the rotation:
 *       - @c R is a 3×3 orthogonal rotation matrix
 *       - @c t is the zero vector
 *
 * @throws std::invalid_argument
 *     If the input axis has zero length.
 *
 * @note
 *     This routine is the fundamental building block for generating
 *     all proper rotational elements in axial, polyhedral, and
 *     icosahedral point groups.
 */
SymOp PointGroupGenerators::rotation_about_axis(const double axis[3], double angle_rad)
{
    SymOp op{};

    // normalize axis
    double x = axis[0];
    double y = axis[1];
    double z = axis[2];

    double norm = std::sqrt(x*x + y*y + z*z);
    if (norm <= 0.0)
        throw std::invalid_argument("rotation_about_axis: zero-length axis");

    x /= norm;
    y /= norm;
    z /= norm;

    double c = std::cos(angle_rad);
    double s = std::sin(angle_rad);
    double t = 1.0 - c;

    // Rodrigues' rotation formula
    op.R[0][0] = t*x*x + c;
    op.R[0][1] = t*x*y - s*z;
    op.R[0][2] = t*x*z + s*y;

    op.R[1][0] = t*x*y + s*z;
    op.R[1][1] = t*y*y + c;
    op.R[1][2] = t*y*z - s*x;

    op.R[2][0] = t*x*z - s*y;
    op.R[2][1] = t*y*z + s*x;
    op.R[2][2] = t*z*z + c;

    // point group → no translation
    op.t[0] = 0.0;
    op.t[1] = 0.0;
    op.t[2] = 0.0;

    return op;
}

/*******************************************
 *                                         *
 *    PointGroupGenerators::identity       *
 *                                         *
 *******************************************/
/**
 * @brief Construct the identity symmetry operation.
 *
 * Generates the identity element of a point group, corresponding to
 * no rotation and no translation.  The rotation matrix is the 3×3
 * identity matrix and the translation vector is zero.
 *
 * This operation acts as the neutral element under composition of
 * symmetry operations and is included explicitly in all generated
 * point groups.
 *
 * @return
 *     A @c SymOp representing the identity operation:
 *       - @c R = I (3×3 identity matrix)
 *       - @c t = (0, 0, 0)
 */
SymOp PointGroupGenerators::identity()
{
    SymOp op{};

    // R = I
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j)
            op.R[i][j] = (i == j ? 1.0 : 0.0);
        op.t[i] = 0.0;
    }

    return op;
}



/*******************************************
 *                                         *
 *    PointGroupGenerators:inversion       *
 *                                         *
 *******************************************/
/**
 * @brief Construct the inversion symmetry operation.
 *
 * Generates the inversion operation through the origin, corresponding
 * to the mapping @f$ \mathbf{r} \rightarrow -\mathbf{r} @f$.
 *
 * The rotation matrix is the negative identity matrix (−I), yielding
 * a symmetry operation with determinant −1.  The translation vector
 * is zero, appropriate for point-group symmetry acting about the origin.
 *
 * Inversion serves as the generator of centrosymmetric point groups
 * such as @f$ C_i @f$, @f$ T_h @f$, @f$ O_h @f$, and @f$ I_h @f$ when
 * combined with proper rotational symmetries.
 *
 * @return
 *     A @c SymOp representing inversion:
 *       - @c R = −I (3×3 negative identity matrix)
 *       - @c t = (0, 0, 0)
 */
SymOp PointGroupGenerators::inversion()
{
    SymOp op{};

    // R = -I
    op.R[0][0] = -1.0;  op.R[0][1] =  0.0;  op.R[0][2] =  0.0;
    op.R[1][0] =  0.0;  op.R[1][1] = -1.0;  op.R[1][2] =  0.0;
    op.R[2][0] =  0.0;  op.R[2][1] =  0.0;  op.R[2][2] = -1.0;

    // point-group symmetry → zero translation
    op.t[0] = 0.0;
    op.t[1] = 0.0;
    op.t[2] = 0.0;

    return op;
}

/*******************************************
 *                                         *
 *    PointGroupGenerators:mirror          *
 *                                         *
 *******************************************/
/**
 * @brief Construct a mirror (reflection) symmetry operation.
 *
 * Generates a reflection through a plane passing through the origin.
 * The mirror plane is specified by its outward normal vector.
 *
 * The reflection matrix is constructed using the formula
 * @f[
 *   R = I - 2 \, \mathbf{n}\mathbf{n}^T
 * @f]
 * where @f$ \mathbf{n} @f$ is the unit normal vector to the mirror plane.
 *
 * This operation is an improper symmetry with determinant −1 and
 * zero translation, appropriate for point-group symmetry acting
 * about the origin.
 *
 * @param normal
 *     A 3-element array specifying the normal vector of the mirror
 *     plane. The vector need not be normalized but must be nonzero.
 *
 * @return
 *     A @c SymOp representing the mirror operation:
 *       - @c R is a 3×3 reflection matrix
 *       - @c t is the zero vector
 *
 * @throws std::invalid_argument
 *     If the input normal vector has zero length.
 *
 * @note
 *     This routine is used to generate mirror planes in point groups
 *     such as @f$ C_s @f$, @f$ C_{nv} @f$, @f$ D_{nh} @f$, @f$ T_d @f$,
 *     and related improper symmetry groups.
 */
SymOp PointGroupGenerators::mirror(const double normal[3])
{
   SymOp op{};

   // normalize the normal vector
   double nx = normal[0];
   double ny = normal[1];
   double nz = normal[2];

   double norm = std::sqrt(nx*nx + ny*ny + nz*nz);
   if (norm <= 0.0)
       throw std::invalid_argument("mirror: zero-length normal vector");

   nx /= norm;
   ny /= norm;
   nz /= norm;

   // R = I - 2 n n^T
   op.R[0][0] = 1.0 - 2.0 * nx * nx;
   op.R[0][1] =      - 2.0 * nx * ny;
   op.R[0][2] =      - 2.0 * nx * nz;

   op.R[1][0] =      - 2.0 * ny * nx;
   op.R[1][1] = 1.0 - 2.0 * ny * ny;
   op.R[1][2] =      - 2.0 * ny * nz;

   op.R[2][0] =      - 2.0 * nz * nx;
   op.R[2][1] =      - 2.0 * nz * ny;
   op.R[2][2] = 1.0 - 2.0 * nz * nz;

   // point-group symmetry → zero translation
   op.t[0] = 0.0;
   op.t[1] = 0.0;
   op.t[2] = 0.0;

   return op;
}


}





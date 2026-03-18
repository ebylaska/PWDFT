#pragma once
#include <string>
#include <vector>
#include <stdexcept>

namespace pwdft {

/**
 * @brief Minimal character-table support for point-group irrep decomposition.
 *
 * Stores conjugacy-class labels/sizes and irreducible-representation (irrep)
 * characters for finite point groups (molecular symmetry).
 *
 * @note This is intended for point groups (finite). Space-group irreps and
 *       k-dependent little-group representations are out of scope.
 */
class PointGroupCharacterTable {
public:
    struct Irrep {
        std::string name;           // e.g. "A1"
        int dim = 1;                // irrep dimension
        std::vector<double> chi;    // characters per conjugacy class (same order as classes())
    };

    const std::string& group() const { return group_; }
    int order() const { return order_; } // |G|
    const std::vector<std::string>& classes() const { return class_names_; }
    const std::vector<int>& class_sizes() const { return class_sizes_; }
    const std::vector<Irrep>& irreps() const { return irreps_; }

    int irrep_index(const std::string& name) const;

    /** Factory for supported groups by Schoenflies symbol (e.g. "C2v", "Td"). */
    static PointGroupCharacterTable from_symbol(const std::string& symbol);

    // Internal builder so implementation can populate private fields.
    static PointGroupCharacterTable build(std::string group,
                                         int order,
                                         std::vector<std::string> class_names,
                                         std::vector<int> class_sizes,
                                         std::vector<Irrep> irreps);
private:
    std::string group_;
    int order_ = 0;

    std::vector<std::string> class_names_;
    std::vector<int> class_sizes_;
    std::vector<Irrep> irreps_;
};

} // namespace pwdft

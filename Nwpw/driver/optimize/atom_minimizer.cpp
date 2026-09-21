// atom_minimizer.cpp
//
// Atom-coordinate optimizer for PWDFT.
//
// This file currently implements only the permutation computation: for
// each symmetry operation and each atom, which atom does the op map it
// to. Gradient projection and position symmetrization build on top of
// this, and both need it precomputed once per optimization run.

#include "atom_minimizer.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>

namespace pwdft {

AtomPermutation compute_atom_permutation(
    const std::vector<FracCoord>& coords_frac,
    const std::vector<FracSymOp>& ops,
    double tol)
{
    const int n_atoms = static_cast<int>(coords_frac.size());
    const int n_ops   = static_cast<int>(ops.size());

    if (n_atoms == 0)
        throw std::runtime_error(
            "compute_atom_permutation: empty coordinate list");
    if (n_ops == 0)
        throw std::runtime_error(
            "compute_atom_permutation: empty symmetry op list");

    // Wrap reference positions once; every comparison uses the wrapped form.
    std::vector<FracCoord> f_ref(n_atoms);
    for (int i = 0; i < n_atoms; ++i)
        f_ref[i] = wrap_frac(coords_frac[i]);

    const double tol_sq = tol * tol;

    AtomPermutation result;
    result.perm.assign(n_ops, std::vector<int>(n_atoms, -1));

    for (int g = 0; g < n_ops; ++g)
    {
        const FracSymOp& op = ops[g];

        for (int i = 0; i < n_atoms; ++i)
        {
            const FracCoord f_image = wrap_frac(op.apply(f_ref[i]));

            // Find the atom closest to f_image under the minimum-image metric.
            int    best_j   = -1;
            double best_d2  = std::numeric_limits<double>::max();

            for (int j = 0; j < n_atoms; ++j)
            {
                const double d2 = frac_distance_sq(f_image, f_ref[j]);
                if (d2 < best_d2) {
                    best_d2 = d2;
                    best_j  = j;
                }
            }

            if (best_d2 > tol_sq)
            {
                std::ostringstream oss;
                oss << "compute_atom_permutation: op " << g
                    << " maps atom " << i
                    << " to no atom within tolerance " << tol
                    << " (nearest distance = " << std::sqrt(best_d2) << ")";
                throw std::runtime_error(oss.str());
            }

            result.perm[g][i] = best_j;
        }
    }

    return result;
}

} // namespace pwdft

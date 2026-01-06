#include <string>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <array>
#include <iomanip>
#include <cmath>

#include "json.hpp"
using json = nlohmann::json;

#include "Symmetry.hpp"

#include "apply_kpoint_symmetry_pruning.hpp"


namespace pwdft {


/**************************************************
 *                                                *
 *               wrap_min_image                   *
 *                                                *
 **************************************************/
/**
 * @brief Apply minimum-image wrapping for a periodic coordinate.
 *
 * Maps a real-valued coordinate difference onto the interval
 * (-0.5, 0.5], corresponding to the minimum-image convention
 * on a unit-period torus.
 *
 * This operation is used for distance calculations in fractional
 * coordinates, ensuring that periodic images differing by an
 * integer lattice translation are treated as equivalent.
 *
 * In the context of Brillouin-zone and k-point symmetry pruning,
 * this guarantees that distances are measured along the shortest
 * path on the reciprocal-space torus.
 *
 * @param[in] x
 *   Coordinate difference (typically a[i] - b[i]) in fractional units.
 *
 * @return
 *   Wrapped coordinate difference in the range (-0.5, 0.5].
 *
 * @note
 *   This implementation uses std::round(), which yields symmetric
 *   wrapping around zero and is robust for symmetry-equivalence tests
 *   with small numerical tolerances.
 */
static inline double wrap_min_image(double x)
{
    // map to (-0.5, 0.5]
    return x - std::round(x);
}


/**************************************************
 *                                                *
 *                  dist_torus                    *
 *                                                *
 **************************************************/
/**
 * @brief Squared distance between two points on a 3D periodic torus.
 *
 * Computes the squared Euclidean distance between two points defined in
 * fractional coordinates, accounting for periodic boundary conditions
 * by applying the minimum-image convention in each dimension.
 *
 * This metric is appropriate for comparing k-points in the Brillouin zone,
 * where points differing by a reciprocal lattice vector are physically
 * equivalent. Distances are evaluated on the unit torus [0,1)^3.
 *
 * The function returns the squared distance to avoid unnecessary square
 * roots during tight symmetry-equivalence tests.
 *
 * @param[in] a
 *   First point in fractional coordinates.
 *
 * @param[in] b
 *   Second point in fractional coordinates.
 *
 * @return
 *   Squared minimum-image distance between @p a and @p b on the 3D torus.
 *
 * @note
 *   This routine assumes coordinates have already been wrapped into a
 *   consistent fractional convention (e.g. Monkhorst–Pack mapping).
 */
static inline double dist2_torus(const std::array<double,3>& a,
                                 const std::array<double,3>& b)
{
    double dx = wrap_min_image(a[0] - b[0]);
    double dy = wrap_min_image(a[1] - b[1]);
    double dz = wrap_min_image(a[2] - b[2]);
    return dx*dx + dy*dy + dz*dz;
}


/**************************************************
 *                                                *
 *      monkhorst_pack_pointgroup_prune           *
 *                                                *
 **************************************************/
/**
 * @brief Prune Monkhorst–Pack k-points using point-group and time-reversal symmetry.
 *
 * This routine reduces a Monkhorst–Pack k-point set by identifying symmetry-
 * equivalent points under the rotational part of the effective space-group
 * symmetry and time-reversal symmetry. K-points related by symmetry are merged
 * and their weights accumulated, producing an irreducible Brillouin zone.
 *
 * Only rotational components of the symmetry operations are used; fractional
 * translations are ignored, as they do not act in reciprocal space.
 *
 * K-points are compared using the Monkhorst–Pack convention, internally wrapped
 * into the interval [0,1) via a (-0.5,0.5) → [0,1) mapping before symmetry tests.
 * Equivalence is determined on the reciprocal-space torus with a tight numerical
 * tolerance.
 *
 * PWDFT uses strict group-theoretic symmetry pruning based on the supplied
 * symmetry operators. This may result in a larger irreducible k-point set than
 * legacy NWChem, but preserves full k-star structure required for modern
 * electronic-structure methods (band structure, Berry phases, magnetism,
 * Wannier constructions, etc.).
 *
 * @param[in,out] ks
 *   Monkhorst–Pack k-point list.
 *   Each entry is of the form {kx, ky, kz, weight}, where k-components are in
 *   fractional reciprocal coordinates. On exit, symmetry-equivalent points are
 *   merged and the list is compacted.
 *
 * @param[in] sym
 *   Effective symmetry object providing the rotational symmetry operators.
 *   If the symmetry is trivial (identity), no pruning is performed.
 *
 * @note
 *   Time-reversal symmetry is always applied explicitly by testing k ↔ −k
 *   equivalence.
 *
 * @warning
 *   This routine assumes a *full* Monkhorst–Pack grid as input. Applying it to
 *   already reduced or hand-selected k-point sets may lead to incorrect weights.
 */
static void monkhorst_pack_pointgroup_prune(std::vector<std::vector<double>>& ks, const pwdft::Symmetry& sym)
{
    constexpr double tol  = 1.0e-9;
    constexpr double tol2 = tol * tol;

    // Identity (P1 / C1): nothing to do
    if (sym.is_trivial())
        return;

    const size_t nk = ks.size();

    // --------------------------------------------------
    // Extract unique rotation matrices (ignore translations)
    // --------------------------------------------------
    std::vector<std::array<double,9>> rotations;

    auto same_rotation = [&](const std::array<double,9>& A,
                             const std::array<double,9>& B)
    {
        for (int i=0;i<9;i++)
            if (std::abs(A[i] - B[i]) > 1e-12)
                return false;
        return true;
    };

    for (const auto& op : sym.operators())
    {
        std::array<double,9> R;
        for (int i=0;i<3;i++)
            for (int j=0;j<3;j++)
                R[3*i + j] = op.R[i][j];

        bool found = false;
        for (const auto& Q : rotations)
            if (same_rotation(R,Q)) { found = true; break; }

        if (!found)
            rotations.push_back(R);
    }

    // --------------------------------------------------
    // Wrap all ks into [0,1) using Monkhorst–Pack convention
    // --------------------------------------------------

    std::vector<std::array<double,3>> k01(nk);
    for (size_t i=0;i<nk;i++)
    {
       double x = ks[i][0];
       double y = ks[i][1];
       double z = ks[i][2];
 
       // MP convention: stored in [-0.5,0.5)
       //x = x + 0.5; x -= std::floor(x);
       //y = y + 0.5; y -= std::floor(y);
       //z = z + 0.5; z -= std::floor(z);

       auto wrap_to_01_mp = [](double& x, double& y, double& z)
       {
          x = x + 0.5; x -= std::floor(x);
          y = y + 0.5; y -= std::floor(y);
          z = z + 0.5; z -= std::floor(z);
       };
       wrap_to_01_mp(x,y,z);
 
       k01[i] = {x,y,z};
    }

//    // --------------------------------------------------
//    // Wrap all ks into [0,1)
//    // --------------------------------------------------
//    std::vector<std::array<double,3>> k01(nk);
//    for (size_t i=0;i<nk;i++)
//    {
//        double x = ks[i][0];
//        double y = ks[i][1];
//        double z = ks[i][2];
//
//        x -= std::floor(x);
//        y -= std::floor(y);
//        z -= std::floor(z);
//
//        k01[i] = {x,y,z};
//    }

    // --------------------------------------------------
    // Symmetry pruning
    // --------------------------------------------------
    for (size_t i=0;i<nk;i++)
    {
        if (ks[i][3] <= 0.0) continue;

        const auto& ki = k01[i];

        for (size_t j=i+1;j<nk;j++)
        {
            if (ks[j][3] <= 0.0) continue;

            const auto& kj = k01[j];
            bool equivalent = false;

            for (const auto& R : rotations)
            {
                // R * k_j
                std::array<double,3> rk {
                    R[0]*kj[0] + R[1]*kj[1] + R[2]*kj[2],
                    R[3]*kj[0] + R[4]*kj[1] + R[5]*kj[2],
                    R[6]*kj[0] + R[7]*kj[1] + R[8]*kj[2]
                };

                // Wrap into [0,1)
                for (int a=0;a<3;a++)
                    rk[a] -= std::floor(rk[a]);

                // Rk ≈ k
                if (dist2_torus(rk, ki) < tol2)
                {
                    equivalent = true;
                }
                else
                {
                    // Rk ≈ -k  (time reversal)
                    std::array<double,3> mk {
                        -rk[0], -rk[1], -rk[2]
                    };
                    for (int a=0;a<3;a++)
                        mk[a] -= std::floor(mk[a]);

                    if (dist2_torus(mk, ki) < tol2)
                        equivalent = true;
                }

                if (equivalent)
                {
                    ks[i][3] += ks[j][3];
                    ks[j][3]  = 0.0;
                    break;
                }
            }

            if (equivalent)
                break;
        }
    }

    // --------------------------------------------------
    // Compact
    // --------------------------------------------------
    std::vector<std::vector<double>> reduced;
    reduced.reserve(ks.size());
    for (const auto& k : ks)
        if (k[3] > tol)
            reduced.push_back(k);

    ks.swap(reduced);
}





/**************************************************
 *                                                *
 *        apply_kpoint_symmetry_pruning           *
 *                                                *
 **************************************************/
/**
 * @brief Apply crystallographic symmetry pruning to Monkhorst–Pack k-points.
 *
 * This routine reduces a Monkhorst–Pack Brillouin-zone k-point set using
 * the effective symmetry of the system. When a space-group symmetry is
 * present, k-points related by symmetry operations are identified and
 * merged, with their weights accumulated, producing an irreducible
 * Brillouin zone consistent with the full space group.
 *
 * Unlike legacy NWChem k-point reduction, this function performs a
 * strictly group-theoretic pruning based on the supplied symmetry
 * operations and does not apply heuristic or energy-preserving
 * collapses beyond the symmetry action. As a result, the irreducible
 * k-point set may contain more points, but preserves the full k-star
 * structure required for band structure calculations, Berry-phase
 * methods, magnetism, and Wannier-based analyses.
 *
 * On exit, the RTDB is updated with:
 *   - pruned k-point list
 *   - number of original k-points
 *   - number of symmetry-reduced k-points
 *
 * @param rtdbstr JSON-encoded RTDB string containing the Brillouin-zone
 *                and effective symmetry information.
 *
 * @return Updated RTDB string with symmetry-pruned k-points.
 */
std::string apply_kpoint_symmetry_pruning(std::string rtdbstr)
{
   auto rtdb = json::parse(rtdbstr);

   if (!rtdb.contains("effective_symmetry"))
      return rtdbstr;

   auto& bz = rtdb["nwpw"]["brillouin_zone"];
   if (!bz.contains("kvectors"))
      return rtdbstr;

   std::vector<std::vector<double>> ks =
       bz["kvectors"].get<std::vector<std::vector<double>>>();

   const int nk_original = ks.size();
   const auto& es = rtdb["effective_symmetry"];
   if (es["type"] == "space_group")
   {
      Symmetry sym(es["name"].get<std::string>());
      monkhorst_pack_pointgroup_prune(ks, sym);
   }
   const int nk_effective = ks.size();

   bz["kvectors"] = ks;
   bz["num_kpoints_original"]   = nk_original;
   bz["num_kpoints_effective"]  = nk_effective;

   return rtdb.dump();
}


}


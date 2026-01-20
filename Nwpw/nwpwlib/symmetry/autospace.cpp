/* autospace.cpp - 
*/

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <array>
#include <cmath>
#include <algorithm>

#include "spglib.h"

#include "autospace.hpp"

namespace pwdft {

static void invert_3x3(const double A[9], double Ainv[9])
{
    const double det =
          A[0]*(A[4]*A[8] - A[5]*A[7])
        - A[1]*(A[3]*A[8] - A[5]*A[6])
        + A[2]*(A[3]*A[7] - A[4]*A[6]);

    if (std::abs(det) < 1e-14)
        throw std::runtime_error("invert_3x3: singular matrix");

    const double invdet = 1.0 / det;

    Ainv[0] =  (A[4]*A[8] - A[5]*A[7]) * invdet;
    Ainv[1] = -(A[1]*A[8] - A[2]*A[7]) * invdet;
    Ainv[2] =  (A[1]*A[5] - A[2]*A[4]) * invdet;

    Ainv[3] = -(A[3]*A[8] - A[5]*A[6]) * invdet;
    Ainv[4] =  (A[0]*A[8] - A[2]*A[6]) * invdet;
    Ainv[5] = -(A[0]*A[5] - A[2]*A[3]) * invdet;

    Ainv[6] =  (A[3]*A[7] - A[4]*A[6]) * invdet;
    Ainv[7] = -(A[0]*A[7] - A[1]*A[6]) * invdet;
    Ainv[8] =  (A[0]*A[4] - A[1]*A[3]) * invdet;
}

static inline double wrap01(double x)
{
    x -= std::floor(x);
    return (x < 0.0) ? x + 1.0 : x;
}


/*******************************************
 *                                         *
 *          detect_space_group             *
 *                                         *
 *******************************************/
AutoSpaceResult detect_space_group(const double unita[9],
                                   const std::vector<std::string>& symbols,
                                   const std::vector<double>& coords_xyz,
                                   double tolerance)
{
    AutoSpaceResult result;
    result.method = "spglib";
    result.tolerance_used = tolerance;

    const int natom = static_cast<int>(symbols.size());

    // -------------------------------
    // Sanity checks (you already had these right)
    // -------------------------------
    if (natom == 0)
    {
        result.message = "No atoms provided";
        return result;
    }

    if (coords_xyz.size() != 3 * natom)
    {
        result.message = "coords_xyz size does not match symbols";
        return result;
    }

    // -------------------------------
    // Build lattice for spglib (row-major)
    // -------------------------------
    double lattice[3][3];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            lattice[i][j] = unita[3*i + j];

    // -------------------------------
    // Build inverse lattice
    // -------------------------------
    double unita_inv[9];
    invert_3x3(unita, unita_inv);

    // -------------------------------
    // Convert symbols → integer types
    // -------------------------------
    std::vector<int> types(natom);
    std::vector<std::string> unique_symbols;

    for (int i = 0; i < natom; ++i)
    {
        auto it = std::find(unique_symbols.begin(),
                            unique_symbols.end(),
                            symbols[i]);
        if (it == unique_symbols.end())
        {
            unique_symbols.push_back(symbols[i]);
            types[i] = static_cast<int>(unique_symbols.size());
        }
        else
        {
            types[i] = static_cast<int>(it - unique_symbols.begin()) + 1;
        }
    }

    // -------------------------------
    // Cartesian → fractional positions
    // -------------------------------
    std::vector<std::array<double,3>> positions(natom);

    for (int i = 0; i < natom; ++i)
    {
        const double x = coords_xyz[3*i+0];
        const double y = coords_xyz[3*i+1];
        const double z = coords_xyz[3*i+2];

        double fx = unita_inv[0]*x + unita_inv[1]*y + unita_inv[2]*z;
        double fy = unita_inv[3]*x + unita_inv[4]*y + unita_inv[5]*z;
        double fz = unita_inv[6]*x + unita_inv[7]*y + unita_inv[8]*z;

        positions[i] = { wrap01(fx), wrap01(fy), wrap01(fz) };
    }

    // -------------------------------
    // Call spglib
    // -------------------------------
    SpglibDataset *dataset =
        spg_get_dataset(lattice,
                        reinterpret_cast<double (*)[3]>(positions.data()),
                        types.data(),
                        natom,
                        tolerance);

    if (!dataset)
    {
        result.message = "spglib failed to identify space group";
        return result;
    }

    // -------------------------------
    // Populate result
    // -------------------------------
    result.success = true;
    result.group_number = dataset->spacegroup_number;
    result.group_name   = dataset->international_symbol;
    //result.primitive_suggested = (dataset->primitive_lattice != nullptr);
    result.primitive_suggested = false;

    result.message = "OK";

    // Cleanup
    spg_free_dataset(dataset);

    return result;
}



} // namespace pwdft


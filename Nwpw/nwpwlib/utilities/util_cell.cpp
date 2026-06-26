
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>

/** @ingroup nwpw_utilities */
namespace pwdft {

/**
 * Helper function to calculate the determinant of a 3x3 matrix.
 * Equivalent to the 'deter3' external function in the Fortran code.
 */
double determinant3x3(const std::array<std::array<double, 3>, 3>& m) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

/**
 * Calculates the lattice_a matrix from lattice parameters.
 * 
 * @param lattice: Array of 6 doubles [a, b, c, alpha, beta, gamma]
 *                 where angles are in radians.
 * @param lattice_a: Output 3x3 matrix representing the lattice vectors.
 */
void util_cell_abc_abg_unita(const double *lattice, double *lattice_a)
{
    
    // Mapping inputs to descriptive names (0-based indexing)
    // cdist: [0]=a, [1]=b, [2]=c
    // cang:  [0]=alpha, [1]=beta, [2]=gamma
    std::array<double, 3> cdist = { lattice[0], lattice[1], lattice[2] };
    std::array<double, 3> cang  = { lattice[3], lattice[4], lattice[5] };

    std::array<std::array<double, 3>, 3> gmat = {{{0.0}}};

    // 1. Build the metric matrix (gmat)
    // Diagonal elements: g_ii = d_i^2
    for (int i = 0; i < 3; ++i) {
        gmat[i][i] = cdist[i] * cdist[i];
    }

    // Off-diagonal elements: g_ij = d_i * d_j * cos(angle)
    // We follow the Fortran loop logic which descends from gamma -> beta -> alpha
    int iang = 2; // Index for angles (0:alpha, 1:beta, 2:gamma)
    for (int i = 0; i < 3; ++i) {
        for (int j = i + 1; j < 3; ++j) {
            double val = cdist[i] * cdist[j] * std::cos(cang[iang]);
            gmat[i][j] = val;
            gmat[j][i] = val;
            iang--; 
        }
    }

    // 2. Get Volume: vol = sqrt(det(gmat))
    double vol = std::sqrt(std::abs(determinant3x3(gmat)));

    // 3. Generate lattice_a matrix
    double c1 = std::cos(cang[0]); // cos(alpha)
    double c2 = std::cos(cang[1]); // cos(beta)
    double c3 = std::cos(cang[2]); // cos(gamma)
    double s3 = std::sin(cang[2]); // sin(gamma)

    // Row 0
    lattice_a[0][0] = cdist[0] * s3;
    lattice_a[0][1] = 0.0;
    lattice_a[0][2] = (cdist[2] * (c2 - c1 * c3) / s3);

    // Row 1
    lattice_a[1][0] = cdist[0] * c3;
    lattice_a[1][1] = cdist[1];
    lattice_a[1][2] = cdist[2] * c1;

    // Row 2
    lattice_a[2][0] = 0.0;
    lattice_a[2][1] = 0.0;
    lattice_a[2][2] = (vol / (cdist[0] * cdist[1] * s3));
}


}

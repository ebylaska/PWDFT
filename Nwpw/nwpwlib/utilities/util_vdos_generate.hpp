#ifndef UTIL_VDOS_GENERATE_HPP
#define UTIL_VDOS_GENERATE_HPP

#include <algorithm>
#include <cmath>
#include <limits>

namespace pwdft {

/*
**************** Definitions of Cubes and Tetrahedrons *****************
*                                                                      *
*  Cube vertices are labeled:                                          *
*                                                                      *
*      (011)------------(111)                    (3)------------(7)    *
*        +                +                        +              +     *
*       /.               /|                       /.             /|     *
*      / .              / |                      / .            / |     *
*     /  .             /  |                     /  .           /  |     *
*    /   .            /   |                    /   .          /   |     *
* (001)------------(101) |      <====>      (1)------------(5)  |     *
*   |   .            |   |                    |   .          |   |     *
*   | (010)..........|.(110)                  | (2)..........|.(6)     *
*   |   .            |   /                    |   .          |   /     *
*   |  .             |  /                     |  .           |  /      *
*   | .              | /                      | .            | /       *
*   |.               |/                       |.             |/        *
*   +                +                        +              +         *
* (000)------------(100)                    (0)------------(4)         *
*                                                                      *
*  The cube is divided into 6 tetrahedra using the shortest body       *
*  diagonal of the reciprocal cell.                                    *
*                                                                      *
**************** Definitions of Cubes and Tetrahedrons *****************
*/

/******************************************
 *                                        *
 *          _util_vdos_generate_          *
 *                                        *
 ******************************************/
/*
 * Generate the vibrational density of states g(w) and integrated number
 * of states N(w) using the tetrahedron method on a uniform reciprocal grid.
 *
 * Inputs:
 *   idx,idy,idz - dimensions of the q-point grid
 *   eigs        - eigenvalues stored band-major, then k-point:
 *
 *                   eigs[ib * nk + ik]
 *
 *                 where:
 *                   nk = idx * idy * idz
 *                   ik = i0 + idx * (i1 + idy * i2)
 *                   ib = 0..neigs-1
 *
 *   neigs       - number of eigenvalues/bands per q-point
 *   npoints     - number of output frequency grid points
 *   emin,emax   - lower and upper bounds of the DOS grid
 *   l_unitg     - reciprocal lattice vectors (3x3, column-major)
 *
 * Output:
 *   efn         - flat array of size 3*npoints
 *
 *                   efn[3*k + 0] = e_k
 *                   efn[3*k + 1] = g(e_k)
 *                   efn[3*k + 2] = N(e_k)
 *
 * Notes:
 *   - Periodic wrap-around is applied on the grid.
 *   - The output frequencies are in the same units as eigs.
 */

inline double vib_Dstates_Tetra(double e, const double ee[4])
{
    double ds;
    double e1, e2, e4;
    double e21, e31, e41, e32, e42, e43;

    if ((ee[0] <= e) && (e < ee[1])) {
        e1  = e - ee[0];
        e21 = ee[1] - ee[0];
        e31 = ee[2] - ee[0];
        e41 = ee[3] - ee[0];
        ds  = 3.0 * e1 * e1 / (e21 * e31 * e41);
    } else if ((ee[1] <= e) && (e < ee[2])) {
        e2  = e - ee[1];
        e21 = ee[1] - ee[0];
        e31 = ee[2] - ee[0];
        e41 = ee[3] - ee[0];
        e32 = ee[2] - ee[1];
        e42 = ee[3] - ee[1];
        ds  = (3.0 * e21 + 6.0 * e2
             - 3.0 * (e31 + e42) * e2 * e2 / (e32 * e42))
             / (e31 * e41);
    } else if ((ee[2] <= e) && (e < ee[3])) {
        e4  = ee[3] - e;
        e41 = ee[3] - ee[0];
        e42 = ee[3] - ee[1];
        e43 = ee[3] - ee[2];
        ds  = 3.0 * e4 * e4 / (e41 * e42 * e43);
    } else {
        ds = 0.0;
    }

    return ds;
}

inline double vib_Nstates_Tetra(double e, const double ee[4])
{
    double ds;
    double e1, e2, e4;
    double e21, e31, e41, e32, e42, e43;

    if ((ee[0] <= e) && (e < ee[1])) {
        e1  = e - ee[0];
        e21 = ee[1] - ee[0];
        e31 = ee[2] - ee[0];
        e41 = ee[3] - ee[0];
        ds  = e1 * e1 * e1 / (e21 * e31 * e41);
    } else if ((ee[1] <= e) && (e < ee[2])) {
        e2  = e - ee[1];
        e21 = ee[1] - ee[0];
        e31 = ee[2] - ee[0];
        e41 = ee[3] - ee[0];
        e32 = ee[2] - ee[1];
        e42 = ee[3] - ee[1];
        ds  = (e21 * e21
             + 3.0 * e21 * e2
             + 3.0 * e2 * e2
             - (e31 + e42) * e2 * e2 * e2 / (e32 * e42))
             / (e31 * e41);
    } else if ((ee[2] <= e) && (e < ee[3])) {
        e4  = ee[3] - e;
        e41 = ee[3] - ee[0];
        e42 = ee[3] - ee[1];
        e43 = ee[3] - ee[2];
        ds  = 1.0 - e4 * e4 * e4 / (e41 * e42 * e43);
    } else if (e >= ee[3]) {
        ds = 1.0;
    } else {
        ds = 0.0;
    }

    return ds;
}

inline double vib_Dstates_Cube(double e,
                               const int itetra[4][6],
                               const double ecube[8])
{
    double ds = 0.0;

    for (int k = 0; k < 6; ++k) {
        double etetra[4] = {
            ecube[itetra[0][k]],
            ecube[itetra[1][k]],
            ecube[itetra[2][k]],
            ecube[itetra[3][k]]
        };

        std::sort(etetra, etetra + 4);
        ds += vib_Dstates_Tetra(e, etetra);
    }

    return ds;
}

inline double vib_Nstates_Cube(double e,
                               const int itetra[4][6],
                               const double ecube[8])
{
    double ds = 0.0;

    for (int k = 0; k < 6; ++k) {
        double etetra[4] = {
            ecube[itetra[0][k]],
            ecube[itetra[1][k]],
            ecube[itetra[2][k]],
            ecube[itetra[3][k]]
        };

        std::sort(etetra, etetra + 4);
        ds += vib_Nstates_Tetra(e, etetra);
    }

    return ds;
}

inline void util_vdos_generate(int idx, int idy, int idz,
                               const double* eigs,
                               int neigs,
                               int npoints,
                               double emin,
                               double emax,
                               const double* l_unitg,
                               double* efn)
{
    auto B = [l_unitg](int r, int c) -> double {
        return l_unitg[r + 3 * c];
    };

    const int nk = idx * idy * idz;

    auto IK = [=](int i, int j, int k) -> int {
        return i + idx * (j + idy * k);
    };

    auto EIG = [=](int i, int j, int k, int ib) -> double {
        return eigs[ib * nk + IK(i, j, k)];
    };

    auto EFN = [efn](int ip, int col) -> double& {
        return efn[3 * ip + col];
    };

    int dosgrid[3] = {idx, idy, idz};

    double unitg[3][3];
    double VG, VT;

    /* reciprocal cell volume */
    unitg[0][0] = B(1,1) * B(2,2) - B(2,1) * B(1,2);
    unitg[1][0] = B(2,1) * B(0,2) - B(0,1) * B(2,2);
    unitg[2][0] = B(0,1) * B(1,2) - B(1,1) * B(0,2);

    unitg[0][1] = B(1,2) * B(2,0) - B(2,2) * B(1,0);
    unitg[1][1] = B(2,2) * B(0,0) - B(0,2) * B(2,0);
    unitg[2][1] = B(0,2) * B(1,0) - B(1,2) * B(0,0);

    unitg[0][2] = B(1,0) * B(2,1) - B(2,0) * B(1,1);
    unitg[1][2] = B(2,0) * B(0,1) - B(0,0) * B(2,1);
    unitg[2][2] = B(0,0) * B(1,1) - B(1,0) * B(0,1);

    VG = B(0,0) * unitg[0][0]
       + B(1,0) * unitg[1][0]
       + B(2,0) * unitg[2][0];

    const int ncubes = dosgrid[0] * dosgrid[1] * dosgrid[2];
    const int ntetra = ncubes * 6;
    VT = VG / static_cast<double>(ntetra);

    /* find shortest diagonal */
    int k1_d[4]  = {0, 1, 0, 1};
    int k2_d[4]  = {0, 0, 1, 1};
    int k3_d[4]  = {0, 0, 0, 0};
    int k1_dd[4] = {1, 0, 1, 0};
    int k2_dd[4] = {1, 1, 0, 0};
    int k3_dd[4] = {1, 1, 1, 1};

    int id = 0;
    double rmax = std::numeric_limits<double>::max();

    for (int i = 0; i < 4; ++i) {
        const double kx  = k1_d[i]  * B(0,0) + k2_d[i]  * B(0,1) + k3_d[i]  * B(0,2);
        const double ky  = k1_d[i]  * B(1,0) + k2_d[i]  * B(1,1) + k3_d[i]  * B(1,2);
        const double kz  = k1_d[i]  * B(2,0) + k2_d[i]  * B(2,1) + k3_d[i]  * B(2,2);

        const double kkx = k1_dd[i] * B(0,0) + k2_dd[i] * B(0,1) + k3_dd[i] * B(0,2);
        const double kky = k1_dd[i] * B(1,0) + k2_dd[i] * B(1,1) + k3_dd[i] * B(1,2);
        const double kkz = k1_dd[i] * B(2,0) + k2_dd[i] * B(2,1) + k3_dd[i] * B(2,2);

        const double r = (kx - kkx) * (kx - kkx)
                       + (ky - kky) * (ky - kky)
                       + (kz - kkz) * (kz - kkz);

        if (r < rmax) {
            rmax = r;
            id = i;
        }
    }

    int itetra[4][6];

    if (id == 0) {
        int tmp[4][6] = {
            {0,0,0,0,0,0},
            {7,7,7,7,7,7},
            {1,1,2,2,4,4},
            {3,5,3,6,5,6}
        };
        std::copy(&tmp[0][0], &tmp[0][0] + 24, &itetra[0][0]);
    } else if (id == 1) {
        int tmp[4][6] = {
            {1,1,1,1,1,1},
            {6,6,6,6,6,6},
            {0,0,3,3,5,5},
            {2,4,2,7,4,7}
        };
        std::copy(&tmp[0][0], &tmp[0][0] + 24, &itetra[0][0]);
    } else if (id == 2) {
        int tmp[4][6] = {
            {2,2,2,2,2,2},
            {5,5,5,5,5,5},
            {3,3,0,0,6,6},
            {1,7,1,4,7,4}
        };
        std::copy(&tmp[0][0], &tmp[0][0] + 24, &itetra[0][0]);
    } else {
        int tmp[4][6] = {
            {3,3,3,3,3,3},
            {4,4,4,4,4,4},
            {2,2,1,1,7,7},
            {0,6,0,5,6,5}
        };
        std::copy(&tmp[0][0], &tmp[0][0] + 24, &itetra[0][0]);
    }

    const double de = (emax - emin) / static_cast<double>(npoints - 1);

    for (int ip = 0; ip < npoints; ++ip) {
        const double e = emin + ip * de;

        double f = 0.0;
        double g = 0.0;

        for (int kk = 0; kk < dosgrid[2]; ++kk) {
            for (int jj = 0; jj < dosgrid[1]; ++jj) {
                for (int ii = 0; ii < dosgrid[0]; ++ii) {
                    int ishft = ii + 1;
                    int jshft = jj + 1;
                    int kshft = kk + 1;

                    if (ishft >= dosgrid[0]) ishft = 0;
                    if (jshft >= dosgrid[1]) jshft = 0;
                    if (kshft >= dosgrid[2]) kshft = 0;

                    for (int ib = 0; ib < neigs; ++ib) {
                        double ecube[8];
                        ecube[0] = EIG(ii,    jj,    kk,    ib);
                        ecube[1] = EIG(ishft, jj,    kk,    ib);
                        ecube[2] = EIG(ii,    jshft, kk,    ib);
                        ecube[3] = EIG(ishft, jshft, kk,    ib);
                        ecube[4] = EIG(ii,    jj,    kshft, ib);
                        ecube[5] = EIG(ishft, jj,    kshft, ib);
                        ecube[6] = EIG(ii,    jshft, kshft, ib);
                        ecube[7] = EIG(ishft, jshft, kshft, ib);

                        f += vib_Dstates_Cube(e, itetra, ecube);
                        g += vib_Nstates_Cube(e, itetra, ecube);
                    }
                }
            }
        }

        f *= (VT / VG);
        g *= (VT / VG);

        EFN(ip, 0) = e;
        EFN(ip, 1) = f;
        EFN(ip, 2) = g;
    }
}

} // namespace pwdft

#endif

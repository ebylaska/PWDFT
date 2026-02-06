/* Psp1d_Hamaan.cpp -
   Author - Eric Bylaska
*/

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
//#include        "blas.h"

#include "Parallel.hpp"
#include "Psp1d_Hamann.hpp"
#include "util.hpp"

namespace pwdft {

#define FMT1 "%lf"
#define FMT2 " %lf %lf"
#define FMT10 "%10.3lf %10.3lf %10.3lf"

/*******************************************
 *                                         *
 *              util_matinvert             *
 *                                         *
 *******************************************/
/* Calculates matrix inverse based on Gauss-Jordan elimination
  method with partial pivoting.
*/
static void util_matinvert(int n, int nmax, double *a) {
  int irow;
  double big, tmp, pivinv;
  int *indx = new int[nmax];
  for (auto i = 0; i < n; ++i) {
    big = 0.0;
    for (auto j = i; j < n; ++j)
      if (std::fabs(a[j + i * nmax]) >= big) {
        big = fabs(a[j + j * nmax]);
        irow = j;
      }
    if (big <= 1.0e-9) {
      printf("Failed to invert matix\n");
      exit(99);
    }
    indx[i] = irow;

    if (irow != i) {
      for (auto j = 0; j < n; ++j) {
        tmp = a[irow + j * nmax];
        a[irow + j * nmax] = a[i + j * nmax];
        a[i + j * nmax] = tmp;
      }
    }

    pivinv = 1.0 / a[i + i * nmax];
    a[i + i * nmax] = 1.0;

    for (auto j = 0; j < n; ++j)
      a[i + j * nmax] *= pivinv;

    for (auto l = 0; l < n; ++l)
      if (l != i) {
        tmp = a[l + i * nmax];
        a[l + i * nmax] = 0.0;
        for (auto j = 0; j < n; ++j)
          a[l + j * nmax] -= a[i + j * nmax] * tmp;
      }
  }

  delete[] indx;
}

/*******************************************
 *                                         *
 *              util_simpson               *
 *                                         *
 *******************************************/
/**
 * @brief Composite Simpson integration on a uniform 1D grid.
 *
 * Computes the integral of y(x) sampled on an evenly spaced grid
 * using Simpson’s rule:
 *
 *   ∫ f(x) dx ≈ h/3 [ f0 + fN
 *                     + 4 Σ f_odd
 *                     + 2 Σ f_even ]
 *
 * Requirements:
 *   - n must be odd (even number of intervals)
 *   - grid spacing h must be uniform
 *
 * Notes:
 *   - No internal validation of n is performed.
 *   - Endpoint values y[0] and y[n-1] are included with weight 1.
 */
static double util_simpson(int n, double *y, double h) {
  double s = -y[0] - y[n - 1];
  for (auto i = 0; i < n; i += 2)
    s += 2.0 * y[i];
  for (auto i = 1; i < n; i += 2)
    s += 4.0 * y[i];

  return s * h / 3.0;
}


/*******************************************
 *                                         *
 *          convert_psp_type               *
 *                                         *
 *******************************************/
/**
 * @brief Decode the pseudopotential type identifier from the PSP header.
 *
 * The pseudopotential type is encoded as the first character of the PSP
 * header record and is expected to be a single digit character '0'–'9'.
 * This routine converts that character to its corresponding integer value.
 *
 * If the header does not contain a valid digit, the type defaults to 0,
 * matching legacy NWPW behavior and ensuring backward compatibility with
 * older or malformed pseudopotential files.
 *
 * This function intentionally performs no validation beyond the first
 * character and avoids ASCII arithmetic to preserve explicit, auditable
 * behavior.
 *
 * @param test  Pointer to the PSP header string.
 * @return      Integer pseudopotential type (0–9).
 */
static int convert_psp_type(char *test) {
  int psp_type = 0;
  if (test[0] == '0')
    psp_type = 0;
  if (test[0] == '1')
    psp_type = 1;
  if (test[0] == '2')
    psp_type = 2;
  if (test[0] == '3')
    psp_type = 3;
  if (test[0] == '4')
    psp_type = 4;
  if (test[0] == '5')
    psp_type = 5;
  if (test[0] == '6')
    psp_type = 6;
  if (test[0] == '7')
    psp_type = 7;
  if (test[0] == '8')
    psp_type = 8;
  if (test[0] == '9')
    psp_type = 9;

  return psp_type;
}


/*******************************************
 *                                         *
 *               read_vpwpup               *
 *                                         *
 *******************************************/
/**
 * @brief Read radial pseudopotential components (v_p, w_p, u_p) from file.
 *
 * This routine reads the radial grids and corresponding angular-momentum–
 * resolved pseudopotential components from a Hamann-style PSP file.
 * The data are read on the master task and broadcast to all MPI ranks.
 *
 * The input file is expected to contain:
 *
 *   (1) Radial grid and local / nonlocal potential components:
 *       - rho(i)
 *       - v_p(rho, l)
 *
 *   (2) Radial grid followed by projector and auxiliary functions:
 *       - rho(i)
 *       - w_p(rho, l)
 *       - u_p(rho, l)
 *
 * where l = 0..lmax0 in the file, but only components with l <= lmax
 * are retained and stored.
 *
 * The arrays are stored in column-major layout with the radial index
 * varying fastest:
 *
 *   vp[i + l*nrho], wp[i + l*nrho], up[i + l*nrho]
 *
 * This routine is used for Hamann-type pseudopotentials that require
 * auxiliary u_p functions (e.g., psp_type == 9).
 *
 * @param myparall Parallel context for MPI broadcast.
 * @param fp       Open file pointer to the PSP file.
 * @param nrho     Number of radial grid points.
 * @param lmax0    Maximum angular momentum present in the file.
 * @param lmax     Maximum angular momentum retained by the calculation.
 * @param rho      Radial grid (size nrho).
 * @param vp       Radial potential components v_p(r,l).
 * @param wp       Radial projector components w_p(r,l).
 * @param up       Radial auxiliary components u_p(r,l).
 */
static void read_vpwpup(Parallel *myparall, FILE *fp, int nrho, int lmax0,
                        int lmax, double *rho, double *vp, double *wp,
                        double *up) 
{
  int i, l;
  double xx;

  if (myparall->is_master()) {
    for (i = 0; i < nrho; ++i) {
      std::fscanf(fp, FMT1, &rho[i]);
      for (l = 0; l <= lmax0; ++l) {
        std::fscanf(fp, FMT1, &xx);
        if (l <= lmax)
          vp[i + l * nrho] = xx;
      }
    }
    for (i = 0; i < nrho; ++i) {
      std::fscanf(fp, FMT1, &rho[i]);
      for (l = 0; l <= lmax0; ++l) {
        std::fscanf(fp, FMT1, &xx);
        if (l <= lmax)
          wp[i + l * nrho] = xx;
      }
      for (l = 0; l <= lmax0; ++l) {
        std::fscanf(fp, FMT1, &xx);
        if (l <= lmax)
          up[i + l * nrho] = xx;
      }
    }
  }
  myparall->Brdcst_Values(0, 0, nrho, rho);
  myparall->Brdcst_Values(0, 0, nrho * (lmax + 1), vp);
  myparall->Brdcst_Values(0, 0, nrho * (lmax + 1), wp);
  myparall->Brdcst_Values(0, 0, nrho * (lmax + 1), up);
}



/*******************************************
 *                                         *
 *               read_vpwp                 *
 *                                         *
 *******************************************/
/**
 * @brief Read radial pseudopotential components (v_p and w_p) from file.
 *
 * This routine reads the radial grid and angular-momentum–resolved
 * pseudopotential data from a Hamann-style PSP file on the master task
 * and broadcasts the results to all MPI ranks.
 *
 * The input file is expected to contain:
 *
 *   (1) Radial grid and local / nonlocal potential components:
 *       - rho(i)
 *       - v_p(rho, l),  l = 0..lmax0
 *
 *   (2) Radial grid and projector components:
 *       - rho(i)
 *       - w_p(rho, l),  l = 0..(lmax0 + n_extra)
 *
 * Only components with l <= lmax are retained for v_p.  For the projector
 * functions w_p, any additional channels (l > lmax0) are stored in the
 * extended portion of the wp array, offset to maintain a contiguous
 * indexing scheme compatible with expanded KB projector sets.
 *
 * Arrays are stored in column-major layout with the radial index
 * varying fastest:
 *
 *   vp[i + l*nrho]
 *   wp[i + l*nrho]
 *
 * This routine is used for standard Hamann/Kleinman–Bylander
 * pseudopotentials that do not require auxiliary u_p functions.
 *
 * @param myparall Parallel context for MPI broadcast.
 * @param fp       Open file pointer to the PSP file.
 * @param nrho     Number of radial grid points.
 * @param lmax0    Maximum angular momentum present in the file.
 * @param lmax     Maximum angular momentum retained by the calculation.
 * @param n_extra  Number of additional projector channels beyond lmax0.
 * @param rho      Radial grid (size nrho).
 * @param vp       Radial potential components v_p(r,l).
 * @param wp       Radial projector components w_p(r,l), including extra channels.
 */
static void read_vpwp(Parallel *myparall, FILE *fp, int nrho, int lmax0,
                      int lmax, int n_extra, double *rho, double *vp,
                      double *wp) 
{
  int i, l;
  double xx;

  if (myparall->is_master()) {
    for (i = 0; i < nrho; ++i) {
      std::fscanf(fp, FMT1, &rho[i]);
      for (l = 0; l <= lmax0; ++l) {
        std::fscanf(fp, FMT1, &xx);
        if (l <= lmax)
          vp[i + l * nrho] = xx;
      }
    }
    for (i = 0; i < nrho; ++i) {
      std::fscanf(fp, FMT1, &rho[i]);
      for (l = 0; l <= (lmax0 + n_extra); ++l) {
        std::fscanf(fp, FMT1, &xx);
        if (l <= lmax)
          wp[i + l * nrho] = xx;
        if (l > lmax0)
          wp[i + (l + lmax - lmax0) * nrho] = xx;
      }
    }
  }
  myparall->Brdcst_Values(0, 0, nrho, rho);
  myparall->Brdcst_Values(0, 0, nrho * (lmax + 1), vp);
  myparall->Brdcst_Values(0, 0, nrho * (lmax + 1 + n_extra), wp);
}



/*******************************************
 *                                         *
 *               read_semicore             *
 *                                         *
 *******************************************/
/**
 * @brief Read optional semicore charge density data from a PSP file.
 *
 * This routine attempts to read an optional semicore charge-density section
 * from the pseudopotential file.  The presence of semicore data is detected
 * implicitly by attempting to read an additional scalar value beyond the
 * standard PSP record; reaching EOF indicates that no semicore information
 * is present.
 *
 * When present, the semicore section consists of:
 *
 *   - r_core           : semicore cutoff radius
 *   - 2*nrho entries   : (r, rho_sc(r)) pairs
 *
 * Only the semicore density values rho_sc(r) are retained; the radial
 * coordinate is assumed to match the primary rho grid.
 *
 * All file I/O is performed on the master task, and the resulting data
 * are broadcast to all MPI ranks.  If no semicore section is present,
 * isemicore is set to zero and no additional data are broadcast.
 *
 * This behavior mirrors the legacy NWPW PSP format and preserves
 * backward compatibility with older pseudopotential files.
 *
 * @param myparall     Parallel context for MPI broadcast.
 * @param fp           Open file pointer to the PSP file.
 * @param isemicore    Output flag: 1 if semicore data present, 0 otherwise.
 * @param rcore        Semicore cutoff radius (valid only if isemicore == 1).
 * @param nrho         Number of radial grid points.
 * @param semicore_r   Semicore charge density array (size 2*nrho).
 */
static void read_semicore(Parallel *myparall, FILE *fp, int *isemicore,
                          double *rcore, int nrho, double *semicore_r) {
  int i, isemicore0;
  double xx, yy, rcore0;

  isemicore0 = 0;
  rcore0 = 0.0;
  if (myparall->is_master()) {
    if (std::fscanf(fp, FMT1, &xx) != EOF) {
      rcore0 = xx;
      isemicore0 = 1;
      for (i = 0; i < (2 * nrho); ++i) {
        std::fscanf(fp, FMT2, &xx, &yy);
        semicore_r[i] = yy;
      }
    }
  }

  myparall->Brdcst_iValue(0, 0, &isemicore0);
  myparall->Brdcst_Values(0, 0, 1, &rcore0);
  if (isemicore0 > 0)
    myparall->Brdcst_Values(0, 0, 2 * nrho, semicore_r);
  *rcore = rcore0;
  *isemicore = isemicore0;
}

// static double dsum(int n, double *x, int incx)
//{
//    double stemp = 0.0;
//    for (int i=0; i<(n*incx); i+=incx) stemp += x[i];
//    return stemp;
// }

// static double simpson(int n, double *y, double h)
//{
//    int ne = n/2;
//    int no = ne+1;
//
//    double s = 3.0*dsum(no,y,2) + 4.0*dsum(ne,&y[1],2) - y[0] - y[n-1];
//    return (s*h/3.0);
// }



/*******************************************
 *                                         *
 *             generate_r3_matrix          *
 *                                         *
 *******************************************/
/*
   Computes the matrix elements for EFG tensor calculations

   r3_matrix(li,lj) = <uli|1/r3|ulj> - <wli|1/r3|wlj>
*/
/**
 * @brief Construct the radial ⟨1/r³⟩ matrix used in EFG tensor calculations.
 *
 * This routine computes matrix elements of the form
 *
 *   r3_matrix(l_i, l_j) =
 *        ⟨ u_{l_i} | 1/r³ | u_{l_j} ⟩
 *      − ⟨ w_{l_i} | 1/r³ | w_{l_j} ⟩
 *
 * where u_l(r) and w_l(r) are the auxiliary and projector radial functions
 * associated with a Hamann-type pseudopotential.  The subtraction removes
 * the corresponding projector contribution, yielding the effective
 * semicore-corrected operator required for electric field gradient (EFG)
 * evaluations.
 *
 * The radial integrals are evaluated numerically using Simpson integration
 * on the logarithmic radial grid, with the angular integration contributing
 * a factor of 4π.
 *
 * Only (l_i + l_j) > 0 terms are evaluated, consistent with the singular
 * behavior of the 1/r³ operator at the origin.
 *
 * The resulting matrix is symmetric in (l_i, l_j) and stored in a
 * flattened (lmax+1) × (lmax+1) layout.
 *
 * @param nrho       Number of radial grid points.
 * @param lmax       Maximum angular momentum channel.
 * @param drho       Radial grid spacing.
 * @param rho        Radial grid values.
 * @param wp         Projector radial functions w_l(r).
 * @param up         Auxiliary radial functions u_l(r).
 * @param r3_matrix  Output matrix of ⟨1/r³⟩ elements.
 */
static void generate_r3_matrix(int nrho, int lmax, double drho, double *rho,
                               double *wp, double *up, double *r3_matrix) {
  int i, li, lj;
  int lmax2 = (lmax + 1) * (lmax + 1);
  double coeff;
  double fourpi = 16.0 * atan(1.0);
  double *f = (double *)new double[nrho];

  for (i = 0; i < lmax2; ++i)
    r3_matrix[i] = 0.0;

  for (lj = 0; lj <= lmax; ++lj)
    for (li = lj; li <= lmax; ++li)
      if ((li + lj) > 0) {
        for (i = 0; i < nrho; ++i)
          f[i] = (up[i + li * nrho] * up[i + lj * nrho] -
                  wp[i + li * nrho] * wp[i + lj * nrho]) /
                 (rho[i] * rho[i] * rho[i]);
        coeff = fourpi * util_simpson(nrho, f, drho);
        r3_matrix[li + lj * lmax] = coeff;
        if (li != lj)
          r3_matrix[lj + li * lmax] = coeff;
      }

  delete[] f;
}

/* Constructors */

/*******************************************
 *                                         *
 *     Psp1d_Hamann::Psp1d_Hamann          *
 *                                         *
 *******************************************/
/**
 * @brief Construct and initialize a 1D Hamann-type pseudopotential.
 *
 * This constructor reads a Hamann-format pseudopotential file and initializes
 * all radial data, projector metadata, normalization constants, and optional
 * semicore and EFG-related quantities required for plane-wave DFT calculations.
 *
 * The following steps are performed:
 *
 * 1. **File parsing (master task only)**
 *    - Read atom label, pseudopotential type, valence charge, atomic mass
 *    - Read angular momentum limits (lmax0, lmax), local channel (locp)
 *    - Read local cutoff radii and projector expansion counts
 *    - Read radial grid definition (nrho, drho) and descriptive comment
 *
 * 2. **MPI broadcast**
 *    - All scalar metadata and radial dimensions are broadcast to all tasks
 *    - Ensures identical pseudopotential state across the parallel communicator
 *
 * 3. **Radial data allocation and reading**
 *    - Allocate rho(r), vp_l(r), wp_{l,n}(r)
 *    - For type-9 pseudopotentials, also read up_l(r) and build the ⟨1/r³⟩ matrix
 *      used in electric field gradient (EFG) calculations
 *
 * 4. **Semicore handling**
 *    - Optionally read semicore density data and mark semicore availability
 *
 * 5. **Nonlocal pseudopotential construction**
 *    - Convert vp_l(r) to nonlocal form by subtracting the local channel
 *
 * 6. **Projector normalization**
 *    - Compute normalization matrices for each angular momentum channel using
 *      Simpson integration
 *    - Invert normalization blocks analytically (n=1,2) or numerically (n>2)
 *
 * 7. **Projector bookkeeping**
 *    - Build flattened lists of projector quantum numbers:
 *        (n_prj, l_prj, m_prj, b_prj)
 *    - Ordering matches packed G-space layouts used elsewhere in PWDFT
 *
 * The resulting object contains a complete, self-consistent representation of
 * the pseudopotential suitable for local, nonlocal, stress, force, and EFG
 * calculations in reciprocal space.
 *
 * @param myparall     Parallel communication object.
 * @param psp_name     Path to the pseudopotential file.
 * @param psp_version  Internal version identifier for compatibility control.
 */
Psp1d_Hamann::Psp1d_Hamann(Parallel *myparall, const char *psp_name,
                           const int psp_version) {
  double xx;
  FILE *fp;
  int nn;

  if (myparall->is_master()) {
    fp = std::fopen(psp_name, "r");
    std::fscanf(fp, "%s", atom);
    psp_type = convert_psp_type(atom);
    if (psp_type > 0)
      std::fscanf(fp, "%s", atom);

    std::fscanf(fp, FMT1, &zv);
    std::fscanf(fp, FMT1, &amass);
    std::fscanf(fp, "%d", &lmax0);
    std::fscanf(fp, "%d", &lmax);
    std::fscanf(fp, "%d", &locp);
    std::fscanf(fp, FMT1, &rlocal);
    for (int l = 0; l <= lmax0; ++l)
      std::fscanf(fp, FMT1, &rc[l]);

    n_extra = 0;
    for (int l = 0; l <= lmax0; ++l)
      n_expansion[l] = 1;
    if (psp_type == 2) {
      std::fscanf(fp, "%d", &n_extra);
      for (int l = 0; l <= lmax0; ++l) {
        std::fscanf(fp, "%d", &nn);
        n_expansion[l] = nn;
      }
    }

    std::fscanf(fp, "%d", &nrho);
    std::fscanf(fp, FMT1, &drho);
    std::fscanf(fp, " %79[^\n]", comment);
    if (lmax > lmax0)
      lmax = lmax0;
    if (lmax < 0)
      lmax = lmax0;
    if (locp > lmax)
      locp = lmax;
    if (locp < 0)
      locp = lmax;
  }
  myparall->Brdcst_cValues(0, 0, 2, atom);
  // myparall->Brdcst_iValue(0,0,&ihasae);
  myparall->Brdcst_iValue(0, 0, &psp_type);
  myparall->Brdcst_Values(0, 0, 1, &zv);
  myparall->Brdcst_Values(0, 0, 1, &amass);
  myparall->Brdcst_iValue(0, 0, &lmax0);
  myparall->Brdcst_iValue(0, 0, &lmax);
  myparall->Brdcst_iValue(0, 0, &locp);
  myparall->Brdcst_Values(0, 0, 1, &rlocal);
  myparall->Brdcst_Values(0, 0, (lmax0 + 1), rc);
  myparall->Brdcst_iValue(0, 0, &nrho);
  myparall->Brdcst_Values(0, 0, 1, &drho);
  myparall->Brdcst_cValues(0, 0, 80, comment);

  myparall->Brdcst_iValue(0, 0, &n_extra);
  myparall->Brdcst_iValues(0, 0, (lmax0 + 1), n_expansion);
  nprj = 0;
  nmax = 1;
  for (int l = 0; l <= (lmax); ++l)
    if (l != locp) {
      nprj += n_expansion[l] * (2 * l + 1);
      if (n_expansion[l] > nmax)
        nmax = n_expansion[l];
    }

  rho = (double *)new double[nrho];
  vp = (double *)new double[(lmax + 1) * nrho];
  wp = (double *)new double[(lmax + 1 + n_extra) * nrho];
  vnlnrm = (double *)new double[nmax * nmax * (lmax + 1)];
  rho_sc_r = (double *)new double[2 * nrho];
  if (psp_type == 9) {
    up = (double *)new double[(lmax + 1) * nrho];
    r3_matrix = (double *)new double[(lmax + 1) * (lmax + 1)];
    read_vpwpup(myparall, fp, nrho, lmax0, lmax, rho, vp, wp, up);
    generate_r3_matrix(nrho, lmax, drho, rho, wp, up, r3_matrix);
  } else
    read_vpwp(myparall, fp, nrho, lmax0, lmax, n_extra, rho, vp, wp);
  read_semicore(myparall, fp, &isemicore, &rcore, nrho, rho_sc_r);
  semicore = (isemicore == 1);

  /* Define non-local pseudopotential  */
  for (auto l = 0; l <= lmax; ++l)
    if (l != locp)
      for (auto i = 0; i < nrho; ++i)
        vp[i + l * nrho] = vp[i + l * nrho] - vp[i + locp * nrho];

  /* set up indx(n,l) --> to wp */
  int indx[5 * 4];
  int nb = lmax + 1;
  for (auto l = 0; l <= lmax; ++l) {
    indx[l * 5] = l;
    for (auto n1 = 1; n1 < n_expansion[l]; ++n1) {
      indx[n1 + l * 5] = nb;
      ++nb;
    }
  }

  // version = 3;
  version = psp_version;
  /* Normarization constants */
  double a;
  double *f = new double[nrho];
  for (auto l = 0; l <= lmax; ++l) {
    if (l != locp) {
      for (auto n2 = 0; n2 < n_expansion[l]; ++n2) {
        for (auto i = 0; i < nrho; ++i)
          f[i] = vp[i + l * nrho] * wp[i + indx[n2 + l * 5] * nrho] *
                 wp[i + indx[n2 + l * 5] * nrho];

        a = util_simpson(nrho, f, drho);
        vnlnrm[n2 + n2 * nmax + l * nmax * nmax] = a;

        for (auto n1 = n2 + 1; n1 < n_expansion[l]; ++n1) {
          for (auto i = 0; i < nrho; ++i)
            f[i] = vp[i + l * nrho] * wp[i + indx[n1 + l * 5] * nrho] *
                   wp[i + indx[n2 + l * 5] * nrho];
          a = util_simpson(nrho, f, drho);
          vnlnrm[n1 + n2 * nmax + l * nmax * nmax] = a;
          vnlnrm[n2 + n1 * nmax + l * nmax * nmax] = a;
        }
      }
      if (n_expansion[l] == 1) {
        vnlnrm[l * nmax * nmax] = 1 / a;
      } else if (n_expansion[l] == 2) {
        double d =
            vnlnrm[l * nmax * nmax] * vnlnrm[1 + 1 * nmax + l * nmax * nmax] -
            vnlnrm[0 + 1 * nmax + l * nmax * nmax] *
                vnlnrm[1 + 0 * nmax + l * nmax * nmax];
        double q = vnlnrm[l * nmax * nmax];
        vnlnrm[l * nmax * nmax] = vnlnrm[1 + 1 * nmax + l * nmax * nmax] / d;
        vnlnrm[1 + 1 * nmax + l * nmax * nmax] = q / d;
        vnlnrm[0 + 1 * nmax + l * nmax * nmax] =
            -vnlnrm[0 + 1 * nmax + l * nmax * nmax] / d;
        vnlnrm[1 + 0 * nmax + l * nmax * nmax] =
            -vnlnrm[1 + 0 * nmax + l * nmax * nmax] / d;
      } else {
        util_matinvert(n_expansion[l], nmax, &(vnlnrm[l * nmax * nmax]));
      }
    } else {
      for (int n2 = 0; n2 < nmax; ++n2)
        for (int n1 = n2; n1 < nmax; ++n1) {
          vnlnrm[n1 + n2 * nmax + l * nmax * nmax] = 0.0;
          vnlnrm[n2 + n1 * nmax + l * nmax * nmax] = 0.0;
        }
    }
  }
  delete[] f;

  /* define n_prj, l_prj, m_prj, b_prj */
  if (nprj > 0) {
    n_prj = new int[nprj];
    l_prj = new int[nprj];
    m_prj = new int[nprj];
    b_prj = new int[nprj];

    int lcount = nprj;
    /* f projectors */
    if ((locp != 3) && (lmax > 2))
      for (auto n = 0; n < n_expansion[3]; ++n) {
        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 3;
        m_prj[lcount] = -3;
        b_prj[lcount] = indx[n + 3 * 5];

        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 3;
        m_prj[lcount] = -2;
        b_prj[lcount] = indx[n + 3 * 5];

        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 3;
        m_prj[lcount] = -1;
        b_prj[lcount] = indx[n + 3 * 5];

        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 3;
        m_prj[lcount] = 0;
        b_prj[lcount] = indx[n + 3 * 5];

        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 3;
        m_prj[lcount] = 1;
        b_prj[lcount] = indx[n + 3 * 5];

        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 3;
        m_prj[lcount] = 2;
        b_prj[lcount] = indx[n + 3 * 5];

        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 3;
        m_prj[lcount] = 3;
        b_prj[lcount] = indx[n + 3 * 5];
      }

    /* d projectors */
    if ((locp != 2) && (lmax > 1))
      for (auto n = 0; n < n_expansion[2]; ++n) {
        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 2;
        m_prj[lcount] = -2;
        b_prj[lcount] = indx[n + 2 * 5];

        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 2;
        m_prj[lcount] = -1;
        b_prj[lcount] = indx[n + 2 * 5];

        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 2;
        m_prj[lcount] = 0;
        b_prj[lcount] = indx[n + 2 * 5];

        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 2;
        m_prj[lcount] = 1;
        b_prj[lcount] = indx[n + 2 * 5];

        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 2;
        m_prj[lcount] = 2;
        b_prj[lcount] = indx[n + 2 * 5];
      }

    /* p projectors */
    if ((locp != 1) && (lmax > 0))
      for (auto n = 0; n < n_expansion[1]; ++n) {
        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 1;
        m_prj[lcount] = -1;
        b_prj[lcount] = indx[n + 1 * 5];

        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 1;
        m_prj[lcount] = 0;
        b_prj[lcount] = indx[n + 1 * 5];

        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 1;
        m_prj[lcount] = 1;
        b_prj[lcount] = indx[n + 1 * 5];
      }

    /* s projectors */
    if (locp != 0)
      for (auto n = 0; n < n_expansion[0]; ++n) {
        --lcount;
        n_prj[lcount] = n + 1;
        l_prj[lcount] = 0;
        m_prj[lcount] = 0;
        b_prj[lcount] = indx[n + 0 * 5];
      }
  }
}



/*******************************************
 *                                         *
 *     Psp1d_Hamann::vpp_generate_ray      *
 *                                         *
 *******************************************/
/**
 * @brief Generate radial (ray-based) pseudopotential kernels in reciprocal space.
 *
 * This routine computes 1D reciprocal-space representations of the local and
 * nonlocal pseudopotential components as functions of |G|, suitable for later
 * interpolation onto the full 3D plane-wave grid.
 *
 * Specifically, it generates:
 *
 *  (1) Local pseudopotential kernel
 *      ------------------------------------------------------------
 *      vl_ray(|G|) = 4π / |G| ∫ dr r v_loc(r) sin(|G| r)
 *                    + analytic G→0 correction (version-dependent)
 *
 *      Two supported forms are handled:
 *        - version 3 : Coulomb tail treated analytically
 *        - version 4 : Error-function screened Coulomb form
 *
 *  (2) Nonlocal Kleinman–Bylander projector kernels
 *      ------------------------------------------------------------
 *      vnl_ray_p(|G|) for each projector channel p = (n,l),
 *      computed via analytic Fourier–Bessel transforms of the
 *      radial projector functions wp(r) and channel potentials vp_l(r).
 *
 *      The resulting array is indexed as:
 *
 *        vnl_ray[ p*nray + k ] = radial kernel for projector p at |G_k|
 *
 *      where p is the flattened (n,l) index defined by indx(n,l).
 *
 *  (3) Semicore density kernels (optional)
 *      ------------------------------------------------------------
 *      If semicore charge density is present, two isotropic radial kernels
 *      are generated:
 *
 *        rho_sc_k_ray[0*nray + k] : density-like term
 *        rho_sc_k_ray[1*nray + k] : derivative-related term
 *
 *      These are later expanded into scalar and vector components on the
 *      3D reciprocal-space grid.
 *
 * Implementation details:
 *  - Radial integrals are evaluated using Simpson integration.
 *  - Angular dependence is handled analytically using spherical Bessel
 *    function identities for each angular momentum channel (s, p, d, f).
 *  - Parallelization is over ray index k; results are summed using
 *    Vector_SumAll to assemble complete ray arrays.
 *  - G = 0 limits are handled explicitly to enforce regular behavior and
 *    translational invariance.
 *
 * The ray-formatted outputs produced here are consumed by:
 *   - cpp_generate_local_spline()
 *   - cpp_generate_nonlocal_spline()
 *   - stress-generation variants (vpp2 / cpp2 paths)
 *
 * @param myparall        Parallel communication object.
 * @param nray            Number of |G| samples in the ray grid.
 * @param G_ray           Array of ray magnitudes |G|.
 * @param vl_ray          Output local pseudopotential kernel vl(|G|).
 * @param vnl_ray         Output nonlocal projector kernels vnl_p(|G|).
 * @param rho_sc_k_ray    Output semicore density kernels (optional).
 */
void Psp1d_Hamann::vpp_generate_ray(Parallel *myparall, int nray, double *G_ray,
                                    double *vl_ray, double *vnl_ray,
                                    double *rho_sc_k_ray) 
{
   /* set up indx(n,l) --> to wp */
   int indx[5 * 4];
   int nb = lmax + 1;
   for (auto l = 0; l <= lmax; ++l) {
     indx[l * 5] = l;
     for (auto n1 = 1; n1 < n_expansion[l]; ++n1) {
       indx[n1 + l * 5] = nb;
       ++nb;
     }
   }
 
   double pi = 4.00 * atan(1.0);
   double twopi = 2.0 * pi;
   double forpi = 4.0 * pi;
 
   double P0 = sqrt(forpi);
   double P1 = sqrt(3.0 * forpi);
   double P2 = sqrt(15.0 * forpi);
   double P3 = sqrt(105.0 * forpi);
 
   double zero = 0.0;
   int izero = 0;
   int ione = 1;
   int nray2 = 2 * nray;
   int lmaxnray = (lmax + 1 + n_extra) * nray;
 
   double q;
   double *cs = new double[nrho];
   double *sn = new double[nrho];
   double *f = new double[nrho];
   double a, xx;
 
   memset(vl_ray, 0, nray * sizeof(double));
   memset(vnl_ray, 0, lmaxnray * sizeof(double));
   memset(rho_sc_k_ray, 0, nray2 * sizeof(double));
 
   for (auto k1 = (1 + myparall->taskid()); k1 < nray; k1 += myparall->np()) {
     q = G_ray[k1];
     for (auto i = 0; i < nrho; ++i) {
       cs[i] = cos(q * rho[i]);
       sn[i] = sin(q * rho[i]);
     }
 
     /* h projectors */
     /* h projectors */
     /* f projectors */
     if ((locp != 3) && (lmax > 2)) {
       for (auto n = 0; n < n_expansion[3]; ++n) {
         f[0] = 0.0;
         for (auto i = 1; i < nrho; ++i) {
           xx = q * rho[i];
           a = sn[i] / xx;
           a = 15.0 * (a - cs[i]) / (xx * xx) - 6 * a + cs[i];
           f[i] = a * wp[i + indx[n + 3 * 5] * nrho] * vp[i + 3 * nrho];
         }
         vnl_ray[k1 + indx[n + 3 * 5] * nray] =
             P3 * util_simpson(nrho, f, drho) / q;
       }
     }
 
     /* d projectors */
     if ((locp != 2) && (lmax > 1)) {
       for (auto n = 0; n < n_expansion[2]; ++n) {
         f[0] = 0.0;
         for (auto i = 1; i < nrho; ++i) {
           a = 3.0 * (sn[i] / (q * rho[i]) - cs[i]) / (q * rho[i]) - sn[i];
           f[i] = a * wp[i + indx[n + 2 * 5] * nrho] * vp[i + 2 * nrho];
         }
         vnl_ray[k1 + indx[n + 2 * 5] * nray] =
             P2 * util_simpson(nrho, f, drho) / q;
       }
     }
 
     /* p projectors */
     if ((locp != 1) && (lmax > 0)) {
       for (auto n = 0; n < n_expansion[1]; ++n) {
         f[0] = 0.0;
         for (auto i = 1; i < nrho; ++i) {
           a = (sn[i] / (q * rho[i]) - cs[i]);
           f[i] = a * wp[i + indx[n + 1 * 5] * nrho] * vp[i + 1 * nrho];
         }
         vnl_ray[k1 + indx[n + 1 * 5] * nray] =
             P1 * util_simpson(nrho, f, drho) / q;
       }
     }
 
     /* s projectors */
     if (locp != 0) {
       for (auto n = 0; n < n_expansion[0]; ++n) {
         for (auto i = 0; i < nrho; ++i)
           f[i] = sn[i] * wp[i + indx[n + 0 * 5] * nrho] * vp[i + 0 * nrho];
         vnl_ray[k1 + indx[n + 0 * 5] * nray] =
             P0 * util_simpson(nrho, f, drho) / q;
       }
     }
 
     /* local */
     if (version == 3) {
       for (auto i = 0; i < nrho; ++i)
         f[i] = rho[i] * vp[i + locp * nrho] * sn[i];
       vl_ray[k1] = util_simpson(nrho, f, drho) * forpi / q -
                    zv * forpi * cs[nrho - 1] / (q * q);
     } else if (version == 4) {
       for (auto i = 0; i < nrho; ++i)
         f[i] = (rho[i] * vp[i + locp * nrho] + zv * std::erf(rho[i] / rlocal)) *
                sn[i];
       vl_ray[k1] = util_simpson(nrho, f, drho) * forpi / q;
     }
 
     /* semicore density */
     if (semicore) {
       for (auto i = 0; i < nrho; ++i)
         f[i] = rho[i] * sqrt(rho_sc_r[i]) * sn[i];
       rho_sc_k_ray[k1] = util_simpson(nrho, f, drho) * forpi / q;
 
       for (auto i = 0; i < nrho; ++i)
         f[i] = (sn[i] / (q * rho[i]) - cs[i]) * rho_sc_r[i + nrho] * rho[i];
       rho_sc_k_ray[k1 + nray] = util_simpson(nrho, f, drho) * forpi / q;
     }
   }
   myparall->Vector_SumAll(0, 2 * nray, rho_sc_k_ray);
   myparall->Vector_SumAll(0, nray, vl_ray);
   myparall->Vector_SumAll(0, lmaxnray, vnl_ray);
 
   /* G==0 local */
   if (version == 3) {
     for (auto i = 0; i < nrho; ++i) {
       f[i] = vp[i + locp * nrho] * rho[i] * rho[i];
     }
     vl_ray[0] = forpi * util_simpson(nrho, f, drho) +
                 twopi * zv * rho[nrho - 1] * rho[nrho - 1];
   } else if (version == 4) {
     for (auto i = 0; i < nrho; ++i)
       f[i] = (vp[i + locp * nrho] * rho[i] + zv * std::erf(rho[i] / rlocal)) *
              rho[i];
     vl_ray[0] = forpi * util_simpson(nrho, f, drho);
   }
 
   /* G==0 semicore */
   if (semicore) {
     for (auto i = 0; i < nrho; ++i)
       f[i] = sqrt(rho_sc_r[i]) * rho[i] * rho[i];
     rho_sc_k_ray[0] = forpi * util_simpson(nrho, f, drho);
     rho_sc_k_ray[0 + nray] = 0.0;
   }
 
   /* G==0 vnl */
   for (auto l = 0; l <= lmax; ++l)
     for (auto n = 0; n < n_expansion[l]; ++n)
       vnl_ray[0 + indx[n + l * 5] * nray] = 0.0;
 
   /* only j0 is non-zero at zero */
   if (locp != 0)
     for (auto n = 0; n < n_expansion[0]; ++n) {
       for (auto i = 0; i < nrho; ++i)
         f[i] = rho[i] * wp[i + indx[n + 0 * 5] * nrho] * vp[i + 0 * nrho];
       vnl_ray[0 + indx[n + 0 * 5] * nray] = P0 * util_simpson(nrho, f, drho);
     }
 
   delete[] f;
   delete[] sn;
   delete[] cs;
}


/*******************************************
 *                                         *
 *   Psp1d_Hamann::vpp_generate_spline     *
 *                                         *
 *******************************************/
/**
 * @brief Interpolate ray-based pseudopotential kernels onto the packed 3D G-grids.
 *
 * This routine takes radial (|G|-dependent) pseudopotential kernels generated
 * by vpp_generate_ray() and interpolates them onto the full reciprocal-space
 * plane-wave grids used in the calculation.
 *
 * Specifically, it constructs cubic splines in |G| and evaluates them on:
 *
 *  (1) The G = 0 (scalar) grid for the local potential
 *      ------------------------------------------------------------
 *      vl(G) is generated by spline interpolation of vl_ray(|G|),
 *      with explicit handling of the G → 0 limit.
 *
 *  (2) The G ≠ 0 (vector) grid for nonlocal KB projectors
 *      ------------------------------------------------------------
 *      For each nonlocal projector channel p = (n,l), the scalar
 *      radial kernel vnl_ray_p(|G|) is combined with analytic angular
 *      factors to produce the Cartesian projector components stored
 *      in vnl.
 *
 *      The ordering and angular forms match the packed projector
 *      conventions used throughout PWDFT.
 *
 *  (3) Semicore density kernels (optional)
 *      ------------------------------------------------------------
 *      If a semicore charge density is present, both scalar and
 *      vector components are generated:
 *
 *        rho_sc_k[0]              : isotropic component
 *        rho_sc_k[1..3]           : directional components ∝ Ĝ
 *
 * Implementation notes:
 *  - Cubic splines are constructed once per radial channel and reused
 *    across all G-points.
 *  - Angular dependence for s, p, d, and f channels is handled analytically.
 *  - G = 0 values are treated explicitly to avoid singular behavior and
 *    to enforce translational invariance.
 *
 * This routine produces the fully expanded reciprocal-space pseudopotential
 * operators required for applying the local and nonlocal pseudopotential
 * in plane-wave Hamiltonian operations.
 *
 * @param mygrid         Plane-wave reciprocal-space grid.
 * @param nray           Number of radial |G| samples.
 * @param G_ray          Radial grid of |G| values.
 * @param vl_ray         Ray-formatted local pseudopotential kernel.
 * @param vnl_ray        Ray-formatted nonlocal projector kernels.
 * @param rho_sc_k_ray   Ray-formatted semicore density kernels (optional).
 * @param vl             Output packed local potential on the G=0 grid.
 * @param vnl            Output packed nonlocal projector components.
 * @param rho_sc_k       Output packed semicore density components.
 */
void Psp1d_Hamann::vpp_generate_spline(PGrid *mygrid, int nray, double *G_ray,
                                       double *vl_ray, double *vnl_ray,
                                       double *rho_sc_k_ray, double *vl,
                                       double *vnl, double *rho_sc_k) 
{

  /* set up indx(n,l) --> to wp */
  int indx[5 * 4];
  int nb = lmax + 1;
  for (auto l = 0; l <= lmax; ++l) {
    indx[l * 5] = l;
    for (auto n1 = 1; n1 < n_expansion[l]; ++n1) {
      indx[n1 + l * 5] = nb;
      ++nb;
    }
  }

  double pi = 4.00 * atan(1.0);

  /* allocate spline grids */
  double *vl_splineray = new double[nray];
  double *vnl_splineray = new double[(lmax + 1 + n_extra) * nray];
  double *rho_sc_k_splineray = new double[2 * nray];
  double *tmp_splineray = new double[nray];

  /* setup cubic bsplines */
  double dG = G_ray[2] - G_ray[1];

  /* five point formula */
  double yp1 = (-50.00 * vl_ray[1] + 96.00 * vl_ray[2] - 72.00 * vl_ray[3] +
                32.00 * vl_ray[4] - 6.00 * vl_ray[5]) /
               (24.00 * dG);
  util_spline(&(G_ray[1]), &(vl_ray[1]), nray - 1, yp1, 0.00,
              &(vl_splineray[1]), tmp_splineray);

  for (auto l = 0; l <= lmax; ++l)
    if (l != locp)
      for (auto n = 0; n < n_expansion[l]; ++n)
        util_spline(G_ray, &(vnl_ray[indx[n + 5 * l] * nray]), nray, 0.00, 0.00,
                    &(vnl_splineray[indx[n + 5 * l] * nray]), tmp_splineray);

  if (semicore) {
    util_spline(G_ray, rho_sc_k_ray, nray, 0.00, 0.00, rho_sc_k_splineray,
                tmp_splineray);
    util_spline(G_ray, &(rho_sc_k_ray[nray]), nray, 0.00, 0.00,
                &(rho_sc_k_splineray[nray]), tmp_splineray);
  }

  double q, qx, qy, qz, xx;
  double *gx, *gy, *gz;
  int npack0 = mygrid->npack(0);
  int npack1 = mygrid->npack(1);
  int nx, lcount;
  mygrid->t_pack_nzero(0, 1, vl);
  mygrid->t_pack_nzero(1, nprj, vnl);
  if (semicore)
    mygrid->t_pack_nzero(0, 4, rho_sc_k);

  /* generate vl and rho_sc_k */
  gx = mygrid->Gpackxyz(0, 0);
  gy = mygrid->Gpackxyz(0, 1);
  gz = mygrid->Gpackxyz(0, 2);
  for (auto k = 0; k < npack0; ++k) {
    qx = gx[k];
    qy = gy[k];
    qz = gz[k];
    q = sqrt(qx * qx + qy * qy + qz * qz);
    nx = (int)floor(q / dG);

    if (q > 1.0e-9) {
      qx /= q;
      qy /= q;
      qz /= q;
      vl[k] = util_splint(&(G_ray[1]), &(vl_ray[1]), &(vl_splineray[1]),
                          nray - 1, nx, q);
      if (semicore) {
        rho_sc_k[k] =
            util_splint(G_ray, rho_sc_k_ray, rho_sc_k_splineray, nray, nx, q);
        xx = util_splint(G_ray, &(rho_sc_k_ray[nray]),
                         &(rho_sc_k_splineray[nray]), nray, nx, q);
        rho_sc_k[k + npack0] = xx * qx;
        rho_sc_k[k + 2 * npack0] = xx * qy;
        rho_sc_k[k + 3 * npack0] = xx * qz;
      }
    } else {

      vl[k] = vl_ray[0];
      if (semicore) {
        rho_sc_k[k] = rho_sc_k_ray[0];
        rho_sc_k[k + npack0] = 0.0;
        rho_sc_k[k + 2 * npack0] = 0.0;
        rho_sc_k[k + 3 * npack0] = 0.0;
      }
    }
  }

  /* generate vnl */
  gx = mygrid->Gpackxyz(1, 0);
  gy = mygrid->Gpackxyz(1, 1);
  gz = mygrid->Gpackxyz(1, 2);

  for (auto k = 0; k < npack1; ++k) {
    qx = gx[k];
    qy = gy[k];
    qz = gz[k];
    q = sqrt(qx * qx + qy * qy + qz * qz);
    nx = (int)floor(q / dG);

    if (q > 1.0e-9) {
      qx /= q;
      qy /= q;
      qz /= q;
      lcount = nprj;

      /* f projectors */

      if ((locp != 3) && (lmax > 2))
        for (auto n = 0; n < n_expansion[3]; ++n) {
          xx = util_splint(G_ray, &(vnl_ray[indx[n + 3 * 5] * nray]),
                           &(vnl_splineray[indx[n + 3 * 5] * nray]), nray, nx,
                           q);
          --lcount;
          vnl[k + lcount * npack1] =
              xx * qy * (3.00 * (1.00 - qz * qz) - 4.00 * qy * qy) /
              sqrt(24.00);
          --lcount;
          vnl[k + lcount * npack1] = xx * qx * qy * qz;
          --lcount;
          vnl[k + lcount * npack1] =
              xx * qy * (5.00 * qz * qz - 1.00) / sqrt(40.00);
          --lcount;
          vnl[k + lcount * npack1] =
              xx * qz * (5.00 * qz * qz - 3.00) / sqrt(60.00);
          --lcount;
          vnl[k + lcount * npack1] =
              xx * qx * (5.00 * qz * qz - 1.00) / sqrt(40.00);
          --lcount;
          vnl[k + lcount * npack1] = xx * qz * (qx * qx - qy * qy) / 2.00;
          --lcount;
          vnl[k + lcount * npack1] =
              xx * qx * (4.00 * qx * qx - 3.00 * (1.00 - qz * qz)) /
              sqrt(24.00);
        }

      /* d projectors */
      if ((locp != 2) && (lmax > 1))
        for (auto n = 0; n < n_expansion[2]; ++n) {
          xx = util_splint(G_ray, &(vnl_ray[indx[n + 2 * 5] * nray]),
                           &(vnl_splineray[indx[n + 2 * 5] * nray]), nray, nx,
                           q);
          --lcount;
          vnl[k + lcount * npack1] = xx * qx * qy;
          --lcount;
          vnl[k + lcount * npack1] = xx * qy * qz;
          --lcount;
          vnl[k + lcount * npack1] =
              xx * (3.00 * qz * qz - 1.00) / (2.00 * sqrt(3.00));
          --lcount;
          vnl[k + lcount * npack1] = xx * qz * qx;
          --lcount;
          vnl[k + lcount * npack1] = xx * (qx * qx - qy * qy) / (2.00);
        }

      /* p projectors */
      if ((locp != 1) && (lmax > 0))
        for (auto n = 0; n < n_expansion[1]; ++n) {
          xx = util_splint(G_ray, &(vnl_ray[indx[n + 1 * 5] * nray]),
                           &(vnl_splineray[indx[n + 1 * 5] * nray]), nray, nx,
                           q);
          --lcount;
          vnl[k + lcount * npack1] = xx * qy;
          --lcount;
          vnl[k + lcount * npack1] = xx * qz;
          --lcount;
          vnl[k + lcount * npack1] = xx * qx;
        }

      /* s projectors */
      if (locp != 0)
        for (auto n = 0; n < n_expansion[0]; ++n) {
          xx = util_splint(G_ray, &(vnl_ray[indx[n + 0 * 5] * nray]),
                           &(vnl_splineray[indx[n + 0 * 5] * nray]), nray, nx,
                           q);
          --lcount;
          vnl[k + lcount * npack1] = xx;
        }
    } else {
      for (auto l = 0; l < nprj; ++l)
        vnl[k + l * npack1] = 0.0;

      /* only j0 is non-zero at zero */
      if (locp != 0)
        for (auto n = 0; n < n_expansion[0]; ++n)
          vnl[k + indx[n + 0 * 5] * npack1] =
              vnl_ray[0 + indx[n + 0 * 5] * nray];
    }
  }

  /*  deallocate spineray formatted grids */
  delete[] tmp_splineray;
  delete[] rho_sc_k_splineray;
  delete[] vnl_splineray;
  delete[] vl_splineray;
}


/*******************************************
 *                                         *
 * Psp1d_Hamann::vpp2_generate_stress_ray  *
 *                                         *
 *******************************************/
/**
 * @brief Generate ray-formatted pseudopotential stress kernels in G-space.
 *
 * This routine computes the reciprocal-space radial kernels needed to assemble
 * the pseudopotential contribution to the stress tensor.  The kernels are
 * produced as functions of the ray magnitude |G| (a 1D radial grid in
 * reciprocal space) and are later interpolated onto the full 3D plane-wave
 * grid for contraction into the Cartesian stress components.
 *
 * The generated contributions are:
 *
 *  (1) Nonlocal (Kleinman–Bylander) stress kernels
 *      ------------------------------------------------------------
 *      The KB nonlocal stress contribution can be written as:
 *
 *        σ^{NL}_{αβ} = Σ_l Σ_G [ A_l(G) δ_{αβ} + B_l(G) (G_α G_β)/|G| ] .
 *
 *      This routine generates the radial functions A_l(G) and B_l(G) for each
 *      projector channel (n,l) and stores them in dvnl_ray as:
 *
 *        dvnl_ray[ 0*nray + p*nray + k ] -> A_p(|G_k|)   (isotropic / metric term)
 *        dvnl_ray[ 1*nray + p*nray + k ] -> B_p(|G_k|)   (anisotropic term)
 *
 *      where p is the flattened projector-channel index produced by indx(n,l),
 *      and k indexes the ray sample (k=0..nray-1).
 *
 *  (2) Local pseudopotential stress kernel
 *      ------------------------------------------------------------
 *      The local part contributes an isotropic stress kernel:
 *
 *        σ^{L}_{αβ} = Σ_G A_loc(G) δ_{αβ},
 *
 *      where A_loc(G) is generated here as dvl_ray(G).
 *
 *  (3) Semicore density stress kernel (optional)
 *      ------------------------------------------------------------
 *      If a semicore charge density is present, an additional isotropic kernel
 *      is required for consistent stress evaluation:
 *
 *        σ^{SC}_{αβ} = Σ_G A_sc(G) δ_{αβ},
 *
 *      where A_sc(G) is generated here as rho_sc_k_ray(G).
 *
 * Implementation notes:
 *  - Kernels are computed by analytic Fourier–Bessel transforms of radial
 *    quantities (Hamann/Kleinman–Bylander form) and evaluated on the ray grid.
 *  - Parallelization is over ray index k; results are summed with Vector_SumAll.
 *  - G=0 components are explicitly set to zero where required to avoid
 *    singular terms and enforce translational invariance.
 */
void Psp1d_Hamann::vpp2_generate_stress_ray(Parallel *myparall,
                                            int nray, double *G_ray,
                                            double *dvl_ray, double *dvnl_ray,
                                            double *rho_sc_k_ray)
{
   /* set up indx(n,l) --> wp block index */
   int indx[5 * 4];
   int nb = lmax + 1;
   for (int l = 0; l <= lmax; ++l) {
     indx[l * 5] = l;
     for (int n1 = 1; n1 < n_expansion[l]; ++n1) {
       indx[n1 + l * 5] = nb;
       ++nb;
     }
   }
 
   const double pi    = 4.0 * atan(1.0);
   const double twopi = 2.0 * pi;
   const double forpi = 4.0 * pi;
 
   const double P0 = sqrt(forpi);
   const double P1 = sqrt(3.0   * forpi);
   const double P2 = sqrt(15.0  * forpi);
   const double P3 = sqrt(105.0 * forpi);
 
   const int nray2     = 2 * nray;
   const int nproj_blk = (lmax + 1 + n_extra);             // ray blocks per channel
   const int lmaxnray  = nproj_blk * nray;                 // A-block size (or B-block size)
 
   double *cs = new double[nrho];
   double *sn = new double[nrho];
   double *f  = new double[nrho];
 
   std::memset(dvl_ray,      0, nray      * sizeof(double));
   std::memset(dvnl_ray,     0, 2*lmaxnray * sizeof(double));
   std::memset(rho_sc_k_ray, 0, nray2     * sizeof(double));
 
   /* parallel loop over ray points, skipping k=0 (G=0 handled later) */
   for (int k1 = 1 + myparall->taskid(); k1 < nray; k1 += myparall->np()) {
 
     const double q = G_ray[k1];
 
     for (int i = 0; i < nrho; ++i) {
       cs[i] = std::cos(q * rho[i]);
       sn[i] = std::sin(q * rho[i]);
     }
 
     /* f projectors (l=3) */
     if ((locp != 3) && (lmax > 2)) {
       for (int n = 0; n < n_expansion[3]; ++n) {
 
         f[0] = 0.0;
         for (int i = 1; i < nrho; ++i) {
           const double xx = q * rho[i];
           double a = sn[i] / xx;
           a = 15.0 * (a - cs[i]) / (xx * xx) - 6.0 * a + cs[i];
           f[i] = a * wp[i + indx[n + 3*5] * nrho] * vp[i + 3*nrho];
         }
         dvnl_ray[0*lmaxnray + k1 + indx[n + 3*5]*nray] = P3 * util_simpson(nrho, f, drho) / q;
 
         f[0] = 0.0;
         for (int i = 1; i < nrho; ++i) {
           const double xx = q * rho[i];
           double a =
               -60.0 * sn[i] / (xx*xx*xx * q*q)
               +60.0  * cs[i] / (xx*xx     * q*q)
               +27.0  * sn[i] / (xx        * q*q)
               -7.0   * cs[i] / (q*q)
               -rho[i]* sn[i] / q;
           f[i] = a * wp[i + indx[n + 3*5]*nrho] * vp[i + 3*nrho];
         }
         dvnl_ray[1*lmaxnray + k1 + indx[n + 3*5]*nray] = P3 * util_simpson(nrho, f, drho);
       }
     }
 
     /* d projectors (l=2) */
     if ((locp != 2) && (lmax > 1)) {
       for (int n = 0; n < n_expansion[2]; ++n) {
 
         f[0] = 0.0;
         for (int i = 1; i < nrho; ++i) {
           const double xx = q * rho[i];
           const double a =
               3.0 * (sn[i]/xx - cs[i]) / xx - sn[i];
           f[i] = a * wp[i + indx[n + 2*5]*nrho] * vp[i + 2*nrho];
         }
         dvnl_ray[0*lmaxnray + k1 + indx[n + 2*5]*nray] = P2 * util_simpson(nrho, f, drho) / q;
 
         f[0] = 0.0;
         for (int i = 1; i < nrho; ++i) {
           const double xx = q * rho[i];
           const double a =
               -9.0 * sn[i] / (xx*xx * q*q)
               +9.0 * cs[i] / (xx    * q*q)
               +4.0 * sn[i] / (q*q)
               -rho[i]*cs[i] / q;
           f[i] = a * wp[i + indx[n + 2*5]*nrho] * vp[i + 2*nrho];
         }
         dvnl_ray[1*lmaxnray + k1 + indx[n + 2*5]*nray] = P2 * util_simpson(nrho, f, drho);
       }
     }
 
     /* p projectors (l=1) */
     if ((locp != 1) && (lmax > 0)) {
       for (int n = 0; n < n_expansion[1]; ++n) {
 
         f[0] = 0.0;
         for (int i = 1; i < nrho; ++i) {
           const double xx = q * rho[i];
           const double a = (sn[i]/xx - cs[i]);
           f[i] = a * wp[i + indx[n + 1*5]*nrho] * vp[i + 1*nrho];
         }
         dvnl_ray[0*lmaxnray + k1 + indx[n + 1*5]*nray] = P1 * util_simpson(nrho, f, drho) / q;
 
         f[0] = 0.0;
         for (int i = 1; i < nrho; ++i) {
           const double xx = q * rho[i];
           const double a =
               -2.0 * sn[i] / (xx * q*q)
               +2.0 * cs[i] / (q*q)
               +rho[i]*sn[i] / q;
           f[i] = a * wp[i + indx[n + 1*5]*nrho] * vp[i + 1*nrho];
         }
         dvnl_ray[1*lmaxnray + k1 + indx[n + 1*5]*nray] = P1 * util_simpson(nrho, f, drho);
       }
     }
 
     /* s projectors (l=0) */
     if (locp != 0) {
       for (int n = 0; n < n_expansion[0]; ++n) {
 
         for (int i = 0; i < nrho; ++i) {
           const double a = -sn[i]/(q*q) + rho[i]*cs[i]/q;
           f[i] = a * wp[i + indx[n + 0*5]*nrho] * vp[i + 0*nrho];
         }
         dvnl_ray[0*lmaxnray + k1 + indx[n + 0*5]*nray] = P0 * util_simpson(nrho, f, drho) / q;
 
         /* B-term for l=0 is typically zero/unused; keep if your formulation requires it.
            Here we leave dvnl_ray[1*...] as already zeroed. */
       }
     }
 
     /* local stress kernel */
     f[0] = 0.0;
     for (int i = 1; i < nrho; ++i) 
     {
        f[i] = rho[i] * vp[i + locp*nrho] * (rho[i]*cs[i] - sn[i]/q);
     }
     dvl_ray[k1] = util_simpson(nrho, f, drho) * forpi / q + zv * forpi / (q*q) * (2.0*cs[nrho-1]/q + rho[nrho-1]*sn[nrho-1]);
 
     /* semicore stress kernel (isotropic) */
     if (semicore) {
       f[0] = 0.0;
       for (int i = 1; i < nrho; ++i) 
       {
          f[i] = rho[i] * std::sqrt(rho_sc_r[i]) * (rho[i]*cs[i] - sn[i]/q);
       }
       rho_sc_k_ray[k1] = util_simpson(nrho, f, drho) * forpi / q;
     }
   }
 
   /* global reductions */
   myparall->Vector_SumAll(0, 2*nray,     rho_sc_k_ray);
   myparall->Vector_SumAll(0, nray,       dvl_ray);
   myparall->Vector_SumAll(0, 2*lmaxnray, dvnl_ray);
 
   /* G == 0 handling */
   dvl_ray[0]        = 0.0;
   rho_sc_k_ray[0]   = 0.0;
 
   for (int l = 0; l <= lmax; ++l)
     for (int n = 0; n < n_expansion[l]; ++n) 
     {
       dvnl_ray[0*lmaxnray + 0 + indx[n + l*5]*nray] = 0.0;
       dvnl_ray[1*lmaxnray + 0 + indx[n + l*5]*nray] = 0.0;
     }
 
   delete[] f;
   delete[] sn;
   delete[] cs;
}



/*********************************************
 *                                           *
 * Psp1d_Hamann::vpp2_generate_stress_spline *
 *                                           *
 *********************************************/
/**
 * @brief Interpolate ray-based pseudopotential stress kernels onto the packed 3D G-grid.
 *
 * Given ray-formatted stress kernels generated by vpp_generate_stress_ray():
 *   - dvl_ray(|G|)               : local stress kernel (scalar)
 *   - dvnl_ray(|G|)              : nonlocal KB stress kernels (two radial pieces per (n,l))
 *                                 dvnl_ray_A(|G|) and dvnl_ray_B(|G|)
 *   - rho_sc_k_ray(|G|)          : semicore stress kernel (scalar, optional)
 *
 * this routine constructs cubic splines in |G| and evaluates them on the packed
 * reciprocal-space grids to produce:
 *   - dvl(G)                     : packed scalar field (npack(0))
 *   - dvnl(G)                    : packed vector fields (x,y,z) per projector component (npack(1))
 *   - rho_sc_k(G)                : packed scalar field (npack(0), optional)
 *
 * The KB nonlocal stress uses the canonical decomposition
 *
 *   σ^{NL}_{αβ}(G) = Σ_l [ A_l(|G|) δ_{αβ} + B_l(|G|) G_α G_β / |G| ]
 *
 * In the plane-wave implementation this appears as derivatives of the KB projector
 * angular factors. For each nonlocal projector component T(u) (u = G/|G|),
 * the vector kernel stored in dvnl is assembled as
 *
 *   ∂/∂G [ D(|G|) T(u) ]  =  DD(|G|) T(u) u  +  D(|G|) (∂T/∂u) (∂u/∂G)
 *
 * which matches the reference Fortran integrate_stress_new implementation.
 */
void Psp1d_Hamann::vpp2_generate_stress_spline(PGrid *mygrid, int nray, double *G_ray,
                                               double *dvl_ray, double *dvnl_ray,
                                               double *rho_sc_k_ray, double *dvl,
                                               double *dvnl, double *rho_sc_k)
{
  /* set up indx(n,l) --> to wp (same mapping used by vpp_generate_ray/spline) */
  int indx[5 * 4];
  int nb = lmax + 1;
  for (auto l = 0; l <= lmax; ++l) {
    indx[l * 5] = l;
    for (auto n1 = 1; n1 < n_expansion[l]; ++n1) {
      indx[n1 + l * 5] = nb;
      ++nb;
    }
  }

  /* allocate spline second-derivative arrays */
  const int lmaxnray = (lmax + 1 + n_extra) * nray;
  double *dvl_splineray      = new double[nray];
  double *dvnl_splineray     = new double[2 * lmaxnray];
  double *rho_sc_splineray   = new double[nray];
  double *tmp_splineray      = new double[nray];

  double *dvnlD  = dvnl_ray;
  double *dvnlDD = dvnl_ray + lmaxnray;

  double *dvnlD_splineray  = dvnl_splineray;
  double *dvnlDD_splineray = dvnl_splineray + lmaxnray;


  /* setup cubic splines */
  const double dG = G_ray[2] - G_ray[1];

  /* local stress kernel spline:
     use same 5-point slope as your vl spline (avoid q=0 singular behavior) */
  double yp1 = (-50.00 * dvl_ray[1] + 96.00 * dvl_ray[2] - 72.00 * dvl_ray[3] +
                32.00 * dvl_ray[4] - 6.00 * dvl_ray[5]) /
               (24.00 * dG);

  util_spline(&(G_ray[1]), &(dvl_ray[1]), nray - 1, yp1, 0.00,
              &(dvl_splineray[1]), tmp_splineray);


   /* nonlocal stress kernels: two splines per (n,l) channel: D and DD
      Layout in dvnl_ray:
        D  at dvnl_ray[ k + 0*lmaxnray + channel*nray ]
        DD at dvnl_ray[ k + 1*lmaxnray + channel*nray ]  == dvnl_ray[k + lmaxnray + channel*nray]
   */
   for (auto l = 0; l <= lmax; ++l)
      if (l != locp)
         for (auto n = 0; n < n_expansion[l]; ++n) 
         {
            const int ch = indx[n + 5*l];
           
            //util_spline(G_ray, &(dvnl_ray[0*lmaxnray    + ch*nray]), nray, 0.0, 0.0,
            //            &(dvnl_splineray[0*lmaxnray    + ch*nray]), tmp_splineray);
           
            //util_spline(G_ray, &(dvnl_ray[1*lmaxnray    + ch*nray]), nray, 0.0, 0.0,
            //            &(dvnl_splineray[1*lmaxnray    + ch*nray]), tmp_splineray);

            util_spline(G_ray, &(dvnlD [ch*nray]),  nray, 0.0, 0.0, &(dvnlD_splineray [ch*nray]),  tmp_splineray);
            util_spline(G_ray, &(dvnlDD[ch*nray]),  nray, 0.0, 0.0, &(dvnlDD_splineray[ch*nray]),  tmp_splineray);
         }



  /* semicore stress kernel spline (scalar) */
  if (semicore) {
    util_spline(G_ray, rho_sc_k_ray, nray, 0.00, 0.00, rho_sc_splineray, tmp_splineray);
  }

  /* output packing */
  const int npack0 = mygrid->npack(0);
  const int npack1 = mygrid->npack(1);

  mygrid->t_pack_nzero(0, 1, dvl);
  mygrid->t_pack_nzero(1, 3 * nprj, dvnl);   // 3 components per projector component
  if (semicore)
    mygrid->t_pack_nzero(0, 1, rho_sc_k);

  /* ---- local + semicore on pack0 grid ---- */
  {
    double *gx = mygrid->Gpackxyz(0, 0);
    double *gy = mygrid->Gpackxyz(0, 1);
    double *gz = mygrid->Gpackxyz(0, 2);

    for (auto k = 0; k < npack0; ++k) {
      const double qx0 = gx[k];
      const double qy0 = gy[k];
      const double qz0 = gz[k];
      const double q   = std::sqrt(qx0*qx0 + qy0*qy0 + qz0*qz0);
      const int nx     = (int)std::floor(q / dG);

      if (q > 1.0e-9) {
        dvl[k] = util_splint(&(G_ray[1]), &(dvl_ray[1]), &(dvl_splineray[1]),
                             nray - 1, nx, q);
        if (semicore) {
          rho_sc_k[k] = util_splint(G_ray, rho_sc_k_ray, rho_sc_splineray, nray, nx, q);
        }
      } else {
        dvl[k] = 0.0;
        if (semicore) rho_sc_k[k] = 0.0;
      }
    }
  }

  /* ---- nonlocal on pack1 grid (vector kernels) ---- */
  {
    double *gx = mygrid->Gpackxyz(1, 0);
    double *gy = mygrid->Gpackxyz(1, 1);
    double *gz = mygrid->Gpackxyz(1, 2);

    for (auto k = 0; k < npack1; ++k) {

      double Gx = gx[k];
      double Gy = gy[k];
      double Gz = gz[k];
      const double q = std::sqrt(Gx*Gx + Gy*Gy + Gz*Gz);
      const int nx   = (int)std::floor(q / dG);

      if (q <= 1.0e-9) {
        for (auto p = 0; p < nprj; ++p) {
          dvnl[k + (0 + 3*p)*npack1] = 0.0;
          dvnl[k + (1 + 3*p)*npack1] = 0.0;
          dvnl[k + (2 + 3*p)*npack1] = 0.0;
        }
        continue;
      }

      /* unit vector u = G/|G| */
      double ux = Gx / q;
      double uy = Gy / q;
      double uz = Gz / q;

      /* du_i / dG_j (matches Fortran) */
      const double duxdGx = 1.0/q - ux*ux/q;
      const double duxdGy = -ux*uy/q;
      const double duxdGz = -ux*uz/q;

      const double duydGx = -uy*ux/q;
      const double duydGy = 1.0/q - uy*uy/q;
      const double duydGz = -uy*uz/q;

      const double duzdGx = -uz*ux/q;
      const double duzdGy = -uz*uy/q;
      const double duzdGz = 1.0/q - uz*uz/q;

      int lcount = nprj;

      auto emit = [&](double D, double DD, double T, double dTdux, double dTduy, double dTduz) {
        const double sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx;
        const double sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy;
        const double sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz;

        --lcount;
        dvnl[k + (0 + 3*lcount)*npack1] = DD*T*ux + D*sumx;
        dvnl[k + (1 + 3*lcount)*npack1] = DD*T*uy + D*sumy;
        dvnl[k + (2 + 3*lcount)*npack1] = DD*T*uz + D*sumz;
      };

      /* f projectors (l=3): 7 components */
      if ((locp != 3) && (lmax > 2))
         for (auto n = 0; n < n_expansion[3]; ++n) 
         {
            const int ch = indx[n + 3*5];
           
            //const double D  = util_splint(G_ray, &(dvnl_ray[0    + ch*nray]),
            //                              &(dvnl_splineray[0    + ch*nray]), nray, nx, q);
            //const double DD = util_splint(G_ray, &(dvnl_ray[nray + ch*nray]),
            //                              &(dvnl_splineray[nray + ch*nray]), nray, nx, q);
           
            //const double D  = util_splint(G_ray, &(dvnl_ray[0*lmaxnray + ch*nray]),
            //                              &(dvnl_splineray[0*lmaxnray + ch*nray]), nray, nx, q);
            //const double DD = util_splint(G_ray, &(dvnl_ray[1*lmaxnray + ch*nray]),
            //                              &(dvnl_splineray[1*lmaxnray + ch*nray]), nray, nx, q);
           
            const double D  = util_splint(G_ray, &(dvnlD [ch*nray]), &(dvnlD_splineray [ch*nray]), nray, nx, q);
            const double DD = util_splint(G_ray, &(dvnlDD[ch*nray]), &(dvnlDD_splineray[ch*nray]), nray, nx, q);
           
            /* Copying exactly the Fortran T and dT/du blocks */
            {
              double T = uy*(3.0*(1.0-uz*uz) - 4.0*uy*uy)/std::sqrt(24.0);
              double dTdux = 0.0;
              double dTduy = (3.0*(1.0-uz*uz) - 12.0*uy*uy)/std::sqrt(24.0);
              double dTduz = -6.0*uy*uz/std::sqrt(24.0);
              emit(D, DD, T, dTdux, dTduy, dTduz);
            }
            {
              double T = ux*uy*uz;
              double dTdux = uy*uz;
              double dTduy = ux*uz;
              double dTduz = ux*uy;
              emit(D, DD, T, dTdux, dTduy, dTduz);
            }
            {
              double T = uy*(5.0*uz*uz - 1.0)/std::sqrt(40.0);
              double dTdux = 0.0;
              double dTduy = (5.0*uz*uz - 1.0)/std::sqrt(40.0);
              double dTduz = 10.0*uy*uz/std::sqrt(40.0);
              emit(D, DD, T, dTdux, dTduy, dTduz);
            }
            {
              double T = uz*(5.0*uz*uz - 3.0)/std::sqrt(60.0);
              double dTdux = 0.0;
              double dTduy = 0.0;
              double dTduz = (15.0*uz*uz - 3.0)/std::sqrt(60.0);
              emit(D, DD, T, dTdux, dTduy, dTduz);
            }
            {
              double T = ux*(5.0*uz*uz - 1.0)/std::sqrt(40.0);
              double dTdux = (5.0*uz*uz - 1.0)/std::sqrt(40.0);
              double dTduy = 0.0;
              double dTduz = 10.0*ux*uz/std::sqrt(40.0);
              emit(D, DD, T, dTdux, dTduy, dTduz);
            }
            {
              double T = uz*(ux*ux - uy*uy)/2.0;
              double dTdux = ux*uz;
              double dTduy = -uy*uz;
              double dTduz = (ux*ux - uy*uy)/2.0;
              emit(D, DD, T, dTdux, dTduy, dTduz);
            }
            {
              double T = ux*(4.0*ux*ux - 3.0*(1.0-uz*uz))/std::sqrt(24.0);
              double dTdux = (12.0*ux*ux - 3.0*(1.0-uz*uz))/std::sqrt(24.0);
              double dTduy = 0.0;
              double dTduz = 6.0*ux*uz/std::sqrt(24.0);
              emit(D, DD, T, dTdux, dTduy, dTduz);
            }
         }

      /* d projectors (l=2): 5 components */
      if ((locp != 2) && (lmax > 1))
         for (auto n = 0; n < n_expansion[2]; ++n) 
         {
            const int ch = indx[n + 2*5];
           
            //const double D  = util_splint(G_ray, &(dvnl_ray[0    + ch*nray]),
            //                              &(dvnl_splineray[0    + ch*nray]), nray, nx, q);
            //const double DD = util_splint(G_ray, &(dvnl_ray[nray + ch*nray]),
            //                              &(dvnl_splineray[nray + ch*nray]), nray, nx, q);
            //const double D  = util_splint(G_ray, &(dvnl_ray[0*lmaxnray + ch*nray]),
            //                              &(dvnl_splineray[0*lmaxnray + ch*nray]), nray, nx, q);
            //const double DD = util_splint(G_ray, &(dvnl_ray[1*lmaxnray + ch*nray]),
            //                              &(dvnl_splineray[1*lmaxnray + ch*nray]), nray, nx, q);
            const double D  = util_splint(G_ray, &(dvnlD [ch*nray]), &(dvnlD_splineray[ch*nray]), nray, nx, q);
            const double DD = util_splint(G_ray, &(dvnlDD[ch*nray]), &(dvnlDD_splineray[ch*nray]), nray, nx, q);
           
            {
              double T = ux*uy;
              emit(D, DD, T, uy, ux, 0.0);
            }
            {
              double T = uy*uz;
              emit(D, DD, T, 0.0, uz, uy);
            }
            {
              double T = (3.0*uz*uz - 1.0)/(2.0*std::sqrt(3.0));
              emit(D, DD, T, 0.0, 0.0, 6.0*uz/(2.0*std::sqrt(3.0)));
            }
            {
              double T = uz*ux;
              emit(D, DD, T, uz, 0.0, ux);
            }
            {
              double T = (ux*ux - uy*uy)/2.0;
              emit(D, DD, T, ux, -uy, 0.0);
            }
         }

      /* p projectors (l=1): 3 components */
      if ((locp != 1) && (lmax > 0))
         for (auto n = 0; n < n_expansion[1]; ++n) 
         {
            const int ch = indx[n + 1*5];
           
            //const double D  = util_splint(G_ray, &(dvnl_ray[0    + ch*nray]),
            //                              &(dvnl_splineray[0    + ch*nray]), nray, nx, q);
            //const double DD = util_splint(G_ray, &(dvnl_ray[nray + ch*nray]),
            //                              &(dvnl_splineray[nray + ch*nray]), nray, nx, q);
            //const double D  = util_splint(G_ray, &(dvnl_ray[0*lmaxnray + ch*nray]), &(dvnl_splineray[0*lmaxnray + ch*nray]), nray, nx, q);
            //const double DD = util_splint(G_ray, &(dvnl_ray[1*lmaxnray + ch*nray]), &(dvnl_splineray[1*lmaxnray + ch*nray]), nray, nx, q);
           
            const double D  = util_splint(G_ray, &(dvnlD [ch*nray]), &(dvnlD_splineray[ch*nray]), nray, nx, q);
            const double DD = util_splint(G_ray, &(dvnlDD[ch*nray]), &(dvnlDD_splineray[ch*nray]), nray, nx, q);
           
            {
              double T = uy;
              emit(D, DD, T, 0.0, 1.0, 0.0);
            }
            {
              double T = uz;
              emit(D, DD, T, 0.0, 0.0, 1.0);
            }
            {
              double T = ux;
              emit(D, DD, T, 1.0, 0.0, 0.0);
            }
         }

      /* s projectors (l=0): 1 component */
      if (locp != 0)
         for (auto n = 0; n < n_expansion[0]; ++n) 
         {
            const int ch = indx[n + 0*5];
           
            //const double D  = util_splint(G_ray, &(dvnl_ray[0    + ch*nray]),
            //                              &(dvnl_splineray[0    + ch*nray]), nray, nx, q);
            //const double D  = util_splint(G_ray, &(dvnl_ray[0*lmaxnray + ch*nray]),
            //                              &(dvnl_splineray[0*lmaxnray + ch*nray]), nray, nx, q);
            //const double DD = util_splint(G_ray, &(dvnl_ray[1*lmaxnray + ch*nray]),
            //                              &(dvnl_splineray[1*lmaxnray + ch*nray]), nray, nx, q);
            /* For s: DD isn’t used in Fortran (pure direction), but keep the same formula:
               T=1, dT/du = 0 -> dvnl = DD*u */
            //const double DD = util_splint(G_ray, &(dvnl_ray[nray + ch*nray]),
            //                              &(dvnl_splineray[nray + ch*nray]), nray, nx, q);
           
            /* For s: DD isn’t used in Fortran (pure direction), but keep the same formula:
               T=1, dT/du = 0 -> dvnl = DD*u */
            const double D  = util_splint(G_ray, &(dvnlD [ch*nray]), &(dvnlD_splineray[ch*nray]), nray, nx, q);
            const double DD = util_splint(G_ray, &(dvnlDD[ch*nray]), &(dvnlDD_splineray[ch*nray]), nray, nx, q);
           
            emit(D, DD, 1.0, 0.0, 0.0, 0.0);
         }
    }
  }

  /* cleanup */
  delete[] tmp_splineray;
  delete[] rho_sc_splineray;
  delete[] dvnl_splineray;
  delete[] dvl_splineray;
}





/*******************************************
 *                                         *
 *   Psp1d_Hamann::cpp_generate_ray        *
 *                                         *
 *******************************************/
/**
 * @brief Generate ray-formatted reciprocal-space pseudopotential kernels.
 *
 * This routine computes the radial (|G|-dependent) Fourier-space kernels
 * for a norm-conserving Hamann pseudopotential directly from the real-space
 * radial functions.
 *
 * For each sampled |G| value in G_ray, the following quantities are evaluated:
 *
 *   (1) Local potential kernel
 *       ------------------------------------
 *       vl_ray(|G|) = 4π / |G| ∫ dr r V_loc(r) sin(|G| r)
 *
 *       with explicit handling of the G → 0 limit and optional Gaussian
 *       compensation (psp version dependent).
 *
 *   (2) Nonlocal KB projector kernels
 *       ------------------------------------
 *       vnl_ray_p(|G|) = ∫ dr r^2 V_l(r) w_p(r) j_l(|G| r)
 *
 *       where p ≡ (n,l) indexes the KB projectors and j_l are spherical
 *       Bessel functions expanded analytically for s, p, d, and f channels.
 *
 *       The output vnl_ray contains only the radial dependence; angular
 *       factors are applied later when mapping onto the full G-grid.
 *
 *   (3) Semicore density kernels (optional)
 *       ------------------------------------
 *       If a semicore charge density is present, two radial kernels are
 *       generated:
 *
 *         rho_sc_k_ray[0:nray-1]         : isotropic component
 *         rho_sc_k_ray[nray:2*nray-1]    : first-order (∇ρ) component
 *
 * Parallelization:
 *   - The |G| grid is distributed across MPI tasks.
 *   - Each task computes a subset of G_ray entries.
 *   - Results are summed across tasks using collective reductions.
 *
 * Special handling:
 *   - G = 0 values are computed analytically to avoid singular behavior.
 *   - Only the s-channel (j₀) contributes at |G| = 0.
 *
 * This routine produces ray-formatted data intended to be interpolated
 * onto full reciprocal-space grids by cpp_generate_spline() or
 * vpp_generate_spline().
 *
 * @param myparall        Parallel execution context.
 * @param nray            Number of radial |G| samples.
 * @param G_ray           Radial grid of |G| values.
 * @param vl_ray          Output local pseudopotential kernel (ray format).
 * @param vnl_ray         Output nonlocal KB projector kernels (ray format).
 * @param rho_sc_k_ray    Output semicore density kernels (ray format).
 */
void Psp1d_Hamann::cpp_generate_ray(Parallel *myparall, int nray, double *G_ray,
                                          double *vl_ray, double *vnl_ray,
                                          double *rho_sc_k_ray) 
{
   /* set up indx(n,l) --> to wp */
   int indx[5 * 4];
   int nb = lmax + 1;
   for (auto l = 0; l <= lmax; ++l) {
     indx[l * 5] = l;
     for (auto n1 = 1; n1 < n_expansion[l]; ++n1) {
       indx[n1 + l * 5] = nb;
       ++nb;
     }
   }
 
   double pi = 4.00 * atan(1.0);
   double twopi = 2.0 * pi;
   double forpi = 4.0 * pi;
 
   double P0 = sqrt(forpi);
   double P1 = sqrt(3.0 * forpi);
   double P2 = sqrt(15.0 * forpi);
   double P3 = sqrt(105.0 * forpi);
 
   double zero = 0.0;
   int izero = 0;
   int ione = 1;
   int nray2 = 2 * nray;
   int lmaxnray = (lmax + 1 + n_extra) * nray;
 
   double q;
   double *cs = new double[nrho];
   double *sn = new double[nrho];
   double *f = new double[nrho];
   double a, xx;
 
   memset(vl_ray, 0, nray * sizeof(double));
   memset(vnl_ray, 0, lmaxnray * sizeof(double));
   memset(rho_sc_k_ray, 0, nray2 * sizeof(double));
 
   for (auto k1 = (1 + myparall->taskid()); k1 < nray; k1 += myparall->np()) {
     q = G_ray[k1];
     for (auto i = 0; i < nrho; ++i) {
       cs[i] = cos(q * rho[i]);
       sn[i] = sin(q * rho[i]);
     }
 
     /* h projectors */
     /* h projectors */
     /* f projectors */
     if ((locp != 3) && (lmax > 2)) {
       for (auto n = 0; n < n_expansion[3]; ++n) {
         f[0] = 0.0;
         for (auto i = 1; i < nrho; ++i) {
           xx = q * rho[i];
           a = sn[i] / xx;
           a = 15.0 * (a - cs[i]) / (xx * xx) - 6 * a + cs[i];
           f[i] = a * wp[i + indx[n + 3 * 5] * nrho] * vp[i + 3 * nrho];
         }
         vnl_ray[k1 + indx[n + 3 * 5] * nray] =
             P3 * util_simpson(nrho, f, drho) / q;
       }
     }
 
     /* d projectors */
     if ((locp != 2) && (lmax > 1)) {
       for (auto n = 0; n < n_expansion[2]; ++n) {
         f[0] = 0.0;
         for (auto i = 1; i < nrho; ++i) {
           a = 3.0 * (sn[i] / (q * rho[i]) - cs[i]) / (q * rho[i]) - sn[i];
           f[i] = a * wp[i + indx[n + 2 * 5] * nrho] * vp[i + 2 * nrho];
         }
         vnl_ray[k1 + indx[n + 2 * 5] * nray] =
             P2 * util_simpson(nrho, f, drho) / q;
       }
     }
 
     /* p projectors */
     if ((locp != 1) && (lmax > 0)) {
       for (auto n = 0; n < n_expansion[1]; ++n) {
         f[0] = 0.0;
         for (auto i = 1; i < nrho; ++i) {
           a = (sn[i] / (q * rho[i]) - cs[i]);
           f[i] = a * wp[i + indx[n + 1 * 5] * nrho] * vp[i + 1 * nrho];
         }
         vnl_ray[k1 + indx[n + 1 * 5] * nray] =
             P1 * util_simpson(nrho, f, drho) / q;
       }
     }
 
     /* s projectors */
     if (locp != 0) {
       for (auto n = 0; n < n_expansion[0]; ++n) {
         for (auto i = 0; i < nrho; ++i)
           f[i] = sn[i] * wp[i + indx[n + 0 * 5] * nrho] * vp[i + 0 * nrho];
         vnl_ray[k1 + indx[n + 0 * 5] * nray] =
             P0 * util_simpson(nrho, f, drho) / q;
       }
     }
 
     /* local */
     if (version == 3) {
       for (auto i = 0; i < nrho; ++i)
         f[i] = rho[i] * vp[i + locp * nrho] * sn[i];
       vl_ray[k1] = util_simpson(nrho, f, drho) * forpi / q -
                    zv * forpi * cs[nrho - 1] / (q * q);
     } else if (version == 4) {
       for (auto i = 0; i < nrho; ++i)
         f[i] = (rho[i] * vp[i + locp * nrho] + zv * std::erf(rho[i] / rlocal)) *
                sn[i];
       vl_ray[k1] = util_simpson(nrho, f, drho) * forpi / q;
     }
 
     /* semicore density */
     if (semicore) {
       for (auto i = 0; i < nrho; ++i)
         f[i] = rho[i] * sqrt(rho_sc_r[i]) * sn[i];
       rho_sc_k_ray[k1] = util_simpson(nrho, f, drho) * forpi / q;
 
       for (auto i = 0; i < nrho; ++i)
         f[i] = (sn[i] / (q * rho[i]) - cs[i]) * rho_sc_r[i + nrho] * rho[i];
       rho_sc_k_ray[k1 + nray] = util_simpson(nrho, f, drho) * forpi / q;
     }
   }
   myparall->Vector_SumAll(0, 2 * nray, rho_sc_k_ray);
   myparall->Vector_SumAll(0, nray, vl_ray);
   myparall->Vector_SumAll(0, lmaxnray, vnl_ray);
 
   /* G==0 local */
   if (version == 3) {
     for (auto i = 0; i < nrho; ++i) {
       f[i] = vp[i + locp * nrho] * rho[i] * rho[i];
     }
     vl_ray[0] = forpi * util_simpson(nrho, f, drho) +
                 twopi * zv * rho[nrho - 1] * rho[nrho - 1];
   } else if (version == 4) {
     for (auto i = 0; i < nrho; ++i)
       f[i] = (vp[i + locp * nrho] * rho[i] + zv * std::erf(rho[i] / rlocal)) *
              rho[i];
     vl_ray[0] = forpi * util_simpson(nrho, f, drho);
   }
 
   /* G==0 semicore */
   if (semicore) {
     for (auto i = 0; i < nrho; ++i)
       f[i] = sqrt(rho_sc_r[i]) * rho[i] * rho[i];
     rho_sc_k_ray[0] = forpi * util_simpson(nrho, f, drho);
     rho_sc_k_ray[0 + nray] = 0.0;
   }
 
   /* G==0 vnl */
   for (auto l = 0; l <= lmax; ++l)
     for (auto n = 0; n < n_expansion[l]; ++n)
       vnl_ray[0 + indx[n + l * 5] * nray] = 0.0;
 
   /* only j0 is non-zero at zero */
   if (locp != 0)
     for (auto n = 0; n < n_expansion[0]; ++n) {
       for (auto i = 0; i < nrho; ++i)
         f[i] = rho[i] * wp[i + indx[n + 0 * 5] * nrho] * vp[i + 0 * nrho];
       vnl_ray[0 + indx[n + 0 * 5] * nray] = P0 * util_simpson(nrho, f, drho);
     }
 
   delete[] f;
   delete[] sn;
   delete[] cs;
}



/*************************************************
 *                                               *
 *   Psp1d_Hamann::cpp_generate_local_spline     *
 *                                               *
 *************************************************/
/**
 * @brief Interpolate ray-formatted local pseudopotential kernels onto the packed 3D G-grid.
 *
 * This routine takes radial (ray-based) representations of the local
 * pseudopotential and (optionally) the semicore charge-density kernel
 * generated in reciprocal space and interpolates them onto the packed
 * plane-wave G-grid used by PWDFT.
 *
 * Inputs:
 *   - G_ray(|G|)        : radial G-grid
 *   - vl_ray(|G|)       : local pseudopotential kernel V_loc(|G|)
 *   - rho_sc_k_ray(|G|) : semicore kernel in reciprocal space
 *                         [isotropic component + radial derivative]
 *
 * Outputs:
 *   - vl(G)             : local pseudopotential evaluated on packed G-grid
 *   - rho_sc_k(G)       : semicore contributions packed as:
 *                           rho_sc_k[0]           : scalar (isotropic) term
 *                           rho_sc_k[1..3]        : vector components
 *
 * Method:
 *   - Constructs cubic splines of the radial kernels in |G|.
 *   - Evaluates splines at |G| for each packed reciprocal vector.
 *   - Handles the G → 0 limit explicitly to avoid singular behavior.
 *
 * Notes:
 *   - This routine handles only the *local* pseudopotential contribution.
 *   - No angular momentum decomposition is required.
 *   - Semicore contributions are treated consistently with legacy NWPW
 *     implementations, producing both scalar and vector components.
 *
 * Reference:
 *   Corresponds to the local pseudopotential interpolation logic used in
 *   legacy NWPW Fortran routines (non-stress path).
 */

void Psp1d_Hamann::cpp_generate_local_spline(CGrid *mygrid, int nray, double *G_ray,
                                       double *vl_ray, double *rho_sc_k_ray, 
                                       double *vl, double *rho_sc_k) 
{
   /* set up indx(n,l) --> to wp */
   int indx[5 * 4];
   int nb = lmax + 1;
   for (auto l = 0; l <= lmax; ++l) {
     indx[l * 5] = l;
     for (auto n1 = 1; n1 < n_expansion[l]; ++n1) {
       indx[n1 + l * 5] = nb;
       ++nb;
     }
   }
 
   double pi = 4.00 * atan(1.0);
 
   /* allocate spline grids */
   double *vl_splineray = new double[nray];
   double *rho_sc_k_splineray = new double[2 * nray];
   double *tmp_splineray = new double[nray];
 
   /* setup cubic bsplines */
   double dG = G_ray[2] - G_ray[1];
 
   /* five point formula */
   double yp1 = (-50.00 * vl_ray[1] + 96.00 * vl_ray[2] - 72.00 * vl_ray[3] +
                 32.00 * vl_ray[4] - 6.00 * vl_ray[5]) /
                (24.00 * dG);
   util_spline(&(G_ray[1]), &(vl_ray[1]), nray - 1, yp1, 0.00,
               &(vl_splineray[1]), tmp_splineray);
 
   if (semicore) {
     util_spline(G_ray, rho_sc_k_ray, nray, 0.00, 0.00, rho_sc_k_splineray,
                 tmp_splineray);
     util_spline(G_ray, &(rho_sc_k_ray[nray]), nray, 0.00, 0.00,
                 &(rho_sc_k_splineray[nray]), tmp_splineray);
   }
 
   double q, qx, qy, qz, xx;
   double *gx, *gy, *gz;
   int npack0 = mygrid->npack(0);
   int nx, lcount;
   mygrid->t_pack_nzero(0, 1, vl);
   if (semicore)
     mygrid->t_pack_nzero(0, 4, rho_sc_k);
 
   /* generate vl and rho_sc_k */
   gx = mygrid->Gpackxyz(0, 0);
   gy = mygrid->Gpackxyz(0, 1);
   gz = mygrid->Gpackxyz(0, 2);
   for (auto k = 0; k < npack0; ++k) 
   {
      qx = gx[k];
      qy = gy[k];
      qz = gz[k];
      q = sqrt(qx * qx + qy * qy + qz * qz);
      nx = (int)floor(q / dG);
     
      if (q > 1.0e-9) {
        qx /= q;
        qy /= q;
        qz /= q;
        vl[k] = util_splint(&(G_ray[1]), &(vl_ray[1]), &(vl_splineray[1]),
                            nray - 1, nx, q);
        if (semicore) {
          rho_sc_k[k] =
              util_splint(G_ray, rho_sc_k_ray, rho_sc_k_splineray, nray, nx, q);
          xx = util_splint(G_ray, &(rho_sc_k_ray[nray]),
                           &(rho_sc_k_splineray[nray]), nray, nx, q);
          rho_sc_k[k + npack0] = xx * qx;
          rho_sc_k[k + 2 * npack0] = xx * qy;
          rho_sc_k[k + 3 * npack0] = xx * qz;
        }
      } else {
     
        vl[k] = vl_ray[0];
        if (semicore) {
          rho_sc_k[k] = rho_sc_k_ray[0];
          rho_sc_k[k + npack0] = 0.0;
          rho_sc_k[k + 2 * npack0] = 0.0;
          rho_sc_k[k + 3 * npack0] = 0.0;
        }
      }
   }
 
 
 
   /*  deallocate spineray formatted grids */
   delete[] tmp_splineray;
   delete[] rho_sc_k_splineray;
   delete[] vl_splineray;
}



/**************************************************
 *                                                *
 *   Psp1d_Hamann::cpp_generate_nonlocal_spline   *
 *                                                *
 **************************************************/
/**
 * @brief Interpolate nonlocal KB pseudopotential projectors onto the packed 3D G-grid.
 *
 * This routine constructs the nonlocal Kleinman–Bylander (KB) projector
 * values on the packed reciprocal-space plane-wave grid.  Radial projector
 * functions provided on a 1D |G|-ray are spline-interpolated and combined with
 * analytic angular factors to form the full nonlocal projectors in G-space.
 *
 * Input (ray representation):
 *   - vnl_ray(|G|) :
 *       Radial nonlocal projector functions for each projector channel
 *       p = (n,l), defined on a 1D reciprocal-space ray.
 *
 * Input (reciprocal grid):
 *   - G :
 *       Packed reciprocal-space vectors on the FFT grid.
 *   - kvec :
 *       Brillouin-zone k-point associated with the current band.
 *
 * Output:
 *   - vnl(G) :
 *       Packed nonlocal KB projectors evaluated on the 3D G-grid, with one
 *       scalar value per projector component and plane-wave index.
 *
 * Method:
 *   - Construct cubic splines in |G| for each radial projector v_p(|G|).
 *   - For each reciprocal vector G, form G + k and its magnitude |G + k|.
 *   - Evaluate the radial projector at |G + k| via spline interpolation.
 *   - Multiply by the corresponding real spherical-harmonic angular factors
 *     to assemble the full KB projector.
 *
 * Notes:
 *   - Brillouin-zone dependence enters explicitly through (G + k).
 *   - Projector ordering and normalization match the legacy NWPW implementation.
 *   - For G + k = 0, only s-channel (l = 0) projectors are nonzero; higher-l
 *     components are explicitly set to zero.
 */
void Psp1d_Hamann::cpp_generate_nonlocal_spline(CGrid *mygrid, int nbq, double *kvec, int nray, double *G_ray, double *vnl_ray, double *vnl)
{
   /* set up indx(n,l) --> to wp */
   int indx[5*4];
   int nb = lmax + 1;
   for (auto l=0; l<=lmax; ++l) 
   {
      indx[l*5] = l;
      for (auto n1=1; n1<n_expansion[l]; ++n1) 
      {
         indx[n1 + l*5] = nb;
         ++nb;
      }
   }

   double pi = 4.00 * atan(1.0);
   double q, qx, qy, qz, xx;
   double *gx, *gy, *gz;
   int nx, lcount;
   int npack1_max = mygrid->npack1_max();
   int npack1 = mygrid->npack(1+nbq);

   /* allocate spline grids */
   double *vnl_splineray = new double[(lmax + 1 + n_extra)*nray];
   double *tmp_splineray = new double[nray];

   /* setup cubic bsplines */
   double dG = G_ray[2] - G_ray[1];

   for (auto l = 0; l<=lmax; ++l)
      if (l != locp)
         for (auto n = 0; n < n_expansion[l]; ++n)
            util_spline(G_ray, &(vnl_ray[indx[n+5*l]*nray]), nray, 0.00, 0.00,
                               &(vnl_splineray[indx[n+5*l]*nray]), tmp_splineray);

   /* generate vnl */
   //mygrid->t_pack_nzero(nbq,nprj,vnl);
   mygrid->t_pack_max_nzero(nprj,vnl);
   gx = mygrid->Gpackxyz(1+nbq,0);
   gy = mygrid->Gpackxyz(1+nbq,1);
   gz = mygrid->Gpackxyz(1+nbq,2);

   for (auto k=0; k<npack1; ++k) 
   {
      qx = gx[k] + kvec[0];
      qy = gy[k] + kvec[1];
      qz = gz[k] + kvec[2];
      q = sqrt(qx*qx + qy*qy + qz*qz);
      nx = (int)floor(q / dG);
     
      if (q > 1.0e-9) 
      {
         qx /= q;
         qy /= q;
         qz /= q;
         lcount = nprj;
        
         /* f projectors */
        
         if ((locp != 3) && (lmax > 2))
            for (auto n=0; n<n_expansion[3]; ++n) 
            {
               xx = util_splint(G_ray, &(vnl_ray[indx[n + 3*5]*nray]),
                                       &(vnl_splineray[indx[n + 3*5]*nray]), nray, nx, q);
               --lcount;
               vnl[k + lcount*npack1_max] = xx*qy*(3.00*(1.00 - qz*qz) - 4.00*qy*qy) / sqrt(24.00);
               --lcount;
               vnl[k + lcount*npack1_max] = xx*qx*qy*qz;
               --lcount;
               vnl[k + lcount*npack1_max] = xx*qy*(5.00*qz*qz - 1.00) / sqrt(40.00);
               --lcount;
               vnl[k + lcount*npack1_max] = xx*qz*(5.00*qz*qz - 3.00) / sqrt(60.00);
               --lcount;
               vnl[k + lcount*npack1_max] = xx*qx*(5.00*qz*qz - 1.00) / sqrt(40.00);
               --lcount;
               vnl[k + lcount*npack1_max] = xx*qz*(qx*qx - qy*qy) / 2.00;
               --lcount;
               vnl[k + lcount*npack1_max] = xx*qx*(4.00*qx*qx - 3.00*(1.00 - qz*qz)) / sqrt(24.00);
            }
        
         /* d projectors */
         if ((locp != 2) && (lmax > 1))
            for (auto n=0; n<n_expansion[2]; ++n) 
            {
               xx = util_splint(G_ray, &(vnl_ray[indx[n + 2*5]*nray]),
                                       &(vnl_splineray[indx[n + 2*5]*nray]), nray, nx, q);
               --lcount;
               vnl[k + lcount*npack1_max] = xx*qx*qy;
               --lcount;
               vnl[k + lcount*npack1_max] = xx*qy*qz;
               --lcount;
               vnl[k + lcount*npack1_max] = xx*(3.00*qz*qz - 1.00) / (2.00*sqrt(3.00));
               --lcount;
               vnl[k + lcount*npack1_max] = xx*qz*qx;
               --lcount;
               vnl[k + lcount*npack1_max] = xx*(qx*qx - qy*qy) / (2.00);
            }
        
         /* p projectors */
         if ((locp != 1) && (lmax > 0))
           for (auto n = 0; n < n_expansion[1]; ++n) 
           {
              xx = util_splint(G_ray, &(vnl_ray[indx[n + 1*5]*nray]),
                                      &(vnl_splineray[indx[n + 1*5]*nray]), nray, nx, q);
              --lcount;
              vnl[k + lcount*npack1_max] = xx*qx;
              --lcount;
              vnl[k + lcount*npack1_max] = xx*qy;
              --lcount;
              vnl[k + lcount*npack1_max] = xx*qz;
           }
        
         /* s projectors */
         if (locp != 0)
            for (auto n=0; n<n_expansion[0]; ++n) 
            {
               xx = util_splint(G_ray, &(vnl_ray[indx[n + 0*5]*nray]),
                                       &(vnl_splineray[indx[n + 0*5]*nray]), nray, nx, q);
               --lcount;
               vnl[k + lcount * npack1_max] = xx;
            }
      } 
      else 
      {
         for (auto l=0; l<nprj; ++l)
            vnl[k + l*npack1_max] = 0.0;
        
        
         /* only j0 is non-zero at zero */
         if (locp != 0)
            for (auto n = 0; n < n_expansion[0]; ++n)
               vnl[k + indx[n + 0*5]*npack1_max] = vnl_ray[0 + indx[n + 0*5]*nray];
      }
   }

   /*  deallocate spineray formatted grids */
   delete[] tmp_splineray;
   delete[] vnl_splineray;
}



/*******************************************
 *                                         *
 * Psp1d_Hamann::cpp2_generate_stress_ray  *
 *                                         *
 *******************************************/
/**
 * @brief Generate ray-formatted pseudopotential stress kernels in G-space.
 *
 * k-point (BZ) dependence enters later via (G+k) contraction and BZ summation
 *
 * This routine computes the reciprocal-space radial kernels needed to assemble
 * the pseudopotential contribution to the stress tensor.  The kernels are
 * produced as functions of the ray magnitude |G| (a 1D radial grid in
 * reciprocal space) and are later interpolated onto the full 3D plane-wave
 * grid for contraction into the Cartesian stress components.
 *
 * The generated contributions are:
 *
 *  (1) Nonlocal (Kleinman–Bylander) stress kernels
 *      ------------------------------------------------------------
 *      The KB nonlocal stress contribution can be written as:
 *
 *        σ^{NL}_{αβ} = Σ_l Σ_G [ A_l(G) δ_{αβ} + B_l(G) (G_α G_β)/|G| ] .
 *
 *      This routine generates the radial functions A_l(G) and B_l(G) for each
 *      projector channel (n,l) and stores them in dvnl_ray as:
 *
 *        dvnl_ray[ 0*nray + p*nray + k ] -> A_p(|G_k|)   (isotropic / metric term)
 *        dvnl_ray[ 1*nray + p*nray + k ] -> B_p(|G_k|)   (anisotropic term)
 *
 *      where p is the flattened projector-channel index produced by indx(n,l),
 *      and k indexes the ray sample (k=0..nray-1).
 *
 *  (2) Local pseudopotential stress kernel
 *      ------------------------------------------------------------
 *      The local part contributes an isotropic stress kernel:
 *
 *        σ^{L}_{αβ} = Σ_G A_loc(G) δ_{αβ},
 *
 *      where A_loc(G) is generated here as dvl_ray(G).
 *
 *  (3) Semicore density stress kernel (optional)
 *      ------------------------------------------------------------
 *      If a semicore charge density is present, an additional isotropic kernel
 *      is required for consistent stress evaluation:
 *
 *        σ^{SC}_{αβ} = Σ_G A_sc(G) δ_{αβ},
 *
 *      where A_sc(G) is generated here as rho_sc_k_ray(G).
 *
 * Implementation notes:
 *  - Kernels are computed by analytic Fourier–Bessel transforms of radial
 *    quantities (Hamann/Kleinman–Bylander form) and evaluated on the ray grid.
 *  - Parallelization is over ray index k; results are summed with Vector_SumAll.
 *  - G=0 components are explicitly set to zero where required to avoid
 *    singular terms and enforce translational invariance.
 */
void Psp1d_Hamann::cpp2_generate_stress_ray(Parallel *myparall,
                                            int nray, double *G_ray,
                                            double *dvl_ray, double *dvnl_ray,
                                            double *rho_sc_k_ray)
{
   /* set up indx(n,l) --> to wp */
   int indx[5 * 4];
   int nb = lmax + 1;
   for (auto l = 0; l <= lmax; ++l) {
     indx[l * 5] = l;
     for (auto n1 = 1; n1 < n_expansion[l]; ++n1) {
       indx[n1 + l * 5] = nb;
       ++nb;
     }
   }

   double pi = 4.00 * atan(1.0);
   double twopi = 2.0 * pi;
   double forpi = 4.0 * pi;

   double P0 = sqrt(forpi);
   double P1 = sqrt(3.0 * forpi);
   double P2 = sqrt(15.0 * forpi);
   double P3 = sqrt(105.0 * forpi);

   double zero = 0.0;
   int izero = 0;
   int ione = 1;
   int nray2 = 2 * nray;
   int lmaxnray = (lmax + 1 + n_extra) * nray;

   double q;
   double *cs = new double[nrho];
   double *sn = new double[nrho];
   double *f = new double[nrho];
   double a, xx;

   std::memset(dvl_ray, 0, nray * sizeof(double));
   std::memset(dvnl_ray, 0, 2 * lmaxnray * sizeof(double));
   std::memset(rho_sc_k_ray, 0, nray2 * sizeof(double));




   /* parallel loop over ray points, skipping k=0 (G=0 handled later) */
   for (int k1 = 1 + myparall->taskid(); k1 < nray; k1 += myparall->np()) {

     const double q = G_ray[k1];

     for (int i = 0; i < nrho; ++i) {
       cs[i] = std::cos(q * rho[i]);
       sn[i] = std::sin(q * rho[i]);
     }

     /* f projectors (l=3) */
     if ((locp != 3) && (lmax > 2)) {
       for (int n = 0; n < n_expansion[3]; ++n) {

         f[0] = 0.0;
         for (int i = 1; i < nrho; ++i) {
           const double xx = q * rho[i];
           double a = sn[i] / xx;
           a = 15.0 * (a - cs[i]) / (xx * xx) - 6.0 * a + cs[i];
           f[i] = a * wp[i + indx[n + 3*5] * nrho] * vp[i + 3*nrho];
         }
         dvnl_ray[0*lmaxnray + k1 + indx[n + 3*5]*nray] = P3 * util_simpson(nrho, f, drho) / q;

         f[0] = 0.0;
         for (int i = 1; i < nrho; ++i) {
           const double xx = q * rho[i];
           double a =
               -60.0 * sn[i] / (xx*xx*xx * q*q)
               +60.0  * cs[i] / (xx*xx     * q*q)
               +27.0  * sn[i] / (xx        * q*q)
               -7.0   * cs[i] / (q*q)
               -rho[i]* sn[i] / q;
           f[i] = a * wp[i + indx[n + 3*5]*nrho] * vp[i + 3*nrho];
         }
         dvnl_ray[1*lmaxnray + k1 + indx[n + 3*5]*nray] = P3 * util_simpson(nrho, f, drho);
       }
     }


     /* d projectors (l=2) */
     if ((locp != 2) && (lmax > 1)) {
       for (int n = 0; n < n_expansion[2]; ++n) {

         f[0] = 0.0;
         for (int i = 1; i < nrho; ++i) {
           const double xx = q * rho[i];
           const double a =
               3.0 * (sn[i]/xx - cs[i]) / xx - sn[i];
           f[i] = a * wp[i + indx[n + 2*5]*nrho] * vp[i + 2*nrho];
         }
         dvnl_ray[0*lmaxnray + k1 + indx[n + 2*5]*nray] = P2 * util_simpson(nrho, f, drho) / q;

         f[0] = 0.0;
         for (int i = 1; i < nrho; ++i) {
           const double xx = q * rho[i];
           const double a =
               -9.0 * sn[i] / (xx*xx * q*q)
               +9.0 * cs[i] / (xx    * q*q)
               +4.0 * sn[i] / (q*q)
               -rho[i]*cs[i] / q;
           f[i] = a * wp[i + indx[n + 2*5]*nrho] * vp[i + 2*nrho];
         }
         dvnl_ray[1*lmaxnray + k1 + indx[n + 2*5]*nray] = P2 * util_simpson(nrho, f, drho);
       }
     }

     /* p projectors (l=1) */
     if ((locp != 1) && (lmax > 0)) {
       for (int n = 0; n < n_expansion[1]; ++n) {

         f[0] = 0.0;
         for (int i = 1; i < nrho; ++i) {
           const double xx = q * rho[i];
           const double a = (sn[i]/xx - cs[i]);
           f[i] = a * wp[i + indx[n + 1*5]*nrho] * vp[i + 1*nrho];
         }
         dvnl_ray[0*lmaxnray + k1 + indx[n + 1*5]*nray] = P1 * util_simpson(nrho, f, drho) / q;

         f[0] = 0.0;
         for (int i = 1; i < nrho; ++i) {
           const double xx = q * rho[i];
           const double a =
               -2.0 * sn[i] / (xx * q*q)
               +2.0 * cs[i] / (q*q)
               +rho[i]*sn[i] / q;
           f[i] = a * wp[i + indx[n + 1*5]*nrho] * vp[i + 1*nrho];
         }
         dvnl_ray[1*lmaxnray + k1 + indx[n + 1*5]*nray] = P1 * util_simpson(nrho, f, drho);
       }
     }



     /* s projectors (l=0) */
     if (locp != 0) {
       for (int n = 0; n < n_expansion[0]; ++n) {

         for (int i = 0; i < nrho; ++i) {
           const double a = -sn[i]/(q*q) + rho[i]*cs[i]/q;
           f[i] = a * wp[i + indx[n + 0*5]*nrho] * vp[i + 0*nrho];
         }
         dvnl_ray[0*lmaxnray + k1 + indx[n + 0*5]*nray] = P0 * util_simpson(nrho, f, drho) / q;

         /* B-term for l=0 is typically zero/unused; keep if your formulation requires it.
            Here we leave dvnl_ray[1*...] as already zeroed. */
       }
     }

     /* local stress kernel */
     f[0] = 0.0;
     for (int i = 1; i < nrho; ++i)
     {
        f[i] = rho[i] * vp[i + locp*nrho] * (rho[i]*cs[i] - sn[i]/q);
     }
     dvl_ray[k1] = util_simpson(nrho, f, drho) * forpi / q + zv * forpi / (q*q) * (2.0*cs[nrho-1]/q + rho[nrho-1]*sn[nrho-1]);

     /* semicore stress kernel (isotropic) */
     // semicore_small_cell mirrors control_psp_semicore_small() in legacy NWPW
     if (semicore) 
     {
        f[0] = 0.0;
        for (int i=1; i<nrho; ++i)
        {
           const double rsc = rho_sc_r[i + 0*nrho];   // matches Fortran (i,1)
           const double fac = semicore_small_cell ? rsc : std::sqrt(rsc);
 
           //f[i] = rho[i] * std::sqrt(rho_sc_r[i]) * (rho[i]*cs[i] - sn[i]/q);
           //f[i] = rho[i] * std::sqrt(rsc) * (rho[i]*cs[i] - sn[i]/q);
           f[i] = rho[i] * fac * (rho[i]*cs[i] - sn[i]/q);
        }
        rho_sc_k_ray[k1] = util_simpson(nrho, f, drho) * forpi / q;
      }
   }

   /* global reductions */
   myparall->Vector_SumAll(0, 2*nray,     rho_sc_k_ray);
   myparall->Vector_SumAll(0, nray,       dvl_ray);
   myparall->Vector_SumAll(0, 2*lmaxnray, dvnl_ray);

   /* G == 0 handling */
   dvl_ray[0]        = 0.0;
   rho_sc_k_ray[0]   = 0.0;

   for (int l = 0; l <= lmax; ++l)
     for (int n = 0; n < n_expansion[l]; ++n)
     {
       dvnl_ray[0*lmaxnray + 0 + indx[n + l*5]*nray] = 0.0;
       dvnl_ray[1*lmaxnray + 0 + indx[n + l*5]*nray] = 0.0;
     }

   delete[] f;
   delete[] sn;
   delete[] cs;
}


/******************************************************
 *                                                    *
 *   Psp1d_Hamann::cpp2_generate_stress_local_spline  *
 *                                                    *
 ******************************************************/
/**
 * @brief Interpolate local (and semicore) stress kernels from ray form
 *        onto the packed 3D reciprocal-space grid.
 *
 * This routine is the C++ analogue of the legacy Fortran routine
 * integrate_kbpp_band_stress_local_new.  It preserves Fortran spline
 * semantics (1-based indexing) and data layout exactly.
 *
 * Semicore layout (if enabled):
 *   rho_sc_k[0*npack + k] : scalar rho_sc(|G|)
 *   rho_sc_k[1*npack + k] : d/dGx contribution
 *   rho_sc_k[2*npack + k] : d/dGy contribution
 *   rho_sc_k[3*npack + k] : d/dGz contribution
 *
 * Notes:
 *  - G=0 contribution is explicitly zeroed.
 *  - util_spline / util_splint are assumed Fortran-style (Numerical Recipes).
 */
void Psp1d_Hamann::cpp2_generate_stress_local_spline(CGrid *mygrid,
                                                     int     nray,
                                                     double *G_ray, double *dvl_ray, double *rho_sc_k_ray,
                                                     double *dvl, double *rho_sc_k)
{
   /* ---------------- spline workspace ---------------- */

   double *dvl_splineray    = new double[nray];
   double *rho_sc_splineray = semicore ? new double[2*nray] : nullptr;
   double *tmp_splineray    = new double[nray];

   const double dG = G_ray[2] - G_ray[1];

   /* ---------------- local stress spline ----------------
    * Matches:
    *   call nwpw_spline(G_ray,dvl_ray(1,1),nray,0,0,dvl_ray(1,2))
    *
    * with 1-based arrays emulated by pointer shifts.
    */

   double yp1 = (-50.0 * dvl_ray[1] +
                  96.0 * dvl_ray[2] -
                  72.0 * dvl_ray[3] +
                  32.0 * dvl_ray[4] -
                   6.0 * dvl_ray[5]) / (24.0 * dG);

   util_spline(&(G_ray[1]), &(dvl_ray[1]), nray-1, yp1, 0.0, &(dvl_splineray[1]), tmp_splineray);

   /* ---------------- semicore spline ---------------- */

   if (semicore) 
   {
      /* scalar component */
      util_spline(G_ray, rho_sc_k_ray, nray, 0.0, 0.0, rho_sc_splineray, tmp_splineray);

      /* derivative magnitude component */
      util_spline(G_ray, &(rho_sc_k_ray[nray]), nray, 0.0, 0.0, &(rho_sc_splineray[nray]), tmp_splineray);
   }

   /* ---------------- output packing ---------------- */

   const int npack0 = mygrid->npack(0);

   mygrid->t_pack_nzero(0, 1, dvl);
   if (semicore)
       mygrid->t_pack_nzero(0, 4, rho_sc_k);

   /* ---------------- evaluate on grid ---------------- */

   const double *gx = mygrid->Gpackxyz(0, 0);
   const double *gy = mygrid->Gpackxyz(0, 1);
   const double *gz = mygrid->Gpackxyz(0, 2);

   for (int k=0; k<npack0; ++k) 
   {
      const double qx = gx[k];
      const double qy = gy[k];
      const double qz = gz[k];
      const double q  = std::sqrt(qx*qx + qy*qy + qz*qz);

      if (q <= 1.0e-9) {
          dvl[k] = 0.0;
          if (semicore) {
              rho_sc_k[0*npack0 + k] = 0.0;
              rho_sc_k[1*npack0 + k] = 0.0;
              rho_sc_k[2*npack0 + k] = 0.0;
              rho_sc_k[3*npack0 + k] = 0.0;
          }
          continue;
      }

      /* Fortran-style bracket index:
       *   nx = (q/dG) + 1   in Fortran
       */
      const int nloc = nray - 1;
      int nx = static_cast<int>(std::floor(q / dG)) + 1;  // <-- FIX
      if (nx < 1) nx = 1;
      if (nx > nloc - 1) nx = nloc - 1;                   // = nray-2

      int nx_sc  = static_cast<int>(std::floor(q / dG)) + 1;
      if (nx_sc < 1) nx_sc = 1;
      if (nx_sc > nray - 1) nx_sc = nray - 1;             // semicore uses n=nray


      /* local stress */
      dvl[k] = util_splint(&(G_ray[1]), &(dvl_ray[1]), &(dvl_splineray[1]), nray-1, nx, q);

      /* semicore */
      if (semicore) 
      {
         const double rho  = util_splint(G_ray, rho_sc_k_ray, rho_sc_splineray, nray, nx_sc, q);
         const double drho = util_splint(G_ray, &(rho_sc_k_ray[nray]), &(rho_sc_splineray[nray]), nray, nx_sc, q);

         rho_sc_k[0*npack0 + k] = rho;
         rho_sc_k[1*npack0 + k] = drho * qx;
         rho_sc_k[2*npack0 + k] = drho * qy;
         rho_sc_k[3*npack0 + k] = drho * qz;
      }
   }

   /* ---------------- cleanup ---------------- */
   delete[] tmp_splineray;
   delete[] dvl_splineray;
   if (semicore) delete[] rho_sc_splineray;
}




/******************************************************
 *                                                    *
 * Psp1d_Hamann::cpp2_generate_stress_nonlocal_spline *
 *                                                    *
 ******************************************************/
/**
 * @brief Assemble nonlocal KB pseudopotential stress kernels on the packed 3D G-grid.
 *
 * This routine is the C++ / cpp2 analogue of the legacy Fortran routine
 * integrate_kbpp_band_stress_nonlocal_new.  It interpolates ray-based nonlocal
 * stress kernels onto the packed reciprocal-space plane-wave grid and assembles
 * the vector-valued kernels required for stress contraction.
 *
 * Input (ray representation):
 *   - dvnl_ray(|G|) :
 *       Radial nonlocal stress kernels generated by cpp2_generate_stress_ray(),
 *       provided for each projector channel p = (n,l) as two components:
 *
 *         A_p(|G|)  : radial kernel D(|G|) associated with angular projector derivatives
 *         B_p(|G|)  : radial kernel D'(|G|) multiplying the directional term proportional
 *                     to u = (G+k)/|G+k|
 *
 *       corresponding to the Kleinman–Bylander decomposition
 *
 *         σ^{NL}_{αβ}(G) = Σ_p [ A_p(|G|) δ_{αβ}
 *                              + B_p(|G|) G_α G_β / |G| ] .
 *
 * Input (reciprocal grid):
 *   - G :
 *       Packed reciprocal-space vectors on the FFT grid.
 *   - kvec :
 *       Brillouin-zone k-point associated with the current band.
 *
 * Output:
 *   - dvnl(G) :
 *       Packed vector-valued nonlocal stress kernels, one 3-vector per
 *       plane-wave index and projector component.
 *
 *   The output dvnl array is packed as 3 Cartesian components per projector
 *   channel, ordered consistently with the legacy KBPP stress kernel layout.

 *
 * Method:
 *   - Construct cubic splines in |G| for each radial kernel A_p(|G|) and B_p(|G|).
 *   - For each reciprocal vector G, form G + k and its magnitude |G + k|.
 *   - Evaluate the radial kernels and their derivatives at |G + k|.
 *   - Assemble the vector-valued kernel corresponding to derivatives of the
 *     Kleinman–Bylander projector factors with respect to reciprocal vectors G.
 *
 * In the plane-wave formulation, the nonlocal stress contribution arises from
 * derivatives of projector factors of the form
 *
 *     D(|G+k|) T(u),   with  u = (G + k) / |G + k|.
 *
 * The resulting vector kernel is assembled as
 *
 *     ∂/∂G [ D(|G+k|) T(u) ]
 *       = D'(|G+k|) T(u) u
 *         + D(|G+k|) (∂T/∂u) (∂u/∂G),
 *
 * which is algebraically identical to the legacy Fortran implementation.
 *
 * Notes:
 *   - Brillouin-zone dependence enters explicitly through (G + k).
 *   - The routine is valid for Γ-point and general k-point calculations.
 *   - G + k = 0 contributions are explicitly zeroed to avoid singular terms.
 *   - The ordering and normalization of projector components matches the
 *     legacy NWPW KBPP stress implementation.
 */
void Psp1d_Hamann::cpp2_generate_stress_nonlocal_spline(CGrid *mygrid, 
                                                        int nbq, double *kvec,
                                                        int nray, double *G_ray,
                                                        double *dvnl_ray,
                                                        double *dvnl)
{
   /* set up indx(n,l) --> to wp (same mapping used by vpp_generate_ray/spline) */
   int indx[5 * 4];
   int nb = lmax + 1;
   for (auto l = 0; l <= lmax; ++l) {
     indx[l * 5] = l;
     for (auto n1 = 1; n1 < n_expansion[l]; ++n1) {
       indx[n1 + l * 5] = nb;
       ++nb;
     }
   }
 
   /* allocate spline second-derivative arrays */
   const int lmaxnray = (lmax + 1 + n_extra) * nray;
   double *dvnl_splineray     = new double[2 * lmaxnray];
   double *tmp_splineray      = new double[nray];
 
   double *dvnlD  = dvnl_ray;
   double *dvnlDD = dvnl_ray + lmaxnray;
 
   double *dvnlD_splineray  = dvnl_splineray;
   double *dvnlDD_splineray = dvnl_splineray + lmaxnray;
 
 
   /* setup cubic splines */
   const double dG = G_ray[2] - G_ray[1];
 
    /* nonlocal stress kernels: two splines per (n,l) channel: D and DD
       Layout in dvnl_ray:
         D  at dvnl_ray[ k + 0*lmaxnray + channel*nray ]
         DD at dvnl_ray[ k + 1*lmaxnray + channel*nray ]  == dvnl_ray[k + lmaxnray + channel*nray]
    */
    for (auto l = 0; l <= lmax; ++l)
       if (l != locp)
          for (auto n = 0; n < n_expansion[l]; ++n)
          {
             const int ch = indx[n + 5*l];
 
             util_spline(G_ray, &(dvnlD [ch*nray]),  nray, 0.0, 0.0, &(dvnlD_splineray [ch*nray]),  tmp_splineray);
             util_spline(G_ray, &(dvnlDD[ch*nray]),  nray, 0.0, 0.0, &(dvnlDD_splineray[ch*nray]),  tmp_splineray);
          }
 
 
       
   /* output packing */ 
   const int npack1_max = mygrid->npack1_max();
   const int npack1 = mygrid->npack(1+nbq);
       
   mygrid->t_pack_nzero(1, 3 * nprj, dvnl);   // 3 components per projector component
       
 
   /* ---- nonlocal on pack1 grid (vector kernels) ---- */
   {
      double *gx = mygrid->Gpackxyz(1, 0);
      double *gy = mygrid->Gpackxyz(1, 1);
      double *gz = mygrid->Gpackxyz(1, 2);
     
      for (auto k = 0; k < npack1; ++k) 
      {
         double Gx = gx[k] + kvec[0];
         double Gy = gy[k] + kvec[1];
         double Gz = gz[k] + kvec[2];
         const double q = std::sqrt(Gx*Gx + Gy*Gy + Gz*Gz);
         //const int nx   = (int)std::floor(q / dG);
         int nx = (int)std::floor(q / dG);
         if (nx < 0) nx = 0;
         if (nx > nray-2) nx = nray-2;

        
         if (q <= 1.0e-9) {
           for (auto p = 0; p < nprj; ++p) {
             dvnl[k + (0 + 3*p)*npack1_max] = 0.0;
             dvnl[k + (1 + 3*p)*npack1_max] = 0.0;
             dvnl[k + (2 + 3*p)*npack1_max] = 0.0;
           }
           continue;
         }

         /* unit vector u = (G+k)/|G+k| */
         const double ux = Gx / q;
         const double uy = Gy / q;
         const double uz = Gz / q;

         /* du_i / dG_j */
         const double duxdGx = 1.0/q - ux*ux/q;
         const double duxdGy = -ux*uy/q;
         const double duxdGz = -ux*uz/q;

         const double duydGx = -uy*ux/q;
         const double duydGy = 1.0/q - uy*uy/q;
         const double duydGz = -uy*uz/q;

         const double duzdGx = -uz*ux/q;
         const double duzdGy = -uz*uy/q;
         const double duzdGz = 1.0/q - uz*uz/q;

         int lcount = nprj;

         auto emit = [&](double D, double DD, double T,
                         double dTdux, double dTduy, double dTduz)
         {
            const double sumx = dTdux*duxdGx + dTduy*duydGx + dTduz*duzdGx;
            const double sumy = dTdux*duxdGy + dTduy*duydGy + dTduz*duzdGy;
            const double sumz = dTdux*duxdGz + dTduy*duydGz + dTduz*duzdGz;

            --lcount;
            dvnl[k + (0 + 3*lcount)*npack1_max] = DD*T*ux + D*sumx;
            dvnl[k + (1 + 3*lcount)*npack1_max] = DD*T*uy + D*sumy;
            dvnl[k + (2 + 3*lcount)*npack1_max] = DD*T*uz + D*sumz;
         };
        
        
         /* f projectors (l=3): 7 components */
         if ((locp != 3) && (lmax > 2))
            for (auto n = 0; n < n_expansion[3]; ++n)
            {
               const int ch = indx[n + 3*5];
        
               const double D  = util_splint(G_ray, &(dvnlD [ch*nray]), &(dvnlD_splineray [ch*nray]), nray, nx, q);
               const double DD = util_splint(G_ray, &(dvnlDD[ch*nray]), &(dvnlDD_splineray[ch*nray]), nray, nx, q);
        
               /* Copying exactly the Fortran T and dT/du blocks */
               {
                 double T = uy*(3.0*(1.0-uz*uz) - 4.0*uy*uy)/std::sqrt(24.0);
                 double dTdux = 0.0;
                 double dTduy = (3.0*(1.0-uz*uz) - 12.0*uy*uy)/std::sqrt(24.0);
                 double dTduz = -6.0*uy*uz/std::sqrt(24.0);
                 emit(D, DD, T, dTdux, dTduy, dTduz);
               }
               {
                 double T = ux*uy*uz;
                 double dTdux = uy*uz;
                 double dTduy = ux*uz;
                 double dTduz = ux*uy;
                 emit(D, DD, T, dTdux, dTduy, dTduz);
               }
               {
                 double T = uy*(5.0*uz*uz - 1.0)/std::sqrt(40.0);
                 double dTdux = 0.0;
                 double dTduy = (5.0*uz*uz - 1.0)/std::sqrt(40.0);
                 double dTduz = 10.0*uy*uz/std::sqrt(40.0);
                 emit(D, DD, T, dTdux, dTduy, dTduz);
               }
               {
                 double T = uz*(5.0*uz*uz - 3.0)/std::sqrt(60.0);
                 double dTdux = 0.0;
                 double dTduy = 0.0;
                 double dTduz = (15.0*uz*uz - 3.0)/std::sqrt(60.0);
                 emit(D, DD, T, dTdux, dTduy, dTduz);
               }
               {
                 double T = ux*(5.0*uz*uz - 1.0)/std::sqrt(40.0);
                 double dTdux = (5.0*uz*uz - 1.0)/std::sqrt(40.0);
                 double dTduy = 0.0;
                 double dTduz = 10.0*ux*uz/std::sqrt(40.0);
                 emit(D, DD, T, dTdux, dTduy, dTduz);
               }
               {
                 double T = uz*(ux*ux - uy*uy)/2.0;
                 double dTdux = ux*uz;
                 double dTduy = -uy*uz;
                 double dTduz = (ux*ux - uy*uy)/2.0;
                 emit(D, DD, T, dTdux, dTduy, dTduz);
               }
               {
                 double T = ux*(4.0*ux*ux - 3.0*(1.0-uz*uz))/std::sqrt(24.0);
                 double dTdux = (12.0*ux*ux - 3.0*(1.0-uz*uz))/std::sqrt(24.0);
                 double dTduy = 0.0;
                 double dTduz = 6.0*ux*uz/std::sqrt(24.0);
                 emit(D, DD, T, dTdux, dTduy, dTduz);
               }
            }
        
        
         /* d projectors (l=2): 5 components */
         if ((locp != 2) && (lmax > 1))
            for (auto n = 0; n < n_expansion[2]; ++n)
            {
               const int ch = indx[n + 2*5];
        
               const double D  = util_splint(G_ray, &(dvnlD [ch*nray]), &(dvnlD_splineray[ch*nray]), nray, nx, q);
               const double DD = util_splint(G_ray, &(dvnlDD[ch*nray]), &(dvnlDD_splineray[ch*nray]), nray, nx, q);
        
               {
                 double T = ux*uy;
                 emit(D, DD, T, uy, ux, 0.0);
               }
               {
                 double T = uy*uz;
                 emit(D, DD, T, 0.0, uz, uy);
               }
               {
                 double T = (3.0*uz*uz - 1.0)/(2.0*std::sqrt(3.0));
                 emit(D, DD, T, 0.0, 0.0, 6.0*uz/(2.0*std::sqrt(3.0)));
               }
               {
                 double T = uz*ux;
                 emit(D, DD, T, uz, 0.0, ux);
               }
               {
                 double T = (ux*ux - uy*uy)/2.0;
                 emit(D, DD, T, ux, -uy, 0.0);
               }
            }
        
         /* p projectors (l=1): 3 components */
         if ((locp != 1) && (lmax > 0))
            for (auto n = 0; n < n_expansion[1]; ++n)
            {
               const int ch = indx[n + 1*5];
        
               //const double D  = util_splint(G_ray, &(dvnl_ray[0    + ch*nray]),
               //                              &(dvnl_splineray[0    + ch*nray]), nray, nx, q);
               //const double DD = util_splint(G_ray, &(dvnl_ray[nray + ch*nray]),
               //                              &(dvnl_splineray[nray + ch*nray]), nray, nx, q);
               //const double D  = util_splint(G_ray, &(dvnl_ray[0*lmaxnray + ch*nray]), &(dvnl_splineray[0*lmaxnray + ch*nray]), nray, nx, q);
               //const double DD = util_splint(G_ray, &(dvnl_ray[1*lmaxnray + ch*nray]), &(dvnl_splineray[1*lmaxnray + ch*nray]), nray, nx, q);
        
               const double D  = util_splint(G_ray, &(dvnlD [ch*nray]), &(dvnlD_splineray[ch*nray]), nray, nx, q);
               const double DD = util_splint(G_ray, &(dvnlDD[ch*nray]), &(dvnlDD_splineray[ch*nray]), nray, nx, q);
        
               {
                 double T = uy;
                 emit(D, DD, T, 0.0, 1.0, 0.0);
               }
               {
                 double T = uz;
                 emit(D, DD, T, 0.0, 0.0, 1.0);
               }
               {
                 double T = ux;
                 emit(D, DD, T, 1.0, 0.0, 0.0);
               }
            }
        
         /* s projectors (l=0): 1 component */
         if (locp != 0)
            for (auto n = 0; n < n_expansion[0]; ++n)
            {
               const int ch = indx[n + 0*5];
        
               //const double D  = util_splint(G_ray, &(dvnl_ray[0    + ch*nray]),
               //                              &(dvnl_splineray[0    + ch*nray]), nray, nx, q);
               //const double D  = util_splint(G_ray, &(dvnl_ray[0*lmaxnray + ch*nray]),
               //                              &(dvnl_splineray[0*lmaxnray + ch*nray]), nray, nx, q);
               //const double DD = util_splint(G_ray, &(dvnl_ray[1*lmaxnray + ch*nray]),
               //                              &(dvnl_splineray[1*lmaxnray + ch*nray]), nray, nx, q);
               /* For s: DD isn’t used in Fortran (pure direction), but keep the same formula:
                  T=1, dT/du = 0 -> dvnl = DD*u */
               //const double DD = util_splint(G_ray, &(dvnl_ray[nray + ch*nray]),
               //                              &(dvnl_splineray[nray + ch*nray]), nray, nx, q);
        
               /* For s: DD isn’t used in Fortran (pure direction), but keep the same formula:
                  T=1, dT/du = 0 -> dvnl = DD*u */
               const double D  = util_splint(G_ray, &(dvnlD [ch*nray]), &(dvnlD_splineray[ch*nray]), nray, nx, q);
               const double DD = util_splint(G_ray, &(dvnlDD[ch*nray]), &(dvnlDD_splineray[ch*nray]), nray, nx, q);
        
               emit(D, DD, 1.0, 0.0, 0.0, 0.0);
            }
      }
   }
 
   /* cleanup */
   delete[] tmp_splineray;
   delete[] dvnl_splineray;
}












     
         










     
         








} // namespace pwdft

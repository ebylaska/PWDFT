
// extern "C" {
//#include        "compressed_io.h"
// }
#include "compressed_io.hpp"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "Parallel.hpp"
#include "util.hpp"
#include "util_date.hpp"
#include "iofmt.hpp"

namespace pwdft {

/*#include <algorithm>
void c_aindexcopy(const int n, const int *indx, double *A, double *B) {
    for (auto i = 0; i < n; ++i) {
        auto jj = 2 * indx[i];
        std::copy(A + jj, A + jj + 2, B + 2 * i);
    }
}
*/

/******************************************
 *                                        *
 *            transpose2DArray            *
 *                                        *
 ******************************************/
/**
 * @brief Converts a 2D array A(n2ft3d, nffts) to B(nffts, n2ft3d) assuming Fortran ordering.
 *
 * This function transposes the 2D array A from a shape of (n2ft3d, nffts) 
 * to a shape of (nffts, n2ft3d) and stores the result in array B, both 
 * assuming Fortran (column-major) ordering.
 *
 * @param n2ft3d The number of rows in the input array A.
 * @param nffts The number of columns in the input array A.
 * @param A Pointer to the input array.
 * @param B Pointer to the output array.
 */
void transpose2DArray(const int n2ft3d, const int nffts, double *A, double *B)
{
    for (int i = 0; i < n2ft3d; ++i)
    {
        for (int j = 0; j < nffts; ++j)
        {
            B[i * nffts + j] = A[j * n2ft3d + i];
        }
    }
}




/******************************************
 *                                        *
 *              c_aindexcopy              *
 *                                        *
 ******************************************/
/**
 * @brief Copies selected complex numbers from one array to another based on specified indices.
 *
 * This function reads complex numbers from the source array `A` at positions specified by the `indx` array
 * and places them sequentially into the destination array `B`. It is used to extract and reorder complex
 * numbers from one data structure to another, facilitating operations that require specific data alignments
 * or subsets.
 *
 * @param n The number of complex numbers to copy, as specified by the length of the `indx` array.
 * @param indx Array of indices indicating the positions in the source array `A` from which complex numbers
 *             should be extracted.
 * @param A Source array containing complex numbers. Each complex number is composed of two consecutive doubles
 *          (real and imaginary parts).
 * @param B Destination array where the complex numbers are to be stored sequentially.
 *
 * @note The function assumes that `B` is large enough to accommodate `n` complex numbers. It does not perform
 *       bounds checking on array accesses, and it assumes that `indx` contains valid indices into `A`.
 */
void c_aindexcopy(const int n, const int *indx, double *A, double *B) 
{
   int ii = 0;
   for (auto i=0; i<n; ++i) 
   {
      auto jj = 2*indx[i];
      B[ii]   = A[jj];
      B[ii+1] = A[jj+1];
      ii += 2;
   }
}


/******************************************
 *                                        *
 *         c_aindexcopy_stride            *
 *                                        *
 ******************************************/
/**
 * @brief Copies selected complex numbers from one array to another with a specified stride in the destination array.
 *
 * This function copies complex numbers from the source array `A` to the destination array `B` based on specified indices.
 * Each complex number is taken from `A` at positions indexed by `indx` and placed in `B` at intervals defined by `stride`.
 * This allows complex numbers to be spaced out or reorganized in `B`, facilitating different data layouts or alignments.
 *
 * @param stride The gap between successive complex numbers in `B`, measured in terms of the number of real values (each complex number has 2 real values).
 * @param n The number of complex numbers to copy.
 * @param indx Array of indices where each entry specifies the position of a complex number in `A` to be copied to `B`.
 * @param A Source array containing complex numbers, with each complex number composed of two consecutive doubles (real and imaginary parts).
 * @param B Destination array where the complex numbers are to be stored according to the stride.
 *
 * @note It is assumed that `B` is sufficiently large to accommodate the final element being copied into it, considering the `stride`.
 *       The function does not perform bounds checking on array accesses.
 */
void c_aindexcopy_stride(const int stride, const int n, const int *indx, double *A, double *B) 
{
   int ii = 0;
   for (auto i=0; i<n; ++i) 
   {
      auto jj = 2*indx[i];
      B[ii]   = A[jj];
      B[ii+1] = A[jj+1];
      ii += 2*stride;
   }
}


/******************************************
 *                                        *
 *             c_bindexcopy               *
 *                                        *
 ******************************************/
/**
 * @brief Copies complex numbers from one array to specified positions in another array based on given indices.
 *
 * This function reads complex numbers sequentially from the source array `A` and places them into the destination
 * array `B` at positions specified by the `indx` array. Each index in `indx` corresponds to the position in `B` where
 * a complex number from `A` should be placed, making it suitable for reorganizing or redistributing data.
 *
 * @param n The number of complex numbers to copy.
 * @param indx Array of indices indicating the target positions in `B` for the complex numbers. The indices are
 *             expected to be in the range that is valid within the bounds of `B`.
 * @param A Source array containing complex numbers. Each complex number consists of two consecutive doubles
 *          (real and imaginary parts).
 * @param B Destination array where the complex numbers are placed according to the indices in `indx`.
 *
 * @note The function assumes that `B` is sufficiently large to accommodate the indices provided in `indx`. It does not perform
 *       bounds checking on array accesses, and it assumes that `indx` contains valid indices into `B`.
 */
void c_bindexcopy(const int n, const int *indx, double *A, double *B) 
{
   int ii = 0;
   for (auto i=0; i<n; ++i) 
   {
      auto jj = 2*indx[i];
      B[jj]   = A[ii];
      B[jj+1] = A[ii+1];
      ii += 2;
   }
}

/******************************************
 *                                        *
 *         c_bindexcopy_stride            *
 *                                        *
 ******************************************/
 /**
 * @brief Copies complex numbers from one array to specified indexed positions in another array with a source stride.
 *
 * This function reads complex numbers from the source array `A` and places them into the destination array `B` at positions
 * specified by the `indx` array. The copying from `A` is controlled by a stride, allowing for skipping elements in `A`.
 *
 * @param stride The gap between successive reads in the source array `A`, measured in terms of the number of real values (each complex number has 2 real values).
 * @param n The number of complex numbers to copy.
 * @param indx Array of indices where each entry specifies the position in `B` to store a complex number from `A`.
 * @param A Source array containing complex numbers, with each complex number composed of two consecutive doubles (real and imaginary parts).
 * @param B Destination array where the complex numbers are to be stored according to the indices in `indx`.
 *
 * @note It is assumed that `B` is sufficiently large to accommodate the highest index accessed based on `indx` and the number of elements `n`.
 *       The function does not perform bounds checking on array accesses.
 */
void c_bindexcopy_stride(const int stride, const int n, const int *indx, double *A, double *B)
{
   int ii = 0;
   for (auto i=0; i<n; ++i)
   {
      auto jj = 2*indx[i];
      B[jj]   = A[ii];
      B[jj+1] = A[ii+1];
      ii += 2*stride;
   }
}


/******************************************
 *                                        *
 *           c_bindexcopy_conjg           *
 *                                        *
 ******************************************/
/**
 * @brief Copies and conjugates complex numbers from one array to specified positions in another array.
 *
 * This function reads complex numbers sequentially from the source array `A`, conjugates them (negates the imaginary part),
 * and places them into the destination array `B` at positions specified by the `indx` array. This operation is useful for
 * mathematical computations that require the complex conjugate of a set of numbers.
 *
 * @param n The number of complex numbers to process.
 * @param indx Array of indices indicating the target positions in `B` for the conjugated complex numbers. The indices are
 *             assumed to be valid within the bounds of `B`.
 * @param A Source array containing complex numbers. Each complex number consists of two consecutive doubles
 *          (real and imaginary parts).
 * @param B Destination array where the conjugated complex numbers are placed according to the indices in `indx`.
 *
 * @note It is assumed that `B` is sufficiently large to accommodate the highest index provided in `indx`. The function does
 *       not perform bounds checking on array accesses, and it assumes that `indx` contains valid indices into `B`.
 */
void c_bindexcopy_conjg(const int n, const int *indx, double *A, double *B) 
{
   int ii, jj;
   ii = 0;
   for (auto i = 0; i < n; ++i) 
   {
      jj = 2 * indx[i];
      B[jj] = A[ii];
      B[jj + 1] = -A[ii + 1];
      ii += 2;
   }
}

/******************************************
 *                                        *
 *       c_bindexcopy_conjg_stride        *
 *                                        *
 ******************************************/
/**
 * @brief Copies and conjugates complex numbers from one array to indexed positions in another array using a specified stride.
 *
 * This function reads complex numbers from the source array `A`, conjugates them, and places them into the destination
 * array `B` at positions specified by the `indx` array. The stride parameter controls the reading interval in `A`, allowing
 * for skipping elements between reads.
 *
 * @param stride The gap between successive reads in the source array `A`, measured in terms of the number of real values 
 *               (each complex number is composed of two real values).
 * @param n The number of complex numbers to process.
 * @param indx Array of indices where each entry specifies the position in `B` to store the conjugated complex number from `A`.
 * @param A Source array containing complex numbers.
 * @param B Destination array where the conjugated complex numbers are to be stored.
 *
 * @note It is assumed that `B` is sufficiently large to accommodate the highest index accessed based on `indx` and `n`.
 *       The function does not perform bounds checking on array accesses. Ensure that `A` has at least `stride * (n-1) + 1`
 *       complex numbers to avoid reading beyond the end of the array.
 */
void c_bindexcopy_conjg_stride(const int stride, const int n, const int *indx, double *A, double *B)
{
   int ii, jj;
   ii = 0;
   for (auto i = 0; i < n; ++i)
   {
      jj = 2 * indx[i];
      B[jj] = A[ii];
      B[jj + 1] = -A[ii + 1];
      ii += 2*stride;
   }
}


/******************************************
 *                                        *
 *             c_bindexzero               *
 *                                        *
 ******************************************/
/**
 * @brief Sets specified complex numbers in an array to zero.
 *
 * This function zeroes complex numbers in the array `B` at positions specified by the `indx` array. Each index in `indx`
 * corresponds to a complex number in `B` that will be set to zero. This is typically used to selectively reset elements 
 * in an array, useful in various numerical and signal processing applications.
 *
 * @param n The number of complex numbers to zero out, as specified by the length of the `indx` array.
 * @param indx Array of indices indicating the positions of the complex numbers in `B` that should be zeroed.
 * @param B Array containing complex numbers where specified elements will be set to zero. Each complex number is composed
 *          of two consecutive doubles (real and imaginary parts).
 *
 * @note The function assumes that `B` is large enough to contain all indices specified in `indx`. It does not perform
 *       bounds checking on array accesses, and it assumes that `indx` contains valid indices into `B`.
 */
void c_bindexzero(const int n, const int *indx, double *B) 
{
   int jj;
   for (auto i=0; i<n; ++i) 
   {
      jj = 2 * indx[i];
      B[jj]   = 0.0;
      B[jj+1] = 0.0;
   }
}

/******************************************
 *                                        *
 *             t_aindexcopy               *
 *                                        *
 ******************************************/
void t_aindexcopy(const int n, const int *indx, double *A, double *B) 
{
   for (auto i = 0; i < n; ++i)
      B[i] = A[indx[i]];
}

/******************************************
 *                                        *
 *             t_bindexcopy               *
 *                                        *
 ******************************************/
void t_bindexcopy(const int n, const int *indx, double *A, double *B) 
{
   for (auto i = 0; i < n; ++i)
      B[indx[i]] = A[i];
}

/******************************************
 *                                        *
 *             i_aindexcopy               *
 *                                        *
 ******************************************/
void i_aindexcopy(const int n, const int *indx, int *A, int *B) 
{
   for (auto i = 0; i < n; ++i)
      B[i] = A[indx[i]];
}


/**************************************
 *                                    *
 *            eigsrt                  *
 *                                    *
 **************************************/
/**
 * @brief Sorts eigenvalues in descending order and reorders corresponding eigenvectors.
 *
 * This function sorts the eigenvalues contained in array `D` in descending order. It also reorders the
 * corresponding eigenvectors in matrix `V` to match the sorted eigenvalues. Matrix `V` should be stored
 * in column-major format, where each column corresponds to an eigenvector.
 *
 * @param D Array containing the eigenvalues to be sorted.
 * @param V Matrix where each column represents an eigenvector corresponding to the eigenvalue in `D`.
 * @param n The number of eigenvalues and eigenvectors, also the dimension of the square matrix `V`.
 *
 * @note This function modifies both `D` and `V` in-place. Ensure that `D` and `V` are properly initialized
 *       and that `V` has enough space allocated to accommodate the matrix of size `n` by `n`.
 */
void eigsrt(double *D, double *V, int n) 
{
   int i, j, k;
   double p;
 
   for (i = 0; i < (n - 1); ++i) 
   {
      k = i;
      p = D[i];
      for (j = i + 1; j < n; ++j)
         if (D[j] >= p) 
         {
            k = j;
            p = D[j];
         }
     
      if (k != i) 
      {
         D[k] = D[i];
         D[i] = p;
         for (j = 0; j < n; ++j) 
         {
            p = V[j + i * n];
            V[j + i * n] = V[j + k * n];
            V[j + k * n] = p;
         }
      }
   }
}


/**************************************
 *                                    *
 *         util_getfilling            *
 *                                    *
 **************************************/
void util_getfilling(int f, int nfft[], int *filling, double zvalue[]) 
{
   int h = (f % 2);
   int f2 = (f - h) / 2;
   int k = 0;
   while ((k + 1) * (k + 2) * (k + 3) <= (6 * f2))
      ++k;
 
   f2 -= k * (k + 1) * (k + 2) / 6;
   int j = 0;
   while ((j + 1) * (j + 2) <= (2 * f2))
      ++j;
 
   int i = f2 - j * (j + 1) / 2;
 
   filling[0] = i;
   if (i == 0) 
   {
      filling[1] = j;
      filling[2] = k;
   } 
   else 
   {
      if (((j + 1) % 2) == 0)
         filling[1] = (j + 1) / 2;
      else
         filling[1] = -(j + 1) / 2;
     
      if (((k + 1) % 2) == 0)
         filling[2] = (k + 1) / 2;
      else
         filling[2] = -(k + 1) / 2;
   }
 
   if ((i == 0) && (j == 0) && (k == 0)) 
   {
      filling[3] = 0;
      zvalue[0] = 1.0;
      zvalue[1] = 0.0;
   } 
   else if (h == 0) 
   {
      filling[3] = 2;
      zvalue[0] = 1.0 / sqrt(2.0);
      zvalue[1] = 0.0;
   } 
   else 
   {
      filling[3] = -2;
      zvalue[0] = 0.0;
      zvalue[1] = 1.0 / sqrt(2.0);
   }
 
   /* modularize the filling */
   int inc2c = (nfft[0] / 2 + 1);
   filling[0] = (filling[0] + inc2c) % inc2c;
   filling[1] = (filling[1] + nfft[1]) % nfft[1];
   filling[2] = (filling[2] + nfft[2]) % nfft[2];
}

/**************************************
 *                                    *
 *           util_random              *
 *                                    *
 **************************************/

/* returns a random number between 0 and 1

   Entry - seed - if zero set a seed with srand
   Exit - returns a random number between 0 and 1
   Uses - rand and srand stdlib functions
*/
double util_random(const int seed) 
{
   if (seed > 0)
      std::srand(((double)seed));
   return ((double)std::rand() / RAND_MAX);
}

/**************************************
 *                                    *
 *           util_filefind            *
 *                                    *
 **************************************/
bool util_filefind(Parallel *myparall, char *fname) 
{
   int ifound;
 
   if (myparall->is_master())
      ifound = cfileexists(fname);
 
   myparall->Brdcst_iValue(0, 0, &ifound);
 
   return (ifound > 0);
}

/**************************************
 *                                    *
 *           util_spline              *
 *                                    *
 **************************************/
void util_spline(const double *x, const double *y, const int n,
                 const double yp1, const double ypn, double *y2, double *utmp) 
{
   double sig, p, qn, un;
 
   if (yp1 > 0.99e30) 
   {
      y2[0] = 0.0;
      utmp[0] = 0.0;
   } 
   else 
   {
      y2[0] = -0.5;
      utmp[0] = 3.0 / (x[1] - x[0]) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
   }
   for (auto i = 1; i < (n - 1); ++i) 
   {
      sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
      p = sig * y2[i - 1] + 2.00;
      y2[i] = (sig - 1.00) / p;
      utmp[i] = (6.00 *
                     ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
                      (y[i] - y[i - 1]) / (x[i] - x[i - 1])) /
                     (x[i + 1] - x[i - 1]) -
                 sig * utmp[i - 1]) / p;
   }
 
   if (ypn > 0.99e30) 
   {
      qn = 0.0;
      un = 0.0;
   } 
   else 
   {
      qn = 0.5;
      un = 3.00 / (x[n - 1] - x[n - 2]) *
           (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
   }
 
   y2[n - 1] = (un - qn * utmp[n - 2]) / (qn * y2[n - 2] + 1.00);
   for (auto k = n - 2; k >= 0; --k)
      y2[k] = y2[k] * y2[k + 1] + utmp[k];
}

/**************************************
 *                                    *
 *           util_splint              *
 *                                    *
 **************************************/
double util_splint(const double *xa, const double *ya, const double *y2a,
                   const int n, const int nx, const double x) 
{
   int khi = nx;
   int klo = nx - 1;
 
   while ((xa[klo] > x) || (xa[khi] < x)) 
   {
      if (xa[klo] > x) 
      {
         klo = klo - 1;
         khi = khi - 1;
      }
      if (xa[khi] < x) 
      {
         klo = klo + 1;
         khi = khi + 1;
      }
   }
 
   double h = xa[khi] - xa[klo];
   double a = (xa[khi] - x) / h;
   double b = (x - xa[klo]) / h;
 
   return (a * ya[klo] + b * ya[khi] +
           ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * h * h / 6.0);
}

/**************************************
 *                                    *
 *           util_filter              *
 *                                    *
 **************************************/
void util_filter(int nray, double *g_ray, double ecut, double *v_ray) {
  int ncut = 15;
  double qmax = sqrt(ecut + ecut);
  double g;

  for (auto i = 0; i < nray; ++i) {
    g = g_ray[i];
    if (g > (qmax - 0.2)) {
      v_ray[i] *= (1.0 - pow((1.0 - exp(-pow((g / qmax), ncut))), ncut));
    }
  }
}


/**************************************
 *                                    *
 *        util_fattebert_dielec       *
 *                                    *
 **************************************/
/* Computes the dielectric function using the local model of the density from
   Jean-Luc Fattebert. Reference: Fattebert JL, Gygi F. First‐principles
   molecular dynamics simulations in a continuum solvent. International journal
   of quantum chemistry. 2003;93(2):139-47. Suggested parameter values: water -
   eps=78.36, rho0=0.0004 au, and beta=1.3 bohr

*/
void util_fattebert_dielec(const int n2ft3d, const double eps,
                           const double beta, const double rho0,
                           const double *rho, double *epsilon) {
  for (auto k = 0; k < n2ft3d; ++k) {
    auto x = std::pow((rho[k] / rho0), 2.0 * beta);
    epsilon[k] = 1.0 + 0.5 * (eps - 1.0) * (1.0 + (1.0 - x) / (1.0 + x));
  }
}

/**************************************
 *                                    *
 *        util_dfattebert_dielec      *
 *                                    *
 **************************************/
/* Computes the derivative of dielectric function wrt rho using the local model
   of the density from Jean-Luc Fattebert. Reference: Fattebert JL, Gygi F.
   First‐principles molecular dynamics simulations in a continuum solvent.
              International journal of quantum chemistry. 2003;93(2):139-47.
   Suggested parameter values:
      water - eps=78.36, rho0=0.0004 au, and beta=1.3 bohr

*/
void util_dfattebert_dielec(const int n2ft3d, const double eps,
                            const double beta, const double rho0,
                            const double *rho, double *depsilon) {
  for (auto k = 0; k < n2ft3d; ++k) {
    auto x = std::pow((rho[k] / rho0), 2.0 * beta);
    auto y = std::pow((rho[k] / rho0), 2.0 * beta - 1.0) / rho0;
    auto z = (x + 1.0) * (x + 1.0);
    // epsilon[k] = 1.0 + 0.5*(eps-1.0)*(1.0 + (1.0-x)/(1.0+x));
    depsilon[k] = -2.0 * beta * (eps - 1.0) * y / z;
  }
}

/**************************************
 *                                    *
 *   util_weighted_fattebert_dielec   *
 *                                    *
 **************************************/
/* Computes the weighted dielectric function using the local model of the
   density from Jean-Luc Fattebert. Reference: Fattebert JL, Gygi F.
   First‐principles molecular dynamics simulations in a continuum solvent.
              International journal of quantum chemistry. 2003;93(2):139-47.
   Suggested parameter values:
      water - eps=78.36, rho0=0.0004 au, and beta=1.3 bohr

*/
void util_weighted_fattebert_dielec(const int n2ft3d, const double eps,
                                    const double beta, const double rho0,
                                    const double *rho, const double *w,
                                    double *epsilon) {
  for (auto k = 0; k < n2ft3d; ++k) {
    auto x = std::pow((rho[k] / rho0), 2 * beta);
    epsilon[k] = 1.0 + w[k] * 0.5 * (eps - 1.0) * (1.0 + (1.0 - x) / (1.0 + x));
  }
}

/**************************************
 *                                    *
 *        util_andreussi_dielec       *
 *                                    *
 **************************************/
/* Computes the dielectric function using the first local model of the density
   from Andreussi, Dabo, and Marzari Reference: Andreussi, Dabo, and Marzari, J.
   Chem. Phys.136, 064102 (2012) Suggested parameter values: water - eps=78.36,
   rhomin=0.0001 au, and rhomax=0.0035 bohr ρmax=0.0035 a.u. and ρmin=0.0001
   a.u.

*/
void util_andreussi_dielec(const int n2ft3d, const double eps,
                           const double rhomin, const double rhomax,
                           const double *rho, double *epsilon) 
{
   double twopi = 8.0*std::atan(1.0);
 
   for (auto k=0; k<n2ft3d; ++k) 
   {
      if (rho[k]>rhomax)
         epsilon[k] = 1.0;
      else if (rho[k]<rhomin)
         epsilon[k] = eps;
      else 
      {
         auto x = twopi*(rhomax-rho[k])/(rhomax-rhomin);
         epsilon[k] = 1.0+((eps-1.0)/twopi)*(x-std::sin(x));
      }
   }
}

/**************************************
 *                                    *
 *        util_dandreussi_dielec      *
 *                                    *
 **************************************/
/* Computes the derivative of dielectric function wrt to rho  using the first
   local model of the density from Andreussi, Dabo, and Marzari Reference:
   Andreussi, Dabo, and Marzari, J. Chem. Phys.136, 064102 (2012) Suggested
   parameter values: water - eps=78.36, rhomin=0.0001 au, and rhomax=0.0035 bohr
      ρmax=0.0035 a.u. and ρmin=0.0001 a.u.

*/
void util_dandreussi_dielec(const int n2ft3d, const double eps,
                            const double rhomin, const double rhomax,
                            const double *rho, double *depsilon) 
{
   double twopi = 8.0*std::atan(1.0);
 
   for (auto k=0; k<n2ft3d; ++k) 
   {
      if ((rho[k] > rhomin) && (rho[k] < rhomax)) 
      {
         auto x = twopi*(rhomax-rho[k])/(rhomax-rhomin);
         auto dxdrho = -twopi/(rhomax-rhomin);
         depsilon[k] = ((eps-1.0)/twopi)*(1.0-std::cos(x))*dxdrho;
      } 
      else
         depsilon[k] = 0.0;
   }
}

/**************************************
 *                                    *
 *        util_ddandreussi_dielec     *
 *                                    *
 **************************************/
/* Computes the second derivative of dielectric function wrt to rho  using the
   first local model of the density from Andreussi, Dabo, and Marzari Reference:
   Andreussi, Dabo, and Marzari, J. Chem. Phys.136, 064102 (2012) Suggested
   parameter values: water - eps=78.36, rhomin=0.0001 au, and rhomax=0.0035 bohr
      ρmax=0.0035 a.u. and ρmin=0.0001 a.u.

*/
void util_ddandreussi_dielec(const int n2ft3d, const double eps,
                             const double rhomin, const double rhomax,
                             const double *rho, double *ddepsilon) 
{
   double twopi = 8.0 * std::atan(1.0);
 
   for (auto k=0; k<n2ft3d; ++k) 
   {
      if ((rho[k]>rhomin) && (rho[k]<rhomax)) 
      {
         auto x = twopi * (rhomax-rho[k])/(rhomax-rhomin);
         auto dxdrho = -twopi / (rhomax - rhomin);
         ddepsilon[k] = ((eps-1.0)/twopi)*(std::sin(x))*dxdrho*dxdrho;
      } 
      else
         ddepsilon[k] = 0.0;
   }
}

/**************************************
 *                                    *
 *        util_andreussi2_dielec      *
 *                                    *
 **************************************/
/* Computes the dielectric function using the second local model of the
   density from Andreussi, Dabo, and Marzari
   Reference: Andreussi, Dabo, and Marzari,
              J. Chem. Phys.136, 064102 (2012)
   Suggested parameter values:
      water - eps=78.36, rhomin=0.0001 au, and rhomax=0.0035 bohr
      ρmax=0.0035 a.u. and ρmin=0.0001 a.u.

*/
void util_andreussi2_dielec(const int n2ft3d, const double eps,
                            const double rhomin, const double rhomax,
                            const double *rho, double *epsilon)
{
   double twopi = 8.0*std::atan(1.0);
 
   for (auto k = 0; k < n2ft3d; ++k)
   {
      if (rho[k]>rhomax)
         epsilon[k] = 1.0;
      else if (rho[k]<rhomin)
         epsilon[k] = eps;
      else 
      {
         auto x = twopi*(std::log(rhomax)-std::log(rho[k]))/(std::log(rhomax)-std::log(rhomin));
         auto t = (std::log(eps)/twopi)*(x-std::sin(x));
         epsilon[k] = std::exp(t);
      }
   }
}

/**************************************
 *                                    *
 *        util_dandreussi2_dielec     *
 *                                    *
 **************************************/
/* Computes the derivative of the dielectric function using the second local
   model of the density from Andreussi, Dabo, and Marzari Reference: Andreussi,
   Dabo, and Marzari, J. Chem. Phys.136, 064102 (2012) Suggested parameter
   values: water - eps=78.36, rhomin=0.0001 au, and rhomax=0.0035 bohr
      ρmax=0.0035 a.u. and ρmin=0.0001 a.u.

*/
void util_dandreussi2_dielec(const int n2ft3d, const double eps,
                             const double rhomin, const double rhomax,
                             const double *rho, double *depsilon)
{
   double twopi = 8.0*std::atan(1.0);
 
   for (auto k=0; k<n2ft3d; ++k) 
   {
      if ((rho[k]>rhomin) && (rho[k]<rhomax)) 
      {
         auto x = twopi*(std::log(rhomax)-std::log(rho[k]))/(std::log(rhomax)-std::log(rhomin));
         auto t = (std::log(eps)/twopi)*(x-std::sin(x));
         auto dtdx = (std::log(eps)/twopi)*(1.0-std::cos(x));
         auto dxdrho = twopi/(rho[k]*(std::log(rhomin)-std::log(rhomax)));
         depsilon[k] = std::exp(t)*dtdx*dxdrho;
      } 
      else
         depsilon[k] = 0.0;
   }   
}

/**************************************
 *                                    *
 *        util_ddandreussi2_dielec    *
 *                                    *
 **************************************/
/* Computes the second derivative dielectric function using the second local
   model of the density from Andreussi, Dabo, and Marzari Reference: Andreussi,
   Dabo, and Marzari, J. Chem. Phys.136, 064102 (2012) Suggested parameter
   values: water - eps=78.36, rhomin=0.0001 au, and rhomax=0.0035 bohr
      ρmax=0.0035 a.u. and ρmin=0.0001 a.u.

*/
void util_ddandreussi2_dielec(const int n2ft3d, const double eps,
                              const double rhomin, const double rhomax,
                              const double *rho, double *ddepsilon) {
  double twopi = 8.0 * std::atan(1.0);

  for (auto k = 0; k < n2ft3d; ++k) {
    if ((rho[k] > rhomin) && (rho[k] < rhomax)) {
      auto x =
          twopi * (log(rhomax) - log(rho[k])) / (log(rhomax) - log(rhomin));
      auto t = (log(eps) / twopi) * (x - sin(x));
      auto dtdx = (log(eps) / twopi) * (1.0 - cos(x));
      auto d2tdx2 = (log(eps) / twopi) * (sin(x));
      auto dxdrho = twopi / (rho[k] * (log(rhomin) - log(rhomax)));
      auto d2xdrho2 = twopi / (rho[k] * rho[k] * (log(rhomax) - log(rhomin)));
      ddepsilon[k] = exp(t) * (d2tdx2 * dxdrho * dxdrho + dtdx * d2xdrho2 +
                               dtdx * dtdx * dxdrho * dxdrho);
    } else
      ddepsilon[k] = 0.0;
  }
}

/**************************************
 *                                    *
 *        util_sphere_dielec          *
 *                                    *
 **************************************/
/* 
*/
void util_sphere_dielec(const int n2ft3d, const double *rgrid,const double *rshift,
                        const double eps, const double R0, const double R1, 
                        double *epsilon)
{
   double s_d   = R0;
   double s_rho = R1-s_d;

 
   for (auto i=0; i<n2ft3d; ++i)
   {
      double dx = rgrid[3*i]   - rshift[0];
      double dy = rgrid[3*i+1] - rshift[1];
      double dz = rgrid[3*i+2] - rshift[2];
      double r  = std::sqrt(dx*dx+ dy*dy + dz*dz);
      epsilon[i] = 1+(eps-1)*util_switching_function(s_d,s_rho,r);
   }
}

/**************************************
 *                                    *
 *    util_sphere_gradient_dielec     *
 *                                    *
 **************************************/
void util_sphere_gradient_dielec(const int n2ft3d, const double *rgrid, const double *rshift,
                                 const double eps, const double R0, const double R1, 
                                 double *epsilon_x, double *epsilon_y, double *epsilon_z)
{
   double s_d   = R0;
   double s_rho = R1-s_d;

   for (auto i=0; i<n2ft3d; ++i)
   {
      double dx = rgrid[3*i]   - rshift[0];
      double dy = rgrid[3*i+1] - rshift[1];
      double dz = rgrid[3*i+2] - rshift[2];
      double r  = std::sqrt(dx*dx + dy*dy + dz*dz);
      double dfdr = (eps-1)*util_dswitching_function(s_d,s_rho,r);
      if (r>1.0e-6)
      {
         epsilon_x[i] = (dfdr)*dx/r;
         epsilon_y[i] = (dfdr)*dy/r;
         epsilon_z[i] = (dfdr)*dz/r;
      }
      else
      {
         epsilon_x[i] = 0.0;
         epsilon_y[i] = 0.0;
         epsilon_z[i] = 0.0;
      }
   }
}


/**************************************
 *                                    *
 *        util_switching_function     *
 *                                    *
 **************************************/
/*  The switching function is defined by
    where s = 0.0                              if r<=s_d
          s = 1.0-(1.0-(r-s_d)**2/s_rho**2)**2 if r>s_d and r<(s_d+s_rho)
          s = 1.0                              if r>(s_d+s_rho)
*/
double util_switching_function(const double s_d, const double s_rho,
                               const double r) {
  // calculate dielectric switching function
  double eps = 1.0;

  if (r <= s_d)
    eps = 0.0;
  else if (r < (s_d + s_rho)) {
    auto x = (r - s_d) / s_rho;
    auto y = 1.0 - x * x;
    eps = 1.0 - y * y;
  }
  return eps;
}

/**************************************
 *                                    *
 *       util_dswitching_function     *
 *                                    *
 **************************************/
/* The dswitching/dr function is defined by
   where  s = 0.0                                               if r<=s_d
          s = 4.0*((r-s_d)/s_rho**2)*(1.0-(r-s_d)**2/s_rho**2)  if r>s_d and
   r<(s_d+s_rho) s = 0.0                                               if
   r>(s_d+s_rho)
*/
double util_dswitching_function(const double s_d, const double s_rho, const double r) 
{
   // calculate the derivative of switching function
   double deps = 0.0;
 
   if (r <= s_d)
     deps = 0.0;
   else if (r < (s_d + s_rho)) {
     auto x = (r - s_d) / s_rho;
     auto xx = 1.0 - x * x;
     deps = 4.0 * (x / s_rho) * xx;
   }
   return deps;
}


/**************************************
 *                                    *
 *    util_occupation_distribution    *
 *                                    *
 **************************************/
/**
 * @brief Computes the occupation probability based on the specified smearing algorithm.
 *
 * This function calculates the occupation distribution for a given energy `e` and
 * smearing type `smeartype`. Different smearing algorithms are implemented to model
 * physical systems with partial occupancies or broadened energy states. Each smearing
 * type corresponds to a specific algorithm, with options for Fermi-Dirac, Gaussian,
 * Hermite, Marzari-Vanderbilt, Methfessel-Paxton, Cold Smearing, and Lorentzian distributions.
 *
 * @param smeartype An integer specifying the smearing type:
 *        - 1: Fermi-Dirac smearing
 *        - 2: Gaussian smearing
 *        - 3: Hermite smearing
 *        - 4: Marzari-Vanderbilt smearing
 *        - 5: Methfessel-Paxton smearing (first order)
 *        - 6: Cold smearing
 *        - 7: Lorentzian smearing
 *        - Other: Step function (0 for e > 0, 1 for e <= 0)
 * @param e The energy value for which the occupation is calculated.
 *          - Large positive and negative values are handled explicitly in some cases
 *            to avoid numerical overflow.
 * @return A double representing the occupation probability based on the specified
 *         smearing type. The result typically falls between 0.0 and 1.0, except for
 *         cases like Lorentzian smearing, where the value may not be normalized.
 *
 * @details
 * - **Fermi-Dirac (1)**: A standard statistical distribution used for electronic occupations.
 * - **Gaussian (2)**: Uses a Gaussian function for smoothing the distribution.
 * - **Hermite (3)**: Incorporates a third-order Hermite polynomial correction.
 * - **Marzari-Vanderbilt (4)**: Combines exponential and complementary error functions.
 * - **Methfessel-Paxton (5)**: Uses Hermite polynomials for higher-order corrections
 *   to Gaussian smearing (first-order implemented here).
 * - **Cold Smearing (6)**: Avoids the Gibbs phenomenon and improves energy convergence.
 * - **Lorentzian (7)**: Uses a Lorentzian function, often for density of states.
 * - **Default (Other)**: A simple step function.
 *
 * @note Ensure proper parameterization (e.g., smearing width for Lorentzian) if additional
 *       customization is required for specific smearing types.
 *
 * Example Usage:
 * @code
 * double occupation = util_occupation_distribution(1, -0.5); // Fermi-Dirac
 * std::cout << "Occupation: " << occupation << std::endl;
 * @endcode
 */
double util_occupation_distribution(const int smeartype, const double e)
{
   const double sqrt_pi = std::sqrt(M_PI);
   const double sqrt_two_pi = std::sqrt(2.0 * M_PI);
   const double sqrt_half = std::sqrt(0.5);

   double f = 0.0;

   if (smeartype == 1) { // Fermi-Dirac
      if (e > 30.0) {
          f = 0.0;
      } else if (e < -30.0) {
          f = 1.0;
      } else {
          f = 1.0 / (1.0 + std::exp(e));
      }
   } else if (smeartype == 2) { // Gaussian
      f = 0.5 * std::erfc(e);
   } else if (smeartype == 3) { // Hermite smearing
      double exp_term = std::exp(-e * e);     // Gaussian term
      double hermite_correction = (2.0 * e * e - 1.0) * exp_term; // Hermite correction term
 
      // Compute the original Hermite-Gaussian function
      double f_original = exp_term + hermite_correction / sqrt_pi;
 
      // Find fmax by evaluating the original function at e=0
      double fmax = exp_term + (2.0 * 0 * 0 - 1.0) * exp_term / sqrt_pi; // f_original at e=0
 
      if (e <= 0.0) {
          // Set value to 1.0 for e <= e_max
          f = 1.0;
      } else {
          // Scale the Hermite-Gaussian function for e > e_max
          f = f_original / fmax;
      }
   } else if (smeartype == 4) { // Marzari-Vanderbilt
      double factor = std::sqrt(0.125 / std::atan(1.0)); // atan(1.0) = pi/4
      f = std::exp(-(e + sqrt_half) * (e + sqrt_half)) * factor + 0.5 * std::erfc(e + sqrt_half);
   } else if (smeartype == 5) { // Methfessel-Paxton
      double exp_term = std::exp(-e * e);
      double hermite_poly = 1.0 - 2.0 * e * e; // First-order Methfessel-Paxton
      f = exp_term * hermite_poly / std::sqrt(M_PI);
   } else if (smeartype == 6) { // Cold Smearing
      double exp_term = std::exp(-0.5 * e * e);
      double erfc_term = 0.5 * std::erfc(-e / std::sqrt(2.0));
      f = erfc_term - e * exp_term / std::sqrt(2.0 * M_PI);
   } else if (smeartype == 7) { // Lorentzian Smearing
      double sigma = 1.0; // Assume a default smearing width
      f = (sigma / M_PI) / (e * e + sigma * sigma);
   } else { // Default: Step function
      if (e > 0.0) {
          f = 0.0;
      } else {
          f = 1.0;
      }
   }

   return f;
}

/**************************************
 *                                    *
 *        util_smearcorrection        *
 *                                    *
 **************************************/
/**
 * @brief Computes the correction term for smearing algorithms.
 *
 * This function calculates the correction term based on the specified smearing type,
 * which accounts for the influence of energy-level broadening in physical systems.
 * The correction is used to adjust energy contributions for improved numerical stability
 * and convergence in simulations.
 *
 * @param smeartype An integer specifying the smearing type:
 *        - 1: Fermi-Dirac smearing
 *        - 2: Gaussian smearing
 *        - 3: Hermite smearing
 *        - 4: Marzari-Vanderbilt smearing
 *        - 5: Methfessel-Paxton smearing (first order)
 *        - 6: Cold smearing
 *        - 7: Lorentzian smearing
 *        - Other: Default (no correction)
 * @param smearkT The smearing parameter, often proportional to temperature,
 *                defining the width of the smearing function.
 * @param smearfermi The Fermi level energy around which the correction is applied.
 * @param occ The occupancy probability for the energy level.
 * @param eig The eigenvalue or energy level under consideration.
 *
 * @return A double representing the correction term for the specified smearing type.
 *
 * @details
 * - The function handles a variety of smearing algorithms, adjusting for numerical precision:
 *   - **Fermi-Dirac (1):** Standard thermodynamic correction term for partial occupancies.
 *   - **Gaussian (2):** Uses a Gaussian function for energy-level broadening.
 *   - **Hermite (3):** Third-order Hermite polynomial correction for improved accuracy.
 *   - **Marzari-Vanderbilt (4):** A hybrid exponential and error function approach.
 *   - **Methfessel-Paxton (5):** First-order correction based on Hermite polynomials.
 *   - **Cold Smearing (6):** Smooth correction avoiding Gibbs oscillations.
 *   - **Lorentzian (7):** Broadening correction based on a Lorentzian function.
 *   - **Default:** Returns zero correction for step-function smearing.
 * - Numerical safeguards are included to handle edge cases, such as small or large energy
 *   differences and extreme occupancy probabilities.
 *
 * @note Ensure that `smearkT` and `smearfermi` are appropriately set based on the
 *       physical system under study. For Lorentzian smearing, `smearkT` serves as the
 *       effective smearing width.
 *
 * Example Usage:
 * @code
 * double correction = util_smearcorrection(1, 0.01, 0.5, 0.8, 0.3); // Fermi-Dirac
 * std::cout << "Correction: " << correction << std::endl;
 * @endcode
 */
double util_smearcorrection(const int smeartype, const double smearkT, const double smearfermi, 
                            const double occ, const double eig)
{
   double smearcorrection = 0.0;
   const double sqrt_pi = sqrt(M_PI);
   const double sqrt_half = sqrt(0.5);
   const double sqrt_two  = sqrt(2.0);

   // Energy difference normalized by smearing temperature
   double x = (eig - smearfermi) / smearkT;
   double y = occ; // Occupancy

   // Calculate corrections based on smearing type
   if (smeartype == 1) { // Fermi-Dirac correction
       if ((y > 1.0e-6) && ((1.0-y) > 1.0e-6)) {
           smearcorrection += smearkT * (y*log(y) + (1.0-y)*log(1.0-y));
       }
   } else if (smeartype == 2) { // Gaussian correction
       smearcorrection -= smearkT * std::exp(-x * x) / (4.0 * sqrt_pi);
   } else if (smeartype == 3) { // Hermite smearing
       double exp_term = std::exp(-x * x);
       double hermite_correction = (2.0 * x * x - 1.0) * exp_term;
       double f_original = exp_term + hermite_correction / sqrt_pi;
       double fmax = std::exp(-0.0) + (2.0 * 0.0 * 0.0 - 1.0) * std::exp(-0.0) / sqrt_pi;
     
       if (x >0.0) { // Smearing correction - No correction needed for x <= 0, as f(x) = 1.0
           // Scale the correction for x > 0
           smearcorrection -= smearkT * f_original / fmax;
       }
   } else if (smeartype == 4) { // Marzari-Vanderbilt correction
       smearcorrection -= smearkT * exp(-(x +sqrt_half)*(x+sqrt_half))*(1.0+sqrt_two*x) / (2.0*sqrt_pi);
   } else if (smeartype == 5) { // Methfessel-Paxton
       double exp_term = std::exp(-x * x);
       double hermite_poly = 1.0 - x * x;
       smearcorrection -= smearkT * exp_term * hermite_poly / sqrt_pi;
   } else if (smeartype == 6) { // Cold smearing correction (example for smearing type 5)
       double exp_term = std::exp(-0.5 * x * x);
       smearcorrection -= smearkT * exp_term * x / sqrt_two;
   } else if (smeartype == 7) { // Lorentzian Smearing correction
       double sigma = smearkT;  // Smearing width (can be parameterized further)
       smearcorrection -= smearkT * (sigma / M_PI) / (x * x + sigma * sigma);
   } else { // Default: Step function (no correction needed)
      smearcorrection = 0.0;
   }
   return smearcorrection;
}



/**************************************
 *                                    *
 *    util_kiril_coulomb_transform    *
 *                                    *
 **************************************/
/*     This function returns the fourier transform of
 *
 *             if flag==1  v_kiril = exp(-(r/rcut)**pp)/r
 *     or      if flag==2  v_kiril = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)/r
 *
 *     Entry - gg: g squared
 *             rcut:
 *             pp:
 *
 *     Exit - returns
 *                              /infty
 *                             |
 *      v_kiril(g) = (4*pi)  * | r**2 * v_kiril(r)* j0(gr) dr
 *                             |
 *                            / 0
*/
/**
 * \brief Computes the Fourier transform of the Kiril Coulomb potential with Bessel function.
 *
 * This function returns the Fourier transform of:
 * - \( v_{kiril} = \frac{e^{-\left(\frac{r}{rcut}\right)^{pp}}}{r} \) if flag == 1
 * - \( v_{kiril} = \frac{1.0 - \left(1.0 - e^{-\left(\frac{r}{rcut}\right)^{pp2}}\right)^{pp}}{r} \) if flag == 2
 *
 * The function performs a numerical integration using the trapezoidal rule and includes the Bessel function of the first kind.
 *
 * \param flag Determines the form of the potential.
 * \param gg   The squared magnitude of the wave vector.
 * \param rcut A parameter for the potential.
 * \param pp   Another parameter for the potential.
 * \return     The Fourier transform value at \( g \).
 */
double util_kiril_coulomb_transform(const int flag, const double gg, const double rcut, const double pp) 
{
   //int nrho = 15000;
   int nrho = 100;
   double pp2    = pp + 2.0;
   double drho   = 2.0*rcut/((double) nrho);
   double q      = std::sqrt(gg);
   double fourpi = 16.0 * std::atan(1.0);
      
   double sum = 0.0;
   double r   = 0.0;

   if (flag == 1) 
   {
       for (auto i=2; i<nrho; ++i)
       {
         r += drho;
         sum += std::sin(q*r)* std::exp(-std::pow(r/rcut,pp));
       }
       r = drho*(nrho-1);
       sum += sum + 0.5*std::sin(q*r)* std::exp(-std::pow(r/rcut,pp));
   }
   else 
   {
      for (auto i=2; i<nrho; ++i)
      {
         r += drho;
          sum += std::sin(q * r) * (1.0 - std::pow(1.0 - std::exp(-std::pow(r / rcut, pp2)), pp));
      }
      r = drho*(nrho-1);
      sum += 0.5*std::sin(q * r) * (1.0 - std::pow(1.0 - std::exp(-std::pow(r / rcut, pp2)), pp));
   }
   return (fourpi/q)*sum*drho;
}


/**************************************
 *                                    *
 *   util_kiril_coulomb_transform0    *
 *                                    *
 **************************************/
/*     This function returns the fourier transform of 
*
*           if flag==1   v_kiril = exp(-(r/rcut)**pp)/r
*     or    if flag==2   v_kiril = (1.0d0-(1.0d0-dexp(-(r/rcut)**pp2))**pp)/r
*
*     Entry - 
*             rcut: 
*             pp:
*                                           
*     Exit - returns 
*                              /infty
*                             | 
*      v_kiril(g=0) = (4*pi)* | r**2 * v_kiril(r) dr
*                             |
*                            / 0
*/
/**
 * \brief Computes the Fourier transform of the Kiril Coulomb potential.
 *
 * This function returns the Fourier transform of:
 * - \( v_{kiril} = \frac{e^{-\left(\frac{r}{rcut}\right)^{pp}}}{r} \) if flag == 1
 * - \( v_{kiril} = \frac{1.0 - \left(1.0 - e^{-\left(\frac{r}{rcut}\right)^{pp2}}\right)^{pp}}{r} \) if flag == 2
 *
 * The function performs a numerical integration using the trapezoidal rule.
 *
 * \param flag Determines the form of the potential.
 * \param rcut A parameter for the potential.
 * \param pp   Another parameter for the potential.
 * \return     The Fourier transform value at \( g = 0 \).
 */
double util_kiril_coulomb_transform0(const int flag, const double rcut, const double pp) 
{
   //int nrho = 15000;
   int nrho = 100;
   double pp2    = pp + 2.0;
   double drho   = 2.0*rcut/((double) nrho);
   double fourpi = 16.0 * std::atan(1.0);

   double sum = 0.0;
   double r   = 0.0;

   if (flag == 1)
   {
       for (auto i=2; i<nrho; ++i)
       {
         r += drho;
         sum += r * std::exp(-std::pow(r/rcut,pp));
       }
       r = drho*(nrho-1);
       sum += sum + 0.5 * r * std::exp(-std::pow(r/rcut,pp));
   }
   else 
   {
      for (auto i=2; i<nrho; ++i)
      {
         r += drho;
          sum += r * (1.0 - std::pow(1.0 - std::exp(-std::pow(r / rcut, pp2)), pp));
      }
      r = drho*(nrho-1);
      sum += 0.5 * r * (1.0 - std::pow(1.0 - std::exp(-std::pow(r / rcut, pp2)), pp));
   }

   return fourpi*sum*drho;
}


/**************************************
 *                                    *
 *       util_print_elapsed_time      *
 *                                    *
 **************************************/
/**
 * \brief  Prints the elapsed time of a simulation in different units.
 *
 * This function prints the elapsed time of a simulation in different units 
 * (femtoseconds, picoseconds, or nanoseconds) based on its magnitude.
 *
 * \param  autime The elapsed time in atomic units.
 */
void util_print_elapsed_time(const double autime) {
  double sectime = autime * 2.41889e-17;

  std::cout << std::endl << std::endl;
  if (sectime < 1.0e-12)
    std::cout << " Elapsed time of simulation was" << std::right << std::fixed
              << std::setw(8) << std::setprecision(3) << (sectime / 1.0e-15)
              << " fs" << std::endl;
  else if (sectime < 1.0e-9)
    std::cout << " Elapsed time of simulation was" << std::right << std::fixed
              << std::setw(8) << std::setprecision(3) << (sectime / 1.0e-12)
              << " ps" << std::endl;
  else 
    std::cout << " Elapsed time of simulation was" << std::right << std::fixed
              << std::setw(8) << std::setprecision(3) << (sectime / 1.0e-9)
              << " ns" << std::endl;
}


#ifdef _MATPRINT_
void util_matprint(std::string matlabel, int n, double *A) {
  std::cout << "util_matprint: " << matlabel << std::endl;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j)
      std::cout << A[i + j * n] << " ";
    std::cout << std::endl;
  }
  std::cout << std::endl;
}
#endif

} // namespace pwdft

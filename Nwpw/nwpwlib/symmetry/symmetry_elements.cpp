
#include <iostream>
#include <string>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <cstring> 
#include <map>
#include <array>
#include "blas.h"


#define mindex(i,j) (3*i+j)
#define m2index(i,j) (3*i+j)
#define m3index(i,j,k) (9*i+3*j+k)

namespace pwdft {

enum MirrorPlaneType {
    SIGMA_H,  // Horizontal mirror plane
    SIGMA_V,  // Vertical mirror plane
    SIGMA_D   // Diagonal mirror plane
};

// Function to calculate the GCD of two integers using Euclidean algorithm
int gcd(int a, int b) {
    while (b != 0) {
        int temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}


/*******************************************
 *                                         *
 *         calculate_center                *
 *                                         *
 *******************************************/
/**
 * @brief Calculate the geometric center of a set of 3D points.
 *
 * This function computes the geometric center (centroid) of a collection of 3D points.
 * The center is calculated as the arithmetic mean of the coordinates of all the points.
 *
 * @param points A vector containing 3D points as arrays of doubles, where each point is
 *               represented as {x, y, z}.
 *
 * @return An array containing the x, y, and z coordinates of the calculated center.
 */
static std::array<double, 3> calculate_center(const std::vector<std::array<double, 3>> &points) {
    std::array<double, 3> center = {0, 0, 0};
    for (const auto &point : points) {
        center[0] += point[0];
        center[1] += point[1];
        center[2] += point[2];
    }
    center[0] /= points.size();
    center[1] /= points.size();
    center[2] /= points.size();
    return center;
} 

/*******************************************
 *                                         *
 *         sort_point_spherical            *
 *                                         *
 *******************************************/
/**
 * @brief Sort a vector of 3D points based on their spherical coordinates.
 *
 * This function sorts a vector of 3D points in increasing order of their spherical
 * coordinates with respect to a specified center point. The sorting is performed in
 * spherical coordinates (theta and phi).
 *
 * @param points A vector of 3D points represented as arrays of doubles, where each point
 *               is {x, y, z}.
 * @param center The center point with respect to which the spherical coordinates are
 *               calculated.
 */
static void sort_points_spherical(std::vector<std::array<double, 3>> &points, std::array<double, 3> center) 
{
   std::sort(points.begin(), points.end(), [center](const std::array<double, 3>& a, const std::array<double, 3>& b) -> bool {
       double ax = a[0] - center[0];
       double ay = a[1] - center[1];
       double az = a[2] - center[2];

       double bx = b[0] - center[0];
       double by = b[1] - center[1];
       double bz = b[2] - center[2];

       double theta_a = std::atan2(ay, ax);
       double phi_a = std::atan2(std::sqrt(ax * ax + ay * ay), az);

       double theta_b = std::atan2(by, bx);
       double phi_b = std::atan2(std::sqrt(bx * bx + by * by), bz);

       if (theta_a == theta_b) {
           return phi_a < phi_b;
       }
       return theta_a < theta_b;
   });
}


/*******************************************
 *                                         *
 *             is_regular_triangle         *
 *                                         *
 *******************************************/
/**
 * @brief Check if three points form an equilateral triangle in 3D space.
 *
 * This function checks whether three 3D points form an equilateral triangle within the
 * specified tolerance. An equilateral triangle has three sides of equal length.
 *
 * @param points A vector containing three 3D points represented as arrays of doubles,
 *               where each point is {x, y, z}.
 * @param tol The tolerance within which the side lengths are considered equal.
 * @return True if the points form an equilateral triangle within the specified tolerance;
 *         otherwise, false.
 */
static bool is_regular_triangle(std::vector<std::array<double, 3>> &points, const double tol)
{
   if (points.size() != 3) {
       return false;
   }

   // Calculate distances between every pair of points
   std::vector<double> distances;
   for (int i = 0; i < 3; ++i) {
       for (int j = i + 1; j < 3; ++j) {
           double dx = points[i][0] - points[j][0];
           double dy = points[i][1] - points[j][1];
           double dz = points[i][2] - points[j][2];
           double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
           distances.push_back(dist);
       }
   }

   // Check if the 3 distances are the same
   bool has_equal_distances = ((std::abs(distances[0]-distances[1]) < tol) &&
                               (std::abs(distances[1]-distances[2]) < tol) &&
                               (std::abs(distances[2]-distances[1]) < tol) );

   return has_equal_distances;
}

// overloaded function
/**
 * @brief Check if three points, given by their indices, form an equilateral triangle in 3D space.
 *
 * This overloaded function checks whether three points specified by their indices (i, j, k)
 * from an array of 3D coordinates (rion) form an equilateral triangle within the specified tolerance.
 * An equilateral triangle has equal side lengths and equal angles between all sides.
 *
 * @param i Index of the first point in the coordinate array.
 * @param j Index of the second point in the coordinate array.
 * @param k Index of the third point in the coordinate array.
 * @param rion An array containing the 3D coordinates of atoms. The coordinates are expected to be in the format [x1, y1, z1, x2, y2, z2, ...].
 * @param nion The total number of atoms in the coordinate array.
 * @param tol The tolerance within which the side lengths are considered equal.
 * @return True if the specified points form an equilateral triangle within the specified tolerance;
 *         otherwise, false.
 */
static bool is_regular_triangle(const int i, const int j, const int k,
                 const double *rion, const int nion, double tol)
{
   if (nion < 3) {
       return false;
   }

   // Extract the points based on the provided indices
   std::vector<std::array<double, 3>> points;
   points.push_back({rion[3 * i], rion[3 * i + 1], rion[3 * i + 2]});
   points.push_back({rion[3 * j], rion[3 * j + 1], rion[3 * j + 2]});
   points.push_back({rion[3 * k], rion[3 * k + 1], rion[3 * k + 2]});

   return is_regular_triangle(points, tol);
}


/*******************************************
 *                                         *
 *               is_square                 *
 *                                         *
 *******************************************/
/**
 * @brief Check if four points form a regular square in 3D space.
 *
 * This function checks whether four 3D points form a regular square within the specified tolerance.
 * A regular square has equal side lengths and right angles between adjacent sides.
 *
 * @param points A vector containing four 3D points represented as arrays of doubles,
 *               where each point is {x, y, z}.
 * @param tol The tolerance within which the side lengths are considered equal.
 * @return True if the points form a regular square within the specified tolerance;
 *         otherwise, false.
 */
static bool is_square(std::vector<std::array<double, 3>> &points, const double tol)
{
   if (points.size() != 4) {
       return false;
   }

   // Calculate distances between every pair of points
   std::vector<double> distances;
   for (int i = 0; i < 4; ++i) {
       for (int j = i + 1; j < 4; ++j) {
           double dx = points[i][0] - points[j][0];
           double dy = points[i][1] - points[j][1];
           double dz = points[i][2] - points[j][2];
           double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
           distances.push_back(dist);
       }
   }

   std::sort(distances.begin(), distances.end());

   // Check if there are 4 smaller and 2 larger distances
   double small_dist = distances[0];
   double large_dist = distances[5];
   int small_count = 0, large_count = 0;
   for (const auto& dist : distances) {
       if (std::abs(dist - small_dist) < tol) {
           ++small_count;
       }
       else if (std::abs(dist - large_dist) < tol) {
           ++large_count;
       }
   }

   return small_count == 4 && large_count == 2;
}

// overloaded function
/**
 * @brief Check if four points, given by their indices, form a regular square in 3D space.
 *
 * This overloaded function checks whether four points specified by their indices (i, j, k, l)
 * from an array of 3D coordinates (rion) form a regular square within the specified tolerance.
 * A regular square has equal side lengths and equal angles between adjacent sides.
 *
 * @param i Index of the first point in the coordinate array.
 * @param j Index of the second point in the coordinate array.
 * @param k Index of the third point in the coordinate array.
 * @param l Index of the fourth point in the coordinate array.
 * @param rion An array containing the 3D coordinates of atoms. The coordinates are expected to be in the format [x1, y1, z1, x2, y2, z2, ...].
 * @param nion The total number of atoms in the coordinate array.
 * @param tol The tolerance within which the side lengths are considered equal.
 * @return True if the specified points form a regular square within the specified tolerance;
 *         otherwise, false.
 */
static bool is_square(const int i, const int j, const int k, const int l, 
                 const double *rion, const int nion, double tol)
{
   if (nion < 4) {
       return false;
   }

   // Extract the points based on the provided indices
   std::vector<std::array<double, 3>> points;
   points.push_back({rion[3 * i], rion[3 * i + 1], rion[3 * i + 2]});
   points.push_back({rion[3 * j], rion[3 * j + 1], rion[3 * j + 2]});
   points.push_back({rion[3 * k], rion[3 * k + 1], rion[3 * k + 2]});
   points.push_back({rion[3 * l], rion[3 * l + 1], rion[3 * l + 2]});

   return is_square(points, tol);
}


/*******************************************
 *                                         *
 *             is_pentagon                 *
 *                                         *
 *******************************************/
/**
 * @brief Check if five points form a regular pentagon in 3D space.
 *
 * This function checks whether five 3D points form a regular pentagon within the specified tolerance.
 * A regular pentagon has equal side lengths and equal angles between adjacent sides.
 *
 * @param points A vector containing five 3D points represented as arrays of doubles,
 *               where each point is {x, y, z}.
 * @param tol The tolerance within which the side lengths are considered equal.
 * @return True if the points form a regular pentagon within the specified tolerance;
 *         otherwise, false.
 */
static bool is_pentagon(std::vector<std::array<double, 3>> &points, const double tol) 
{
   if (points.size() != 5) {
       return false;
   }

   // Calculate distances between every pair of points
   std::vector<double> distances;
   for (int i = 0; i < 5; ++i) {
       for (int j = i + 1; j < 5; ++j) {
           double dx = points[i][0] - points[j][0];
           double dy = points[i][1] - points[j][1];
           double dz = points[i][2] - points[j][2];
           double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
           distances.push_back(dist);
       }
   }

   std::sort(distances.begin(), distances.end());

   // Check if there are 5 smaller and 5 larger distances
   double small_dist = distances[0];
   double large_dist = distances[9];
   int small_count = 0, large_count = 0;
   for (const auto& dist : distances) {
       if (std::abs(dist - small_dist) < tol) {
           ++small_count;
       }
       else if (std::abs(dist - large_dist) < tol) {
           ++large_count;
       }
   }

   return small_count == 5 && large_count == 5;
}

// overloaded function
/**
 * @brief Check if five points, given by their indices, form a regular pentagon in 3D space.
 *
 * This overloaded function checks whether five points specified by their indices (i, j, k, l, m)
 * from an array of 3D coordinates (rion) form a regular pentagon within the specified tolerance.
 * A regular pentagon has equal side lengths and equal angles between adjacent sides.
 *
 * @param i Index of the first point in the coordinate array.
 * @param j Index of the second point in the coordinate array.
 * @param k Index of the third point in the coordinate array.
 * @param l Index of the fourth point in the coordinate array.
 * @param m Index of the fifth point in the coordinate array.
 * @param rion An array containing the 3D coordinates of atoms. The coordinates are expected to be in the format [x1, y1, z1, x2, y2, z2, ...].
 * @param nion The total number of atoms in the coordinate array.
 * @param tol The tolerance within which the side lengths are considered equal.
 * @return True if the specified points form a regular pentagon within the specified tolerance;
 *         otherwise, false.
 */
static bool is_pentagon(const int i, const int j, const int k, const int l, const int m,
                 const double *rion, const int nion, double tol) 
{
   if (nion < 5) {
       return false;
   }

   // Extract the points based on the provided indices
   std::vector<std::array<double, 3>> points;
   points.push_back({rion[3 * i], rion[3 * i + 1], rion[3 * i + 2]});
   points.push_back({rion[3 * j], rion[3 * j + 1], rion[3 * j + 2]});
   points.push_back({rion[3 * k], rion[3 * k + 1], rion[3 * k + 2]});
   points.push_back({rion[3 * l], rion[3 * l + 1], rion[3 * l + 2]});
   points.push_back({rion[3 * m], rion[3 * m + 1], rion[3 * m + 2]});

   return is_pentagon(points, tol);
}






/*******************************************
 *                                         *
 *           eigsrt_hitolow                *
 *                                         *
 *******************************************/
/* Description:
     Sorts eigenvalues and corresponding eigenvectors in descending order
     (from highest to lowest eigenvalue).

   Arguments:
     - D: Array of eigenvalues to be sorted.
     - V: Matrix of eigenvectors corresponding to the eigenvalues.
     - n: Number of eigenvalues and eigenvectors.

   Note:
     This function modifies both the D and V arrays.
*/
static void eigsrt_hitolow(double *D, double *V, int n) {
  int i, j, k;
  double p;

  for (i = 0; i < (n - 1); ++i) {
    k = i;
    p = D[i];
    for (j = i + 1; j < n; ++j)
      if (D[j] >= p) {
        k = j;
        p = D[j];
      }

    if (k != i) {
      D[k] = D[i];
      D[i] = p;
      for (j = 0; j < n; ++j) {
        p = V[j + i * n];
        V[j + i * n] = V[j + k * n];
        V[j + k * n] = p;
      }
    }
  }
}


/*******************************************
 *                                         *
 *           jacobi_eigenvalue             *
 *                                         *
 *******************************************/
/* Description:
   Compute the eigenvalues and eigenvectors of a real symmetric 3x3 matrix
   using the Jacobi rotation method with threshold pivoting.

   Arguments:
     - a: Input matrix of size 3x3 (must be real, symmetric).
     - v: Output matrix containing the eigenvectors (size 3x3).
     - d: Output array containing the eigenvalues (size 3).
     - it_num: Output variable storing the total number of iterations performed.
     - rot_num: Output variable storing the total number of rotations applied.

   Note:
   The function assumes that the input matrix 'a' is real and symmetric.

   Returns:
   The eigenvalues are stored in 'd' in descending order. The eigenvectors
   are stored as columns in the matrix 'v'.

   Example Usage:
   double a[9] = { ... }; // Input matrix
   double v[9], d[3];    // Output matrices
   int it_num, rot_num;   // Variables to store iteration and rotation counts
   jacobi_eigenvalue(a, v, d, it_num, rot_num);

   References:
   - Gene Golub, Charles Van Loan, "Matrix Computations," Third Edition,
     Johns Hopkins, 1996, ISBN: 0-8018-4513-X, LC: QA188.G65.
*/

static void jacobi_eigenvalue(double a[], double v[], double d[], int& it_num, int& rot_num) 
{
    // Initialize the eigenvectors matrix v as the identity matrix
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            v[i * 3 + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    double bw[3];
    double zw[3];

    for (int i = 0; i < 3; i++) {
        bw[i] = d[i] = a[i * 3 + i];
        zw[i] = 0.0;
    }

    it_num = 0;
    rot_num = 0;

    while (it_num < 5000) {  // You can adjust the maximum number of iterations
        it_num++;

        // Calculate the threshold for convergence
        double thresh = 0.0;
        for (int j = 0; j < 3; j++) {
            for (int i = 0; i < j; i++) {
                thresh += a[i + j * 3] * a[i + j * 3];
            }
        }
        thresh = std::sqrt(thresh) / 12.0;

        if (thresh == 0.0) {
            break;
        }

        for (int p = 0; p < 3; p++) {
            for (int q = p + 1; q < 3; q++) {
                double gapq = 10.0 * std::abs(a[p + q * 3]);
                double termp = gapq + std::abs(d[p]);
                double termq = gapq + std::abs(d[q]);

                // Annihilate tiny off-diagonal elements
                if (it_num > 4 && termp == std::abs(d[p]) && termq == std::abs(d[q])) {
                    a[p + q * 3] = 0.0;
                }
                // Apply a rotation if the threshold is met
                else if (std::abs(a[p + q * 3]) >= thresh) {
                    double h = d[q] - d[p];
                    double term = std::abs(h) + gapq;

                    double t;
                    if (term == std::abs(h)) {
                        t = a[p + q * 3] / h;
                    } else {
                        double theta = 0.5 * h / a[p + q * 3];
                        t = 1.0 / (std::abs(theta) + std::sqrt(1.0 + theta * theta));
                        if (theta < 0.0) {
                            t = -t;
                        }
                    }

                    double c = 1.0 / std::sqrt(1.0 + t * t);
                    double s = t * c;
                    double tau = s / (1.0 + c);

                    h = t * a[p + q * 3];
                    zw[p] -= h;
                    zw[q] += h;
                    d[p] -= h;
                    d[q] += h;
                    a[p + q * 3] = 0.0;

                    for (int j = 0; j < p; j++) {
                        double g = a[j + p * 3];
                        h = a[j + q * 3];
                        a[j + p * 3] = g - s * (h + g * tau);
                        a[j + q * 3] = h + s * (g - h * tau);
                    }

                    for (int j = p + 1; j < q; j++) {
                        double g = a[p + j * 3];
                        h = a[j + q * 3];
                        a[p + j * 3] = g - s * (h + g * tau);
                        a[j + q * 3] = h + s * (g - h * tau);
                    }

                    for (int j = q + 1; j < 3; j++) {
                        double g = a[p + j * 3];
                        h = a[q + j * 3];
                        a[p + j * 3] = g - s * (h + g * tau);
                        a[q + j * 3] = h + s * (g - h * tau);
                    }

                    for (int j = 0; j < 3; j++) {
                        double g = v[j + p * 3];
                        h = v[j + q * 3];
                        v[j + p * 3] = g - s * (h + g * tau);
                        v[j + q * 3] = h + s * (g - h * tau);
                    }

                    rot_num++;
                }
            }
        }

        for (int i = 0; i < 3; i++) {
            bw[i] += zw[i];
            d[i] = bw[i];
            zw[i] = 0.0;
        }
    }
}


/*******************************************
 *                                         *
 *             cos_dot_product             *
 *                                         *
 *******************************************/
/**
 * @brief Calculate the cosine of the angle between two vectors represented by their coordinates.
 *
 * This function calculates the cosine of the angle between two vectors represented by their 3D coordinates.
 * It computes the dot product of the two vectors and divides it by the product of their magnitudes
 * to find the cosine of the angle.
 *
 * @param r1 The 3D coordinates of the first vector in the format [x1, y1, z1].
 * @param r2 The 3D coordinates of the second vector in the format [x2, y2, z2].
 * @return The cosine of the angle between the two vectors.
 */
static double cos_dot_product(const double *r1, const double *r2)
{
   return (r1[0]*r2[0]+r1[1]*r2[1]+r1[2]*r2[2])/std::sqrt((r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2])*(r2[0]*r2[0]+r2[1]*r2[1]+r2[2]*r2[2]));
}

/*******************************************
 *                                         *
 *         normalized_cross_product        *
 *                                         *
 *******************************************/
/**
 * @brief Calculate the normalized cross product of two vectors.
 *
 * This function computes the cross product of two 3D vectors represented by their coordinates,
 * normalizes the resulting vector to have unit magnitude, and stores it in the provided array `n`.
 * If the cross product results in a zero vector, a warning is issued, and the `n` vector will also be zero.
 *
 * @param r1 The 3D coordinates of the first vector in the format [x1, y1, z1].
 * @param r2 The 3D coordinates of the second vector in the format [x2, y2, z2].
 * @param n  A pointer to an array where the normalized cross product will be stored in the format [nx, ny, nz].
 */
static void normalized_cross_product(const double *r1, const double *r2, double *n)
{
    // Calculate the cross product
    n[0] = r1[1]*r2[2] - r2[1]*r1[2];
    n[1] = r1[2]*r2[0] - r2[2]*r1[0];
    n[2] = r1[0]*r2[1] - r2[0]*r1[1];

    // Calculate the magnitude of the resulting vector
    double magnitude = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);

    // Check for zero magnitude to prevent division by zero when normalizing
    if (magnitude == 0) 
    {
       std::cerr << "Warning: Cannot normalize a zero vector. Returning zero vector." << std::endl;
       n[0]=0.0; n[1]=0.0; n[2]=0.0;
       return;
    }

    // Normalize the vector
    n[0] /= magnitude;
    n[1] /= magnitude;
    n[2] /= magnitude;
}


/*******************************************
 *                                         *
 *           shift_to_center_mass          *
 *                                         *
 *******************************************/
/**
 * @brief Shifts ion coordinates so that they are centered around the system's center of mass.
 * 
 * This function takes the 3D coordinates of ions and shifts them such that their center
 * of mass is at the origin (0,0,0). This is useful for various calculations including
 * the computation of multipole moments like the quadrupole tensor.
 * 
 * @param rion           Pointer to the array containing the original ion coordinates.
 *                       The array is structured as [x1, y1, z1, x2, y2, z2, ..., xn, yn, zn].
 * @param ion_mass       Pointer to the array containing the masses of the ions.
 *                       The array is structured as [m1, m2, ..., mn].
 * @param nion           The number of ions in the system.
 * @param rion2          Output array where the shifted ion coordinates will be stored.
 *                       The array is structured similarly to `rion`.
 * 
 * Preconditions:
 * - `rion` and `ion_mass` must have at least `nion` elements.
 * - `rion2` must be pre-allocated with at least 3 * `nion` elements.
 * 
 * Postconditions:
 * - `rion2` will contain the coordinates of ions shifted such that the center of mass
 *   of the system is at the origin.
 */
static void shift_to_center_mass(const double* rion, const double *ion_mass, const int nion, 
                                 double *rion2)
{
   double mtotal = 0.0;
   double center_mass_x = 0.0;
   double center_mass_y = 0.0;
   double center_mass_z = 0.0;

   for (int ii=0; ii<nion; ++ii)
   {
      double mass = ion_mass[ii];
      double x = rion[3*ii];
      double y = rion[3*ii+1];
      double z = rion[3*ii+2];
      mtotal += mass;
      center_mass_x += mass*x;
      center_mass_y += mass*y;
      center_mass_z += mass*z;
   }
   center_mass_x /= mtotal;
   center_mass_y /= mtotal;
   center_mass_z /= mtotal;

   for (int ii=0; ii<nion; ++ii)
   {
       rion2[3*ii]   = rion[3*ii]   - center_mass_x;
       rion2[3*ii+1] = rion[3*ii+1] - center_mass_y;
       rion2[3*ii+2] = rion[3*ii+2] - center_mass_z;
   }
}

/*******************************************
 *                                         *
 *           align_to_axes                 *
 *                                         *
 *******************************************/
/**
 * @brief Rotates the atomic coordinates to align them with the given axes.
 *
 * This function takes an array of atomic coordinates `rion` for `nion` atoms and
 * aligns them along the new axes defined by the `axes` matrix.
 *
 * @param rion Pointer to the array holding the atomic coordinates. The array is modified in place.
 * @param nion The number of atoms.
 * @param axes The 3x3 rotation matrix defining the new axes.
 */
static void align_to_axes(double* rion, const int nion, const double axes[9])
{
   // Loop over all atoms to rotate coordinates
   for (int ii=0; ii<nion; ++ii)
   {
       // Extract original coordinates
       double xold = rion[3*ii];
       double yold = rion[3*ii+1];
       double zold = rion[3*ii+2];

       // Perform matrix multiplication to find new coordinates
       double xnew = axes[0]*xold + axes[1]*yold + axes[2]*zold;
       double ynew = axes[3]*xold + axes[4]*yold + axes[5]*zold;
       double znew = axes[6]*xold + axes[7]*yold + axes[8]*zold;

       // Store new coordinates
       rion[3*ii]   = xnew;
       rion[3*ii+1] = ynew;
       rion[3*ii+2] = znew;
   }
}




/*******************************************
 *                                         *
 *           generate_quadrupole           *
 *                                         *
 *******************************************/
/**
 * @brief Calculate the quadrupole moments of a distribution of point masses in 3D space.
 *
 * @param rion Pointer to the array containing the positions of the ions (atoms).
 * @param ion_mass Pointer to the array containing the masses of the ions (atoms).
 * @param nion The number of ions (atoms) in the system.
 * @param Qxx,Qyy,Qzz,Qxy,Qxz,Qyz Output variables for each component of the quadrupole tensor.
 *
 * @note
 * The function calculates the center of mass and then uses it to compute the quadrupole moments.
 *
 * @code
 * Example Usage:
 * double rion[3 * nion];      // Point mass positions
 * double ion_mass[nion];      // Masses of point masses
 * int nion;                   // Number of point masses
 * double Qxx, Qyy, Qzz, Qxy, Qxz, Qyz;  // Quadrupole moments
 * generate_quadrupole(rion, ion_mass, nion, Qxx, Qyy, Qzz, Qxy, Qxz, Qyz);
 * @endcode
 */
static void generate_quadrupole(const double* rion, const double *ion_mass, const int nion, 
                                double& Qxx, double& Qyy, double& Qzz,
                                double& Qxy, double& Qxz, double& Qyz)
{
   double mtotal = 0.0;
   double center_mass_x = 0.0;
   double center_mass_y = 0.0;
   double center_mass_z = 0.0;
   
   for (auto ii=0; ii<nion; ++ii)
   {
      auto mass = ion_mass[ii];
      auto x = rion[3*ii];
      auto y = rion[3*ii+1];
      auto z = rion[3*ii+2];
      mtotal += mass;
      center_mass_x += mass*x;
      center_mass_y += mass*y;
      center_mass_z += mass*z;
   }
   center_mass_x /= mtotal;
   center_mass_y /= mtotal;
   center_mass_z /= mtotal;

   Qxx = 0.0; Qyy = 0.0; Qzz = 0.0; 
   Qxy = 0.0; Qxz = 0.0; Qyz = 0.0; 
   for (auto ii=0; ii<nion; ++ii)
   {
      auto mass = ion_mass[ii];
      auto x = rion[3*ii]   - center_mass_x;
      auto y = rion[3*ii+1] - center_mass_y;
      auto z = rion[3*ii+2] - center_mass_z;
      auto r2 = x*x + y*y + z*z;
      Qxx += mass*(3*x*x - r2);
      Qyy += mass*(3*y*y - r2);
      Qzz += mass*(3*z*z - r2);
      Qxy += mass*(3*x*y);
      Qxz += mass*(3*x*z);
      Qyz += mass*(3*y*z);
   }
}

/*******************************************
 *                                         *
 *       generate_quadrupole_tensor        *
 *                                         *
 *******************************************/
/**
 * @brief Generates the quadrupole tensor for a system of ions in 3D space.
 * 
 * This function calculates the 3x3 quadrupole tensor based on the positions and masses
 * of ions. The tensor is calculated with respect to the center of mass of the ion system.
 * The tensor is stored as a flattened 1D array with 9 elements, mapped from the 3x3 tensor.
 *
 * @param rion           Pointer to the array containing the ion coordinates.
 *                       The array is structured as [x1, y1, z1, x2, y2, z2, ..., xn, yn, zn].
 * @param ion_mass       Pointer to the array containing the masses of the ions.
 *                       The array is structured as [m1, m2, ..., mn].
 * @param nion           The number of ions in the system.
 * @param quadrupole_tensor Output array where the quadrupole tensor will be stored.
 *                          The tensor is stored in row-major order.
 * 
 * Preconditions:
 * - `rion` and `ion_mass` must have at least `nion` elements.
 * - `quadrupole_tensor` must be pre-allocated with at least 9 elements.
 * 
 * Postconditions:
 * - `quadrupole_tensor` will contain the calculated tensor elements.
 */
static void generate_quadrupole_tensor(const double* rion, const double *ion_mass, const int nion, double quadrupole_tensor[9])
{
   // Initialize variables for center of mass and total mass
   double mtotal = 0.0;
   double center_mass_x = 0.0;
   double center_mass_y = 0.0;
   double center_mass_z = 0.0;

    // Calculate center of mass
   for (int ii=0; ii<nion; ++ii)
   {
      double mass = ion_mass[ii];
      double x = rion[3*ii];
      double y = rion[3*ii+1];
      double z = rion[3*ii+2];
      mtotal += mass;
      center_mass_x += mass*x;
      center_mass_y += mass*y;
      center_mass_z += mass*z;
   }
   center_mass_x /= mtotal;
   center_mass_y /= mtotal;
   center_mass_z /= mtotal;

   // Initialize quadrupole tensor to zeros
   std::fill(quadrupole_tensor, quadrupole_tensor + 9, 0.0);

   // Calculate quadrupole tensor
   for (int ii=0; ii<nion; ++ii)
   {
      double mass = ion_mass[ii];
      double x[3] = {rion[3*ii]-center_mass_x,rion[3*ii+1]-center_mass_y,rion[3*ii+2]-center_mass_z};
      double r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];

      // Populate quadrupole tensor
      for (int j=0; j<3; ++j)
         for (int i=0; i<3; ++i)
            quadrupole_tensor[m2index(i,j)] += mass*(3*x[j]*x[i] - r2*(i==j));
   }
}


/*******************************************
 *                                         *
 *           generate_octapole             *
 *                                         *
 *******************************************/
/**
 * @brief Calculates the components of the octapole tensor for a system of point masses.
 *
 * This function computes the octapole tensor using a prefactor of 5 instead of the usual 3.
 * It also adjusts the positions of the atoms relative to the center of mass of the system.
 * The octapole tensor components are returned via output variables.
 *
 * @param rion Pointer to the array containing the positions of the ions (atoms).
 * @param ion_mass Pointer to the array containing the masses of the ions (atoms).
 * @param nion The number of ions (atoms) in the system.
 * @param Oxxx to Ozzz Output variables for each component of the octapole tensor.
 */
static void generate_octapole(const double* rion, const double *ion_mass, const int nion,
                              double& Oxxx, double& Oxxy, double& Oxxz,
                              double& Oxyx, double& Oxyy, double& Oxyz,
                              double& Oxzx, double& Oxzy, double& Oxzz,
                              double& Oyxx, double& Oyxy, double& Oyxz,
                              double& Oyyx, double& Oyyy, double& Oyyz,
                              double& Oyzx, double& Oyzy, double& Oyzz,
                              double& Ozxx, double& Ozxy, double& Ozxz,
                              double& Ozyx, double& Ozyy, double& Ozyz,
                              double& Ozzx, double& Ozzy, double& Ozzz)

{
   double mtotal = 0.0;
   double center_mass_x = 0.0;
   double center_mass_y = 0.0;
   double center_mass_z = 0.0;
   
   for (auto ii=0; ii<nion; ++ii)
   {
      auto mass = ion_mass[ii];
      auto x = rion[3*ii];
      auto y = rion[3*ii+1];
      auto z = rion[3*ii+2];
      mtotal += mass;
      center_mass_x += mass*x;
      center_mass_y += mass*y;
      center_mass_z += mass*z;
   }
   center_mass_x /= mtotal;
   center_mass_y /= mtotal;
   center_mass_z /= mtotal;

   Oxxx = 0.0; Oxxy = 0.0; Oxxz = 0.0;
   Oxyx = 0.0; Oxyy = 0.0; Oxyz = 0.0;
   Oxzx = 0.0; Oxzy = 0.0; Oxzz = 0.0;

   Oyxx = 0.0; Oyxy = 0.0; Oyxz = 0.0;
   Oyyx = 0.0; Oyyy = 0.0; Oyyz = 0.0;
   Oyzx = 0.0; Oyzy = 0.0; Oyzz = 0.0;

   Ozxx = 0.0; Ozxy = 0.0; Ozxz = 0.0;
   Ozyx = 0.0; Ozyy = 0.0; Ozyz = 0.0;
   Ozzx = 0.0; Ozzy = 0.0; Ozzz = 0.0;
   for (auto ii=0; ii<nion; ++ii)
   {
      auto mass = ion_mass[ii];
      auto x = rion[3*ii]   - center_mass_x;
      auto y = rion[3*ii+1] - center_mass_y;
      auto z = rion[3*ii+2] - center_mass_z;
      auto r2 = x*x + y*y + z*z;

      // Oxxx...Oxzz
      Oxxx += mass*5*x*x*x - r2*x;
      Oxxy += mass*5*x*x*y - r2*x;
      Oxxz += mass*5*x*x*z - r2*x;
      Oxyx += mass*5*x*y*x - r2*x;
      Oxyy += mass*5*x*y*y - r2*x;
      Oxyz += mass*5*x*y*z - r2*x;
      Oxzx += mass*5*x*z*x - r2*x;
      Oxzy += mass*5*x*z*y - r2*x;
      Oxzz += mass*5*x*z*z - r2*x;

      // Oyxx...Oyzz
      Oyxx += mass*5*y*x*x - r2*y;
      Oyxy += mass*5*y*x*y - r2*y;
      Oyxz += mass*5*y*x*z - r2*y;
      Oyyx += mass*5*y*y*x - r2*y;
      Oyyy += mass*5*y*y*y - r2*y;
      Oyyz += mass*5*y*y*z - r2*y;
      Oyzx += mass*5*y*z*x - r2*y;
      Oyzy += mass*5*y*z*y - r2*y;
      Oyzz += mass*5*y*z*z - r2*y;

      // Ozxx...Ozzz
      Ozxx += mass*5*z*x*x - r2*z;
      Ozxy += mass*5*z*x*y - r2*z;
      Ozxz += mass*5*z*x*z - r2*z;
      Ozyx += mass*5*z*y*x - r2*z;
      Ozyy += mass*5*z*y*y - r2*z;
      Ozyz += mass*5*z*y*z - r2*z;
      Ozzx += mass*5*z*z*x - r2*z;
      Ozzy += mass*5*z*z*y - r2*z;
      Ozzz += mass*5*z*z*z - r2*z;
   }
}

/*******************************************
 *                                         *
 *         generate_octapole_tensor        *
 *                                         *
 *******************************************/
/**
 * @brief Calculates the components of the octapole tensor for a system of point masses.
 *
 * This function computes the octapole tensor using a prefactor of 5 instead of the usual 3.
 * It also adjusts the positions of the atoms relative to the center of mass of the system.
 * The octapole tensor components are returned in the output array.
 *
 * @param rion Pointer to the array containing the positions of the ions (atoms).
 * @param ion_mass Pointer to the array containing the masses of the ions (atoms).
 * @param nion The number of ions (atoms) in the system.
 * @param octapole_tensor Flat array for storing the components of the octapole tensor.
 */
static void generate_octapole_tensor(const double* rion, const double *ion_mass, const int nion, double octapole_tensor[27])
{
   // Your code for calculating the center of mass here
   double mtotal = 0.0;
   double center_mass_x = 0.0;
   double center_mass_y = 0.0;
   double center_mass_z = 0.0;

   for (int ii=0; ii<nion; ++ii)
   {
      double mass = ion_mass[ii];
      double x = rion[3*ii];
      double y = rion[3*ii+1];
      double z = rion[3*ii+2];
      mtotal += mass;
      center_mass_x += mass*x;
      center_mass_y += mass*y;
      center_mass_z += mass*z;
   }
   center_mass_x /= mtotal;
   center_mass_y /= mtotal;
   center_mass_z /= mtotal;

   // Initialize the octapole tensor
   std::fill(octapole_tensor, octapole_tensor + 27, 0.0);

   // Calculate the octapole tensor components
   for (int ii=0; ii<nion; ++ii)
   {
      double mass = ion_mass[ii];
      double x[3] = {rion[3*ii]-center_mass_x, rion[3*ii+1]-center_mass_y, rion[3*ii+2]-center_mass_z};
      double r2 = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];

      // Loop over 3x3x3 octapole tensor components
      // Oxxx...Ozzz
      for (int i=0; i<3; ++i)
         for (int j=0; j<3; ++j)
            for (int k=0; k<3; ++k)
               octapole_tensor[m3index(i,j,k)] += mass*5*x[i]*x[j]*x[k] - r2*x[i];
   }
}




/*******************************************
 *                                         *
 *      eccentricity_and_asphericity       *
 *                                         *
 *******************************************/
/* Description:
   Calculate the eccentricity and asphericity parameters based on the quadrupole moments
   of a distribution of point masses in 3D space.

   Arguments:
     - rion: Array of size 3 * nion containing the positions (x, y, z) of nion point masses.
     - ion_mass: Array of size nion containing the masses of the nion point masses.
     - nion: The number of point masses.

   Output:
     - eccentricity: Variable to store the calculated eccentricity parameter.
     - asphericity: Variable to store the calculated asphericity parameter.

   Note:
   The function calculates the eccentricity and asphericity based on the quadrupole moments
   obtained using the generate_quadrupole function.

   Example Usage:
   double rion[3 * nion];      // Point mass positions
   double ion_mass[nion];      // Masses of point masses
   int nion;                   // Number of point masses
   double eccentricity, asphericity;  // Eccentricity and asphericity parameters
   eccentricity_and_asphericity(rion, ion_mass, nion, eccentricity, asphericity);
*/
static void eccentricity_and_asphericity(const double* rion, const double *ion_mass, const int nion, 
                                         double& eccentricity, double& asphericity)
{
   double Qxx,Qyy,Qzz,Qxy,Qxz,Qyz;
   generate_quadrupole(rion,ion_mass,nion,Qxx,Qyy,Qzz,Qxy,Qxz,Qyz);

   eccentricity = (Qxx-Qyy)/Qzz;

   asphericity  = 1.5*pow((Qxx-Qyy)/Qzz,2) 
                + 0.5*pow(2*Qxy/Qzz,2) 
                + 0.5*pow(2*Qxz/Qyy,2) 
                + 0.5*pow(2*Qyz/Qxx,2);
}

/*******************************************
 *                                         *
 *       is_spherical_multipoles2          *
 *                                         *
 *******************************************/

static bool is_spherical_multipoles2(const double* rion, const double *ion_mass, const int nion, const double sym_tolerance)
{
    double pole_mat[nion*nion];
    double pole_eigs[nion];

    double mtotal   = 0.0;
    double center[3] = {0.0,0.0,0.0};

    for (auto ii=0; ii<nion; ++ii)
    {
        auto mass = ion_mass[ii];
        auto x = rion[3*ii]; 
        auto y = rion[3*ii+1]; 
        auto z = rion[3*ii+2];
        center[0] += mass*x; center[1] += mass*y; center[2] += mass*z;
        mtotal += mass;
    }
    center[0] /= mtotal;
    center[1] /= mtotal;
    center[2] /= mtotal;

    for (auto ii=0; ii<nion; ++ii)
    {
       auto massii = ion_mass[ii];
       auto xii = rion[3*ii]   - center[0]; 
       auto yii = rion[3*ii+1] - center[1]; 
       auto zii = rion[3*ii+2] - center[2];
       pole_mat[ii*nion+ii] = (xii*xii + yii*yii + zii*zii)*massii;
    }

    for (auto ii=0; ii<nion; ++ii)
    {
       auto massii = ion_mass[ii];
       auto xii = rion[3*ii]   - center[0]; 
       auto yii = rion[3*ii+1] - center[1]; 
       auto zii = rion[3*ii+2] - center[2];
       for (auto jj=ii+1; jj<nion; ++jj)
       {
           auto massjj = ion_mass[jj];
           auto xjj = rion[3*jj]   - center[0]; 
           auto yjj = rion[3*jj+1] - center[1]; 
           auto zjj = rion[3*jj+2] - center[2];
           pole_mat[jj*nion+ii] = pole_mat[ii*nion+jj] 
                                = (xii*xjj + yii*yjj + zii*zjj) *std::sqrt(massii*massjj);
       }
    }

   int ierr;
   int n = nion;
   int nn = 3*nion+3;
   double xtmp1[nn];

   EIGEN_PWDFT(n,pole_mat,pole_eigs,xtmp1,nn,ierr);

   //sort the eigevalues and eigenvectors
   //eigsrt_hitolow(pole_eigs, pole_mat, nion);

   //std::vector<double> multipoles(5);
   //multipoles[0] = pole_eigs[0]; // Set multipole(0)
   //multipoles[1] = pole_eigs[1]; // Set multipole(1)
   //multipoles[2] = pole_eigs[2]; // Set multipole(-1)
   //multipoles[3] = pole_eigs[3]; // Set multipole(2)
   //multipoles[4] = pole_eigs[4]; // Set multipole(-2)

   //   # 5 degeneracy of I and Ih group
   //     if (moment2[1] - moment2[-1]) / m_tot < global_var.tol or (moment2[0] - moment2[-2]) / m_tot < global_var.tol:

   //return( ((pole_eigs[1]-pole_eigs[2])/mtotal < sym_tolerance) || 
   //        ((pole_eigs[0]-pole_eigs[4])/mtotal < sym_tolerance) );

   return( ((pole_eigs[3]-pole_eigs[1])/mtotal < sym_tolerance) || 
           ((pole_eigs[2]-pole_eigs[0])/mtotal < sym_tolerance) );
}



/*******************************************
 *                                         *
 *           power_commplex                *
 *                                         *
 *******************************************/
 /* Description:
     Calculate the power of a complex number raised to an integer exponent,
     with consideration for large values.

   Arguments:
     - c: Real part of the complex number.
     - s: Imaginary part of the complex number.
     - n: Exponent to raise the complex number to.
     - big: Threshold value for determining large values.

   Returns:
     Result of raising the complex number to the exponent.
*/
std::complex<double> power_complex(double c, double s, int n, double big) 
{
   const double toll = 1.0e-4;

   auto itisbig = [](double x, int in, double y) {
       return in * std::log(std::abs(x)) > std::log(y);
   };

   bool foundbig = false;

   if (std::abs(c) > toll) {
       foundbig = itisbig(c, n, big);
   }
   if (std::abs(s) > toll) {
       foundbig = foundbig || itisbig(s, n, big);
   }
   if (std::abs(c) > toll && std::abs(s) > toll) {
       int n2 = n / 2;
       if (n % 2 != 0) {
           n2 = n2 + 1;
       }
       foundbig = foundbig || itisbig(c * c + s * s, n2, big);
   }

   if (foundbig) {
       return std::complex<double>(big, 0.0);
   }

   if (std::abs(c) <= toll && std::abs(s) <= toll) {
       return std::complex<double>(0.0, 0.0);
   } else {
       if (n != 1) {
           return std::pow(std::complex<double>(c, s), n);
       } else {
           return std::complex<double>(c, s);
       }
   }
}


/*******************************************
 *                                         *
 *         generate_inertia_tensor         *
 *                                         *
 *******************************************/
/**
 * @brief Calculate the inertia tensor and center of mass for a molecule given ion positions and masses.
 * 
 * @param rion Pointer to the array containing the coordinates (x, y, z) of the ions.
 * @param ion_mass Pointer to the array containing the masses of the ions.
 * @param nion The number of ions in the molecule.
 * @param inertia_tensor Pointer to the output array where the 3x3 inertia tensor will be stored.
 * 
 * @note
 * This function first calculates the center of mass of the molecule.
 * Then, it calculates the inertia tensor elements based on the position and mass of each ion relative to the center of mass.
 * 
 * @code
 * Example Usage:
 * double rion[3 * nion];  // Array to store ion positions
 * double ion_mass[nion];  // Array to store ion masses
 * int nion = 5;           // Number of ions
 * double inertia_tensor[9]; // 3x3 matrix to store inertia tensor
 * 
 * generate_inertia_tensor(rion, ion_mass, nion, inertia_tensor);
 * @endcode
 */
static void generate_inertia_tensor(const double* rion, const double* ion_mass, const int nion, double* inertia_tensor)
{
   // Calculate the center of mass
   double total_mass = 0.0;
   double center_mass[3] = {0.0, 0.0, 0.0};
   for (int ii = 0; ii < nion; ++ii)
   {
      double mass = ion_mass[ii];
      total_mass += mass;
      center_mass[0] += mass * rion[3 * ii];
      center_mass[1] += mass * rion[3 * ii + 1];
      center_mass[2] += mass * rion[3 * ii + 2];
   }
   center_mass[0] /= total_mass;
   center_mass[1] /= total_mass;
   center_mass[2] /= total_mass;

    // Initialize the elements of the inertia tensor
    double Ixx = 0.0, Iyy = 0.0, Izz = 0.0;
    double Ixy = 0.0, Ixz = 0.0, Iyz = 0.0;

   // Calculate the inertia tensor
   for (int ii = 0; ii < nion; ++ii)
   {
      double mass = ion_mass[ii];
      double dx = rion[3*ii]     - center_mass[0];
      double dy = rion[3*ii + 1] - center_mass[1];
      double dz = rion[3*ii + 2] - center_mass[2];

      Ixx += mass*(dy*dy + dz*dz);
      Iyy += mass*(dx*dx + dz*dz);
      Izz += mass*(dx*dx + dy*dy);
      Ixy -= mass*dx*dy;
      Ixz -= mass*dx*dz;
      Iyz -= mass*dy*dz;
   }

   // Fill the inertia tensor
   inertia_tensor[0] = Ixx;
   inertia_tensor[4] = Iyy;
   inertia_tensor[8] = Izz;
   inertia_tensor[3] = inertia_tensor[1] = Ixy;
   inertia_tensor[6] = inertia_tensor[2] = Ixz;
   inertia_tensor[7] = inertia_tensor[5] = Iyz;
}

/*******************************************
 *                                         *
 *          canonicalize_axes              *
 *                                         *
 *******************************************/
/* Canonicalize the moments of inertial axes

   The new axes for the re-oriented molecule are determined as the
   eigenvectors of the moments of inertia. These eigenvectors have
   arbitrary signs and degenerate eigenvectors may appear in arbitrary
   order. This leads to a certain amount of arbitrariness in the
   choice of the final orientation which is a pain for testing the code.

   This subroutine fixes the sign of the axes as well as the order
   in which degenerate axes come out. For the sign, we enforce that

      \sum_{i=1}^{N} (N-i+1) axis_i \ge 0

   where 'axis' is the new axis in terms of the old coordinate
   system.

   For the ordering of the degenerate axes, we count the number of nodes
   as the number of sign changes, and re-order the axes such that the
   ones with the fewest number of nodes come first. A node is counted every time
      axis_i * axis_{i+1} < 0

   Arguments:
      - thresh: [Input] The threshold for degeneracy
      - eig:    [Input] The moments of inertia
      - axis:   [In/Output] The axes
*/


static void canonicalize_axes(const double thresh, double *eig, double *axis)
{
   int nodes[3] = {0}; // Initialize nodes array with zeros

   // Canonicalize the signs of the axes
   for (int j = 0; j < 3; ++j)
   {
      if (axis[mindex(0, j)] * axis[mindex(1, j)] < 0.0) nodes[j] += 1;
      if (axis[mindex(1, j)] * axis[mindex(2, j)] < 0.0) nodes[j] += 1;
      double sum = 0.0;
      for (int i = 0; i < 3; ++i)
         sum += (3 - i) * axis[mindex(i, j)];
      if (sum < 0.0)
         for (int i = 0; i < 3; ++i)
            axis[mindex(i, j)] = -axis[mindex(i, j)];
   }

   // Canonicalize the order of the degenerate axes
   for (int i = 0; i < 2; ++i)
   {
      for (int j = i + 1; j < 3; ++j)
      {
         if ((std::abs(eig[j]) - std::abs(eig[i])) < thresh)
         {
            if (nodes[j] < nodes[i])
            {
               std::swap(eig[i], eig[j]);
               for (int k = 0; k < 3; ++k)
                  std::swap(axis[mindex(k, j)], axis[mindex(k, i)]);
            }
         }
      }
   }
}







/*******************************************
 *                                         *
 *         generate_principle_axes         *
 *                                         *
 *******************************************/
static void generate_principle_axes(const double *rion, const double *ion_mass, const int nion,
                                    double *inertia_tensor, double *principle_values, double *principle_axes)
{
   double tolerance = 0.0001;  // Adjust as needed based on the precision of your data

   generate_inertia_tensor(rion, ion_mass, nion, inertia_tensor);


   //Calculate the eigenvalues and eigenvectors of the inertial_tensor
   //int it_num,rot_num;
   //jacobi_eigenvalue(inertia_tensor, principle_axes, principle_values, it_num, rot_num);
   //eigsrt_hitolow(principle_values, principle_axes, 3);
   //std::cout << "iterations=" << it_num << " " << rot_num << std::endl;

   //Calculate the eigenvalues and eigenvectors of the inertial_tensor
   int ierr;
   int n = 3;
   int nn = 30;
   double xtmp1[nn];

   std::memcpy(principle_axes,inertia_tensor,9*sizeof(double));
   EIGEN_PWDFT(n,principle_axes,principle_values,xtmp1,nn,ierr);

   //sort the eigevalues and eigenvectors
   eigsrt_hitolow(principle_values, principle_axes, 3);


   //canonicalize_axes(tolerance, principle_values, principle_axes);

   //Check for right handedness, correct if not 
   /*
   double test = principle_axes[mindex(0,2)]*(principle_axes[mindex(1,0)]*principle_axes[mindex(2,1)]-principle_axes[mindex(2,0)]*principle_axes[mindex(1,1)]) 
               + principle_axes[mindex(1,2)]*(principle_axes[mindex(2,0)]*principle_axes[mindex(0,1)]-principle_axes[mindex(0,0)]*principle_axes[mindex(2,1)])
               + principle_axes[mindex(2,2)]*(principle_axes[mindex(0,0)]*principle_axes[mindex(1,1)]-principle_axes[mindex(1,0)]*principle_axes[mindex(0,1)]);
   if (test>0.0) return;
   if (std::abs(principle_values[0] - principle_values[1]) <= tolerance)
   {
      std::swap(principle_values[0],principle_values[1]);
      for (auto i=0; i<3; ++i)
         std::swap(principle_axes[mindex(i,0)],principle_axes[mindex(i,1)]);
      return;
   }
   if (std::abs(principle_values[1] - principle_values[2]) <= tolerance)
   {
      std::swap(principle_values[1],principle_values[2]);
      for (auto i=0; i<3; ++i)
         std::swap(principle_axes[mindex(i,1)],principle_axes[mindex(i,2)]);
      return;
   }
   for (auto i=0; i<3; ++i)
      principle_axes[mindex(i,2)] = -principle_axes[mindex(i,2)];
      */


   return;
}

/*******************************************
 *                                         *
 *          has_reflection_axis            *
 *                                         *
 *******************************************/
/**
 * @brief Check for the presence of a reflection axis in a molecule's geometry.
 *
 * This function determines whether a molecule's geometry exhibits reflection symmetry
 * along a specified axis defined by the vector u. It checks if ions in the molecule
 * have symmetric positions with respect to the reflection operation defined by the
 * reflection matrix.
 *
 * @param rion          An array of ion coordinates (x, y, z) for all ions in the molecule.
 * @param ion_mass      An array of ion masses corresponding to each ion in rion.
 * @param nion          The number of ions in the molecule.
 * @param u             The axis vector along which reflection symmetry is checked.
 * @param sym_tolerance A tolerance value for considering coordinates as nearly equal.
 *
 * @return True if the molecule exhibits reflection symmetry along the specified axis,
 *         false otherwise.
 */
static bool has_reflection_axis(const double* rion, const double* ion_mass, const int nion, const double u[3], const double sym_tolerance)
{
   std::vector<bool> ion_checked(nion, false);

   // Define the reflection matrix for reflection along the u axis - reflect_matrix = I - 2*u*u'
   double reflection_matrix[9] = {1.0 -2.0*u[0]*u[0],     -2.0*u[1]*u[0],      -2.0*u[2]*u[0],
                                      -2.0*u[0]*u[1], 1.0 -2.0*u[1]*u[1],      -2.0*u[2]*u[1],
                                      -2.0*u[0]*u[2],     -2.0*u[1]*u[2], 1.0 - 2.0*u[2]*u[2] };

   //Define lambda function
   auto check_coord2_reflection = [sym_tolerance](const double *reflect, const double *rion1, const double *rion2) {
   double x2 = reflect[0]*rion1[0] + reflect[3]*rion1[1] + reflect[6]*rion1[2];
   double y2 = reflect[1]*rion1[0] + reflect[4]*rion1[1] + reflect[7]*rion1[2];
   double z2 = reflect[2]*rion1[0] + reflect[5]*rion1[1] + reflect[8]*rion1[2];
   return ((std::abs(x2 - rion2[0]) < sym_tolerance) &&
           (std::abs(y2 - rion2[1]) < sym_tolerance) &&
           (std::abs(z2 - rion2[2]) < sym_tolerance));
   };
  //Check for atom at origin
   for (int ii = 0; ii<nion; ++ii)
      if ((std::abs(rion[3*ii])   < sym_tolerance)   &&
          (std::abs(rion[3*ii+1]) < sym_tolerance) &&
          (std::abs(rion[3*ii+2]) < sym_tolerance))
         ion_checked[ii] = true;

   for (int ii = 0; ii<nion; ++ii)
      for (int jj = 0; jj<nion; ++jj)
      if (ion_mass[ii]==ion_mass[jj])
      {
         if (check_coord2_reflection(reflection_matrix,rion+3*ii,rion+3*jj))
         {
            ion_checked[ii] = true;
         }
      }

   return  std::all_of(ion_checked.begin(), ion_checked.end(), [](bool val) { return val; });
}


/*******************************************
 *                                         *
 *       has_improper_rotation_axis        *
 *                                         *
 *******************************************/
/**
 * @brief Determine if a set of ions can be aligned to form an improper rotation axis.
 *
 * This function checks if a given set of ions with the same mass can be improperly rotated
 * to align with each other, within a specified tolerance.
 *
 * @param rion           An array of ion coordinates (3D) for all ions.
 * @param ion_mass       An array of ion masses.
 * @param nion           The number of ions in the system.
 * @param u              A 3D vector representing the rotation axis.
 * @param theta          The rotation angle in radians.
 * @param sym_tolerance  A tolerance value for coordinate comparison.
 *
 * @return true if an improper rotation axis can be formed with the given ions and parameters, false otherwise.
 */
static bool has_improper_rotation_axis(const double* rion, const double* ion_mass, const int nion, const double u[3], const double theta, const double sym_tolerance)
{
   std::vector<bool> ion_checked(nion, false);

   //Determine the rotation matrix (rot_matrix) from u and theta
   double wx = std::cos(theta);
   double wy = std::sin(theta);
   double mwx = 1.0-wx;
   double rot_matrix[9] = {wx+u[0]*u[0]*mwx,            u[1]*u[0]*mwx+u[2]*wy,    u[2]*u[0]*mwx-u[1]*wy,
                              u[0]*u[1]*mwx-u[2]*wy, wx+u[1]*u[1]*mwx,            u[2]*u[1]*mwx+u[0]*wy,
                              u[0]*u[2]*mwx+u[1]*wy,    u[1]*u[2]*mwx-u[0]*wy, wx+u[2]*u[2]*mwx};

   // Define the reflection matrix for reflection along the u axis - reflect_matrix = I - 2*u*u'
   double reflection_matrix[9] = {1.0 -2.0*u[0]*u[0],     -2.0*u[1]*u[0],      -2.0*u[2]*u[0],
                                      -2.0*u[0]*u[1], 1.0 -2.0*u[1]*u[1],      -2.0*u[2]*u[1],
                                      -2.0*u[0]*u[2],     -2.0*u[1]*u[2], 1.0 - 2.0*u[2]*u[2] };

   // Multiply the reflection matrix by the rotation matrix to combine reflection and rotation
   // Multiply the reflection matrix by the rotation matrix to combine reflection and rotation
   double imp_rot_matrix[9];
   for (int i = 0; i < 3; ++i) 
   {
      for (int j = 0; j < 3; ++j) 
      {
         imp_rot_matrix[3*i+j] = 0.0;
         for (int k = 0; k < 3; ++k) 
            imp_rot_matrix[3*i+j] += reflection_matrix[3*i+k] * rot_matrix[3*k+j]; 
      }
   }

   //Define lambda function
   auto check_coord2_improper_rotate = [sym_tolerance](const double *rot, const double *rion1, const double *rion2) {
   double x2 = rot[0]*rion1[0] + rot[3]*rion1[1] + rot[6]*rion1[2];
   double y2 = rot[1]*rion1[0] + rot[4]*rion1[1] + rot[7]*rion1[2];
   double z2 = rot[2]*rion1[0] + rot[5]*rion1[1] + rot[8]*rion1[2];
   return ((std::abs(x2 - rion2[0]) < sym_tolerance) &&
           (std::abs(y2 - rion2[1]) < sym_tolerance) &&
           (std::abs(z2 - rion2[2]) < sym_tolerance));
   };


   //Check for atom at origin
   for (int ii = 0; ii<nion; ++ii)
      if ((std::abs(rion[3*ii])   < sym_tolerance)   &&
          (std::abs(rion[3*ii+1]) < sym_tolerance) &&
          (std::abs(rion[3*ii+2]) < sym_tolerance))
         ion_checked[ii] = true;

   for (int ii = 0; ii<nion; ++ii)
      for (int jj = 0; jj<nion; ++jj)
      if (ion_mass[ii]==ion_mass[jj])
      {
         if (check_coord2_improper_rotate(imp_rot_matrix,rion+3*ii,rion+3*jj))
         {
            ion_checked[ii] = true;
         }
      }

   return  std::all_of(ion_checked.begin(), ion_checked.end(), [](bool val) { return val; });
}

    

/*******************************************
 *                                         *
 *          has_rotation_axis              *
 *                                         *
 *******************************************/
 /**
 * @brief Determine if a set of ions can be aligned to form a rotation axis.
 *
 * This function checks if a given set of ions with the same mass can be rotated
 * to align with a specified rotation axis and angle, within a specified tolerance.
 *
 * @param rion         An array of ion coordinates (3D) for all ions.
 * @param ion_mass     An array of ion masses.
 * @param nion         The number of ions in the system.
 * @param u            A 3D vector representing the rotation axis.
 * @param theta        The rotation angle in radians.
 * @param sym_tolerance A tolerance value for coordinate comparison.
 *
 * @return true if a rotation axis can be formed with the given ions and parameters, false otherwise.
 *
 * @note The function checks if ions with the same mass can be aligned to form a rotation axis.
 *       It uses a rotation matrix to transform ion coordinates and compares them within the specified tolerance.
 *       The result indicates whether all ions can be aligned to form a rotation axis.
 */
static bool has_rotation_axis(const double* rion, const double* ion_mass, const int nion, const double u[3], const double theta, const double sym_tolerance)
{
   std::vector<bool> ion_checked(nion, false);

   // Determine the rotation matrix (rot_matrix) from u and theta
   double wx = std::cos(theta);
   double wy = std::sin(theta);
   double mwx = 1.0-wx;
   double rot_matrix[9] = {wx+u[0]*u[0]*mwx,            u[1]*u[0]*mwx+u[2]*wy,    u[2]*u[0]*mwx-u[1]*wy,
                              u[0]*u[1]*mwx-u[2]*wy, wx+u[1]*u[1]*mwx,            u[2]*u[1]*mwx+u[0]*wy,
                              u[0]*u[2]*mwx+u[1]*wy,    u[1]*u[2]*mwx-u[0]*wy, wx+u[2]*u[2]*mwx};

   //Define lambda function
   auto check_coord2_rotate = [sym_tolerance](const double *rot, const double *rion1, const double *rion2) {
      double x2 = rot[0]*rion1[0] + rot[3]*rion1[1] + rot[6]*rion1[2];
      double y2 = rot[1]*rion1[0] + rot[4]*rion1[1] + rot[7]*rion1[2];
      double z2 = rot[2]*rion1[0] + rot[5]*rion1[1] + rot[8]*rion1[2];
      return ((std::abs(x2 - rion2[0]) < sym_tolerance) &&
              (std::abs(y2 - rion2[1]) < sym_tolerance) &&
              (std::abs(z2 - rion2[2]) < sym_tolerance));
   };

   //Check for atom at origin
   for (int ii = 0; ii<nion; ++ii)
      if ((std::abs(rion[3*ii])   < sym_tolerance)   &&
          (std::abs(rion[3*ii+1]) < sym_tolerance) &&
          (std::abs(rion[3*ii+2]) < sym_tolerance))
         ion_checked[ii] = true;

   for (int ii = 0; ii<nion; ++ii)
      for (int jj = 0; jj<nion; ++jj)
      if (ion_mass[ii]==ion_mass[jj])
      {
         if (check_coord2_rotate(rot_matrix,rion+3*ii,rion+3*jj))
         {
            ion_checked[ii] = true;
         }
      }

   return  std::all_of(ion_checked.begin(), ion_checked.end(), [](bool val) { return val; });
}

/*******************************************
 *                                         *
 *             has_C2_axes                 *
 *                                         *
 *******************************************/
static bool has_C2_axes(const double* rion, const double* ion_mass, const int nion, const double sym_tolerance)
{
   bool has_C2_axis = false;

  //mid-point bond - the direction is passes thru a bond mid-point - rc = (ri+rj)/2 - rcom
   for (int i=0;   i<nion; ++i)
   for (int j=i+1; j<nion; ++j)
   if (ion_mass[i] == ion_mass[j])
   {
      double dx = 0.5*(rion[3*i]  +rion[3*j]);
      double dy = 0.5*(rion[3*i+1]+rion[3*j+1]);
      double dz = 0.5*(rion[3*i+2]+rion[3*j+2]);
      double dij = std::sqrt(dx*dx + dy*dy + dz*dz);
      double ij_axis[3] = {dx/dij, dy/dij, dz/dij};

      has_C2_axis = has_rotation_axis(rion,ion_mass,nion,ij_axis,M_PI,sym_tolerance);
      if (has_C2_axis) return has_C2_axis;
   }

  //atom - the direction is passes thru an atom - rc = ri - rcom
   for (int i=0; i<nion; ++i)
   {
      double dx = rion[3*i];
      double dy = rion[3*i+1];
      double dz = rion[3*i+2];
      double dij = std::sqrt(dx*dx + dy*dy + dz*dz);
      double ij_axis[3] = {dx/dij, dy/dij, dz/dij};

      has_C2_axis = has_rotation_axis(rion,ion_mass,nion,ij_axis,M_PI,sym_tolerance);
      if (has_C2_axis) return has_C2_axis;
   }

   //cross bonds - direction is rc = rij x rkl
   for (int i=0;   i<nion; ++i)
   for (int j=i+1; j<nion; ++j)
   if (ion_mass[i] == ion_mass[j])
   for (int k=j+1; k<nion; ++k)
   if (ion_mass[k] == ion_mass[i])
   for (int l=k+1; l<nion; ++l)
   if (ion_mass[l] == ion_mass[i])
   {
      double try_axis[3];
      double rij[3] = {(rion[3*i]-rion[3*j]), (rion[3*i+1]-rion[3*j+1]), (rion[3*i+2]-rion[3*j+2])};
      double rkl[3] = {(rion[3*k]-rion[3*l]), (rion[3*k+1]-rion[3*l+1]), (rion[3*k+2]-rion[3*l+2])};

      normalized_cross_product(rij,rkl,try_axis);
      has_C2_axis = has_rotation_axis(rion,ion_mass,nion,try_axis,M_PI,sym_tolerance);

      if (has_C2_axis) return has_C2_axis;
   }

   return has_C2_axis;
}

/*******************************************
 *                                         *
 *             has_sigmav_axis             *
 *                                         *
 *******************************************/
static bool has_sigma_v_axis(const double* rion, const double* ion_mass, const int nion, const double sym_tolerance)
{
   bool has_sigma_v = false;

   //bonds - the direction is along a bond - rc = (ri-rj)/2
   for (int i=0;   i<nion; ++i)
   for (int j=i+1; j<nion; ++j)
   if (ion_mass[i] == ion_mass[j])
   {
      double dx = (rion[3*i]  -rion[3*j]);
      double dy = (rion[3*i+1]-rion[3*j+1]);
      double dz = (rion[3*i+2]-rion[3*j+2]);
      double dij = std::sqrt(dx*dx + dy*dy + dz*dz);
      double ij_axis[3] = {dx/dij, dy/dij, dz/dij};

      //Only look for sigma_v, skip sigma_h
      if ((std::abs(ij_axis[0]) > sym_tolerance) ||
          (std::abs(ij_axis[1]) > sym_tolerance) ||
          ((std::abs(ij_axis[2])-1.0) > sym_tolerance))
      {
         has_sigma_v = has_reflection_axis(rion,ion_mass,nion,ij_axis,sym_tolerance);
      }

      if (has_sigma_v) 
      {
         return has_sigma_v;
      }
   }
   return has_sigma_v;
}


/*******************************************
 *                                         *
 *         has_proper_rotation_C5          *
 *                                         *
 *******************************************/
static bool has_proper_rotation_C5(const double* rion, const double* ion_mass, const int nion, const double sym_tolerance)
{
   int count=0;
   double try_axis[3],try_axisc[3], r1[3],r2[3];
   double theta5     = 2*M_PI/5.0;
   double theta10    = 2*M_PI/10.0;
   double cos_theta5 = std::cos(theta5);

   bool found = false;

   //std::map<double, std::array<double, 3>> dist_vec;
   std::vector<std::array<double, 3>> c5_axes;
   std::vector<std::array<double, 3>> s10_axes;

   for (int i=0;   i<nion; ++i)
   for (int j=i+1; j<nion; ++j)
   if (ion_mass[i] == ion_mass[j])
   for (int k=j+1; k<nion; ++k)
   if (ion_mass[k] == ion_mass[i])
   for (int l=k+1; l<nion; ++l)
   if (ion_mass[l] == ion_mass[i])
   for (int m=l+1; m<nion; ++m)
   if (ion_mass[m] == ion_mass[i])
   {
      // Check if the 5 lowest keys are approximately equal
      if (is_pentagon(i,j,k,l,m,rion,nion, sym_tolerance))
      {
          ++count;

         //Extract the points based on the i,j,k,l,m indices
         std::vector<std::array<double, 3>> points;
         points.push_back({rion[3 * i], rion[3 * i + 1], rion[3 * i + 2]});
         points.push_back({rion[3 * j], rion[3 * j + 1], rion[3 * j + 2]});
         points.push_back({rion[3 * k], rion[3 * k + 1], rion[3 * k + 2]});
         points.push_back({rion[3 * l], rion[3 * l + 1], rion[3 * l + 2]});
         points.push_back({rion[3 * m], rion[3 * m + 1], rion[3 * m + 2]});

         //Define the center for sorting
         std::array<double, 3> center = calculate_center(points);

         //Sort points in increasing order of their angles with respect to the center
         sort_points_spherical(points, center);
         double magc = std::sqrt(center[0]*center[0] + center[1]*center[1] + center[2]*center[2]);

         r1[0] = points[1][0] - points[0][0];
         r1[1] = points[1][1] - points[0][1];
         r1[2] = points[1][2] - points[0][2];
         r2[0] = points[2][0] - points[1][0];
         r2[1] = points[2][1] - points[1][1];
         r2[2] = points[2][2] - points[1][2];

         //Generate the rotation axis
         normalized_cross_product(r2,r1,try_axis);
         bool has_c5_rotation  = has_rotation_axis(rion,ion_mass,nion,try_axis,theta5,sym_tolerance);
         bool has_s10_rotation = has_improper_rotation_axis(rion,ion_mass,nion,try_axis,theta10,sym_tolerance);

         if (has_c5_rotation)
            c5_axes.push_back({try_axis[0],try_axis[1],try_axis[2]});

         if (has_s10_rotation)
            s10_axes.push_back({try_axis[0],try_axis[1],try_axis[2]});
      }
   }

   double tolerance = sym_tolerance;

   // Sort and count unique points within the specified tolerance
   // Replace auto with explicit parameter types (e.g., const std::array<double, 3>&)
   std::sort(c5_axes.begin(), c5_axes.end(), [tolerance](const std::array<double, 3>& a, const std::array<double, 3>& b) {
    for (int i = 0; i < 3; ++i) {
        if (std::abs(a[i] - b[i]) > tolerance) {
            return a[i] < b[i];
        }
    }
    return false; // a and b are considered equal within tolerance
   });

   c5_axes.erase(std::unique(c5_axes.begin(), c5_axes.end(), [tolerance](const std::array<double, 3>& a, const std::array<double, 3>& b) {
    for (int i = 0; i < 3; ++i) {
        if (std::abs(a[i] - b[i]) > tolerance) {
            return true; // Points are different
        }
    }
    return false; // Points are considered equal within tolerance
   }), c5_axes.end());


   //Count the number of unique points
   int c5_uniqueCount = c5_axes.size();

   //std::cout << "There are " << c5_uniqueCount << " unique C5 axes!" << std::endl;
   found = (c5_uniqueCount == 12);


   // Sort and count unique points within the specified tolerance
   // Replace auto with explicit parameter types (e.g., const std::array<double, 3>&)
   std::sort(s10_axes.begin(), s10_axes.end(), [tolerance](const std::array<double, 3>& a, const std::array<double, 3>& b) {
    for (int i = 0; i < 3; ++i) {
        if (std::abs(a[i] - b[i]) > tolerance) {
            return a[i] < b[i];
        }
    }
    return false; // a and b are considered equal within tolerance
   });

   s10_axes.erase(std::unique(s10_axes.begin(), s10_axes.end(), [tolerance](const std::array<double, 3>& a, const std::array<double, 3>& b) {
    for (int i = 0; i < 3; ++i) {
        if (std::abs(a[i] - b[i]) > tolerance) {
            return true; // Points are different
        }
    }
    return false; // Points are considered equal within tolerance
   }), s10_axes.end());

   //Count the number of unique points
   int s10_uniqueCount = s10_axes.size();

   //std::cout << "There are " << s10_uniqueCount << " unique S10 axes!" << std::endl;

   return found;
}


/*******************************************
 *                                         *
 *         has_inversion_center            *
 *                                         *
 *******************************************/
/**
 * @brief Determine whether a molecule has an inversion center based on ion coordinates and ion masses.
 *
 * @param rion Pointer to the array containing the coordinates (x, y, z) of the ions.
 * @param ion_mass Pointer to the array containing the masses of the ions.
 * @param nion The number of ions in the molecule.
 * @param sym_tolerance The tolerance level for coordinate inversion check.
 *
 * @return True if an inversion center is present; False otherwise.
 *
 * @note
 * The function checks each ion to see if there is another ion with the same mass but at an inverted coordinate.
 * All ions must have such a pair for the molecule to have an inversion center.
 *
 * @code
 * Example Usage:
 * double rion[3 * nion];      // Ion positions
 * double ion_mass[nion];      // Masses of ions
 * int nion;                   // Number of ions
 * double sym_tolerance = 1e-5; // Symmetry tolerance
 * bool result = has_inversion_center(rion, ion_mass, nion, sym_tolerance);
 * @endcode
 */
static bool has_inversion_center(const double* rion, const double* ion_mass, const int nion, const double sym_tolerance)
{
   std::vector<bool> ion_checked(nion, false);

   //Define lambda function
   auto coord2_inverse = [sym_tolerance](const double *rion1, const double *rion2) {
        return ((std::abs(rion1[0] + rion2[0]) < sym_tolerance) &&
                (std::abs(rion1[1] + rion2[1]) < sym_tolerance) &&
                (std::abs(rion1[2] + rion2[2]) < sym_tolerance));
   };

   //Check for atom at origin
   for (int ii = 0; ii<nion; ++ii)
      if ((std::abs(rion[3*ii])   < sym_tolerance)   &&
           (std::abs(rion[3*ii+1]) < sym_tolerance) &&
           (std::abs(rion[3*ii+2]) < sym_tolerance))
         ion_checked[ii] = true;

   for (int ii = 0; ii<nion; ++ii)
   for (int jj = ii+1; jj<nion; ++jj)
   if (((!ion_checked[ii]) && (!ion_checked[jj])) && (ion_mass[ii] == ion_mass[jj]))
   {
      if (coord2_inverse(rion+3*ii,rion+3*jj))
      {
         ion_checked[ii] = true;
         ion_checked[jj] = true;
      }
   }

   //for (int ii = 0; ii<nion; ++ii)
   //    std::cout << "ion ii=" << ii << " checked=" << ion_checked[ii] << std::endl;

   bool all_ions_checked = std::all_of(ion_checked.begin(), ion_checked.end(), [](bool val) { return val; });

   return all_ions_checked;
}


/*******************************************
 *                                         *
 *              is_linear                  *
 *                                         *
 *******************************************/
/**
 * @brief Determine whether a molecule is linear based on ion coordinates.
 * 
 * @param rion Pointer to the array containing the coordinates (x, y, z) of the ions.
 * @param nion The number of ions in the molecule.
 * @param sym_tolerance The tolerance level for checking linearity.
 * 
 * @return True if the molecule is linear; False otherwise.
 * 
 * @note
 * For a molecule to be considered linear, all angles between ions must be nearly zero degrees.
 * The function calculates the cosine of the angle between the first ion and each subsequent ion.
 * If all cosines are within `sym_tolerance` of 1, the molecule is considered linear.
 * 
 * @code
 * Example Usage:
 * double rion[3 * nion];      // Ion positions
 * int nion;                   // Number of ions
 * double sym_tolerance = 1e-5; // Linearity tolerance
 * bool result = is_linear(rion, nion, sym_tolerance);
 * @endcode
 */
static bool is_linear(const double* rion, const int nion, const double sym_tolerance)
{
   double cos_theta = 1.0;     // Cosine of the angle for a linear molecule (0 degrees)

   if (nion < 3)
      return true;
   else
   {
      double r10[3] = {rion[3] - rion[0], rion[4] - rion[1], rion[5] - rion[2]};
      double d10 = std::sqrt(r10[0] * r10[0] + r10[1] * r10[1] + r10[2] * r10[2]);

      for (auto i = 2; i < nion; ++i)
      {
         double ri0[3] = {rion[3 * i] - rion[0], rion[3 * i + 1] - rion[1], rion[3 * i + 2] - rion[2]};
         double di0 = std::sqrt(ri0[0] * ri0[0] + ri0[1] * ri0[1] + ri0[2] * ri0[2]);
         double cos_theta_i = (rion[3 * i] * r10[0] + rion[3 * i + 1] * r10[1] + rion[3 * i + 2] * r10[2]) / (d10 * di0);

         if (std::abs(cos_theta_i) - cos_theta > sym_tolerance)
            return false;
      }
   }
   return true;
}



static bool has_proper_rotation_Cnv(const double* rion, const double* ion_mass, const int nion, const double sym_tolerance, const int nfold)
{
    double cos_theta;

    if (nfold == 2) {
        cos_theta = 0.0;  // 180 degrees for C2 axis
    } else if (nfold == 3) {
        cos_theta = -0.5; // 120 degrees for C3 axis
    } else if (nfold == 5) {
        cos_theta = 0.3090; // 72 degrees for C5 axis
    } else {
        return false;     // Not supported
    }

    for (int i = 0; i < nion; ++i) {
        for (int j = 0; j < nion; ++j) {
            if (i != j && ion_mass[i] == ion_mass[j]) {  // Adding mass comparison
                int count_matches = 0;

                double dot_product = 0.0;
                for (int l = 0; l < 3; ++l) {
                    dot_product += rion[3*i + l] * rion[3*j + l];
                }
                dot_product /= sqrt(rion[3*i] * rion[3*i] + rion[3*i+1] * rion[3*i+1] + rion[3*i+2] * rion[3*i+2]) *
                                sqrt(rion[3*j] * rion[3*j] + rion[3*j+1] * rion[3*j+1] + rion[3*j+2] * rion[3*j+2]);

                if (std::abs(dot_product - cos_theta) < sym_tolerance) {
                    count_matches++;
                }

                if (count_matches == nfold) {
                    bool is_centered = true;
                    for (int k = 0; k < 3; ++k) {
                        if (std::abs(rion[3*i + k] + rion[3*j + k]) > sym_tolerance) {
                            is_centered = false;
                            break;
                        }
                    }
                    if (is_centered) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}





/*******************************************
 *                                         *
 *       determine_spherical_group         *
 *                                         *
 *******************************************/
 /* only call if nion>3
 */
static void determine_spherical_group(const double *rion, const double *ion_mass, const int nion,
                           const double sym_tolerance,
                           std::string&  group_name, int& group_rank, std::string& rotation_type,
                           double *inertia_tensor,  double *inertia_moments,  double *inertia_axes)
{
    //There are enough atoms to compute a multipole2
   double eccentricity, asphericity;

   eccentricity_and_asphericity(rion,ion_mass,nion,eccentricity,asphericity);

   double Qxx,Qyy,Qzz,Qxy,Qxz,Qyz;
   generate_quadrupole(rion,ion_mass,nion,Qxx,Qyy,Qzz,Qxy,Qxz,Qyz);
   double abssum = std::abs(Qxx) + std::abs(Qyy) + std::abs(Qzz) + std::abs(Qxy) + std::abs(Qxz) + std::abs(Qyz);
   double diagonal_sum = Qxx + Qyy + Qzz;

   double Oxxx,Oxxy,Oxxz, Oxyx,Oxyy,Oxyz, Oxzx,Oxzy,Oxzz;
   double Oyxx,Oyxy,Oyxz, Oyyx,Oyyy,Oyyz, Oyzx,Oyzy,Oyzz;
   double Ozxx,Ozxy,Ozxz, Ozyx,Ozyy,Ozyz, Ozzx,Ozzy,Ozzz;
   generate_octapole(rion,ion_mass,nion,
                              Oxxx, Oxxy, Oxxz,
                              Oxyx, Oxyy, Oxyz,
                              Oxzx, Oxzy, Oxzz,
                              Oyxx, Oyxy, Oyxz,
                              Oyyx, Oyyy, Oyyz,
                              Oyzx, Oyzy, Oyzz,
                              Ozxx, Ozxy, Ozxz,
                              Ozyx, Ozyy, Ozyz,
                              Ozzx, Ozzy, Ozzz);

   double abs_octapole = std::abs(Oxxx) + std::abs(Oxxy) + std::abs(Oxxz) 
                       + std::abs(Oxyx) + std::abs(Oxyy) + std::abs(Oxyz) 
                       + std::abs(Oxzx) + std::abs(Oxzy) + std::abs(Oxzz) 
                       + std::abs(Oyxx) + std::abs(Oyxy) + std::abs(Oyxz) 
                       + std::abs(Oyyx) + std::abs(Oyyy) + std::abs(Oyyz) 
                       + std::abs(Oyzx) + std::abs(Oyzy) + std::abs(Oyzz) 
                       + std::abs(Ozxx) + std::abs(Ozxy) + std::abs(Ozxz) 
                       + std::abs(Ozyx) + std::abs(Ozyy) + std::abs(Ozyz) 
                       + std::abs(Ozzx) + std::abs(Ozzy) + std::abs(Ozzz);

   bool has_c5 = has_proper_rotation_C5(rion,ion_mass,nion,sym_tolerance);

   //bool is_spherical = (diagonal_sum < sym_tolerance) && (abssum > sym_tolerance);
   bool is_spherical = has_c5;
   if (is_spherical)
   {
      bool has_inversion = has_inversion_center(rion,ion_mass,nion,sym_tolerance);
      if (has_inversion)
      {
         group_name = "I_h";
         group_rank = 120;
      }
      else 
      {
         group_name = "I";
         group_rank = 60;
      }
   }
   else
   {
       //group_name = "T_d, T, T_h, O_h, O";
       // still have to look for T and O groups!

       bool has_inversion = has_inversion_center(rion,ion_mass,nion,sym_tolerance);

       if (has_inversion)
       {
          if (abs_octapole<sym_tolerance)
          {
             group_name = "O_h";
             group_rank = 48;
          }
          else
          {
             group_name = "T_h";
             group_rank = 24;
          }
       }
       else
       {
          group_name = "T_d";
          group_rank = 24;
       }
       //group_name = "O"; group_rank = 24;
       //group_name = "T"; group_rank = 12;

       //if ((abs_octapole/27.0)<sym_tolerance)
       //   group_name = "O_h O, Th";
       //else
       //   group_name = "T_d, T";
       
   }
}

/*******************************************
 *                                         *
 *       determine_symmetric_group         *
 *                                         *
 *******************************************/
static void determine_symmetric_group(const double *rion, const double *ion_mass, const int nion,
                           const double sym_tolerance,
                           std::string&  group_name, int& group_rank, std::string& rotation_type,
                           double *inertia_tensor,  double *inertia_moments,  double *inertia_axes)
{
   bool has_Cn=false;
   bool has_Cx=false;
   bool has_Cy=false;
   bool has_Cz=false;
   double x_axis[3] = {1.0, 0.0, 0.0};
   double y_axis[3] = {0.0, 1.0, 0.0};
   double z_axis[3] = {0.0, 0.0, 1.0};

   // Create a map to store the counts of each unique element
   std::map<double, int> elementCounts;

   //Iterate through the ion_mass array located off the z-axis and count the occurrences of each element
   for (int ii=0; ii<nion; ++ii) 
      if ((std::abs(rion[3*ii])>sym_tolerance) ||
          (std::abs(rion[3*ii+1])>sym_tolerance))
      {
         elementCounts[ion_mass[ii]]++;
      }

   //Initialize nmax with the count of the first element
   int nmax = elementCounts.begin()->second;

   //Calculate the GCD of all element counts
   for (const auto& pair : elementCounts) 
      nmax = gcd(nmax, pair.second);

   //Decide the power of main axis
   int n = nmax+1;
   do
   {
      --n;
      has_Cz = has_rotation_axis(rion,ion_mass,nion,z_axis, 2*M_PI/((double) n),sym_tolerance);
   } while ((n>2) && (!has_Cz));
   bool has_S2nz = has_improper_rotation_axis(rion,ion_mass,nion,z_axis, M_PI/((double) n),sym_tolerance);

   auto has_sigma_yz = has_reflection_axis(rion,ion_mass,nion,x_axis,sym_tolerance);
   auto has_sigma_xz = has_reflection_axis(rion,ion_mass,nion,y_axis,sym_tolerance);
   auto has_sigma_xy = has_reflection_axis(rion,ion_mass,nion,z_axis,sym_tolerance);

   bool has_sigma_v = has_sigma_v_axis(rion,ion_mass,nion,sym_tolerance);

   bool has_C2_axis = false;
   if ((n==2) && has_Cz)
       has_C2_axis = true;
    else
       has_C2_axis = has_C2_axes(rion,ion_mass,nion,sym_tolerance);

   bool has_inversion = has_inversion_center(rion,ion_mass,nion,sym_tolerance);

   if (has_C2_axis)
   {
      if (has_sigma_v)
      {
         // Check if n is odd or  has_i 
         //if ((n%2==1) ^ has_inversion)
         if ((n%2==1) != has_inversion)
         {
            group_name = "D_" + std::to_string(n) + "h";
            group_rank = 4*n;
         }
         else
         {
            group_name = "D_" + std::to_string(n) + "d";
            group_rank = 4*n;
         }
      }
      else
      {
         group_name = "D_" + std::to_string(n);
         group_rank = 2*n;
      }
   }
   else
   {
      if (has_sigma_v)
      {
         group_name = "C_" + std::to_string(n) + "v";
         group_rank = 2*n;
      }
      else if (has_sigma_xy)
      {
         group_name = "C_" + std::to_string(n) + "h";
         group_rank = 2*n;
      }
      else if (has_S2nz)
      {
         group_name = "S_" + std::to_string(2*n);
         group_rank = 2*n;
      }
      else
      {
         group_name = "C_" + std::to_string(n);
         group_rank = n;
      }
   }
}

/*******************************************
 *                                         *
 *       determine_asymmetric_group        *
 *                                         *
 *******************************************/

static void determine_asymmetric_group(const double *rion, const double *ion_mass, const int nion,
                           const double sym_tolerance,
                           std::string&  group_name, int& group_rank,  std::string& rotation_type,
                           double *inertia_tensor,  double *inertia_moments,  double *inertia_axes)
{
   group_name = "D_2h, D_2, C_2v, C_2h, C2, Ci, Cs, C1";
   double x_axis[3] = {1.0, 0.0, 0.0};
   double y_axis[3] = {0.0, 1.0, 0.0};
   double z_axis[3] = {0.0, 0.0, 1.0};


   bool has_C2_x = has_rotation_axis(rion,ion_mass,nion,x_axis, M_PI,sym_tolerance);
   bool has_C2_y = has_rotation_axis(rion,ion_mass,nion,y_axis, M_PI,sym_tolerance);
   bool has_C2_z = has_rotation_axis(rion,ion_mass,nion,z_axis, M_PI,sym_tolerance);
   int c2_count = has_C2_x + has_C2_y + has_C2_z;

   auto has_sigma_yz = has_reflection_axis(rion,ion_mass,nion,x_axis,sym_tolerance);
   auto has_sigma_xz = has_reflection_axis(rion,ion_mass,nion,y_axis,sym_tolerance);
   auto has_sigma_xy = has_reflection_axis(rion,ion_mass,nion,z_axis,sym_tolerance);

   bool has_inversion = has_inversion_center(rion,ion_mass,nion,sym_tolerance);

   if (c2_count==3)
   {
      if (has_sigma_xy)
      {
         group_name = "D_2h";
         group_rank = 8;
      }
      else
      {
         group_name = "D_2";
         group_rank = 4;
      }
   }
   else if (c2_count==1)
   {
      if ((has_sigma_yz) && (has_sigma_xz))
      {
         group_name = "C_2v";
         group_rank = 4;
      }
      else if (has_sigma_xy)
      {
         group_name = "C_2h";
         group_rank = 4;
      }
      else
      {
         group_name = "C_2";
         group_rank = 2;
      }
   }
   else
   {
      if (has_inversion)
      {
         group_name = "Ci";
         group_rank = 2;
      }
      else if (has_sigma_yz || has_sigma_xz || has_sigma_xy)
      {
         group_name = "Cs";
         group_rank = 2;
      }
      else
      {
         group_name = "C1";
         group_rank = 1;
      }
   }
}


/*******************************************
 *                                         *
 *         determine_point_group           *
 *                                         *
 *******************************************/
void determine_point_group(const double *rion, const double *ion_mass, const int nion,
                           const double sym_tolerance,
                           std::string&  group_name, int &group_rank, std::string& rotation_type,
                           double *inertia_tensor,  double *inertia_moments,  double *inertia_axes, 
                           double *rion2)
{
   //rion2 is rion shift to have a center of mass==0, it will be used thruout
   shift_to_center_mass(rion,ion_mass,nion,rion2);

   double m_total = 0.0;
   for (auto ii=0; ii<nion; ++ii) 
       m_total += ion_mass[ii];

   group_name = "unknown";
   rotation_type="unknown";

   generate_principle_axes(rion2,ion_mass,nion,inertia_tensor,inertia_moments,inertia_axes);

   if ( ( ( inertia_moments[0]*inertia_moments[0] 
          + inertia_moments[1]*inertia_moments[1] 
          + inertia_moments[2]*inertia_moments[2])
        /(m_total*m_total))
       < (sym_tolerance*sym_tolerance))
   {
      rotation_type = "point";
      group_name = "SO(3)";
      group_rank = -1;
      inertia_axes[0] = 1.0; inertia_axes[1] = 0.0; inertia_axes[2] = 0.0;
      inertia_axes[3] = 0.0; inertia_axes[4] = 1.0; inertia_axes[5] = 0.0;
      inertia_axes[6] = 0.0; inertia_axes[7] = 0.0; inertia_axes[8] = 1.0;
   }
   else if (((inertia_moments[2]*inertia_moments[2])/(m_total*m_total)) < (sym_tolerance*sym_tolerance))
   {
      rotation_type = "linear";

      // align the molecular axes along the inertia_axes
      align_to_axes(rion2,nion, inertia_axes);

      bool has_inversion = has_inversion_center(rion2,ion_mass,nion,sym_tolerance);
      if (has_inversion)
         group_name = "D∞h";
      else 
         group_name = "C∞v";
      group_rank = -1;
   }
   else if ((((inertia_moments[0]-inertia_moments[2])/m_total) < sym_tolerance) && (nion>3))
   {
      rotation_type = "spherical";

      // align the molecular axes along the inertia_axes
      align_to_axes(rion2,nion,inertia_axes);

      determine_spherical_group(rion2,ion_mass,nion,
                                sym_tolerance,
                                group_name,group_rank,rotation_type,
                                inertia_tensor,inertia_moments,inertia_axes);
   }
   else if (((inertia_moments[0]-inertia_moments[1])/m_total) < sym_tolerance)
   {
      rotation_type = "prolate top";

      // align the molecular axes along the inertia_axes
      align_to_axes(rion2,nion,inertia_axes);
      
      determine_symmetric_group(rion2,ion_mass,nion,
                                sym_tolerance,
                                group_name,group_rank,rotation_type,
                                inertia_tensor,inertia_moments,inertia_axes);
   }
   else if (((inertia_moments[1]-inertia_moments[2])/m_total) < sym_tolerance)
   {
      rotation_type = "oblate top";

      //rotate inertia_moments and inertia_axes
      double m2 = inertia_moments[2];
      double a2 = inertia_axes[6];
      double b2 = inertia_axes[7];
      double c2 = inertia_axes[8];

      inertia_moments[2] = inertia_moments[0];
      inertia_axes[6] = inertia_axes[0];
      inertia_axes[7] = inertia_axes[1];
      inertia_axes[8] = inertia_axes[2];

      inertia_moments[0] = inertia_moments[1];
      inertia_axes[0] = inertia_axes[3];
      inertia_axes[1] = inertia_axes[4];
      inertia_axes[2] = inertia_axes[5];

      inertia_moments[1] = m2;
      inertia_axes[3] = a2;
      inertia_axes[4] = b2;
      inertia_axes[5] = c2;

      //align the molecular axes along the inertia_axes
      align_to_axes(rion2,nion,inertia_axes);
      
      determine_symmetric_group(rion2,ion_mass,nion,
                                sym_tolerance,
                                group_name,group_rank,rotation_type,
                                inertia_tensor,inertia_moments,inertia_axes);
   }
   else
   {
      rotation_type = "asymmetric";

      //rotate inertia_moments and inertia_axes
      double m2 = inertia_moments[2];
      double a2 = inertia_axes[6];
      double b2 = inertia_axes[7];
      double c2 = inertia_axes[8];

      inertia_moments[2] = inertia_moments[1];
      inertia_axes[6] = inertia_axes[3];
      inertia_axes[7] = inertia_axes[4];
      inertia_axes[8] = inertia_axes[5];

      inertia_moments[1] = inertia_moments[0];
      inertia_axes[3] = inertia_axes[0];
      inertia_axes[4] = inertia_axes[1];
      inertia_axes[5] = inertia_axes[2];

      inertia_moments[0] = m2;
      inertia_axes[0] = a2;
      inertia_axes[1] = b2;
      inertia_axes[2] = c2;
      

      // align the molecular axes along the inertia_axes
      align_to_axes(rion2,nion,inertia_axes);

      determine_asymmetric_group(rion2,ion_mass,nion,
                                 sym_tolerance,
                                 group_name,group_rank,rotation_type,
                                 inertia_tensor,inertia_moments,inertia_axes);
   }


   return;
}


}



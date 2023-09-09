
#include <iostream>
#include <string>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <cstring> 
#include "blas.h"


#define mindex(i,j) (3*i+j)

namespace pwdft {

enum MirrorPlaneType {
    SIGMA_H,  // Horizontal mirror plane
    SIGMA_V,  // Vertical mirror plane
    SIGMA_D   // Diagonal mirror plane
};



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
 *           generate_quadrupole           *
 *                                         *
 *******************************************/
/* Description:
   Calculate the quadrupole moments of a distribution of point masses in 3D space.

   Arguments:
     - rion: Array of size 3 * nion containing the positions (x, y, z) of nion point masses.
     - ion_mass: Array of size nion containing the masses of the nion point masses.
     - nion: The number of point masses.

   Output:
     - Qxx, Qyy, Qzz, Qxy, Qxz, Qyz: Variables to store the calculated quadrupole moments.

   Note:
   The function calculates the center of mass and then uses it to compute the quadrupole moments.

   Example Usage:
   double rion[3 * nion];      // Point mass positions
   double ion_mass[nion];      // Masses of point masses
   int nion;                   // Number of point masses
   double Qxx, Qyy, Qzz, Qxy, Qxz, Qyz;  // Quadrupole moments
   generate_quadrupole(rion, ion_mass, nion, Qxx, Qyy, Qzz, Qxy, Qxz, Qyz);
*/
static void generate_quadrupole(double* rion, double *ion_mass, int nion, 
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
static void eccentricity_and_asphericity(double* rion, double *ion_mass, int nion, 
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
/* Description:
     Calculate the inertia tensor and center of mass for a molecule given ion
     positions and masses.
  
   Arguments:
     - rion: Array of ion positions
     - ion_mass: Array of ion masses
     - nion: Number of ions
     - inertia_tensor: Output inertia tensor (3x3 matrix)
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

   std::cout << "tensor = " << std::endl;
   std::cout << inertia_tensor[0] << " " << inertia_tensor[3] << " " << inertia_tensor[6] << std::endl;
   std::cout << inertia_tensor[1] << " " << inertia_tensor[4] << " " << inertia_tensor[7] << std::endl;
   std::cout << inertia_tensor[2] << " " << inertia_tensor[5] << " " << inertia_tensor[8] << std::endl;

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
 *         has_inversion_center            *
 *                                         *
 *******************************************/
/* Description:
     Determine whether a molecule has an inversion center based on ion coordinates and ion masses.

   Arguments:
     - rion: Pointer to the array of ion coordinates (x, y, z)
     - ion_mass: Pointer to the array of ion masses
     - nion: Number of ions in the molecule

   Returns:
     True if an inversion center is present, False otherwise
*/
static bool has_inversion_center(const double* rion, const double* ion_mass, const int nion)
{
   double tolerance = 0.0001;  // Adjust as needed based on the precision of your data
   std::vector<bool> ion_checked(nion, false);

   auto coord2_inverse = [tolerance](const double *rion1, const double *rion2) {
        return ((std::abs(rion1[0] + rion2[0]) < tolerance) &&
                (std::abs(rion1[1] + rion2[1]) < tolerance) &&
                (std::abs(rion1[2] + rion2[2]) < tolerance));
    };

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

   bool all_ions_checked = std::all_of(ion_checked.begin(), ion_checked.end(), [](bool val) { return val; });

   return all_ions_checked;
}


/*******************************************
 *                                         *
 *              is_linear                  *
 *                                         *
 *******************************************/
/* Description:
     Determine whether a molecule is linear based on ion coordinates.

   Arguments:
     - rion: Pointer to the array of ion coordinates (x, y, z)
     - nion: Number of ions in the molecule

   Returns:
     True if the molecule is linear, False otherwise
*/
static bool is_linear(const double* rion, const int nion)
{
   double tolerance = 0.0001;  // Adjust as needed based on the precision of your data
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

         if (std::abs(cos_theta_i) - cos_theta > tolerance)
            return false;
      }
   }
   return true;
}


/*******************************************
 *                                         *
 *         determine_point_group           *
 *                                         *
 *******************************************/
void determine_point_group(const double *rion, const double *ion_mass, const int nion,
                           const double sym_tolerance,
                           std::string&  group_name, std::string& rotation_type,
                           double *inertia_tensor,  double *inertia_moments,  double *inertia_axes) 
{
   double m_total = 0.0;
   for (auto ii=0; ii<nion; ++ii) 
       m_total += ion_mass[ii];

   bool has_inversion = has_inversion_center(rion,ion_mass,nion);

   group_name = "unknown";
   rotation_type="unknown";

   generate_principle_axes(rion,ion_mass,nion,inertia_tensor,inertia_moments,inertia_axes);

   std::cout << "m_total^2=" << m_total << std::endl;
   std::cout << "point       A=" << (inertia_moments[0]*inertia_moments[0] +inertia_moments[1]*inertia_moments[1] +inertia_moments[2]*inertia_moments[2])/(m_total*m_total) << std::endl;
   std::cout << "I[0]=I[1], I[2]=0 linear     B=" << (inertia_moments[2]*inertia_moments[2])/(m_total*m_total) << std::endl;
   std::cout << "I[0]=I[1]=I[2] spherical     C=" << (inertia_moments[0]-inertia_moments[2])/m_total 
                                                  << " " << ((inertia_moments[0]-inertia_moments[2])/m_total < sym_tolerance) << std::endl;
   std::cout << "I[0]=I[1] > I[2] prolate top D=" << (inertia_moments[0]-inertia_moments[1])/m_total << std::endl;
   std::cout << "I[0] > I[1]=I[2] oblate top  E=" << (inertia_moments[1]-inertia_moments[2])/m_total << std::endl;
   std::cout << "has_inversion = " << has_inversion << std::endl;
   std::cout << "sym_tolerance = " << sym_tolerance << std::endl;

   if ( ( ( inertia_moments[0]*inertia_moments[0] 
          + inertia_moments[1]*inertia_moments[1] 
          + inertia_moments[2]*inertia_moments[2])
        /(m_total*m_total))
       < (sym_tolerance*sym_tolerance))
   {
      rotation_type = "point";
      group_name = "SO(3)";
      inertia_axes[0] = 1.0; inertia_axes[1] = 0.0; inertia_axes[2] = 0.0;
      inertia_axes[3] = 0.0; inertia_axes[4] = 1.0; inertia_axes[5] = 0.0;
      inertia_axes[6] = 0.0; inertia_axes[7] = 0.0; inertia_axes[8] = 1.0;
   }
   else if (((inertia_moments[2]*inertia_moments[2])/(m_total*m_total)) < (sym_tolerance*sym_tolerance))
   {
      rotation_type = "linear";
      if (has_inversion)
         group_name = "D∞v";
      else 
         group_name = "C∞v";
   }
   else if (((inertia_moments[0]-inertia_moments[2])/m_total) < sym_tolerance)
   {
      rotation_type = "spherical";
      //determine_spherical_group(rion,ion_mass,nion,
      //                          sym_tolerance,
      //                          group_name,rotation_type,
      //                          inertia_tensor,inertia_moments,inertia_axes);
   }
   else if (((inertia_moments[0]-inertia_moments[1])/m_total) < sym_tolerance)
   {
      rotation_type = "prolate top";
      //determine_prolate_top_group(rion,ion_mass,nion,
      //                            sym_tolerance,
      //                            group_name,rotation_type,
      //                            inertia_tensor,inertia_moments,inertia_axes);
   }
   else if (((inertia_moments[1]-inertia_moments[2])/m_total) < sym_tolerance)
   {
      rotation_type = "oblate top";
      //determine_olate_top_group(rion,ion_mass,nion,
      //                          sym_tolerance,
      //                          group_name,rotation_type,
      //                          inertia_tensor,inertia_moments,inertia_axes);
   }
   else
   {
      rotation_type = "asymmetric";
      //determine_asymmetric_group(rion,ion_mass,nion,
      //                           sym_tolerance,
      //                           group_name,rotation_type,
      //                           inertia_tensor,inertia_moments,inertia_axes);
   }


   return;
}


}



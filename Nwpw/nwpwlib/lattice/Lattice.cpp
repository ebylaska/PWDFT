/* Lattice.cpp
 *  Author - Eric Bylaska
 *
 * This class represents a lattice in three dimensions. It stores the unit cell vectors
 * and related parameters, such as the reciprocal lattice vectors, volume, and energy
 * cutoffs. The class provides methods to convert between lattice parameters and the
 * unit cell vectors, and to compute various lattice properties, such as the distance
 * between two points in the lattice and the nearest periodic image of a point.
 *
 * Author: Eric Bylaska
 * Last modified: 5/9/2023
 */

/**
 * @class Lattice
 * @brief Represents a lattice in a DFT simulation.
 */

#include "Lattice.hpp"
#include "Control2.hpp"
#include <cmath>

namespace pwdft {

/********************************
 *                              *
 *         get_cube             *
 *                              *
 ********************************/
/* This static function calculates the volume and the reciprocal lattice vectors
 * of a unit cell described by a 3x3 matrix of lattice vectors 'unita'.
 * The reciprocal lattice vectors are returned in 'unitg' as a 3x3 matrix,
 * and the volume of the unit cell is returned in 'omega'.
 *
 * The function uses the formula for the reciprocal lattice vectors derived
 * from the cross products of the lattice vectors, and scales the vectors by
 * 2*pi/volume to get the correct units.
 *
 * Parameters:
 *     unita : pointer to a 3x3 matrix of lattice vectors (input)
 *     unitg : pointer to a 3x3 matrix to store the reciprocal lattice vectors (output)
 *     omega : pointer to a variable to store the volume of the unit cell (output)
 */
/**
 * @brief Calculate the volume and reciprocal lattice vectors of a unit cell described by lattice vectors.
 *
 * This static function calculates the volume and the reciprocal lattice vectors of a unit cell described by a 3x3 matrix
 * of lattice vectors 'unita'. The reciprocal lattice vectors are returned in 'unitg' as a 3x3 matrix,
 * and the volume of the unit cell is returned in 'omega'.
 *
 * The function uses the formula for the reciprocal lattice vectors derived from the cross products of the lattice vectors
 * and scales the vectors by 2*pi/volume to get the correct units.
 *
 * @param unita Pointer to a 3x3 matrix of lattice vectors (input).
 * @param unitg Pointer to a 3x3 matrix to store the reciprocal lattice vectors (output).
 * @param omega Pointer to a variable to store the volume of the unit cell (output).
 */
static void get_cube(double *unita, double *unitg, double *omega) {
   double twopi = 8.0*std::atan(1.0);
   unitg[0] = unita[4]*unita[8] - unita[5]*unita[7];
   unitg[1] = unita[5]*unita[6] - unita[3]*unita[8];
   unitg[2] = unita[3]*unita[7] - unita[4]*unita[6];
   unitg[3] = unita[7]*unita[2] - unita[8]*unita[1];
   unitg[4] = unita[8]*unita[0] - unita[6]*unita[2];
   unitg[5] = unita[6]*unita[1] - unita[7]*unita[0];
   unitg[6] = unita[1]*unita[5] - unita[2]*unita[4];
   unitg[7] = unita[2]*unita[3] - unita[0]*unita[5];
   unitg[8] = unita[0]*unita[4] - unita[1]*unita[3];
   double volume = unita[0]*unitg[0] + unita[1]*unitg[1] + unita[2]*unitg[2];
   for (int i = 0; i < 9; ++i)
      unitg[i] *= twopi / volume;
   *omega = fabs(volume);
}

 
/********************************
 *                              *
 *            get_ub            *
 *                              *
 ********************************/
/**
 * @brief Calculate the reciprocal lattice vectors (ub) from the direct lattice vectors (unita).
 *
 * This static function calculates the reciprocal lattice vectors (ub) from the direct lattice vectors (unita).
 * The calculated reciprocal lattice vectors are stored in the `ub` array. The function also calculates
 * the volume of the unit cell and normalizes the reciprocal lattice vectors.
 *
 * @param unita Pointer to the array of direct lattice vectors (unita).
 * @param ub Pointer to the array to store the calculated reciprocal lattice vectors (ub).
 */
static void get_ub(double *unita, double *ub)
{
   ub[0] = unita[4]*unita[8] - unita[5]*unita[7];
   ub[1] = unita[5]*unita[6] - unita[3]*unita[8];
   ub[2] = unita[3]*unita[7] - unita[4]*unita[6];
   ub[3] = unita[7]*unita[2] - unita[8]*unita[1];
   ub[4] = unita[8]*unita[0] - unita[6]*unita[2];
   ub[5] = unita[6]*unita[1] - unita[7]*unita[0];
   ub[6] = unita[1]*unita[5] - unita[2]*unita[4];
   ub[7] = unita[2]*unita[3] - unita[0]*unita[5];
   ub[8] = unita[0]*unita[4] - unita[1]*unita[3];
   double volume = unita[0]*ub[0] + unita[1]*ub[1] + unita[2]*ub[2];
   for (auto i=0; i<9; ++i) 
      ub[i] /= volume;
}


/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
/* 
 * This constructor initializes the Lattice class instance using a Control2 
 * object. It sets the lattice vectors, the reciprocal lattice vectors 
 * and the unit cell volume, based on the input parameters of the Control2 
 * object. It also calculates some useful values, such as the energy cutoff, 
 * the width cutoff, and whether the calculation is periodic or not.            
 */                                                 
/**
 * @brief Constructor for the Lattice class.
 *
 * This constructor initializes the Lattice class instance using a Control2 object. It sets the lattice vectors,
 * the reciprocal lattice vectors, and the unit cell volume, based on the input parameters of the Control2 object.
 * It also calculates some useful values, such as the energy cutoff, the width cutoff, and whether the calculation is
 * periodic or not.
 *
 * @param control A Control2 object containing input parameters.
 */
Lattice::Lattice(Control2 &control) 
{
   int nx, ny, nz, nxh, nyh, nzh;
   double gx, gy, gz, gg, gg1, gg2, gg3, ecut0, wcut0;
 
   ecut0 = control.ecut();
   wcut0 = control.wcut();
   punita[0] = control.unita(0,0);
   punita[1] = control.unita(1,0);
   punita[2] = control.unita(2,0);
 
   punita[3] = control.unita(0,1);
   punita[4] = control.unita(1,1);
   punita[5] = control.unita(2,1);
 
   punita[6] = control.unita(0,2);
   punita[7] = control.unita(1,2);
   punita[8] = control.unita(2,2);
   get_cube(punita,punitg,&pomega);
   get_ub(punita,pub);
 
   nx = control.ngrid(0);
   ny = control.ngrid(1);
   nz = control.ngrid(2);
   nxh = nx/2;
   nyh = ny/2;
   nzh = nz/2;
 
   gx = punitg[0]*((double)nxh);
   gy = punitg[1]*((double)nxh);
   gz = punitg[2]*((double)nxh);
   gg1 = gx*gx + gy*gy + gz*gz;
 
   gx = punitg[3]*((double)nyh);
   gy = punitg[4]*((double)nyh);
   gz = punitg[5]*((double)nyh);
   gg2 = gx*gx + gy*gy + gz*gz;
 
   gx = punitg[6]*((double)nzh);
   gy = punitg[7]*((double)nzh);
   gz = punitg[8]*((double)nzh);
   gg3 = gx*gx + gy*gy + gz*gz;
 
   gg = gg1;
   if (gg2<gg) gg = gg2;
   if (gg3<gg) gg = gg3;
 
   pecut = 0.50 * gg;
   if (ecut0<pecut) pecut = ecut0;
   pwcut = pecut;
   if (wcut0<pwcut) pwcut = wcut0;
 
   pfast_erf = control.fast_erf();
   paperiodic = (control.version == 4);
}


/********************************
 *                              *
 *       Lattice::abc_abg       *
 *                              *
 ********************************/
/* This function calculates the lattice parameters (a, b, c) and angles (alpha, beta, gamma) 
 * from the unit cell vectors. It takes pointers to the output variables and sets their 
 * values accordingly. The function uses trigonometric equations to calculate the angles and 
 * the Euclidean distance formula to calculate the lengths of the unit cell vectors.
 * The angles are returned in degrees.
 */
/**
 * @brief Calculate lattice parameters and angles.
 *
 * This function calculates the lattice parameters (a, b, c) and angles (alpha, beta, gamma) 
 * from the unit cell vectors. The calculated angles are returned in degrees.
 *
 * @param a1 Pointer to store the lattice parameter 'a'.
 * @param b1 Pointer to store the lattice parameter 'b'.
 * @param c1 Pointer to store the lattice parameter 'c'.
 * @param alpha1 Pointer to store the angle 'alpha' in degrees.
 * @param beta1 Pointer to store the angle 'beta' in degrees.
 * @param gamma1 Pointer to store the angle 'gamma' in degrees.
 */
void Lattice::abc_abg(double *a1, double *b1, double *c1,
                      double *alpha1, double *beta1, double *gamma1) 
{
   double pi = 4.00*std::atan(1.00);
   double a = std::sqrt(std::pow(unita(0,0),2.0) + std::pow(unita(1,0),2.0) + std::pow(unita(2,0),2.0));
   double b = std::sqrt(std::pow(unita(0,1),2.0) + std::pow(unita(1,1),2.0) + std::pow(unita(2,1),2.0));
   double c = std::sqrt(std::pow(unita(0,2),2.0) + std::pow(unita(1,2),2.0) + std::pow(unita(2,2),2.0));
 
   double d2 = std::pow((unita(0,1)-unita(0,2)),2.0)
             + std::pow((unita(1,1)-unita(1,2)),2.0)
             + std::pow((unita(2,1)-unita(2,2)),2.0);
   double alpha = (b*b+c*c-d2)/(2.00*b*c);
   alpha = std::acos(alpha)*180.00/pi;
 
   d2 = std::pow((unita(0,2)-unita(0,0)),2.0)
      + std::pow((unita(1,2)-unita(1,0)),2.0)
      + std::pow((unita(2,2)-unita(2,0)),2.0);
   double beta = (c*c+a*a-d2) / (2.00*c*a);
   beta = std::acos(beta) * 180.00 / pi;
 
   d2 = std::pow((unita(0,0)-unita(0,1)),2.0)
      + std::pow((unita(1,0)-unita(1,1)),2.0)
      + std::pow((unita(2,0)-unita(2,1)),2.0);
   double gamma = (a*a+b*b-d2)/(2.00*a*b);
   gamma = std::acos(gamma)*180.00 / pi;
 
   /* return values */
   *a1 = a;
   *b1 = b;
   *c1 = c;
   *alpha1 = alpha;
   *beta1 = beta;
   *gamma1 = gamma;
}


/********************************
 *                              *
 *     Lattice::min_diff_xyz    *
 *                              *
 ********************************/
void Lattice::min_diff_xyz(double *x, double *y, double *z)
{
   if (!paperiodic)
   {
      double c0 = (*x)*ub(0,0) + (*y)*ub(1,0) + (*z)*ub(2,0);
      double c1 = (*x)*ub(0,1) + (*y)*ub(1,1) + (*z)*ub(2,1);
      double c2 = (*x)*ub(0,2) + (*y)*ub(1,2) + (*z)*ub(2,2);
      *x = unita(0,0)*c0 + unita(0,1)*c1 + unita(0,2)*c2;
      *y = unita(1,0)*c0 + unita(1,1)*c1 + unita(1,2)*c2;
      *z = unita(2,0)*c0 + unita(2,1)*c1 + unita(2,2)*c2;
   }
}

/********************************
 *                              *
 *       Lattice::min_diff      *
 *                              *
 ********************************/
/**
 * @brief Computes minimum differences in coordinates.
 * 
 * This function computes the minimum differences in coordinates using the given
 * periodic boundary conditions (if enabled).
 * 
 * @param rxyz Array representing the coordinates.
 */
void Lattice::min_diff(double *rxyz)
{
   if (!paperiodic)
   {
      double c0 = rxyz[0]*ub(0,0) + rxyz[1]*ub(1,0) + rxyz[2]*ub(2,0);
      double c1 = rxyz[0]*ub(0,1) + rxyz[1]*ub(1,1) + rxyz[2]*ub(2,1);
      double c2 = rxyz[0]*ub(0,2) + rxyz[1]*ub(1,2) + rxyz[2]*ub(2,2);
      rxyz[0] = unita(0,0)*c0 + unita(0,1)*c1 + unita(0,2)*c2;
      rxyz[1] = unita(1,0)*c0 + unita(1,1)*c1 + unita(1,2)*c2;
      rxyz[2] = unita(2,0)*c0 + unita(2,1)*c1 + unita(2,2)*c2;
   }

}

} // namespace pwdft


#include <cmath>
#include <cstring> //memset
#include <iostream>
#include "compressed_io.hpp"

#include "Parallel.hpp"
#include "Pneb.hpp"


namespace pwdft {

/**********************************************
 *                                            *
 *              nwpw_pochhammer               *
 *                                            *
 **********************************************/
/**
 * @brief Computes the Pochhammer symbol (rising factorial) for given parameters.
 * 
 * The Pochhammer symbol, also known as the rising factorial, is defined as:
 * \f$(a)_n = a(a+1)(a+2)...(a+n-1)\f$.
 * This function calculates the Pochhammer symbol for a given value of 'a' and 'n'.
 *
 * @param a The initial value for the Pochhammer symbol computation.
 * @param n The number of terms in the rising factorial.
 * @return The computed value of the Pochhammer symbol.
 */
static double nwpw_pochhammer(const double a, const int n)
{
   double fac = 1.0;
   for (auto i=0; i<n; ++i)
       fac *= (a+((double) i));

   return fac;
}


/**********************************************
 *                                            *
 *                nwpw_Qfac                   *
 *                                            *
 **********************************************/
/**
 * @brief Computes the Q-factor for given parameters.
 * 
 * This function calculates the Q-factor based on the provided parameters 'k', 'm', and 'n'.
 * The Q-factor is determined using a combination of factorials and the Pochhammer symbol.
 *
 * @param k Parameter used in the Q-factor computation.
 * @param m Parameter used in the Q-factor computation.
 * @param n Parameter used in the Q-factor computation.
 * @return The computed value of the Q-factor.
 */
static double  nwpw_Qfac(const int k, const int m, const int n)
{
   double nfac = 1.0;
   if (std::abs(m) > 0)
      nfac = 2.0*nwpw_pochhammer(1.0,n-2*k-std::abs(m)+1)/nwpw_pochhammer(1.0,n-2*k+std::abs(m)+1);
    
   double num = nwpw_pochhammer(0.5,k)*nwpw_pochhammer(1.0,n-k)*(1.0+2.0*n - 4.0*k);
   double den = nwpw_pochhammer(1.0,k)*nwpw_pochhammer(1.5,n-k);
 
   return nfac*num/den;
}


/***************************************************
 *                                                 *
 *                 nwpw_plgndrs                    *
 *                                                 *
 ***************************************************/
/**  This routine calculates all the associated Legendre
 *  polynomials thru lmax at the value x.
 *
 *  P_lm(x) = (1-x^2)^(m/2) * d^m/dx^m Pl(x)
 *
 *  and
 *
 *  P_l(x) = 1/(2^l*l!) * d^l/dx^l (x^2-1)^l
 *
 *  In this routine the following relationships are used to calculate P_lm(x)
 *
 *   P_mm(x)    = (2*m-1)!! (1-x^2)^(m/2)
 *   P_m+1,m(x) = x*(2*m+1)*P_mm(x)
 *
 * and
 *
 *    (l-m)*P_lm(x) = x*(2*l-1)*P_(l-1),m(x) - (l+m-1)*P_(l-2),m(x)
 *
 * @brief Computes the associated Legendre polynomials.
 *
 * This function calculates the associated Legendre polynomials for a given value of x and up to a specified maximum degree (lmax).
 * The results are stored in the provided Plm vector.
 *
 * @param lmax The maximum degree of the polynomial.
 * @param x The point at which the polynomial is evaluated.
 * @param Plm A vector to store the computed polynomial values.
 *
 * @note The size of the Plm vector should be appropriately initialized before calling this function.
 */
static void nwpw_plgndrs(const int lmax, const double x, double *Plm)
{
   // Local variables
   int i, l, m, lm, lm0, lm1, lm2, lmi, lmi0, lmi1, lmi2;
   double fact, pmm, somx2;

   // Compute P_mm(x) thru lmax
   somx2 = std::sqrt((1.0 - x) * (1.0 + x));
   pmm = 1.0;
   fact = 1.0;
   lm = 2;
   Plm[0] = 1.0;
   for (l = 1; l <= lmax; ++l)
   {
       pmm = pmm * fact * somx2;
       fact += 2.0;
       Plm[lm] = pmm;
       lm += (l + 2);
   }

   // Compute P_(m+1),m(x) thru lmax
   lm = 0;
   lm2 = 1;
   for (m = 0; m <= lmax - 1; ++m)
   {
       Plm[lm2] = x * (2 * m + 1) * Plm[lm];
       lm += (m + 2);
       lm2 += (m + 3);
   }

   // Compute the rest of the P_lm's
   lmi0 = 0;
   lmi1 = 1;
   lmi = 3;
   for (m = 0; m <= lmax - 2; ++m)
   {
       lm0 = lmi0;
       lm1 = lmi1;
       lm = lmi;
       for (l = m + 2; l <= lmax; ++l)
       {
           Plm[lm] = (x * static_cast<double>(2 * l - 1) * Plm[lm1]
                      - static_cast<double>(l + m - 1) * Plm[lm0])
                     / static_cast<double>(l - m);
           lm0 = lm1;
           lm1 = lm;
           lm += (l + 1);
       }
       lmi0 += (m + 2);
       lmi1 += (m + 3);
       lmi += (m + 4);
   }
}


/***************************************************
 *                                                 *
 *                 nwpw_Tesserals                  *
 *                                                 *
 ***************************************************/
/**  this routine calculates Tesseral harmonics thru lmax, where
 *
 *  Tlm(theta,phi) =                            Pl|m|(cos(theta))            for m =0
 *  Tlm(theta,phi) = sqrt(2*(l-|m|)!/(l+|m|)) * Pl|m|(cos(theta))*cos(m*phi) for m> 0
 *  Tlm(theta,phi) = sqrt(2*(l-|m|)!/(l+|m|)) * Pl|m|(cos(theta))*sin(m*phi) for m<0
 *
 *  Tlm(theta,phi) = Pl|m|(cos(theta))            for m =0
 *  Tlm(theta,phi) = Pl|m|(cos(theta))*cos(m*phi) for m> 0
 *  Tlm(theta,phi) = Pl|m|(cos(theta))*sin(m*phi) for m<0
 *
 *  Plm(x) = (1-x^2)^(m/2) d^m/dx^m Pl(x)
 *
 *  where -l<= m <= l  and -1 <= x <= 1
 *
 *  Note that the (-1)^m factor is not included
 *
 * @brief Computes the Tesseral harmonics up to a specified maximum degree (lmax).
 *
 * This function calculates the Tesseral harmonics for given values of theta and phi, up to a specified maximum degree (lmax).
 * The Tesseral harmonics are defined as:
 * - \( T_{lm}(\theta,\phi) = P_{l|m|}(\cos(\theta)) \) for \( m = 0 \)
 * - \( T_{lm}(\theta,\phi) = \sqrt{\frac{2(l-|m|)!}{(l+|m|)!}} * P_{l|m|}(\cos(\theta))\cos(m\phi) \) for \( m > 0 \)
 * - \( T_{lm}(\theta,\phi) = \sqrt{\frac{2(l-|m|)!}{(l+|m|)!}} * P_{l|m|}(\cos(\theta))\sin(m\phi) \) for \( m < 0 \)
 *
 * Where:
 * \( P_{lm}(x) \) is the associated Legendre polynomial defined as:
 * \( P_{lm}(x) = (1-x^2)^{\frac{m}{2}} \frac{d^m}{dx^m} P_l(x) \)
 *
 * @param lmax The maximum degree of the polynomial.
 * @param theta The polar angle in spherical coordinates.
 * @param phi The azimuthal angle in spherical coordinates.
 * @param Tlm A matrix to store the computed Tesseral harmonics values.
 *
 * @note The factor \((-1)^m\) is not included in the calculations.
 */
static void nwpw_Tesserals(const int lmax, const double costheta, const double cosphi, const double sinphi, 
                           double *Plm, double *Tlm) 
{
   int l, m, lm, lmc, lms, lmi, lmci, lmsi;
   double cs0, cs1, cs2, sn0, sn1, sn2, fact;

   nwpw_plgndrs(lmax, costheta, Plm);

   // m==0 copy
   lm = 0;
   lmc = 0;
   for (l = 0; l <= lmax; l++) 
   {
      Tlm[lmc] = Plm[lm];
      lm += (l + 1);
      lmc += (2 * l + 2);
   }

   // m>0 copy
   cs0 = 0.0;
   sn0 = 0.0;
   cs1 = 1.0;
   sn1 = 0.0;
   cs2 = cosphi;
   sn2 = sinphi;
   lmi = 2;
   lmci = 3;
   lmsi = 1;
   for (m = 1; m <= l; m++) 
   {
      lm = lmi;
      lmc = lmci;
      lms = lmsi;
      for (l = m; l <= lmax; l++) 
      {
         Tlm[lmc] = cs2 * Plm[lm];
         Tlm[lms] = sn2 * Plm[lm];
         lm += (l + 1);
         lmc += (2 * l + 2);
         lms += (2 * l + 2);
      }
      cs0 = cs1;
      cs1 = cs2;
      cs2 = (2.0 * cosphi) * cs1 - cs0;
      sn0 = sn1;
      sn1 = sn2;
      sn2 = (2.0 * cosphi) * sn1 - sn0;
      lmi += (m + 2);
      lmci += (2 * m + 3);
      lmsi += (2 * m + 1);
   }
}


/**********************************************
 *                                            *
 *            f3_to_theta_min                 *
 *                                            *
 **********************************************/
/**
 * @brief Calculates the value of theta_min based on the input value f3.
 * 
 * This function determines the minimum theta value based on the provided f3 value.
 * 
 * @param f3 A double value representing the input for which theta_min is calculated.
 * @return A double value representing the minimum theta based on f3.
 */
static double f3_to_theta_min(const double f3) 
{
   const double pi = 4.0 * std::atan(1.0);
   double theta_min = 0.0;

   if (std::abs(f3) < 1.0e-6) {
       theta_min = 0.0;
   } else if (f3 > 0.0) {
       theta_min = 0.0;
   } else {
       theta_min = 0.5 * pi;
   }

   return theta_min;
}


/**********************************************
 *                                            *
 *            f3_to_theta_max                 *
 *                                            *
 **********************************************/
/**
 * @brief Calculates the value of theta_max based on the input value f3.
 * 
 * This function determines the maximum theta value based on the provided f3 value.
 * 
 * @param f3 A double value representing the input for which theta_max is calculated.
 * @return A double value representing the maximum theta based on f3.
 */
static double f3_to_theta_max(const double f3) 
{
   const double pi = 4.0 * std::atan(1.0);
   double theta_max = 0.0;

   if (std::abs(f3) < 1.0e-6) {
       theta_max = pi;
   } else if (f3 > 0.0) {
       theta_max = 0.5 * pi;
   } else {
       theta_max = pi;
   }

   return theta_max;
}


/**********************************************
 *                                            *
 *            f1f2_to_phi_min                 *
 *                                            *
 **********************************************/
/**
 * @brief Calculates the value of phi_min based on the input values f1 and f2.
 * 
 * This function determines the minimum phi value based on the provided f1 and f2 values.
 * 
 * @param f1 A double value representing the first input for which phi_min is calculated.
 * @param f2 A double value representing the second input for which phi_min is calculated.
 * @return A double value representing the minimum phi based on f1 and f2.
 */
static double f1f2_to_phi_min(double f1, double f2) 
{
   const double pi = 4.0 * std::atan(1.0);
   double phi_min = 0.0;

   if (std::abs(f1) < 1.0e-6 && std::abs(f2) < 1.0e-6) {
      phi_min = 0.0;
   } else if (f1 > 0.0 && std::abs(f2) < 1.0e-6) {
      phi_min = -0.5 * pi;
   } else if (f1 < 0.0 && std::abs(f2) < 1.0e-6) {
      phi_min = 0.5 * pi;
   } else if (std::abs(f1) < 1.0e-6 && f2 > 0.0) {
      phi_min = 0.0;
   } else if (f1 > 0.0 && f2 > 0.0) {
      phi_min = 0.0;
   } else if (f1 < 0.0 && f2 > 0.0) {
      phi_min = 0.5 * pi;
   } else if (std::abs(f1) < 1.0e-6 && f2 < 0.0) {
      phi_min = pi;
   } else if (f1 > 0.0 && f2 < 0.0) {
      phi_min = 1.5 * pi;
   } else if (f1 < 0.0 && f2 < 0.0) {
      phi_min = pi;
   }

   return phi_min;
}


/**********************************************
 *                                            *
 *             f1f2_to_phi_max                *
 *                                            *
 **********************************************/
/**
 * @brief Calculates the value of phi_max based on the input values f1 and f2.
 * 
 * This function determines the maximum phi value based on the provided f1 and f2 values.
 * 
 * @param f1 A double value representing the first input for which phi_max is calculated.
 * @param f2 A double value representing the second input for which phi_max is calculated.
 * @return A double value representing the maximum phi based on f1 and f2.
 */
static double f1f2_to_phi_max(const double f1, const double f2) 
{
   const double pi = 4.0 * std::atan(1.0);
   double phi_max = 0.0;

   if (std::abs(f1) < 1.0e-6 && std::abs(f2) < 1.0e-6) {
       phi_max = 2.0 * pi;
   } else if (f1 > 0.0 && std::abs(f2) < 1.0e-6) {
       phi_max = 0.5 * pi;
   } else if (f1 < 0.0 && std::abs(f2) < 1.0e-6) {
       phi_max = 1.5 * pi;
   } else if (std::abs(f1) < 1.0e-6 && f2 > 0.0) {
       phi_max = pi;
   } else if (f1 > 0.0 && f2 > 0.0) {
       phi_max = 0.5 * pi;
   } else if (f1 < 0.0 && f2 > 0.0) {
       phi_max = pi;
   } else if (std::abs(f1) < 1.0e-6 && f2 < 0.0) {
       phi_max = 2.0 * pi;
   } else if (f1 > 0.0 && f2 < 0.0) {
       phi_max = 2.0 * pi;
   } else if (f1 < 0.0 && f2 < 0.0) {
       phi_max = 1.5 * pi;
   }

   return phi_max;
}


/**********************************************
 *                                            *
 * setup_integrate_oneoverG2_BZ_Gengenbauer2  *
 *                                            *
 **********************************************/
/**
 * @brief Sets up the integration for one over G^2 in the Brillouin Zone using Gegenbauer polynomials of the second kind.
 * 
 * This function prepares the necessary data structures and performs initial computations for the integration
 * of one over G^2 in the Brillouin Zone using Gegenbauer polynomials of the second kind.
 * 
 * @param Nt3 The total number of grid points in the 3D grid.
 * @param gammasdxdydz Pointer to an array containing the gamma values for each grid point.
 * @param lk Pointer to an array containing the lk values for each grid point.
 * @param nmax Maximum value of n for the Gegenbauer polynomials.
 * @param Mpack Pointer to an array where the Mpack values will be stored.
 * @param QQ Pointer to an array where the QQ values will be stored.
 * @param Plm Pointer to an array where the Plm values will be stored.
 * @param Tlm Pointer to an array where the Tlm values will be stored.
 */
static void setup_integrate_oneoverG2_BZ_Gengenbauer2(const int Nt3, const double *gammasdxdydz, const double *lk, 
                                                      const int nmax, double *Mpack, double *QQ, double *Plm, double *Tlm) 
{
   int indx = 0;
   for (int n = 0; n <= nmax; ++n) 
   {
      for (int k = 0; k <= n/2; ++k) 
      {
         for (int m = -n + 2*k; m <= n - 2*k; ++m) 
         {
            Mpack[indx] = 0.0;
            QQ[indx] = nwpw_Qfac(k, m, n);
            indx++;
         }
      }
   }

   for (int i = 0; i < Nt3; ++i) 
   {
      double gx = lk[3*i];
      double gy = lk[3*i + 1];
      double gz = lk[3*i + 2];
      double g = std::sqrt(gx*gx + gy*gy + gz*gz);
      double gxy = std::sqrt(gx*gx + gy*gy);
      double costheta, cosphi, sinphi;
      if (g > 1.0e-9) 
      {
         costheta = gz/g;
         if (gxy < 1.0e-9) 
         {
            cosphi = 1.0;
            sinphi = 0.0;
         } 
         else 
         {
            cosphi = gx/gxy;
            sinphi = gy/gxy;
         }

         nwpw_Tesserals(nmax, costheta, cosphi, sinphi, Plm, Tlm);

         indx = 0;
         double gn = 1.0;
         for (int n = 0; n <= nmax; ++n) 
         {
            for (int k = 0; k <= n/2; ++k) 
            {
               int l = n - 2*k;
               for (int m = -n + 2*k; m <= n - 2*k; ++m) 
               {
                  int lm = (l+1)*(l+1) - l + m - 1;
                  Mpack[indx] += QQ[indx] * gn * Tlm[lm] * gammasdxdydz[i];
                  indx++;
               }
            }
            gn *= g;
         }
      }
   }
}


/**********************************************
 *                                            *
 *    integrate_oneoverG2_BZ_Gengenbauer2     *
 *                                            *
 **********************************************/
/**
 * @brief Integrates one over G^2 in the Brillouin Zone using Gegenbauer polynomials of the second kind.
 * 
 * This function performs the integration of one over G^2 in the Brillouin Zone using Gegenbauer 
 * polynomials of the second kind. It uses precomputed values and the provided G values to compute the integral.
 * 
 * @param nmax Maximum value of n for the Gegenbauer polynomials.
 * @param Mpack Pointer to an array containing the precomputed Mpack values.
 * @param gx G value in the x direction.
 * @param gy G value in the y direction.
 * @param gz G value in the z direction.
 * @param Plm Pointer to an array where the Plm values will be stored.
 * @param Tlm Pointer to an array where the Tlm values will be stored.
 * @return The result of the integration.
 */
static double integrate_oneoverG2_BZ_Gengenbauer2(const int nmax, const double *Mpack, 
                                                  const double gx, const double gy, const double gz, 
                                                  double *Plm, double *Tlm) 
{
   double g = std::sqrt(gx*gx + gy*gy + gz*gz);
   double gxy = std::sqrt(gx*gx + gy*gy);
   double gn = g * g;

   double costheta = gz / g;
   double cosphi, sinphi;
   if (gxy < 1.0e-9) 
   {
      cosphi = 1.0;
      sinphi = 0.0;
   } 
   else 
   {
      cosphi = gx / gxy;
      sinphi = gy / gxy;
   }

   nwpw_Tesserals(nmax, costheta, cosphi, sinphi, Plm, Tlm);

   int indx = 0;
   double sum = 0.0;
   for (int n = 0; n <= nmax; ++n) 
   {
      double sum0 = 0.0;
      for (int k = 0; k <= n/2; ++k) 
      {
         int l = n - 2*k;
         for (int m = -n + 2*k; m <= n - 2*k; ++m) 
         {
            int lm = (l+1)*(l+1) - l + m - 1;
            sum0 += Mpack[indx] * Tlm[lm];
            indx++;
         }
      }
      sum += sum0 / gn;
      gn *= g;
   }

   return sum;
}


/**********************************************
 *                                            *
 *    integrate_spherical_oneoveralphaG2_BZ   *
 *                                            *
 **********************************************/
/**
 * @brief Integrates one over alpha*G^2 in the Brillouin Zone using spherical coordinates.
 * 
 * This function performs the integration of one over alpha*G^2 in the Brillouin Zone using 
 * spherical coordinates. It takes into account the factors af1, af2, and af3 for the gamma 
 * directions and uses the provided parameters to compute the integral.
 * 
 * @param af1 Factor for the gamma1 direction.
 * @param af2 Factor for the gamma2 direction.
 * @param af3 Factor for the gamma3 direction.
 * @param Nr Number of radial divisions.
 * @param r0 Starting value of radial coordinate.
 * @param r1 Ending value of radial coordinate.
 * @param Ntheta Number of theta divisions.
 * @param theta0 Starting value of theta angle.
 * @param theta1 Ending value of theta angle.
 * @param Nphi Number of phi divisions.
 * @param phi0 Starting value of phi angle.
 * @param phi1 Ending value of phi angle.
 * @param unitg Pointer to an array containing the unit cell vectors.
 * @param alpha Scaling factor for the G^2 term.
 * @return The result of the integration.
 */
static double integrate_spherical_oneoveralphaG2_BZ(const double af1, const double af2, const double af3,
                                                    const int Nr, const double r0, const double r1,
                                                    const int Ntheta, const double theta0, const double theta1,
                                                    const int Nphi, const double phi0, const double phi1,
                                                    const double *unitg, const double alpha) 
{
   const double fourpi = 16.0 * std::atan(1.0);

   double dr = (r1 - r0) / static_cast<double>(Nr);
   double dtheta = (theta1 - theta0) / static_cast<double>(Ntheta);
   double dphi = (phi1 - phi0) / static_cast<double>(Nphi);

   double asum = 0.0;
   for (int iphi = 0; iphi <= Nphi; ++iphi) 
   {
      for (int itheta = 0; itheta <= Ntheta; ++itheta) 
      {
         for (int ir = 0; ir <= Nr; ++ir) 
         {
            double r = r0 + ir * dr;
            double theta = theta0 + itheta * dtheta;
            double phi = phi0 + iphi * dphi;
            double gamma1 = r * std::cos(phi) * std::sin(theta);
            double gamma2 = r * std::sin(phi) * std::sin(theta);
            double gamma3 = r * std::cos(theta);

            double tgamma1 = std::cos(phi) * std::sin(theta);
            double tgamma2 = std::sin(phi) * std::sin(theta);
            double tgamma3 = std::cos(theta);
            double tlkx = tgamma1 * unitg[0] + tgamma2 * unitg[3] + tgamma3 * unitg[6];
            double tlky = tgamma1 * unitg[1] + tgamma2 * unitg[4] + tgamma3 * unitg[7];
            double tlkz = tgamma1 * unitg[2] + tgamma2 * unitg[5] + tgamma3 * unitg[8];
            double tR2 = tlkx * tlkx + tlky * tlky + tlkz * tlkz;

            double dx = (ir % 2 == 0) ? 2.0 * dr / 3.0 : 4.0 * dr / 3.0;
            double dy = (itheta % 2 == 0) ? 2.0 * dtheta / 3.0 : 4.0 * dtheta / 3.0;
            double dz = (iphi % 2 == 0) ? 2.0 * dphi / 3.0 : 4.0 * dphi / 3.0;

            if (ir == 0 || ir == Nr) dx = dr / 3.0;
            if (itheta == 0 || itheta == Ntheta) dy = dtheta / 3.0;
            if (iphi == 0 || iphi == Nphi) dz = dphi / 3.0;

            double w = ((1.0 - af1) * (1.0 - std::abs(gamma1)) + af1 * std::abs(gamma1))
                     * ((1.0 - af2) * (1.0 - std::abs(gamma2)) + af2 * std::abs(gamma2))
                     * ((1.0 - af3) * (1.0 - std::abs(gamma3)) + af3 * std::abs(gamma3));

            asum += fourpi * (std::exp(-alpha * r * r * tR2) / tR2) * w * std::sin(theta) * dx * dy * dz;
         }
      }
   }

   return asum;
}


/**********************************************
 *                                            *
 *        integrate_oneoveralphaG2_BZ         *
 *                                            *
 **********************************************/
/**
 * @brief Integrates one over alpha*G^2 in the Brillouin Zone.
 * 
 * This function performs the integration of one over alpha*G^2 in the Brillouin Zone. 
 * It computes the sum over the given gamma points and lattice vectors, taking into account 
 * the provided alpha scaling factor and the G vector components gx, gy, and gz.
 * 
 * @param Nt3 Total number of gamma points.
 * @param gammasdxdydz Pointer to an array containing the gamma point values.
 * @param lk Pointer to an array containing the lattice vectors.
 * @param alpha Scaling factor for the G^2 term.
 * @param gx G vector x-component.
 * @param gy G vector y-component.
 * @param gz G vector z-component.
 * @return The result of the integration.
 */
static double integrate_oneoveralphaG2_BZ(const int Nt3, const double *gammasdxdydz, const double *lk, 
                                          const double alpha, const double gx, const double gy, const double gz) 
{
   double sum = 0.0;
   for (int i = 0; i < Nt3; ++i) 
   {
      double glkx = gx + lk[3*i];
      double glky = gy + lk[3*i + 1];
      double glkz = gz + lk[3*i + 2];
      double gg = glkx * glkx + glky * glky + glkz * glkz;
      if (gg < 1.0e-9) 
      {
         sum += gammasdxdydz[i] * alpha;
      } 
      else 
      {
         sum += gammasdxdydz[i] * (1.0 - std::exp(-alpha * gg)) / gg;
      }
   }
   return sum;
}


/**********************************************
 *                                            *
 *           coulomb_filter_setup             *
 *                                            *
 **********************************************/
/**
 * @brief Sets up the Coulomb filter for a given set of parameters.
 * 
 * This function prepares the Coulomb filter by computing the gamma points and lattice vectors 
 * based on the provided parameters. The gamma points and lattice vectors are stored in the 
 * gammasdxdydz and lk arrays, respectively.
 * 
 * @param N1 Number of divisions along the first dimension.
 * @param a1 Lower bound for the first dimension.
 * @param b1 Upper bound for the first dimension.
 * @param N2 Number of divisions along the second dimension.
 * @param a2 Lower bound for the second dimension.
 * @param b2 Upper bound for the second dimension.
 * @param N3 Number of divisions along the third dimension.
 * @param a3 Lower bound for the third dimension.
 * @param b3 Upper bound for the third dimension.
 * @param unitg Pointer to an array containing the unit vectors.
 * @param gammasdxdydz Pointer to an array where the gamma point values will be stored.
 * @param lk Pointer to an array where the lattice vectors will be stored.
 */
static void coulomb_filter_setup(const int N1, const double a1, const double b1,
                                 const int N2, const double a2, const double b2,
                                 const int N3, const double a3, const double b3,
                                 const double *unitg,
                                 double *gammasdxdydz,
                                 double *lk) 
{
   const double fourpi = 16.0 * std::atan(1.0);

   double dgamma1 = (b1 - a1) / static_cast<double>(N1);
   double dgamma2 = (b2 - a2) / static_cast<double>(N2);
   double dgamma3 = (b3 - a3) / static_cast<double>(N3);

   int index = 0;
   for (int i3 = 0; i3 <= N3; ++i3) 
   {
      for (int i2 = 0; i2 <= N2; ++i2) 
      {
         for (int i1 = 0; i1 <= N1; ++i1) 
         {
            double gamma1 = a1 + i1 * dgamma1;
            double gamma2 = a2 + i2 * dgamma2;
            double gamma3 = a3 + i3 * dgamma3;

            double dx = (i1 % 2 == 0) ? 2.0 * dgamma1 / 3.0 : 4.0 * dgamma1 / 3.0;
            double dy = (i2 % 2 == 0) ? 2.0 * dgamma2 / 3.0 : 4.0 * dgamma2 / 3.0;
            double dz = (i3 % 2 == 0) ? 2.0 * dgamma3 / 3.0 : 4.0 * dgamma3 / 3.0;

            if (i1 == 0 || i1 == N1) dx = dgamma1 / 3.0;
            if (i2 == 0 || i2 == N2) dy = dgamma2 / 3.0;
            if (i3 == 0 || i3 == N3) dz = dgamma3 / 3.0;

            gammasdxdydz[index] = fourpi * (1.0 - std::abs(gamma1))
                                          * (1.0 - std::abs(gamma2))
                                          * (1.0 - std::abs(gamma3))
                                          * dx * dy * dz;

            lk[3*index]   = gamma1 * unitg[0] + gamma2 * unitg[3] + gamma3 * unitg[6];
            lk[3*index+1] = gamma1 * unitg[1] + gamma2 * unitg[4] + gamma3 * unitg[7];
            lk[3*index+2] = gamma1 * unitg[2] + gamma2 * unitg[5] + gamma3 * unitg[8];

            ++index;
         }
      }
   }
}


/**********************************************
 *                                            *
 *              coulomb_filter                *
 *                                            *
 **********************************************/
/**
 *  This code computes the filtered coulomb potential using
 *  filon integration strategy given in "An Accurate Integration 
 *  Strategy for Calculating Exact Exchange in Periodic Boundary 
 *  Conditions: A Plane-Wave DFT Implementation."  
 *
 *  The I1,I2,and I3 integrals are computed using Simpson integration.
 *
 *  Warning: While this routine this code was lifted from NWChem and 
 * it should be good enough for testing, it probably should  be 
 * restructured to improve parallelization and the use of memory. 
 *
 * @brief Computes the filtered coulomb potential using the filon integration strategy.
 * 
 * This function computes the filtered coulomb potential using the filon 
 * integration strategy as described in "An Accurate Integration Strategy 
 * for Calculating Exact Exchange in Periodic Boundary Conditions: A Plane-Wave 
 * DFT Implementation." The I1, I2, and I3 integrals are computed using 
 * Simpson integration.
 *
 * @warning This code was lifted from NWChem. While it should be good enough for testing, 
 * it might need restructuring to improve parallelization and memory usage.
 *
 * @param mygrid Pointer to the grid data structure.
 * @param vfilter Pointer to the filtered potential array.
 * @param filename Name of the file containing kernel data.
 */
void coulomb_filter(Pneb *mygrid, double *vfilter, const std::string filename)
{
   Parallel *myparall = mygrid->d3db::parall;
   int npack0 = mygrid->npack(0);

   int Nt1 = 228, Nt2 = 228, Nt3 = 228;
   int Ntt3 = (Nt1 + 1) * (Nt2 + 1) * (Nt3 + 1);
   int Ns1 = 432, Ns2 = 432, Ns3 = 432;
   int ifound=0;
   bool found_filename;

   int tngrid[3];
   double tunita[9],tunitg[9];
   double tol = 1.0e-9;
   double *tmp2 = new (std::nothrow) double[mygrid->nfft3d]();


   // check for kernel file 
   if (myparall->is_master()) {ifound = cfileexists(const_cast<char *>(filename.data())); }
   myparall->Brdcst_iValue(0, 0, &ifound);
   found_filename = (ifound!=0);

   if (found_filename)
   {
      ifound = 0;
      if (myparall->is_master()) 
      {
         openfile(5, const_cast<char *>(filename.data()),"r");
         iread(5,tngrid,3);
         dread(5,tunita,9);
         dread(5,tunitg,9);
         closefile(5);
         found_filename = (mygrid->nx == tngrid[0]) &&
                          (mygrid->ny == tngrid[1]) &&
                          (mygrid->nz == tngrid[2]);
         for (auto i=0; i<9; ++i)
            found_filename = found_filename && (std::fabs(mygrid->lattice->unita1d(i) - tunita[i]) < tol) 
                                            && (std::fabs(mygrid->lattice->unitg1d(i) - tunitg[i]) < tol);
      }
      if (found_filename) ifound = 1;
      myparall->Brdcst_iValue(0, 0, &ifound);
      found_filename = (ifound!=0);
   }

   // read kernel from file
   if (found_filename)
   {
      if (myparall->is_master()) 
      {
         openfile(5, const_cast<char *>(filename.data()),"r");
         iread(5,tngrid,3);
         dread(5,tunita,9);
         dread(5,tunitg,9);
      }

      // read in vfilter 3d block
      mygrid->t_read(5, tmp2, -1);
      mygrid->t_pack(0, tmp2);
      mygrid->tt_pack_copy(0, tmp2, vfilter);

      if (myparall->is_master()) 
         closefile(5);

   }

   // generating kernel     
   else
   {
      double twopi = 8.0*std::atan(1.0);
      double oneovertwopi = 1.0/twopi;
      double scut,scut1,rcut,salpha,alpha;
      double *Gx = mygrid->Gpackxyz(0,0);
      double *Gy = mygrid->Gpackxyz(0,1);
      double *Gz = mygrid->Gpackxyz(0,2);

      tngrid[0] = mygrid->nx;
      tngrid[1] = mygrid->ny;
      tngrid[2] = mygrid->nz;
      for (auto i=0; i<9; ++i)
      {
         tunita[i] = mygrid->lattice->unita1d(i);
         tunitg[i] = mygrid->lattice->unitg1d(i);
      }

      int nmax = 10;
      double *gammasdxdydz = new (std::nothrow) double[3*Ntt3];
      double *lk           = new (std::nothrow) double[3*Ntt3];
      double *Mpack        = new (std::nothrow) double[5*(nmax+1)*(nmax+1)];
      double *QQ           = new (std::nothrow) double[5*(nmax+1)*(nmax+1)];
      double *Plm          = new (std::nothrow) double[5*(nmax+1)*(nmax+1)];
      double *Tlm          = new (std::nothrow) double[5*(nmax+1)*(nmax+1)];


      coulomb_filter_setup(Nt1,-1.0,1.0,
                           Nt2,-1.0,1.0,
                           Nt3,-1.0,1.0,
                           tunitg,gammasdxdydz,lk);
                           
      setup_integrate_oneoverG2_BZ_Gengenbauer2(Ntt3,gammasdxdydz,lk,
                                                nmax,Mpack,QQ,Plm,Tlm);

      scut  = std::sqrt(tunitg[0]*tunitg[0]+tunitg[1]*tunitg[1]+tunitg[2]*tunitg[2]);
      scut1 = std::sqrt(tunitg[3]*tunitg[3]+tunitg[4]*tunitg[4]+tunitg[5]*tunitg[5]);
      if (scut1 < scut) scut = scut1;
      scut1 = std::sqrt(tunitg[6]*tunitg[6]+tunitg[7]*tunitg[7]+tunitg[8]*tunitg[8]);
      if (scut1 < scut) scut = scut1;
      rcut = 0.3*scut;
      salpha = 1.0/rcut;
      alpha  = salpha*salpha;

      for (auto k=0; k<npack0; ++k)
      {
         double f1 = oneovertwopi*(tunita[0]*Gx[k] + tunita[3]*Gy[k] + tunita[6]*Gz[k]);
         double f2 = oneovertwopi*(tunita[1]*Gx[k] + tunita[4]*Gy[k] + tunita[7]*Gz[k]);
         double f3 = oneovertwopi*(tunita[2]*Gx[k] + tunita[5]*Gy[k] + tunita[8]*Gz[k]);

         if ((std::abs(f1) < 2.0) && (std::abs(f2) < 2.0) && (std::abs(f2) < 2.0))
         {
            vfilter[k]  = integrate_oneoveralphaG2_BZ(Ntt3,gammasdxdydz,lk,alpha,Gx[k],Gy[k],Gz[k]);
            vfilter[k] += integrate_spherical_oneoveralphaG2_BZ(
                                   std::abs(f1),std::abs(f2),std::abs(f3),
                                   Ns1,0.0,1.0,
                                   Ns2,f3_to_theta_min(f3),f3_to_theta_max(f3),
                                   Ns3,f1f2_to_phi_min(f1,f2),
                                   f1f2_to_phi_max(f1,f2),
                                   tunitg,alpha);
         }
         else
            vfilter[k] = integrate_oneoverG2_BZ_Gengenbauer2(nmax,Mpack,Gx[k],Gy[k],Gz[k],Plm,Tlm);
      }
      delete [] Tlm;
      delete [] Plm;
      delete [] QQ;
      delete [] Mpack;
      delete [] lk;
      delete [] gammasdxdydz;


      // save kernel to file
      tngrid[0] = mygrid->nx;
      tngrid[1] = mygrid->ny;
      tngrid[2] = mygrid->nz;
      for (auto i=0; i<9; ++i)
      {
         tunita[i] = mygrid->lattice->unita1d(i);
         tunitg[i] = mygrid->lattice->unitg1d(i);
      }
      if (myparall->is_master()) 
      {
         openfile(6, const_cast<char *>(filename.data()),"w");
         iwrite(5,tngrid,3);
         dwrite(5,tunita,9);
         dwrite(5,tunitg,9);
      }

      // write vfilter 3d block
      mygrid->tt_pack_copy(0, vfilter,tmp2);
      mygrid->t_unpack(0, tmp2);
      mygrid->t_write(6, tmp2, 0);

      if (myparall->is_master()) 
         closefile(5);
   }
   delete[] tmp2;
}


} // namespace pwdft

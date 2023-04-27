
#include "util_log_integrate.hpp"
#include "Parallel.hpp"
#include "util_gamma.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>

namespace pwdft {

/*************************************
 *                                    *
 *        util_double_factorial       *
 *                                    *
 **************************************/
double util_double_factorial(const int n) {
  int n11[18] = {1,     1,     1,      2,      3,       8,
                 15,    48,    105,    384,    945,     3840,
                 10395, 46080, 135135, 645120, 2027025, 10321920};

  if ((n >= -1) && (n <= 16)) {
    return ((double)n11[n + 1]);
  } else {
    std::cout << "too big parameter in nwpw_double_factorial" << std::endl;
    return -999999.0;
  }
}

/*************************************************
 *                                               *
 *         util_compcharge_gen_rgaussian         *
 *                                               *
 *************************************************/
void util_compcharge_gen_rgaussian(const int l, const double sigma,
                                   const int nr, const double r[],
                                   double gl[]) {
  double pi = 4.00 * atan(1.00);
  double c =
      pow(2.00, (l + 2)) /
      (sqrt(pi) * (util_double_factorial(2 * l + 1)) * pow(sigma, (2 * l + 3)));

  // *** this fixes possible underflow error ***
  for (auto i = 0; i < nr; ++i)
    gl[i] = 0.00;

  for (auto i = 0; i < nr; ++i)
    if (std::abs(r[i]) < (8.00 * sigma))
      gl[i] = c * pow(r[i], l) * exp(-(r[i] / sigma) * (r[i] / sigma));
}

/*******************************************
 *                                         *
 *         util_log_integrate_def          *
 *                                         *
 *******************************************/
double util_log_integrate_def(const int power_f, const double f[],
                              const int power_r, const double r[],
                              const double log_amesh, const int nrange) {
  double integrand[nrange];
  // double *integrand = new double[nrange];

  for (auto k = 0; k < nrange; ++k)
    integrand[k] = f[k] * pow(r[k], power_r + 1);

  // *** integrate from the origin to the first point ***
  double sum_f = integrand[0] / ((double)(power_r + power_f + 1));

  // *** the rest via trapesoidal rule ***
  double tmp_sum = 0.0;
  for (auto k = 0; k < nrange; ++k)
    tmp_sum += integrand[k];

  // *** the rest via trapesoidal rule ***
  sum_f += log_amesh * tmp_sum -
           0.50 * log_amesh * (integrand[0] + integrand[nrange - 1]);

  // delete [] integrand;

  return sum_f;
}

/************************************************
 *                                              *
 *            log_integrate_def0                *
 *                                              *
 ************************************************/
/*
     Computes the following integral


          /               /
          | f(r)/r dr =   | f(i) * log(a) di
          /               /

     where r(i) = r0*a**i   and f(r-->0) = r**power_f
*/

double util_log_integrate_def0(const int power_f, const double f[],
                               const double r[], const double log_amesh,
                               const int nrange) {
  // *** integrate from the origin to the first point ***
  double sum_f = f[0] / ((double)power_f);

  // *** the rest via trapesoidal rule ***
  double tmp_sum = 0.0;
  for (auto i = 0; i < nrange; ++i)
    tmp_sum += f[i];

  // *** the rest via trapesoidal rule ***
  sum_f += log_amesh * tmp_sum - 0.50 * log_amesh * (f[0] + f[nrange - 1]);

  return sum_f;
}

/************************************************
 *                                              *
 *            log_integrate_indef               *
 *                                              *
 ************************************************/
void util_log_integrate_indef(const int power_f, const double f[],
                              const int power_r, const double r[],
                              const double log_amesh, const int nrange,
                              double sum_f[]) {
  double integrand[nrange];
  // double *integrand = new double[nrange];

  for (auto k = 0; k < nrange; ++k)
    integrand[k] = f[k] * pow(r[k], (power_r + 1));

  if (nrange <= 5) {
    for (auto k = 0; k < nrange; ++k)
      sum_f[k] = integrand[k] / ((double)(power_r + power_f + 1));
  } else {
    for (auto k = 0; k < 5; ++k)
      sum_f[k] = integrand[k] / ((double)(power_r + power_f + 1));

    for (auto k = 5; k < nrange; ++k)
      sum_f[k] =
          sum_f[k - 1] + log_amesh * 0.50 * (integrand[k - 1] + integrand[k]);
  }
  // delete [] integrand;
}

/*******************************************
 *                                         *
 *         util_log_multipole_energy       *
 *                                         *
 *******************************************/
double util_log_multipole_energy(const int l, const int nrange,
                                 const double g_r[], const int power_q1,
                                 const double q1[], const int power_q2,
                                 const double q2[], const double log_amesh) {
  double fourpi = 16.0 * atan(1.0);
  double power_f = power_q1 + power_q2;
  double q1_l[nrange], q2_l[nrange], ftmp[nrange];

  util_log_integrate_indef(power_q1, q1, l, g_r, log_amesh, nrange, q1_l);
  util_log_integrate_indef(power_q1, q2, l, g_r, log_amesh, nrange, q2_l);

  for (auto k = 0; k < nrange; ++k)
    ftmp[k] = (q1[k] * q2_l[k] + q1_l[k] * q2[k]) / pow(g_r[k], (l + 1));

  return (util_log_integrate_def(power_f, ftmp, 0, g_r, log_amesh, nrange) *
          fourpi / ((double)(2 * l + 1)));
}

/************************************************
 *                                              *
 *           util_log_r2integrate_eric          *
 *                                              *
 ************************************************/
double util_log_r2integrate_eric(const int ngrid, const double log_amesh,
                                 const double rgrid[], const double f[]) {
  double mysum =
      (9.00 * f[0] * pow(rgrid[0], 3) + 28.00 * f[1] * pow(rgrid[1], 3) +
       23.00 * f[2] * pow(rgrid[2], 3)) /
      24.00;

  for (auto i = 3; i < ngrid; ++i)
    mysum += f[i] * pow(rgrid[i], 3);

  return (log_amesh * mysum + f[0] * pow(rgrid[0], 3) / 3.00);
}

/************************************************
 *                                              *
 *            log_corrector_iF                  *
 *                                              *
 ************************************************/
/*
     Computes y(i) <-- y(i+1) + F(f(i),f(i+1),f(i+2),f(i+3))
    where F is a 5 point corrector                           */

double util_log_corrector_iF(const int i, const double f[]) {
  double oneo24 = 1.0 / 24.0;

  return (-oneo24 *
          (9.00 * f[i] + 19.00 * f[i + 1] - 5.00 * f[i + 2] + 1.00 * f[i + 3]));
}

/************************************************
 *                                              *
 *            util_log_coulomb0_energy          *
 *                                              *
 ************************************************/
double util_log_coulomb0_energy(const double rho[], const double charge,
                                const double rgrid[], const int ngrid,
                                const double log_amesh, const double zion) {
  double fourpi = 16.00 * atan(1.0);
  double tmp[ngrid];

  for (auto i = 0; i < ngrid; ++i)
    tmp[i] = fourpi * rho[i] / rgrid[i];

  double E = util_log_r2integrate_eric(ngrid, log_amesh, rgrid, tmp);

  E = -E * zion;

  return E;
}

/************************************************
 *                                              *
 *            util_log_coulomb_energy           *
 *                                              *
 ************************************************/
double util_log_coulomb_energy(const double rho[], const double charge,
                               const double rgrid[], const int ngrid,
                               const double log_amesh) {
  double fourpi = 16.0 * atan(1.0);
  double tmp[ngrid], vh[ngrid];

  for (auto i = 0; i < ngrid; ++i)
    tmp[i] = fourpi * rho[i] * log_amesh * pow(rgrid[i], 3);

  vh[ngrid - 1] = charge;
  vh[ngrid - 2] = charge;
  vh[ngrid - 3] = charge;

  for (auto i = (ngrid - 4); i >= 0; --i)
    vh[i] = vh[i + 1] + util_log_corrector_iF(i, tmp);

  for (auto i = 0; i < ngrid; ++i)
    tmp[i] = fourpi * rho[i] * log_amesh * pow(rgrid[i], 2);

  double tt = 0.0;

  for (auto i = (ngrid - 4); i >= 0; --i) {
    tt += util_log_corrector_iF(i, tmp);
    vh[i] -= rgrid[i] * tt;
  }

  for (auto i = 0; i < ngrid; ++i)
    vh[i] /= rgrid[i];

  for (auto i = 0; i < ngrid; ++i)
    tmp[i] = fourpi * rho[i] * vh[i];

  double E = util_log_r2integrate_eric(ngrid, log_amesh, rgrid, tmp);

  return (E * 0.50);
}

/******************************************************
 *                                                    *
 *             util_SpecialKummer                     *
 *                                                    *
 ******************************************************/
/*
 *     Calculates a special case of the Kummer confluent hypergeometric
 *     function, M(n+1/2,l+3/2,z) for z .LE. 0
 *
 *     This function was created by  Marat Valiev, and  modified by Eric
 * Bylaska. See Abramowitz and Stegun for the formulas used in this function.
 */
double util_SpecialKummer(const int n, const int l, const double z) {
  double eps = 1.0e-16;
  double result = 0.0;

  //*** cannot handle positive z ***
  if (z > 0.0)
    std::cout << "util_SpecialKummer:invalid parameter, z>0" << std::endl;

  //*** solution for z==0 ***
  if (z == 0.0)
    return (0.0);

  //***** M(a,a+1,z) = a * (-z)**(-a) * igamma(a,-z) = a * (-z)**(-a) * P(a,-z)
  //*Gamma(a)  where z is real and a = (n+0.5)  ****
  if (n == l) {
    result = util_gammp(n + 0.50, (-z)) * (n + 0.50) *
             pow((-z), ((-n) - 0.50)) * util_gamma(n + 0.50);
    return result;
  }

  //***** M(a,a,z) = exp(z)  where a = (n+0.5)  ****
  else if (n == (l + 1)) {
    result = exp(z);
    return result;
  }

  //*** do inifinite series for small z
  if (std::abs(z) <= 1.0) {
    result = 1.0;
    double s = 1.0;
    double a = n + 0.5;
    double b = l + 1.5;
    int i = 1;
    while ((i < 10000) && (std::abs(s) < eps)) {
      s *= (a + i - 1) * z / ((b + i - 1) * i);
      result += s;
      ++i;
    }
    if (i > 10000)
      std::cout << "util_SpecialKummer:cannot converge" << std::endl;
    return result;
  }

  if (n < l) {
    //*** starting point n=l or b=a+1***
    double a = n + 0.5;
    double b = n + 1.5;

    //*** m1 = M(a,b-1) ***
    //*** m2 = M(a,b,z) ***
    double m1 = exp(z);
    result = util_gammp(a, (-z)) * a / pow((-z), a) * util_gamma(a);

    // using recursion formula
    // z(a-b)M(a,b+1,z)=b(b-1)M(a,b-1,z)+b(1-b-z)M(a,b,z)
    // obtain M(1/2,3/2+l  ,z) --> m2
    //        M(1/2,3/2+l-1,z) --> m2
    for (auto i = 1; i <= (l - n); ++i) {
      double m3 =
          (b * (b - 1.0) * m1 + b * (1.0 - b - z) * result) / (z * (a - b));
      b += 1.0;
      m1 = result;
      result = m3;
    }
  } else if (n > (l + 1)) {
    //*** starting point n=l+1 or b=a ***
    double a = l + 1.5;
    double b = l + 1.5;

    //*** m1 = M(a-1,b) ***
    //*** m2 = M(a,a,z) ***
    double m1 = util_gammp(a - 1.0, (-z)) * (a - 1.0) / pow((-z), (a - 1.0)) *
                util_gamma(a - 1.0);
    result = exp(z);

    // using recursion formula
    // aM(a+1,b,z)=(b-a)M(a-1,b,z)+(2a-b+z)M(a,b,z)
    // obtain M(n+1/2-1,3/2,z)   --> m1
    //        M(n+1/2  ,3/2,z)   --> m2
    for (auto i = 1; i < (n - l - 1); ++i) {
      double m3 = ((b - a) * m1 + (2 * a - b + z) * result) / a;
      m1 = result;
      result = m3;
      a += 1.0;
    }
  }

  return result;
}

/*****************************************************
 *                                                   *
 *               util_gauss_weights                  *
 *                                                   *
 *****************************************************/
void util_gauss_weights(const double x1, const double x2, double x[],
                        double w[], const int n) {
  double pp, z1, z;
  double eps = 3.0e-14;
  double pi = 4.0 * atan(1.0);
  double xm = 0.5 * (x2 + x1);
  double xl = 0.5 * (x2 - x1);
  int m = (n + 1) / 2;

  for (auto i = 0; i < m; ++i) {
    // z    = cos(pi*(0-0.250)/(n+0.50));
    z = cos(pi * (i - 1.250) / (n + 0.50));

    int niter = 0;
    bool repeat = true;

    while (repeat && (niter < 1000000)) {
      ++niter;

      double p1 = 1.00;
      double p2 = 0.00;

      for (auto j = 1; j <= n; ++j) {
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.00) * p3) / ((double)j);
      }

      pp = n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;

      if (std::abs(z - z1) <= eps)
        repeat = false;
    }
    if (niter >= 1000000)
      std::cout << "cannot converge in gauss_weights" << std::endl;

    x[i] = xm - xl * z;
    x[n - i - 1] = xm + xl * z;
    w[i] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w[n - i - 1] = w[i];
  }
}

/******************************************************
 *                                                    *
 *               util_GaussBessel                     *
 *                                                    *
 ******************************************************/
/*
     Calculates the Gaussian Bessel function,

                                 /infinity
    GausBessel(n,l,alpha,R)  =   | k**n exp(-alpha**2 * k**2) j_l(R*k) dk
                                 /0

     This function uses the SpecialKummer function. Note it is assumed that
     (n+l) is an even integer.
*/
double util_GaussBessel(const int n, const int l, const double alpha,
                        const double R) {
  if ((n + l % 2) == 1)
    std::cout << "util_GaussBessel: n+l is not even" << std::endl;

  double pi = 4.0 * atan(1.0);
  double c = (sqrt(pi) / (pow(2.0, (l + 2)) * pow(alpha, (n + l + 1)))) *
             exp(util_ln_gamma((n + l + 1) / 2.0) - util_ln_gamma(l + 1.5));

  return (c * pow(R, l) *
          util_SpecialKummer((n + l) / 2, l, pow(-(0.5 * R / alpha), 2)));
}

/******************************************************
 *                                                    *
 *                util_dGaussBessel                   *
 *                                                    *
 ******************************************************/
/*
     Calculates the derivative of the Gaussian Bessel function wrt R,

                                       /infinity
   dGausBessel(n,l,alpha,R)  = (d/dR)  | k**n exp(-alpha**2 * k**2) j_l(R*k) dk
                                       /0

     This function uses the SpecialKummer function. Note it is assumed that
     (n+l) is an even integer.
*/
double util_dGaussBessel(const int n, const int l, const double alpha,
                         const double R) {

  if ((n + l % 2) == 1)
    std::cout << "util_dGaussBessel: n+l is not even" << std::endl;
  double pi = 4.0 * atan(1.0);
  double c = (sqrt(pi) / (pow(2.0, (l + 2)) * pow(alpha, (n + l + 1)))) *
             exp(util_ln_gamma((n + l + 1) / 2.0) - util_ln_gamma(l + 1.5));

  return (c *
          (pow(l * R, (l - 1)) *
               util_SpecialKummer((n + l) / 2, l, pow(-(0.5 * R / alpha), 2)) -
           (0.5 * pow(R, (l + 1)) / pow(alpha, 2)) * ((n + l) / 2.0 + 0.5) /
               (l + 1.5) *
               util_SpecialKummer((n + l) / 2 + 1, l + 1,
                                  pow(-(0.5 * R / alpha), 2))));
}

} // namespace pwdft

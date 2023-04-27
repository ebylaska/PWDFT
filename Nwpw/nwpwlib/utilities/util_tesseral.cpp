
//#include 	<iomanip>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include "Parallel.hpp"
#include "util_gamma.hpp"
#include "util_legendre.hpp"
#include "util_tesseral.hpp"

namespace pwdft {

/**************************************
 *                                    *
 *      util_rSphericalHarmonic3      *
 *                                    *
 **************************************/
/*
   Calculates the r*Spherical Harmonic for x,y,z such that

                                                                                  {exp(i*|m|*phi)  m>0
      rY_lm(cos_theta,phi)=(r**l)*sqrt((l-|m|)!/(l+|m|)!)*rlegendre_lm(cos_theta)*{1
   m==0 {exp(-i*|m|*phi) m<0 where   cos_theta = z/r and phi = atan2(y,x)
*/

std::complex<double> util_rSphericalHarmonic3(const int l, const int m,
                                              const double x, const double y,
                                              const double z) {
  std::complex<double> rSphericalHarmonic3;
  std::complex<double> ctmp;
  int mod_m = abs(m);

  if (mod_m > l)
    std::cout << "Parameter out of order in function rSphericalHarmonic3"
              << std::endl;

  double r = sqrt(x * x + y * y + z * z);
  if (r > 1.0e-9) {
    double coeff;
    double cos_theta = z / r;
    double phi = atan2(y, x);

    // *** find coefficient ***
    if (mod_m == 0) {
      coeff = 1.0;
      ctmp = std::complex<double>(1.0, 0.0);
    } else {
      coeff = 1.0;
      for (auto i = 1; i <= (2 * mod_m); ++i)
        coeff /= ((double)(l - mod_m + i));
      ctmp = std::complex<double>(cos(m * phi), sin(m * phi));
    }
    coeff = sqrt(coeff);
    rSphericalHarmonic3 =
        ctmp * coeff * util_rlegendre_lm(l, mod_m, cos_theta) * std::pow(r, l);
  } else {
    if (l == 0) {
      ctmp = std::complex<double>(1.0, 0.0);
    } else {
      ctmp = std::complex<double>(0.0, 0.0);
    }
    rSphericalHarmonic3 = ctmp * std::pow(r, l);
  }

  return rSphericalHarmonic3;
}

/**************************************
 *                                    *
 *      util_Tesseral3_vector_lm      *
 *                                    *
 **************************************/
/*
  Calculates the tesseral harmonic for x,y,z such that

                                               {cos(|m|*phi)   m>0
      T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
                                               {sin(|m|*phi)   m<0

    where   cos_theta = z/r and phi = atan2(y,x)
*/
void util_Tesseral3_vector_lm(const int l, const int m, const int nfft3d,
                              const double x[], const double y[],
                              const double z[], double Tlm[]) {
  double coeff;
  double twopi = 8.00 * atan(1.0);
  int mod_m = abs(m);

  if (mod_m > l)
    std::cout << "Parameter out of order in function Tesseral3_vector_lm"
              << std::endl;

  // *** find coefficient ***
  if (mod_m == 0) {
    coeff = 0.5;
  } else {
    coeff = 1.0;
    for (auto i = 1; i <= (2 * mod_m); ++i)
      coeff /= ((double)(l - mod_m + i));
  }
  coeff *= (2 * l + 1) / twopi;
  coeff = sqrt(coeff);

  for (auto k = 0; k < nfft3d; ++k) {
    double r = sqrt(x[k] * x[k] + y[k] * y[k] + z[k] * z[k]);
    ;
    if (r > 1.0e-9) {
      double tmp2;
      double cos_theta = z[k] / r;
      double phi = atan2(y[k], x[k]);

      if (m < 0) {
        tmp2 = sin(mod_m * phi);
      } else if (m > 0) {
        tmp2 = cos(mod_m * phi);
      } else {
        tmp2 = 1.0;
      }
      Tlm[k] = tmp2 * coeff * util_rlegendre_lm(l, mod_m, cos_theta);
    } else {
      if (l == 0) {
        Tlm[k] = 1.0 / sqrt(2 * twopi);
      } else {
        Tlm[k] = 0.0;
      }
    }
  }
}

/**************************************
 *                                    *
 *       util_Tesseral3_rgrid_lm      *
 *                                    *
 **************************************/
/*
   Calculates the tesseral harmonic for x,y,z such that

                                                {cos(|m|*phi)   m>0
       T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
                                                {sin(|m|*phi)   m<0

     where   cos_theta = z/r and phi = atan2(y,x)
*/
void util_Tesseral3_rgrid_lm(const int l, const int m, const int nfft3d,
                             const double rgrid[], double Tlm[]) {
  double coeff;
  double twopi = 8.0 * atan(1.0);
  int mod_m = abs(m);

  if (mod_m > l)
    std::cout << "Parameter out of order in function Tesseral3_vector_lm"
              << std::endl;

  // *** find coefficient ***
  if (mod_m == 0) {
    coeff = 0.5;
  } else {
    coeff = 1.0;
    for (auto i = 1; i <= (2 * mod_m); ++i)
      coeff /= ((double)(l - mod_m + i));
  }
  coeff *= (2 * l + 1) / twopi;
  coeff = sqrt(coeff);

  for (auto k = 0; k < nfft3d; ++k) {
    double x = rgrid[3 * k];
    double y = rgrid[3 * k + 1];
    double z = rgrid[3 * k + 2];
    double r = sqrt(x * x + y * y + z * z);
    if (r > 1.0e-9) {
      double tmp2;
      double cos_theta = z / r;
      double phi = atan2(y, x);

      if (m < 0) {
        tmp2 = sin(mod_m * phi);
      } else if (m > 0) {
        tmp2 = cos(mod_m * phi);
      } else {
        tmp2 = 1.0;
      }
      Tlm[k] = tmp2 * coeff * util_rlegendre_lm(l, mod_m, cos_theta);
    } else {
      if (l == 0) {
        Tlm[k] = 1.0 / sqrt(2 * twopi);
      } else {
        Tlm[k] = 0.0;
      }
    }
  }
}

/**************************************
 *                                    *
 *     util_Tesseral3_rgrid_lm_rl     *
 *                                    *
 **************************************/
/*
   Calculates the r**l times the tesseral harmonic for x,y,z such that

                                                          {cos(|m|*phi)   m>0
       r**l*T_lm(cos_theta,phi)=r**l*rtheta_lm(cos_theta)*{1              m==0
                                                          {sin(|m|*phi)   m<0
   where   cos_theta = z/r and phi = atan2(y,x)
*/
void util_Tesseral3_rgrid_lm_rl(const int l, const int m, const int nfft3d,
                                const double rgrid[], double rTlm[]) {
  double coeff;
  double twopi = 8.0 * atan(1.0);
  int mod_m = abs(m);

  if (mod_m > l)
    std::cout << "Parameter out of order in function Tesseral3_vector_lm"
              << std::endl;

  // *** find coefficient ***
  if (mod_m == 0) {
    coeff = 0.5;
  } else {
    coeff = 1.0;
    for (auto i = 1; i <= (2 * mod_m); ++i)
      coeff /= ((double)(l - mod_m + i));
  }
  coeff *= (2 * l + 1) / twopi;
  coeff = sqrt(coeff);

  for (auto k = 0; k < nfft3d; ++k) {
    double tmp2;
    double x = rgrid[3 * k];
    double y = rgrid[3 * k + 1];
    double z = rgrid[3 * k + 2];
    double r = sqrt(x * x + y * y + z * z);
    if (r > 1.0e-9) {
      double cos_theta = z / r;
      double phi = atan2(y, x);

      if (m < 0) {
        tmp2 = sin(mod_m * phi);
      } else if (m > 0) {
        tmp2 = cos(mod_m * phi);
      } else {
        tmp2 = 1.0;
      }
      rTlm[k] = std::pow(r, l) * tmp2 * coeff *
                util_rlegendre_lm(l, mod_m, cos_theta);
    } else {
      rTlm[k] = 0.0;
    }
  }
}

/**************************************
 *                                    *
 *   util_Tesseral3_rgrid_lm_over_rl  *
 *                                    *
 **************************************/
/*
  Calculates the 1/r**(l+1) times the tesseral harmonic for x,y,z such that

                                                                     {cos(|m|*phi)
  m>0 T_lm(cos_theta,phi)/r**(l+1)=(1/r**(l+1))*rtheta_lm(cos_theta)*{1 m==0
                                                                     {sin(|m|*phi)
  m<0 where   cos_theta = z/r and phi = atan2(y,x)
*/
void util_Tesseral3_rgrid_lm_over_rl(const int l, const int m, const int nfft3d,
                                     const double rgrid[], double rTlm[]) {
  double coeff;
  double twopi = 8.0 * atan(1.0);
  int mod_m = abs(m);

  if (mod_m > l)
    std::cout << "parameter out of order in function Tesseral3_vector_lm"
              << std::endl;

  // *** find coefficient ***
  if (mod_m == 0) {
    coeff = 0.5;
  } else {
    coeff = 1.0;
    for (auto i = 1; i <= (2 * mod_m); ++i)
      coeff /= ((double)(l - mod_m + i));
  }
  coeff *= (2 * l + 1) / twopi;
  coeff = sqrt(coeff);

  for (auto k = 0; k < nfft3d; ++k) {
    double x = rgrid[3 * k];
    double y = rgrid[3 * k + 1];
    double z = rgrid[3 * k + 2];
    double r = sqrt(x * x + y * y + z * z);
    if (r > 1.0e-9) {
      double tmp2;
      double cos_theta = z / r;
      double phi = atan2(y, x);

      if (m == 0) {
        tmp2 = sin(mod_m * phi);
      } else if (m > 0) {
        tmp2 = cos(mod_m * phi);
      } else {
        tmp2 = 1.0;
      }
      rTlm[k] = tmp2 * coeff * util_rlegendre_lm(l, mod_m, cos_theta) /
                std::pow(r, (l + 1));
    } else {
      rTlm[k] = 0.0;
    }
  }
}

/**************************************
 *                                    *
 *           util_Tesseral3_lm        *
 *                                    *
 **************************************/
/*
  Calculates the tesseral harmonic for x,y,z such that

                                               {cos(|m|*phi)   m>0
      T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
                                               {sin(|m|*phi)   m<0
    where   cos_theta = z/r and phi = atan2(y,x)
*/
double util_Tesseral3_lm(const int l, const int m, const double x,
                         const double y, const double z) {
  int mod_m = abs(m);
  double tmp2;
  double r = sqrt(x * x + y * y + z * z);
  double cos_theta = z / r;
  double phi = atan2(y, x);

  if (m < 0) {
    tmp2 = sin(mod_m * phi);
  } else if (m > 0) {
    tmp2 = cos(mod_m * phi);
  } else {
    tmp2 = 1.0;
  }

  return (util_rtheta_lm(l, m, cos_theta) * tmp2);
}

/**************************************
 *                                    *
 *      util_ Tesseral3_lm_over_rl    *
 *                                    *
 **************************************/
/*
  Calculates the tesseral harmonic divided by r**(l+1) for x,y,z such that

                                               {cos(|m|*phi)   m>0
      T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
                                               {sin(|m|*phi)   m<0
    where   cos_theta = z/r and phi = atan2(y,x)
*/
double util_Tesseral3_lm_over_rl(const int l, const int m, const double x,
                                 const double y, const double z) {
  int mod_m = abs(m);
  double tmp2;
  double r = sqrt(x * x + y * y + z * z);
  double cos_theta = z / r;
  double phi = atan2(y, x);

  if (m < 0) {
    tmp2 = sin(mod_m * phi);
  } else if (m > 0) {
    tmp2 = cos(mod_m * phi);
  } else {
    tmp2 = 1.0;
  }

  return (util_rtheta_lm(l, m, cos_theta) * tmp2 / std::pow(r, (l + 1)));
}

/**************************************
 *                                    *
 *           util_dTesseral3_lm       *
 *                                    *
 **************************************/
/*
  Calculates the derivative of tesseral harmonic for x,y,z wrt to x,y,z,such
  that

                                               {cos(|m|*phi)   m>0
      T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
                                               {sin(|m|*phi)   m<0
    where   cos_theta = z/r and phi = atan2(y,x)
*/
void util_dTesseral3_lm(const int l, const int m, const double x,
                        const double y, const double z, double &dTx,
                        double &dTy, double &dTz) {

  int mod_m = abs(m);
  double tmp2, f2;
  double r = sqrt(x * x + y * y + z * z);
  double cos_theta = z / r;
  double phi = atan2(y, x);

  if (m < 0) {
    tmp2 = sin(mod_m * phi);
    f2 = mod_m * cos(mod_m * phi) * util_rtheta_lm_div(l, m, cos_theta);
  } else if (m > 0) {
    tmp2 = cos(mod_m * phi);
    f2 = -mod_m * sin(mod_m * phi) * util_rtheta_lm_div(l, m, cos_theta);
  } else {
    tmp2 = 1.0;
    f2 = 0.0;
  }
  double f1 = util_drtheta_lm(l, m, cos_theta) * tmp2;

  dTx = f1 * cos_theta * cos(phi) - f2 * sin(phi);
  dTy = f1 * cos_theta * sin(phi) + f2 * cos(phi);
  dTz = -f1 * sqrt(1.0 - cos_theta * cos_theta);
}

/**************************************
 *                                    *
 *       util_dTesseral3_lm_over_rl   *
 *                                    *
 **************************************/
/*
  Calculates the derivative of tesseral harmonic for x,y,z wrt to x,y,z,such
  that

                                                {cos(|m|*phi)   m>0
      dT_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
                                                {sin(|m|*phi)   m<0
    where   cos_theta = z/r and phi = atan2(y,x)
*/
void util_dTesseral3_lm_over_rl(const int l, const int m, const double x,
                                const double y, const double z, double &dTx,
                                double &dTy, double &dTz) {
  double coeff;
  double twopi = 8.0 * atan(1.0);
  int mod_m = abs(m);

  if (mod_m > l)
    std::cout << "parameter out of order in function dTesseral3_lm_over_rl"
              << std::endl;

  // *** find coefficient ***
  if (mod_m == 0) {
    coeff = 0.5;
  } else {
    coeff = 1.0;
    for (auto i = 1; i <= (2 * mod_m); ++i)
      coeff /= ((double)(l - mod_m + i));
  }
  coeff *= (2 * l + 1) / twopi;
  coeff = sqrt(coeff);

  double r = sqrt(x * x + y * y + z * z);
  if (r > 1.0e-9) {
    double tmp2, f2;
    double rl1 = 1.0 / std::pow(r, (l + 1));
    double drl1 = ((double)(-(l + 1))) / std::pow(r, (l + 2));
    double cos_theta = z / r;
    double phi = atan2(y, x);

    if (m < 0) {
      tmp2 = sin(mod_m * phi);
      f2 = mod_m * cos(mod_m * phi) * util_rtheta_lm_div(l, m, cos_theta);
    } else if (m > 0) {
      tmp2 = cos(mod_m * phi);
      f2 = -mod_m * sin(mod_m * phi) * util_rtheta_lm_div(l, m, cos_theta);
    } else {
      tmp2 = 1.0;
      f2 = 0.0;
    }
    double drTlm = tmp2 * coeff * util_rlegendre_lm(l, mod_m, cos_theta) * drl1;
    double f1 = util_drtheta_lm(l, m, cos_theta) * tmp2;

    dTx = drTlm * x / r + (f1 * cos_theta * cos(phi) - f2 * sin(phi)) * rl1;
    dTy = drTlm * y / r + (f1 * cos_theta * sin(phi) + f2 * cos(phi)) * rl1;
    dTz = drTlm * z / r - f1 * sqrt(1.0 - cos_theta * cos_theta) * rl1;
  } else {
    dTx = 0.0;
    dTy = 0.0;
    dTz = 0.0;
  }
}

/**************************************
 *                                    *
 *      dTesseral3_rgrid_lm_rl        *
 *                                    *
 **************************************/

void dTesseral3_rgrid_lm_rl(const int l, const int m, const int nfft3d,
                            const double rgrid[], double drTlm[]) {
  int mod_m = abs(m);
  double coeff;
  double twopi = 8.0 * atan(1.0);

  if (mod_m > l)
    std::cout << "parameter out of order in function dTesseral3_rgrid_lm_rl"
              << std::endl;

  // *** find coefficient ***
  if (mod_m == 0) {
    coeff = 0.5;
  } else {
    coeff = 1.0;
    for (auto i = 1; i <= (2 * mod_m); ++i)
      coeff /= ((double)(l - mod_m + i));
  }
  coeff *= (2 * l + 1) / twopi;
  coeff = sqrt(coeff);

  for (auto k = 0; k < nfft3d; ++k) {
    double x = rgrid[3 * k];
    double y = rgrid[3 * k + 1];
    double z = rgrid[3 * k + 2];
    double r = sqrt(x * x + y * y + z * z);

    if (r > 1.0e-9) {
      double tmp2, f2, drl;
      double rl = std::pow(r, l);
      if (l == 0) {
        drl = 0.0;
      } else if (l == 1) {
        drl = 1.0;
      } else {
        drl = ((double)l) * std::pow(r, (l - 1));
      }

      double cos_theta = z / r;
      double phi = atan2(y, x);

      if (m < 0) {
        tmp2 = sin(mod_m * phi);
        f2 = mod_m * cos(mod_m * phi) * util_rtheta_lm_div(l, m, cos_theta);
      } else if (m > 0) {
        tmp2 = cos(mod_m * phi);
        f2 = -mod_m * sin(mod_m * phi) * util_rtheta_lm_div(l, m, cos_theta);
      } else {
        tmp2 = 1.0;
        f2 = 0.0;
      }
      double Tlm = drl * tmp2 * coeff * util_rlegendre_lm(l, mod_m, cos_theta);
      double f1 = util_drtheta_lm(l, m, cos_theta) * tmp2;

      drTlm[3 * k] =
          Tlm * x / r + rl * (f1 * cos_theta * cos(phi) - f2 * sin(phi));
      drTlm[3 * k + 1] =
          Tlm * y / r + rl * (f1 * cos_theta * sin(phi) + f2 * cos(phi));
      drTlm[3 * k + 2] =
          Tlm * z / r - rl * f1 * sqrt(1.0 - cos_theta * cos_theta);
    } else {
      drTlm[3 * k] = 0.0;
      drTlm[3 * k + 1] = 0.0;
      drTlm[3 * k + 2] = 0.0;
    }
  }
}

/**************************************
 *                                    *
 *         util_Tesseral_lm           *
 *                                    *
 **************************************/
/*
  Calculates the tesseral harmonic for x,y,z such that

                                               {cos(|m|*phi)   m>0
      T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
                                               {sin(|m|*phi)   m<0
    where   cos_theta = z/r and phi = atan2(y,x)
*/
double util_Tesseral_lm(const int l, const int m, const double cos_theta,
                        const double phi) {
  double tmp2;
  int mod_m = abs(m);

  if (m < 0) {
    tmp2 = sin(mod_m * phi);
  } else if (m > 0) {
    tmp2 = cos(mod_m * phi);
  } else {
    tmp2 = 1.0;
  }

  return (util_rtheta_lm(l, m, cos_theta) * tmp2);
}

/**************************************
 *                                    *
 *          util_YSpherical_lm        *
 *                                    *
 **************************************/
/*
  Calculates the spherical harmonic for x,y,z such that

      Y_lm(cos_theta,phi)=rtheta_lm(cos_theta)* exp(i*m*phi)

    where   cos_theta = z/r and phi = atan2(y,x)
*/
std::complex<double> util_YSpherical_lm(const int l, const int m,
                                        const double cos_theta,
                                        const double phi) {
  return (util_ytheta_lm(l, m, cos_theta) *
          std::complex<double>(cos(m * phi), sin(m * phi)));
}

/**************************************
 *                                    *
 *         util_dTesseral_lm          *
 *                                    *
 **************************************/
/*
  Calculates the derivative of tesseral harmonic for x,y,z wrt to x,y,z,such
  that

                                               {cos(|m|*phi)   m>0
      T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
                                               {sin(|m|*phi)   m<0
    where   cos_theta = z/r and phi = atan2(y,x)
*/
void util_dTesseral_lm(const int l, const int m, const double cos_theta,
                       const double phi, double &dTx, double &dTy,
                       double &dTz) {
  int mod_m = abs(m);
  double tmp2, f2;

  if (m < 0) {
    tmp2 = sin(mod_m * phi);
    f2 = mod_m * cos(mod_m * phi) * util_rtheta_lm_div(l, m, cos_theta);
  } else if (m > 0) {
    tmp2 = cos(mod_m * phi);
    f2 = -mod_m * sin(mod_m * phi) * util_rtheta_lm_div(l, m, cos_theta);
  } else {
    tmp2 = 1.0;
    f2 = 0.0;
  }
  double f1 = util_drtheta_lm(l, m, cos_theta) * tmp2;

  dTx = f1 * cos_theta * cos(phi) - f2 * sin(phi);
  dTy = f1 * cos_theta * sin(phi) + f2 * cos(phi);
  dTz = -f1 * sqrt(1.0 - cos_theta * cos_theta);
}

/**************************************
 *                                    *
 *          util_dYspherical_lm       *
 *                                    *
 **************************************/
/*
  Calculates the derivative of tesseral harmonic for x,y,z wrt to x,y,z,such
  that

      Y_lm(cos_theta,phi)=theta_lm(cos_theta)* exp(i*m*phi)

    where   cos_theta = z/r and phi = atan2(y,x)
*/
void util_dYspherical_lm(const int l, const int m, const double cos_theta,
                         const double phi, std::complex<double> &dYx,
                         std::complex<double> &dYy, std::complex<double> &dYz) {
  std::complex<double> f1, f2;

  // *** derivative with respect to theta ***
  f1 = util_dytheta_lm(l, m, cos_theta) *
       std::complex<double>(cos(m * phi), sin(m * phi));

  // *** derivative with respect to phi ***
  if (m == 0) {
    f2 = std::complex<double>(0.0, 0.0);
  } else {
    f2 = util_ytheta_lm_div(l, m, cos_theta) *
         std::complex<double>(cos(m * phi), sin(m * phi)) *
         std::complex<double>(0.0, 1.0 * m);
  }

  dYx = f1 * cos_theta * cos(phi) - f2 * sin(phi);
  dYy = f1 * cos_theta * sin(phi) + f2 * cos(phi);
  dYz = -f1 * sqrt(1.0 - cos_theta * cos_theta);
}

} // namespace pwdft

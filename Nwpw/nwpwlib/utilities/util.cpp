
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

void t_aindexcopy(const int n, const int *indx, double *A, double *B) 
{
   for (auto i = 0; i < n; ++i)
      B[i] = A[indx[i]];
}

void t_bindexcopy(const int n, const int *indx, double *A, double *B) 
{
   for (auto i = 0; i < n; ++i)
      B[indx[i]] = A[i];
}

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

void eigsrt(double *D, double *V, int n) {
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

/**************************************
 *                                    *
 *         util_getfilling            *
 *                                    *
 **************************************/
void util_getfilling(int f, int nfft[], int *filling, double zvalue[]) {
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
  if (i == 0) {
    filling[1] = j;
    filling[2] = k;
  } else {

    if (((j + 1) % 2) == 0)
      filling[1] = (j + 1) / 2;
    else
      filling[1] = -(j + 1) / 2;

    if (((k + 1) % 2) == 0)
      filling[2] = (k + 1) / 2;
    else
      filling[2] = -(k + 1) / 2;
  }

  if ((i == 0) && (j == 0) && (k == 0)) {
    filling[3] = 0;
    zvalue[0] = 1.0;
    zvalue[1] = 0.0;
  } else if (h == 0) {
    filling[3] = 2;
    zvalue[0] = 1.0 / sqrt(2.0);
    zvalue[1] = 0.0;
  } else {
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
double util_random(const int seed) {
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
                 const double yp1, const double ypn, double *y2, double *utmp) {
  double sig, p, qn, un;

  if (yp1 > 0.99e30) {
    y2[0] = 0.0;
    utmp[0] = 0.0;
  } else {
    y2[0] = -0.5;
    utmp[0] = 3.0 / (x[1] - x[0]) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
  }
  for (auto i = 1; i < (n - 1); ++i) {
    sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
    p = sig * y2[i - 1] + 2.00;
    y2[i] = (sig - 1.00) / p;
    utmp[i] = (6.00 *
                   ((y[i + 1] - y[i]) / (x[i + 1] - x[i]) -
                    (y[i] - y[i - 1]) / (x[i] - x[i - 1])) /
                   (x[i + 1] - x[i - 1]) -
               sig * utmp[i - 1]) /
              p;
  }

  if (ypn > 0.99e30) {
    qn = 0.0;
    un = 0.0;
  } else {
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
                   const int n, const int nx, const double x) {
  int khi = nx;
  int klo = nx - 1;

  while ((xa[klo] > x) || (xa[khi] < x)) {
    if (xa[klo] > x) {
      klo = klo - 1;
      khi = khi - 1;
    }
    if (xa[khi] < x) {
      klo = klo + 1;
      khi = khi + 1;
    }
  }

  double h = xa[khi] - xa[klo];
  double a = (xa[khi] - x) / h;
  double b = (x - xa[klo]) / h;

  return (a * ya[klo] + b * ya[khi] +
          ((a * a * a - a) * y2a[klo] + (b * b * b - b) * y2a[khi]) * h * h /
              6.0);
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
double util_dswitching_function(const double s_d, const double s_rho,
                                const double r) {
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
 *       util_print_elapsed_time      *
 *                                    *
 **************************************/
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
  else if (sectime < 1.0e-9)
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

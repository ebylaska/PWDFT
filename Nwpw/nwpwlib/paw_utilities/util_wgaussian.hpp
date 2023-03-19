#ifndef _UTIL_wgaussian_HPP_
#define _UTIL_wgaussian_HPP_

#pragma once

//#include      <iomanip>
#include        <iostream>
#include        <cstdlib>
#include        <complex>
#include        <cmath>

namespace pwdft {


extern double util_WGaussian(const int, const int, const double, const int, const int, const double, const double *);
extern void util_dWGaussian(const int, const int, const double, const int, const int, const double, const double *, double&, double *);
extern double util_UGaussian(const int, const int, const double, const int, const int, const double);
extern double util_WGaussian3(const int la, const int, const double, const int, const int, const double, const double, const double *);

extern std::complex<double> util_CWGaussian(const int, const int, const double, const int, const int, const double, const double *);
extern void util_dCWGaussian(const int, const int, const double, const int, const int, const double, const double *, std::complex<double>&,  std::complex<double> *);
extern std::complex<double> util_CWGaussian3(const int la, const int ma, const double, const int, const int, const double, const double, const double *);
extern std::complex<double> util_CWGaussian2(const int, const int, const double, const int, const int, const double, const double *);
extern void util_dCWGaussian3(const int, const int, const double, const int, const int, const double, const double, const double *, std::complex<double>&,  std::complex<double> *);

}
#endif

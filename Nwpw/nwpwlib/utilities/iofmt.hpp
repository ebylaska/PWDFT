#ifndef _IOFMT_HPP_
#define _IOFMT_HPP_

#pragma once

#include <cstring>
#include <iomanip>
#include <iostream>

#define Efmt(w, p)                                                             \
  std::right << std::setw(w) << std::setprecision(p) << std::scientific
#define Ffmt(w, p)                                                             \
  std::right << std::setw(w) << std::setprecision(p) << std::fixed
#define Ifmt(w) std::right << std::setw(w)

// static void iofmt_forceprint(std::string forcelabel, double *A) {
//    int nion=2;
//    std::cout << "iofmt_forceprint: " << forcelabel << std::endl;
//    for (int ii=0; ii<nion; ++ii)
//       std::cout << Efmt(20,12) <<  A[3*ii] << " " << Efmt(20,12) << A[3*ii+1]
//       << " " << Efmt(20,12) <<  A[3*ii+2] << std::endl;
//    std::cout << std::endl;

//}

#endif

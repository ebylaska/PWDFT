#ifndef _PSPW_LMBFGS2_HPP_
#define _PSPW_LMBFGS2_HPP_

#pragma once

#include "Geodesic2.hpp"
#include <cmath>

namespace pwdft {

class pspw_lmbfgs2 {

  Geodesic2 *mygeodesic;

  int npack1, neall, nsize, max_m, m;
  int indx[20], size_list;
  double rho[20];
  double *lm_list;

public:
  /* Constructors */
  pspw_lmbfgs2(Geodesic2 *mygeodesic0, const int max_m0) {

    mygeodesic = mygeodesic0;
    max_m = max_m0;
    npack1 = mygeodesic->mygrid->npack(1);
    neall = mygeodesic->mygrid->neq[0] + mygeodesic->mygrid->neq[1];
    nsize = 2 * neall * npack1;
    size_list = 2 * max_m;
    for (auto k = 0; k < size_list; ++k)
      indx[k] = k;

    m = 0;

    lm_list = mygeodesic->mygrid->g_nallocate(1, size_list);

    // std::memcpy(&lm_list[2*m*nsize],g,nsize*sizeof(double));
    // mygeodesic->mygrid->g_Scale(-1.0,&lm_list[2*m*nsize]);
    // std::memcpy(&lm_list[(2*m+1)*nsize],&lm_list[2*m*nsize],nsize*sizeof(double));
  }

  /* destructor */
  ~pspw_lmbfgs2() { mygeodesic->mygrid->g_deallocate(lm_list); }

  void start(const double *g) {
    std::memcpy(&lm_list[2 * m * nsize], g, nsize * sizeof(double));
    mygeodesic->mygrid->g_Scale(-1.0, &lm_list[2 * m * nsize]);
    std::memcpy(&lm_list[(2 * m + 1) * nsize], &lm_list[2 * m * nsize],
                nsize * sizeof(double));
  }

  void fetch(const double tmin, double *g, double *s) {
    double *yy, *ss, sum, sum0, alpha[20], beta;

    mygeodesic->mygrid->g_Scale(-1.0, g);

    std::memcpy(s, g, nsize * sizeof(double));

    yy = &lm_list[indx[2 * m] * nsize];
    ss = &lm_list[indx[2 * m + 1] * nsize];

    // mygeodesic->psi_1Gtransport(tmin,yy);
    mygeodesic->psi_1transport(tmin, ss);

    mygeodesic->mygrid->gg_daxpy(-1.0, g, yy);
    mygeodesic->mygrid->g_Scale(-1.0, yy);

    sum = mygeodesic->mygrid->gg_traceall(yy, ss);

    //*** exit if dividing by small number ***
    if (abs(sum) > 1.0e-15) {
      rho[m] = 1.0 / sum;

      sum = mygeodesic->mygrid->gg_traceall(ss, s);
      alpha[m] = rho[m] * sum;

      mygeodesic->mygrid->gg_daxpy(-alpha[m], yy, s);
      for (auto k = m - 1; k > 0; --k) {
        yy = &lm_list[indx[2 * (k + 1)] * nsize];
        ss = &lm_list[indx[2 * (k + 1) + 1] * nsize];

        // mygeodesic->psi_1Gtransport(tmin,yy);
        // mygeodesic->psi_1Gtransport(tmin,ss);

        sum = mygeodesic->mygrid->gg_traceall(ss, s);
        alpha[k] = rho[k] * sum;

        mygeodesic->mygrid->gg_daxpy(-alpha[k], yy, s);
      }

      //**** preconditioner ****
      // call Grsm_gg_dScale(npack1,neall,h0,s,s)

      for (auto k = 0; k < (m - 1); ++k) {
        sum = mygeodesic->mygrid->gg_traceall(yy, s);
        beta = rho[k] * sum;
        sum0 = alpha[k] - beta;

        mygeodesic->mygrid->gg_daxpy(sum0, ss, s);

        yy = &lm_list[indx[2 * (k + 1)] * nsize];
        ss = &lm_list[indx[2 * (k + 1) + 1] * nsize];
      }

      sum = mygeodesic->mygrid->gg_traceall(yy, s);
      beta = rho[m] * sum;
      sum0 = alpha[m] - beta;

      mygeodesic->mygrid->gg_daxpy(sum0, ss, s);

      if (m < (max_m - 1))
        ++m;
      else {
        int itmp0 = indx[0];
        int itmp1 = indx[1];
        for (auto k = 0; k < (size_list - 2); ++k)
          indx[k] = indx[k + 2];
        indx[size_list - 2] = itmp0;
        indx[size_list - 1] = itmp1;

        for (auto k = 0; k < (m - 1); ++k)
          rho[k] = rho[k + 1];
      }
      mygeodesic->mygrid->g_Scale(-1.0, s);
    }
    std::memcpy(&lm_list[indx[2 * m] * nsize], g, nsize * sizeof(double));
    std::memcpy(&lm_list[indx[2 * m + 1] * nsize], s, nsize * sizeof(double));
    mygeodesic->mygrid->g_Scale(-1.0, g);
  }
};

} // namespace pwdft

#endif

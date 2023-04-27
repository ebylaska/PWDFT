/* Paw_xc.cpp -
   Author - Eric Bylaska
*/

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>

#include "Paw_xc.hpp"
#include "blas.h"
#include "util.hpp"
#include "util_gaunt.hpp"

namespace pwdft {

/*****************************************************
 *                                                   *
 *           Paw_get_spher_grid                      *
 *                                                   *
 *****************************************************/
static void Paw_get_spher_grid(const int ntheta, const int nphi,
                               double angle_phi[], double cos_theta[],
                               double w_theta[], double w_phi[]) {
  double pi = 4.0 * atan(1.0);

  //*** gaussian quadrature angular grid for cos_theta ***
  util_gauss_weights(-1.0, 1.0, cos_theta, w_theta, ntheta);

  if (nphi > 1) {
    //***linear angular grid for angle_phi***
    for (auto i = 0; i < nphi; ++i) {
      angle_phi[i] = 2.0 * pi * (i - 1) / ((double)(nphi - 1));
      w_phi[i] = 2.0 * pi / ((double)(nphi - 1));
    }
    w_phi[0] *= 0.50;
    w_phi[nphi - 1] = w_phi[0];
  } else {
    angle_phi[0] = 0.0;
    w_phi[0] = 2.0 * pi;
  }
}

/*******************************************************
 *                                                     *
 *                 Paw_xc::Paw_xc                      *
 *                                                     *
 *******************************************************/

Paw_xc::Paw_xc(Ion *myion0, Pneb *mypneb0, Control2 &control, const int nprj[],
               const int nbasis[], const int n1dgrid[], const int psp_type[],
               const int lmax0[], int *l_prj[], int *m_prj[], int *b_prj[]) {
  int mtr_size;

  myion = myion0;
  mypneb = mypneb0;

  int nion = myion->nion;
  int nkatm = myion->nkatm;
  gga = control.get_gga();

  nprj_max = 0;
  nbasis_max = 0;
  n1dgrid_max = 0;
  for (int ia = 0; ia < nkatm; ++ia) {
    if (nprj[ia] > nprj_max)
      nprj_max = nprj[ia];
    if (nbasis[ia] > nbasis_max)
      nbasis_max = nbasis[ia];
    if (n1dgrid[ia] > n1dgrid_max)
      n1dgrid_max = n1dgrid[ia];
  }

  mult_l_max = control.lmax_multipole();

  if (mult_l_max < 0) {
    mult_l_max = 0;
    for (auto ia = 0; ia < nkatm; ++ia)
      if (psp_type[ia] == 4)
        if (mult_l_max < (2 * lmax0[ia]))
          mult_l_max = 2 * lmax0[ia];
  }

  //*** paw_xc energies ***
  paw_xc_e = new (std::nothrow) double[nion]();

  //*** xc matrix arrays - nbasis_max**2 * (mult_l_max+1)**2 ***
  mtr_size =
      2 * (nbasis_max * nbasis_max) * ((mult_l_max + 1) * (mult_l_max + 1));
  paw_xc_matr = new (std::nothrow) double[mtr_size]();
  if (gga >= 10)
    paw_xc_dmatr = new (std::nothrow) double[3 * mtr_size]();

  mtr_size = 2 * nprj_max * nprj_max;
  paw_xc_pot = new (std::nothrow) double[mtr_size]();

  //*** spherical grid arrays ***
  if (mult_l_max == 0) {
    paw_xc_nphi = 1;
    paw_xc_ntheta = 1;
  } else {
    paw_xc_nphi = 3 * mult_l_max;
    paw_xc_ntheta = 3 * mult_l_max;
  }

  paw_xc_angle_phi = new (std::nothrow) double[paw_xc_nphi]();
  paw_xc_cos_theta = new (std::nothrow) double[paw_xc_ntheta]();
  paw_xc_w_phi = new (std::nothrow) double[paw_xc_nphi]();
  paw_xc_w_theta = new (std::nothrow) double[paw_xc_ntheta]();
  paw_xc_tlm = new (std::nothrow) double[paw_xc_ntheta * paw_xc_nphi *
                                         (mult_l_max + 1) * (mult_l_max + 1)]();

  //**** used for generating derivatives of tlm's ****
  if (gga >= 10) {
    paw_xc_dtlm_theta =
        new (std::nothrow) double[paw_xc_ntheta * paw_xc_nphi *
                                  (mult_l_max + 1) * (mult_l_max + 1)]();
    paw_xc_dtlm_phi =
        new (std::nothrow) double[paw_xc_ntheta * paw_xc_nphi *
                                  (mult_l_max + 1) * (mult_l_max + 1)]();
  }

  Paw_get_spher_grid(paw_xc_ntheta, paw_xc_nphi, paw_xc_angle_phi,
                     paw_xc_cos_theta, paw_xc_w_theta, paw_xc_w_phi);

  //**** define tlm's ****
  int i_tlm = 0;
  for (auto i_t = 0; i_t < paw_xc_ntheta; ++i_t)
    for (auto i_p = 0; i_p < paw_xc_nphi; ++i_p) {
      for (auto l = 0; l <= mult_l_max; ++l)
        for (auto m = -l; m <= l; ++m) {
          double tmp_phi;
          double tmp_theta = util_rtheta_lm(l, m, paw_xc_cos_theta[i_t]);
          double angle_phi = paw_xc_angle_phi[i_p];
          if (m < 0) {
            tmp_phi = sin(std::abs(m) * angle_phi);
          } else if (m > 0) {
            tmp_phi = cos(std::abs(m) * angle_phi);
          } else {
            tmp_phi = 1.0;
          }
          paw_xc_tlm[i_tlm] = tmp_theta * tmp_phi;
          ++i_tlm;
        }
    }

  if (gga >= 10) {
    //**** define derivative wrt to theta ****
    i_tlm = 0;
    for (auto i_t = 0; i_t < paw_xc_ntheta; ++i_t)
      for (auto i_p = 0; i_p < paw_xc_nphi; ++i_p) {
        for (auto l = 0; l <= mult_l_max; ++l)
          for (auto m = -l; m <= l; ++m) {
            double tmp_phi;
            double tmp_theta = util_drtheta_lm(l, m, paw_xc_cos_theta[i_t - 1]);
            double angle_phi = paw_xc_angle_phi[i_p];
            if (m < 0) {
              tmp_phi = sin(std::abs(m) * angle_phi);
            } else if (m > 0) {
              tmp_phi = cos(std::abs(m) * angle_phi);
            } else {
              tmp_phi = 1.0;
            }
            paw_xc_dtlm_theta[i_tlm] = tmp_theta * tmp_phi;
            ++i_tlm;
          }
      }

    //**** define derivative wrt to phi ****
    i_tlm = 0;
    for (auto i_t = 0; i_t < paw_xc_ntheta; ++i_t)
      for (auto i_p = 0; i_p < paw_xc_nphi; ++i_p) {
        for (auto l = 0; l <= mult_l_max; ++l)
          for (auto m = -l; m <= l; ++m) {
            if (m == 0) {
              paw_xc_dtlm_phi[i_tlm] = 0.0;
            } else {
              double tmp_phi;
              double tmp_theta =
                  util_rtheta_lm_div(l, m, paw_xc_cos_theta[i_t]);
              double angle_phi = paw_xc_angle_phi[i_p];

              if (m < 0) {
                tmp_phi = abs(m) * cos(std::abs(m) * angle_phi);
              } else if (m > 0) {
                tmp_phi = -abs(m) * sin(std::abs(m) * angle_phi);
              } else {
                tmp_phi = 0.0;
              }
              paw_xc_dtlm_phi[i_tlm] = tmp_theta * tmp_phi;
            }
            ++i_tlm;
          }
      }
  }

  //*** temp arrays ***
  int nsize = 2 * n1dgrid_max;
  rho_ae = new (std::nothrow) double[nsize]();
  rho_ps = new (std::nothrow) double[nsize]();

  //**** allocate gradient's and agr's ****
  if (gga >= 10) {
    nsize = 2 * 3 * n1dgrid_max;
    rho_ae_prime = new (std::nothrow) double[nsize]();
    rho_ps_prime = new (std::nothrow) double[nsize]();

    nsize = (2 * 2 - 1) * n1dgrid_max;
    agr_ae = new (std::nothrow) double[nsize]();
    agr_ps = new (std::nothrow) double[nsize]();
    fdn_ae = new (std::nothrow) double[nsize]();
    fdn_ps = new (std::nothrow) double[nsize]();
  }

  nsize = 2 * n1dgrid_max;
  vxc_ae = new (std::nothrow) double[nsize]();
  vxc_ps = new (std::nothrow) double[nsize]();
  exc_ae = new (std::nothrow) double[nsize]();
  exc_ps = new (std::nothrow) double[nsize]();
  xc_temp = new (std::nothrow) double[nsize]();
  xc_temp2 = new (std::nothrow) double[nsize]();

  //*** allocate vxclm multipole expansion  arrays ****
  nsize = 2 * n1dgrid_max * (mult_l_max + 1) * (mult_l_max + 1);

  paw_vxc_ae = new (std::nothrow) double[nsize]();
  paw_vxc_ps = new (std::nothrow) double[nsize]();

  //*** allocate dvxclm multipole expansion  arrays ****
  if (gga >= 10) {
    paw_dvxc_ae = new (std::nothrow) double[3 * nsize]();
    paw_dvxc_ps = new (std::nothrow) double[3 * nsize]();
  }

  //*** allocate rholm multipole expansion arrays ****
  paw_rho2_ae = new (std::nothrow) double[nsize]();
  paw_rho2_ps = new (std::nothrow) double[nsize]();

  //*** allocate rholm_prime multipole expansion arrays - need for GGA's****
  if (gga >= 10) {
    paw_rho2_ae_prime = new (std::nothrow) double[nsize]();
    paw_rho2_ps_prime = new (std::nothrow) double[nsize]();
  }

  //***** set up index arrays for nwpw_xc_density_solve2 *****
  nindx_rholm = new (std::nothrow) int[nkatm]();
  shift_rholm = new (std::nothrow) int[nkatm]();

  nsize = 0;
  for (auto ia = 0; ia < nkatm; ++ia) {
    shift_rholm[ia] = nsize;
    for (auto jprj = 0; jprj < nprj[ia]; ++jprj) {
      int lj = l_prj[ia][jprj];
      int mj = m_prj[ia][jprj];
      for (auto iprj = 0; iprj < nprj[ia]; ++iprj) {
        int li = l_prj[ia][iprj];
        int mi = m_prj[ia][iprj];
        for (auto l = 0; l <= mult_l_max; ++l) {
          for (auto m = -l; m <= l; ++m) {
            if ((l <= (li + lj)) && (l >= abs(li - lj))) {
              double tmp_gaunt = util_gaunt(false, l, m, li, mi, lj, mj);
              if (std::abs(tmp_gaunt) > 1.0e-11)
                ++nsize;
            }
          }
        }
      }
    }
    nindx_rholm[ia] = nsize - shift_rholm[ia];
  }

  bi_rholm = new (std::nothrow) int[nsize]();
  bj_rholm = new (std::nothrow) int[nsize]();
  lm_rholm = new (std::nothrow) int[nsize]();
  iprj_rholm = new (std::nothrow) int[nsize]();
  jprj_rholm = new (std::nothrow) int[nsize]();
  coeff_rholm = new (std::nothrow) double[nsize]();
  coeff_rholm2 = new (std::nothrow) double[nsize]();
  coeff_rholm3 = new (std::nothrow) double[nsize]();

  nsize = 0;
  for (auto ia = 0; ia < nkatm; ++ia) {
    for (auto jprj = 0; jprj < nprj[ia]; ++jprj) {
      int lj = l_prj[ia][jprj];
      int mj = m_prj[ia][jprj];
      int bj = b_prj[ia][jprj];
      for (auto iprj = 0; iprj < nprj[ia]; ++iprj) {
        int li = l_prj[ia][iprj];
        int mi = m_prj[ia][iprj];
        int bi = b_prj[ia][iprj];
        int lm = 0;
        for (auto l = 0; l <= mult_l_max; ++l) {
          for (auto m = -l; m <= l; ++m) {
            if ((l <= (li + lj)) && (l >= abs(li - lj))) {
              double tmp_gaunt2, tmp_gaunt3;
              double tmp_gaunt = util_gaunt(false, l, m, li, mi, lj, mj);
              if ((gga >= 10) && (mult_l_max > 0)) {
                tmp_gaunt2 = util_gaunt2(false, l, m, li, mi, lj, mj);
                tmp_gaunt3 = util_gaunt3(false, l, m, li, mi, lj, mj);
              } else {
                tmp_gaunt2 = 0.0;
                tmp_gaunt3 = 0.0;
              }
              if ((std::abs(tmp_gaunt) > 1.0e-11) ||
                  (std::abs(tmp_gaunt2) > 1.0e-11) ||
                  (std::abs(tmp_gaunt3) > 1.0e-11)) {
                coeff_rholm[nsize] = tmp_gaunt;
                coeff_rholm2[nsize] = tmp_gaunt2;
                coeff_rholm3[nsize] = tmp_gaunt3;
                lm_rholm[nsize] = lm;
                iprj_rholm[nsize] = iprj;
                jprj_rholm[nsize] = jprj;
                bi_rholm[nsize] = bi;
                bj_rholm[nsize] = bj;
                ++nsize;
              }
            }
            ++lm;
          }
        }
      }
    }
  }
}

} // namespace pwdft

/* PGrid.cpp
   Author - Eric Bylaska

        this class is used for defining 3d parallel maps
*/

#include "Control2.hpp"
#include "Lattice.hpp"
#include "fft.h"
#include "gdevice.hpp"
#include "nwpw_timing.hpp"
#include "util.hpp"
#include <cmath>
#include <cstring> //memset
#include <iostream>

#include "PGrid.hpp"

#include "blas.h"
#define mytaskid 1

namespace pwdft {

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

// PGrid::PGrid(Parallel *inparall, Lattice *inlattice, Control2& control) :
// d3db(inparall,control.mapping(),control.ngrid(0),control.ngrid(1),control.ngrid(2))

PGrid::PGrid(Parallel *inparall, Lattice *inlattice, int mapping0, int balance0,
             int nx0, int ny0, int nz0, int pfft3_qsize0)
    : d3db(inparall, mapping0, nx0, ny0, nz0) {
  int nxh, nyh, nzh, p, q, indx, nb;
  int nwave_in[2], nwave_out[2];
  double *G1, *G2, *G3;
  double ggcut, eps, ggmax, ggmin;
  bool *zero_arow3, *zero_arow2;
  bool yzslab, zrow;

  lattice = inlattice;

  eps = 1.0e-12;
  Garray = new (std::nothrow) double[3 * nfft3d]();
  G1 = Garray;
  G2 = (double *)&Garray[nfft3d];
  G3 = (double *)&Garray[2 * nfft3d];
  nxh = nx / 2;
  nyh = ny / 2;
  nzh = nz / 2;
  ggmin = 9.9e9;
  ggmax = 0.0;
  for (auto k3 = (-nzh + 1); k3 <= nzh; ++k3)
    for (auto k2 = (-nyh + 1); k2 <= nyh; ++k2)
      for (auto k1 = 0; k1 <= nxh; ++k1) {
        auto gx = k1 * lattice->unitg(0, 0) + k2 * lattice->unitg(0, 1) +
                  k3 * lattice->unitg(0, 2);
        auto gy = k1 * lattice->unitg(1, 0) + k2 * lattice->unitg(1, 1) +
                  k3 * lattice->unitg(1, 2);
        auto gz = k1 * lattice->unitg(2, 0) + k2 * lattice->unitg(2, 1) +
                  k3 * lattice->unitg(2, 2);
        auto gg = gx * gx + gy * gy + gz * gz;
        if (gg > ggmax)
          ggmax = gg;
        if ((gg < ggmin) && (gg > 1.0e-6))
          ggmin = gg;
        auto i = k1;
        if (i < 0)
          i = i + nx;
        auto j = k2;
        if (j < 0)
          j = j + ny;
        auto k = k3;
        if (k < 0)
          k = k + nz;

        auto indx = ijktoindex(i, j, k);
        auto p = ijktop(i, j, k);
        if (p == parall->taskid_i()) {
          G1[indx] = gx;
          G2[indx] = gy;
          G3[indx] = gz;
        }
      }
  Gmax = sqrt(ggmax);
  Gmin = sqrt(ggmin);
  masker[0] = new (std::nothrow) int[2 * nfft3d]();
  masker[1] = (int *)&(masker[0][nfft3d]);
  parall->Barrier();
  for (int k = 0; k < (nfft3d); ++k) {
    masker[0][k] = 1;
    masker[1][k] = 1;
  }

  for (auto nb = 0; nb <= 1; ++nb) {
    nwave[nb] = 0;
    if (nb == 0)
      ggcut = lattice->eggcut();
    else
      ggcut = lattice->wggcut();

    for (auto k3 = (-nzh + 1); k3 < nzh; ++k3)
      for (auto k2 = (-nyh + 1); k2 < nyh; ++k2)
        for (auto k1 = 0; k1 < nxh; ++k1) {
          auto gx = k1 * lattice->unitg(0, 0) + k2 * lattice->unitg(0, 1) +
                    k3 * lattice->unitg(0, 2);
          auto gy = k1 * lattice->unitg(1, 0) + k2 * lattice->unitg(1, 1) +
                    k3 * lattice->unitg(1, 2);
          auto gz = k1 * lattice->unitg(2, 0) + k2 * lattice->unitg(2, 1) +
                    k3 * lattice->unitg(2, 2);
          auto i = k1;
          if (i < 0)
            i = i + nx;
          auto j = k2;
          if (j < 0)
            j = j + ny;
          auto k = k3;
          if (k < 0)
            k = k + nz;

          auto indx = ijktoindex(i, j, k);
          auto p = ijktop(i, j, k);
          if (p == parall->taskid_i()) {
            auto gg = gx * gx + gy * gy + gz * gz;
            gg = gg - ggcut;
            if (gg < (-eps)) {
              masker[nb][indx] = 0;
              ++nwave[nb];
            }
          }
        }
    nwave_entire[nb] = nwave[nb];
    nwave_entire[nb] = parall->ISumAll(1, nwave_entire[nb]);
  }

  packarray[0] = new (std::nothrow) int[2 * nfft3d]();
  packarray[1] = (int *)&(packarray[0][nfft3d]);
  // packarray[1] = new (std::nothrow) int [nfft3d]();

  for (auto nb = 0; nb <= 1; ++nb) {
    nida[nb] = 0;
    nidb2[nb] = 0;

    /* k=(0,0,0)  */
    auto k1 = 0;
    auto k2 = 0;
    auto k3 = 0;
    auto indx = ijktoindex(k1, k2, k3);
    auto p = ijktop(k1, k2, k3);
    if (p == parall->taskid_i())
      if (!masker[nb][indx]) {
        packarray[nb][nida[nb]] = indx;
        ++nida[nb];
      }

    // k=(0,0,k3) - neglect (0,0,-k3) points
    for (auto k = 1; k <= (nzh - 1); ++k) {
      k1 = 0;
      k2 = 0;
      k3 = k;
      indx = ijktoindex(k1, k2, k3);
      p = ijktop(k1, k2, k3);
      if (p == parall->taskid_i())
        if (!masker[nb][indx]) {
          packarray[nb][nida[nb] + nidb2[nb]] = indx;
          ++nidb2[nb];
        }
    }

    // k=(0,k2,k3) - neglect (0,-k2, -k3) points
    for (auto k = (-nzh + 1); k <= (nzh - 1); ++k)
      for (auto j = 1; j <= (nyh - 1); ++j) {
        k1 = 0;
        k2 = j;
        k3 = k;
        if (k3 < 0)
          k3 = k3 + nz;
        indx = ijktoindex(k1, k2, k3);
        p = ijktop(k1, k2, k3);
        if (p == parall->taskid_i())
          if (!masker[nb][indx]) {
            packarray[nb][nida[nb] + nidb2[nb]] = indx;
            ++nidb2[nb];
          }
      }

    // k=(k1,k2,k3)
    for (auto k = (-nzh + 1); k <= (nzh - 1); ++k)
      for (auto j = (-nyh + 1); j <= (nyh - 1); ++j)
        for (auto i = 1; i <= (nxh - 1); ++i) {
          k1 = i;
          k2 = j;
          k3 = k;
          if (k2 < 0)
            k2 = k2 + ny;
          if (k3 < 0)
            k3 = k3 + nz;
          indx = ijktoindex(k1, k2, k3);
          p = ijktop(k1, k2, k3);
          if (p == parall->taskid_i())
            if (!masker[nb][indx]) {
              packarray[nb][nida[nb] + nidb2[nb]] = indx;
              ++nidb2[nb];
            }
        }
  }

  nwave_in[0] = nida[0] + nidb2[0];
  nwave_in[1] = nida[1] + nidb2[1];

  // if (control.balance())
  if (balance0) {
    balanced = 1;
    mybalance = new Balance(parall, 2, nwave_in, nwave_out);
  } else {
    balanced = 0;
    nwave_out[0] = nwave_in[0];
    nwave_out[1] = nwave_in[1];
  }
  nidb[0] = nidb2[0] + (nwave_out[0] - nwave_in[0]);
  nidb[1] = nidb2[1] + (nwave_out[1] - nwave_in[1]);

  nwave_all[0] = nida[0] + nidb2[0];
  nwave_all[1] = nida[1] + nidb2[1];
  parall->Vector_ISumAll(1, 2, nwave_all);

  if (maptype == 1) {

    zero_row3[0] = new (std::nothrow) bool[(nxh + 1) * nq];
    zero_row3[1] = new (std::nothrow) bool[(nxh + 1) * nq];
    zero_row2[0] = new (std::nothrow) bool[(nxh + 1) * nq];
    zero_row2[1] = new (std::nothrow) bool[(nxh + 1) * nq];
    zero_slab23[0] = new (std::nothrow) bool[nxh + 1];
    zero_slab23[1] = new (std::nothrow) bool[nxh + 1];

    zero_arow3 = new bool[(nxh + 1) * ny];
    for (auto nb = 0; nb <= 1; ++nb) {
      if (nb == 0)
        ggcut = lattice->eggcut();
      else
        ggcut = lattice->wggcut();

      /* find zero_row3 - (i,j,*) rows that are zero */
      for (auto i = 0; i < ((nxh + 1) * nq); ++i)
        zero_row3[nb][i] = true;
      for (auto i = 0; i < ((nxh + 1) * ny); ++i)
        zero_arow3[i] = true;

      for (auto k2 = (-nyh + 1); k2 < nyh; ++k2)
        for (auto k1 = 0; k1 < nxh; ++k1) {
          auto i = k1;
          auto j = k2;
          if (i < 0)
            i = i + nx;
          if (j < 0)
            j = j + ny;
          zrow = true;
          for (auto k3 = (-nzh + 1); k3 < nzh; ++k3) {
            auto gx = k1 * lattice->unitg(0, 0) + k2 * lattice->unitg(0, 1) +
                      k3 * lattice->unitg(0, 2);
            auto gy = k1 * lattice->unitg(1, 0) + k2 * lattice->unitg(1, 1) +
                      k3 * lattice->unitg(1, 2);
            auto gz = k1 * lattice->unitg(2, 0) + k2 * lattice->unitg(2, 1) +
                      k3 * lattice->unitg(2, 2);
            auto gg = gx * gx + gy * gy + gz * gz;
            gg = gg - ggcut;
            if (gg < (-eps))
              zrow = false;
          }
          if (!zrow) {
            zero_arow3[i + (nxh + 1) * j] = false;
            q = ijktoq1(0, j, 0);
            p = ijktop1(0, j, 0);
            if (p == parall->taskid_i()) {
              zero_row3[nb][i + (nxh + 1) * q] = false;
            }
          }
        }
      // call D3dB_c_ptranspose_jk_init(nb,log_mb(zero_arow3(1)))
      d3db::c_ptranspose_jk_init(nb, zero_arow3);

      /* find zero_slab23 - (i,*,*) slabs that are zero */
      for (auto i = 0; i < nxh; ++i)
        zero_slab23[nb][i] = true;

      for (auto k1 = 0; k1 < nxh; ++k1) {
        auto i = k1;
        if (i < 0)
          i = i + nx;
        yzslab = true;
        for (auto k3 = (-nzh + 1); k3 < nzh; ++k3)
          for (auto k2 = (-nyh + 1); k2 < nyh; ++k2) {
            auto gx = k1 * lattice->unitg(0, 0) + k2 * lattice->unitg(0, 1) +
                      k3 * lattice->unitg(0, 2);
            auto gy = k1 * lattice->unitg(1, 0) + k2 * lattice->unitg(1, 1) +
                      k3 * lattice->unitg(1, 2);
            auto gz = k1 * lattice->unitg(2, 0) + k2 * lattice->unitg(2, 1) +
                      k3 * lattice->unitg(2, 2);
            auto gg = gx * gx + gy * gy + gz * gz;
            gg = gg - ggcut;
            if (gg < (-eps))
              yzslab = false;
          }
        if (!yzslab)
          zero_slab23[nb][i] = false;
      }

      /* find zero_row2 - (i,*,k) rows that are zero after fft of (i,j,*) */
      for (auto k = 0; k < nz; ++k)
        for (auto i = 0; i < (nxh + 1); ++i) {
          q = ijktoq(i, 0, k);
          p = ijktop(i, 0, k);
          if (p == parall->taskid_i())
            zero_row2[nb][q] = zero_slab23[nb][i];
        }
    }

    delete[] zero_arow3;
  } else {
    zero_row3[0] = new (std::nothrow) bool[nq3]();
    zero_row3[1] = new (std::nothrow) bool[nq3]();
    zero_row2[0] = new (std::nothrow) bool[nq2]();
    zero_row2[1] = new (std::nothrow) bool[nq2]();
    zero_slab23[0] = new (std::nothrow) bool[nxh + 1]();
    zero_slab23[1] = new (std::nothrow) bool[nxh + 1]();

    zero_arow2 = new (std::nothrow) bool[(nxh + 1) * nz]();
    zero_arow3 = new (std::nothrow) bool[(nxh + 1) * ny]();

    for (auto nb = 0; nb <= 1; ++nb) {
      if (nb == 0)
        ggcut = lattice->eggcut();
      else
        ggcut = lattice->wggcut();

      // find zero_row3 - (i,j,*) rows that are zero
      for (auto q = 0; q < nq3; ++q)
        zero_row3[nb][q] = true;

      for (auto q = 0; q < (nxh + 1) * ny; ++q)
        zero_arow3[q] = true;

      for (auto k2 = (-nyh + 1); k2 < nyh; ++k2)
        for (auto k1 = 0; k1 < nxh; ++k1) {
          auto i = k1;
          auto j = k2;
          if (i < 0)
            i = i + nx;
          if (j < 0)
            j = j + ny;
          zrow = true;
          for (auto k3 = (-nzh + 1); k3 < nzh; ++k3) {
            auto gx = k1 * lattice->unitg(0, 0) + k2 * lattice->unitg(0, 1) +
                      k3 * lattice->unitg(0, 2);
            auto gy = k1 * lattice->unitg(1, 0) + k2 * lattice->unitg(1, 1) +
                      k3 * lattice->unitg(1, 2);
            auto gz = k1 * lattice->unitg(2, 0) + k2 * lattice->unitg(2, 1) +
                      k3 * lattice->unitg(2, 2);
            auto gg = gx * gx + gy * gy + gz * gz;
            gg = gg - ggcut;
            if (gg < (-eps))
              zrow = false;
          }
          if (!zrow) {
            // zero_arow3[i-1+(nxh+1)*(j-1)] = 0;
            zero_arow3[i + (nxh + 1) * j] = false;
            q = ijktoq(i, j, 0);
            p = ijktop(i, j, 0);
            if (p == parall->taskid_i()) {
              zero_row3[nb][q] = false;
            }
          }
        }

      /* find zero_slab23 - (i,*,*) slabs that are zero */
      for (auto i = 0; i < (nxh + 1); ++i)
        zero_slab23[nb][i] = true;

      for (auto k1 = 0; k1 < nxh; ++k1) {
        auto i = k1;
        if (i < 0)
          i = i + nx;
        yzslab = true;
        for (auto k3 = (-nzh + 1); k3 < nzh; ++k3)
          for (auto k2 = (-nyh + 1); k2 < nyh; ++k2) {
            auto gx = k1 * lattice->unitg(0, 0) + k2 * lattice->unitg(0, 1) +
                      k3 * lattice->unitg(0, 2);
            auto gy = k1 * lattice->unitg(1, 0) + k2 * lattice->unitg(1, 1) +
                      k3 * lattice->unitg(1, 2);
            auto gz = k1 * lattice->unitg(2, 0) + k2 * lattice->unitg(2, 1) +
                      k3 * lattice->unitg(2, 2);
            auto gg = gx * gx + gy * gy + gz * gz;
            gg = gg - ggcut;
            if (gg < (-eps))
              yzslab = false;
          }
        if (!yzslab)
          zero_slab23[nb][i] = false;
      }

      // find zero_row2 - (i,*,k) rows that are zero after fft of (i,j,*)
      for (auto k = 0; k < nz; ++k)
        for (auto i = 0; i < (nxh + 1); ++i) {
          q = ijktoq1(i, 0, k);
          p = ijktop1(i, 0, k);
          zero_arow2[i + (nxh + 1) * k] = zero_slab23[nb][i];
          if (p == parall->taskid_i())
            zero_row2[nb][q] = zero_slab23[nb][i];
        }

      // call
      // D3dB_c_ptranspose_ijk_init(nb,log_mb(zero_arow2(1)),log_mb(zero_arow3(1)))
      d3db::c_ptranspose_ijk_init(nb, zero_arow2, zero_arow3);
    }
    delete[] zero_arow3;
    delete[] zero_arow2;
  }

  Gpack[0] = new (
      std::nothrow) double[3 * (nida[0] + nidb[0]) + 3 * (nida[1] + nidb[1])]();
  Gpack[1] = (double *)&(Gpack[0][3 * (nida[0] + nidb[0])]);

  double *Gtmp = new (std::nothrow) double[nfft3d]();
  int one = 1;
  for (auto nb = 0; nb <= 1; ++nb) {
    for (auto i = 0; i < 3; ++i) {
      DCOPY_PWDFT(nfft3d, &(Garray[i * nfft3d]), one, Gtmp, one);

      this->t_pack(nb, Gtmp);
      this->tt_pack_copy(nb, Gtmp, &(Gpack[nb][i * (nida[nb] + nidb[nb])]));
    }
  }

  delete[] Gtmp;

  zplane_tmp1 = new (std::nothrow) double[2 * zplane_size + 8];
  zplane_tmp2 = new (std::nothrow) double[2 * zplane_size + 8];

  /* initialize r_grid */
  has_r_grid = (lattice->aperiodic());
  if (has_r_grid) {
    r_grid = r_nalloc(3);
    r_nzero(3, r_grid);
    double a[9];
    for (auto i = 0; i < 3; ++i) {
      a[i] = lattice->unita1d(0 + i) / ((double)nx);
      a[3 + i] = lattice->unita1d(3 + i) / ((double)ny);
      a[6 + i] = lattice->unita1d(6 + i) / ((double)nz);
    }

    /* grid points in coordination space */
    for (auto k3 = (-nzh); k3 < nzh; ++k3)
      for (auto k2 = (-nyh); k2 < nyh; ++k2)
        for (auto k1 = (-nxh); k1 < nxh; ++k1) {
          int i = k1 + nxh;
          int j = k2 + nyh;
          int k = k3 + nzh;
          int indx = ijktoindex2(i, j, k);
          int p = ijktop2(i, j, k);

          if (p == parall->taskid_i()) {
            r_grid[3 * indx] = a[0] * k1 + a[3] * k2 + a[6] * k3;
            r_grid[3 * indx + 1] = a[1] * k1 + a[4] * k2 + a[7] * k3;
            r_grid[3 * indx + 2] = a[2] * k1 + a[5] * k2 + a[8] * k3;
          }
        }
  }

  /* initialize pfft3 queues */
  aqmax = pfft3_qsize0;
  aqsize = 0;
  alast_index = aqmax - 1;
  aqindx = new (std::nothrow) int[aqmax]();
  aqstatus = new (std::nothrow) int[aqmax]();
  atmp = new (std::nothrow) double[2 * aqmax * n2ft3d]();

  bqmax = pfft3_qsize0;
  bqsize = 0;
  blast_index = bqmax - 1;
  bqindx = new (std::nothrow) int[bqmax]();
  bqstatus = new (std::nothrow) int[bqmax]();
  btmp = new (std::nothrow) double[2 * bqmax * n2ft3d]();

  /* initialize async buffer data for pfft */
  for (auto q = 0; q < aqmax; ++q)
    parall->astart(3 + q, 2 * parall->np_i() + 1);
}

PGrid::PGrid(Parallel *inparall, Lattice *inlattice, Control2 &control)
    : PGrid(inparall, inlattice, control.mapping(), control.balance(),
            control.ngrid(0), control.ngrid(1), control.ngrid(2),
            control.pfft3_qsize()) {}

/********************************
 *                              *
 *       PGrid:c_unpack         *
 *                              *
 ********************************/
void PGrid::c_unpack(const int nb, double *a) {
  int one = 1;
  int nn = 2 * (nida[nb] + nidb2[nb]);
  double *tmp1, *tmp2;
  double *tmp = new (std::nothrow) double[2 * nfft3d];
  if (balanced)
    mybalance->c_unbalance(nb, a);

  // DCOPY_PWDFT(nn,a,one,tmp,one);
  // dcopy_(&n2ft3d,&rzero,&zero,a,&one);
  std::memcpy(tmp, a, nn * sizeof(double));
  std::memset(a, 0, 2 * nfft3d * sizeof(double));
  c_bindexcopy(nida[nb] + nidb2[nb], packarray[nb], tmp, a);
  // c_bindexcopy(nida[nb]+nidb[nb],packarray[nb],tmp,a);

  // tmp1 = new (std::nothrow) double[2*zplane_size+1];
  // tmp2 = new (std::nothrow) double[2*zplane_size+1];
  c_timereverse(a, zplane_tmp1, zplane_tmp2);
  // delete [] tmp2;
  // delete [] tmp1;
  delete[] tmp;
}

/********************************
 *                              *
 *       PGrid:c_pack           *
 *                              *
 ********************************/
void PGrid::c_pack(const int nb, double *a) {
  int one = 1;
  double *tmp = new (std::nothrow) double[2 * nfft3d];

  // DCOPY_PWDFT(n2ft3d,a,one,tmp,one);
  std::memcpy(tmp, a, 2 * nfft3d * sizeof(double));

  std::memset(a, 0, 2 * nfft3d * sizeof(double));

  c_aindexcopy(nida[nb] + nidb2[nb], packarray[nb], tmp, a);

  if (balanced)
    mybalance->c_balance(nb, a);

  delete[] tmp;
  return;
}

/********************************
 *                              *
 *       PGrid:cc_pack_copy     *
 *                              *
 ********************************/
void PGrid::cc_pack_copy(const int nb, const double *a, double *b) 
{
   int one = 1;
   // int ng  = 2*(nida[nb]+nidb[nb]);
   int ng = 2*(nida[nb]+nidb[nb]);

   // DCOPY_PWDFT(ng,a,one,b,one);
   std::memcpy(b,a,ng*sizeof(double));
}

/********************************
 *                              *
 *       PGrid:cc_pack_dot      *
 *                              *
 ********************************/
double PGrid::cc_pack_dot(const int nb, double *a, double *b) {
  int one = 1;
  // int ng  = 2*(nida[nb]+nidb[nb]);
  int ng = 2 * (nida[nb] + nidb[nb]);
  int ng0 = 2 * nida[nb];
  double tsum;

  tsum = 2.0 * DDOT_PWDFT(ng, a, one, b, one);
  tsum -= DDOT_PWDFT(ng0, a, one, b, one);

  return d3db::parall->SumAll(1, tsum);
}

/********************************
 *                              *
 *       PGrid:tt_pack_dot      *
 *                              *
 ********************************/
double PGrid::tt_pack_dot(const int nb, double *a, double *b) {
  int one = 1;
  int ng = (nida[nb] + nidb[nb]);
  int ng0 = nida[nb];
  double tsum;

  tsum = 2.0 * DDOT_PWDFT(ng, a, one, b, one);
  tsum -= DDOT_PWDFT(ng0, a, one, b, one);

  return d3db::parall->SumAll(1, tsum);
}

/********************************
 *                              *
 *       PGrid:cc_pack_idot     *
 *                              *
 ********************************/
double PGrid::cc_pack_idot(const int nb, double *a, double *b) {
  int one = 1;
  // int ng  = 2*(nida[nb]+nidb[nb]);
  int ng = 2 * (nida[nb] + nidb[nb]);
  int ng0 = 2 * nida[nb];
  double tsum = 0.0;

  tsum = 2.0 * DDOT_PWDFT(ng, a, one, b, one);
  tsum -= DDOT_PWDFT(ng0, a, one, b, one);

  return tsum;
}

/********************************
 *                              *
 *       PGrid:tt_pack_idot     *
 *                              *
 ********************************/
double PGrid::tt_pack_idot(const int nb, double *a, double *b) {
  int one = 1;
  // int ng  = 2*(nida[nb]+nidb[nb]);
  int ng = (nida[nb] + nidb[nb]);
  int ng0 = nida[nb];
  double tsum = 0.0;

  tsum = 2.0 * DDOT_PWDFT(ng, a, one, b, one);
  tsum -= DDOT_PWDFT(ng0, a, one, b, one);

  return tsum;
}

/********************************
 *                              *
 *       PGrid:cc_pack_indot    *
 *                              *
 ********************************/
void PGrid::cc_pack_indot(const int nb, const int nn, double *a, double *b,
                          double *sum) {
  int one = 1;
  // int ng  = 2*(nida[nb]+nidb[nb]);
  int ng = 2 * (nida[nb] + nidb[nb]);
  int ng0 = 2 * nida[nb];

  for (int i = 0; i < nn; ++i) {
    sum[i] = 2.0 * DDOT_PWDFT(ng, &(a[i * ng]), one, b, one);
    sum[i] -= DDOT_PWDFT(ng0, &(a[i * ng]), one, b, one);
  }
}

/********************************
 *                              *
 *    PGrid:cc_pack_inprjdot    *
 *                              *
 ********************************/
void PGrid::cc_pack_inprjdot(const int nb, int nn, int nprj, double *a,
                             double *b, double *sum) {
  int ng = 2 * (nida[nb] + nidb[nb]);
  int ng0 = 2 * nida[nb];
  int one = 1;
  double rtwo = 2.0;
  double rone = 1.0;
  double rmone = -1.0;
  double rzero = 0.0;

  // DGEMM_PWDFT((char *) "T",(char *) "N",nn,nprj,ng,
  // 	       rtwo,
  // 	       a,ng,
  // 	       b,ng,
  // 	       rzero,
  // 	       sum,nn);
  gdevice_TN_dgemm(nn, nprj, ng, rtwo, a, b, rzero, sum);

  if (ng0 > 0) {
    DGEMM_PWDFT((char *)"T", (char *)"N", nn, nprj, ng0, rmone, a, ng, b, ng,
                rone, sum, nn);
  }
}

/********************************
 *                              *
 *       PGrid:t_unpack         *
 *                              *
 ********************************/
void PGrid::t_unpack(const int nb, double *a) {
  int one = 1;
  int nn = (nida[nb] + nidb2[nb]);
  double *tmp1, *tmp2;
  double *tmp = new (std::nothrow) double[nfft3d];
  if (balanced)
    mybalance->t_unbalance(nb, a);

  // DCOPY_PWDFT(nn,a,one,tmp,one);
  std::memcpy(tmp, a, nn * sizeof(double));

  // dcopy_(&n2ft3d,&rzero,&zero,a,&one);
  std::memset(a, 0, nfft3d * sizeof(double));

  t_bindexcopy(nida[nb] + nidb2[nb], packarray[nb], tmp, a);

  // tmp1 = new (std::nothrow) double[2*zplane_size+1];
  // tmp2 = new (std::nothrow) double[2*zplane_size+1];
  t_timereverse(a, zplane_tmp1, zplane_tmp2);
  // delete [] tmp2;
  // delete [] tmp1;
  delete[] tmp;
}

/********************************
 *                              *
 *       PGrid:t_pack           *
 *                              *
 ********************************/
void PGrid::t_pack(const int nb, double *a) {
  int one = 1;
  int zero = 0;
  double rzero = 0.0;
  double *tmp = new (std::nothrow) double[nfft3d];

  // DCOPY_PWDFT(nfft3d,a,one,tmp,one);
  std::memcpy(tmp, a, nfft3d * sizeof(double));

  // dcopy_(&nfft3d,&rzero,&zero,a,&one);
  std::memset(a, 0, nfft3d * sizeof(double));

  t_aindexcopy(nida[nb] + nidb2[nb], packarray[nb], tmp, a);

  if (balanced)
    mybalance->t_balance(nb, a);

  delete[] tmp;
}

/********************************
 *                              *
 *       PGrid:tt_pack_copy     *
 *                              *
 ********************************/
void PGrid::tt_pack_copy(const int nb, double *a, double *b) {
  int one = 1;
  int ng = nida[nb] + nidb[nb];

  // DCOPY_PWDFT(ng,a,one,b,one);
  std::memcpy(b, a, ng * sizeof(double));
}

/********************************
 *                              *
 *       PGrid:t_pack_nzero     *
 *                              *
 ********************************/
void PGrid::t_pack_nzero(const int nb, const int n, double *a) {
  // int one   = 1;
  // int zero = 0;
  // double azero = 0.0;
  int ng = n * (nida[nb] + nidb[nb]);
  std::memset(a, 0, ng * sizeof(double));
}

/********************************
 *                              *
 *       PGrid:i_pack           *
 *                              *
 ********************************/
void PGrid::i_pack(const int nb, int *a) {
  int i;
  int *tmp = new (std::nothrow) int[nfft3d];

  for (i = 0; i < nfft3d; ++i)
    tmp[i] = a[i];
  for (i = 0; i < nfft3d; ++i)
    a[i] = 0;

  i_aindexcopy(nida[nb] + nidb2[nb], packarray[nb], tmp, a);

  if (balanced)
    mybalance->i_balance(nb, a);

  delete[] tmp;
}

/********************************
 *                              *
 *       PGrid:ii_pack_copy     *
 *                              *
 ********************************/
void PGrid::ii_pack_copy(const int nb, int *a, int *b) {
  int i;
  int ng = nida[nb] + nidb[nb];
  for (i = 0; i < ng; ++i)
    b[i] = a[i];
}

/********************************
 *                              *
 *    PGrid:cr_pfft3b           *
 *                              *
 ********************************/
/* This routine performs the operation of a three dimensional
  complex to complex inverse fft:
   A(nx,ny(nb),nz(nb)) <- FFT3^(-1)[A(kx,ky,kz)]
*/
void PGrid::cr_pfft3b(const int nb, double *a)
{

  nwpw_timing_function ftime(1);
  int i, j, k, jj, kk, q, indx, indx0, indx2, nxh, nxh2, nxhy, nxhy2, nxhz,
      nxhz2, nn, shift;
  double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
  double *tmp3 = new (std::nothrow) double[2 * nfft3d]();

  nxh = nx / 2 + 1;
  nxhy = nxh * ny;
  nxhz = nxh * nz;
  nxh2 = nx + 2;
  nxhy2 = nxh2 * ny;
  nxhz2 = nxh2 * nz;

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) {
    /***************************************************
     ***     do fft along kz dimension               ***
     ***     A(kx,nz,ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
     ***************************************************/
    indx0 = 0;
    indx2 = 0;
    nn = 0;
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        if (!zero_row3[nb][indx2]) {
          kk = 0;
          indx = 2 * i + indx0;
          shift = 2 * nz * nn;
          for (k = 0; k < nz; ++k) {
            tmp2[kk + shift] = a[indx];
            tmp2[kk + 1 + shift] = a[indx + 1];
            kk += 2;
            indx += nxh2;
          }
          ++nn;
        }
        ++indx2;
      }
      indx0 += nxhz2;
    }

    gdevice_batch_cfftz_tmpz(false, nz, nn, n2ft3d, tmp2, d3db::tmpz);

    indx0 = 0;
    indx2 = 0;
    nn = 0;
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        if (!zero_row3[nb][indx2]) {
          kk = 0;
          indx = 2 * i + indx0;
          shift = 2 * nz * nn;
          for (k = 0; k < nz; ++k) {
            a[indx] = tmp2[kk + shift];
            a[indx + 1] = tmp2[kk + 1 + shift];
            kk += 2;
            indx += nxh2;
          }
          ++nn;
        }
        ++indx2;
      }
      indx0 += nxhz2;
    }

    /*
    indx0=0;
    indx2=0;
    for (q=0; q<nq; ++q)
    {
       for (i=0; i<nxh; ++i)
       {
          if (!zero_row3[nb][indx2])
          {
             kk   = 0;
             indx = 2*i+indx0;
             for (k=0; k<nz; ++k)
             {
                tmp2[kk]   = a[indx];
                tmp2[kk+1] = a[indx+1];
                kk   += 2;
                indx += nxh2;
             }

             dcfftb_(&nz,tmp2,d3db::tmpz);

             kk   = 0;
             indx = 2*i+indx0;
             for (k=0; k<nz; ++k)
             {
                a[indx]   = tmp2[kk];
                a[indx+1] = tmp2[kk+1];
                kk   += 2;
                indx += nxh2;
             }
          }
          ++indx2;
       }
       indx0 += nxhz2;
    }
    */

    /***********************************************
     ***         Do a ptranspose of A            ***
     ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
     ************************************************/
    d3db::c_ptranspose1_jk(nb, a, tmp2, tmp3);

    /*************************************************
     ***        do fft along ky dimension          ***
     ***    A(kx,ny,nz) <- fft1d^(-1)[A(kx,ky,nz)] ***
     *************************************************/
    indx0 = 0;
    indx2 = 0;
    nn = 0;
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        if (!zero_row2[nb][indx2]) {
          jj = 0;
          indx = 2 * i + indx0;
          shift = 2 * ny * nn;
          for (j = 0; j < ny; ++j) {
            tmp2[jj + shift] = a[indx];
            tmp2[jj + 1 + shift] = a[indx + 1];
            jj += 2;
            indx += nxh2;
          }
          ++nn;
        }
        ++indx2;
      }
      indx0 += nxhy2;
    }

    gdevice_batch_cffty_tmpy(false, ny, nn, n2ft3d, tmp2, d3db::tmpy);

    indx0 = 0;
    indx2 = 0;
    nn = 0;
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        if (!zero_row2[nb][indx2]) {
          jj = 0;
          indx = 2 * i + indx0;
          shift = 2 * ny * nn;
          for (j = 0; j < ny; ++j) {
            a[indx] = tmp2[jj + shift];
            a[indx + 1] = tmp2[jj + 1 + shift];
            jj += 2;
            indx += nxh2;
          }
          ++nn;
        }
        ++indx2;
      }
      indx0 += nxhy2;
    }

    /*
    indx0=0;
    indx2=0;
    for (q=0; q<nq; ++q)
    {
       for (i=0; i<nxh; ++i)
       {
          if (!zero_row2[nb][indx2])
          {
             jj   = 0;
             indx = 2*i+indx0;
             for (j=0; j<ny; ++j)
             {
                tmp2[jj]   = a[indx];
                tmp2[jj+1] = a[indx+1];
                jj   += 2;
                indx += nxh2;
             }

             dcfftb_(&ny,tmp2,d3db::tmpy);

             jj   = 0;
             indx = 2*i+indx0;
             for (j=0; j<ny; ++j)
             {
                a[indx]   = tmp2[jj];
                a[indx+1] = tmp2[jj+1];
                jj   += 2;
                indx += nxh2;
             }
          }
          ++indx2;
       }
       indx0 += nxhy2;
    }
    */

    /************************************************
     ***     do fft along kx dimension            ***
     ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
     ************************************************/
    /*
     d3db::cshift1_fftb(nx,ny,nq,1,a);
     indx = 0;
     for (q=0; q<nq; ++q)
     for (j=0; j<ny; ++j)
     {
        drfftb_(&nx,a+indx,d3db::tmpx);
        indx += nxh2;
     }
     d3db::zeroend_fftb(nx,ny,nq,1,a);
    */
    gdevice_batch_cfftx_tmpx(false, nx, ny * nq, n2ft3d, a, d3db::tmpx);
    d3db::zeroend_fftb(nx, ny, nq, 1, a);

  }

  /*************************
   **** hilbert mapping ****
   *************************/
  else {

    /************************************************
     ***     do fft along kz dimension            ***
     ***   A(nz,kx,ky) <- fft1d^(-1)[A(kz,kx,ky)] ***
     ************************************************/
    gdevice_batch_cfftz_tmpz_zero(false, nz, nq3, n2ft3d, a, d3db::tmpz,
                                  zero_row3[nb]);
    /*
#if (defined NWPW_SYCL) || (defined NWPW_CUDA)
    gdevice_batch_cfftz(false,nz,nq3,n2ft3d,a);
#else
    indx = 0;
    for (q=0; q<nq3; ++q)
    {
       if (!zero_row3[nb][q])
          dcfftb_(&nz,a+indx,d3db::tmpz);
       indx += 2*nz;
    }
#endif
    */
    d3db::c_ptranspose_ijk(nb, 2, a, tmp2, tmp3);

    /************************************************
     ***     do fft along ky dimension            ***
     ***   A(ny,nz,kx) <- fft1d^(-1)[A(ky,nz,kx)] ***
     ************************************************/
    gdevice_batch_cffty_tmpy_zero(false, ny, nq2, n2ft3d, a, d3db::tmpy,
                                  zero_row2[nb]);
    /*
#if (defined NWPW_SYCL) || (defined NWPW_CUDA)
    gdevice_batch_cffty(false,ny,nq2,n2ft3d,a);
#else
    indx = 0;
    for (q=0; q<nq2; ++q)
    {
       if (!zero_row2[nb][q])
          dcfftb_(&ny,a+indx,d3db::tmpy);
       indx += (2*ny);
    }
#endif
    */
    d3db::c_ptranspose_ijk(nb, 3, a, tmp2, tmp3);

    /************************************************
     ***     do fft along kx dimension            ***
     ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
     ************************************************/
    gdevice_batch_cfftx_tmpx(false, nx, nq1, n2ft3d, a, d3db::tmpx);
    /*
    #if (defined NWPW_SYCL) || (defined NWPW_CUDA)
           gdevice_batch_cfftx(false,nx,nq1,n2ft3d,a);
    #else
           d3db::cshift1_fftb(nx,nq1,1,1,a);
           indx = 0;
           for (q=0; q<nq1; ++q)
           {
              drfftb_(&nx,a+indx,d3db::tmpx);
              indx += nxh2;
           }
    #endif
     */
    d3db::zeroend_fftb(nx, nq1, 1, 1, a);
    if (n2ft3d_map < n2ft3d)
      std::memset(a + n2ft3d_map, 0, (n2ft3d - n2ft3d_map) * sizeof(double));
  }

  delete[] tmp3;
  delete[] tmp2;
}

/********************************
 *                              *
 *       PGrid::rc_pfft3f       *
 *                              *
 ********************************/
/*
   This routine performs the operation of a three dimensional
   complex to complex fft:
      A(kx,ky,kz) <- FFT3[A(nx(nb),ny(nb),nz(nb))]
*/
void PGrid::rc_pfft3f(const int nb, double *a) 
{
  nwpw_timing_function ftime(1);
  int i, j, k, jj, kk, q, indx, indx0, indx2, nxh, nxh2, nxhy, nxhy2, nxhz,
      nxhz2, nn, shift;
  double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
  double *tmp3 = new (std::nothrow) double[2 * nfft3d]();

  nxh = nx / 2 + 1;
  nxhy = nxh * ny;
  nxhz = nxh * nz;
  nxh2 = nx + 2;
  nxhy2 = nxh2 * ny;
  nxhz2 = nxh2 * nz;

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) {
    /********************************************
     ***     do fft along nx dimension        ***
     ***   A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]  ***
     ********************************************/
    /*
    indx = 0;
    for (q=0; q<nq; ++q)
    for (j=0; j<ny; ++j)
    {
       drfftf_(&nx,a+indx,tmpx);
       indx += nxh2;
    }
    d3db::cshift_fftf(nx,ny,nq,1,a);
    */
    gdevice_batch_cfftx_tmpx(true, nx, ny * nq, n2ft3d, a, d3db::tmpx);

    /********************************************
     ***     do fft along ny dimension        ***
     ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
     ********************************************/
    indx0 = 0;
    indx2 = 0;
    nn = 0;
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        if (!zero_row2[nb][indx2]) {
          jj = 0;
          indx = 2 * i + indx0;
          shift = 2 * ny * nn;
          for (j = 0; j < ny; ++j) {
            tmp2[jj + shift] = a[indx];
            tmp2[jj + 1 + shift] = a[indx + 1];
            jj += 2;
            indx += nxh2;
          }
          ++nn;
        }
        ++indx2;
      }
      indx0 += nxhy2;
    }

    gdevice_batch_cffty_tmpy(true, ny, nn, n2ft3d, tmp2, d3db::tmpy);

    indx0 = 0;
    indx2 = 0;
    nn = 0;
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        if (!zero_row2[nb][indx2]) {
          jj = 0;
          indx = 2 * i + indx0;
          shift = 2 * ny * nn;
          for (j = 0; j < ny; ++j) {
            a[indx] = tmp2[jj + shift];
            a[indx + 1] = tmp2[jj + 1 + shift];
            jj += 2;
            indx += nxh2;
          }
          ++nn;
        }
        ++indx2;
      }
      indx0 += nxhy2;
    }

    /*
    indx0=0;
    indx2=0;
    for (q=0; q<nq; ++q)
    {
       for (i=0; i<nxh; ++i)
       {
          if (!zero_row2[nb][indx2])
          {
             jj   = 0;
             indx = 2*i+indx0;
             for (j=0; j<ny; ++j)
             {
                tmp2[jj]   = a[indx];
                tmp2[jj+1] = a[indx+1];
                jj   += 2;
                indx += nxh2;
             }

             dcfftf_(&ny,tmp2,tmpy);

             jj   = 0;
             indx = 2*i+indx0;
             for (j=0; j<ny; ++j)
             {
                a[indx]   = tmp2[jj];
                a[indx+1] = tmp2[jj+1];
                jj   += 2;
                indx += nxh2;
             }
          }
          ++indx2;
       }
       indx0 += nxhy2;
    }
    */

    /********************************************
     ***         Do a transpose of A          ***
     ***      A(ky,nz,ky) <- A(kx,ky,nz)      ***
     ********************************************/
    d3db::c_ptranspose2_jk(nb, a, tmp2, tmp3);
    // d3db::c_transpose_jk(a,tmp2,tmp3);

    /********************************************
     ***     do fft along nz dimension        ***
     ***   A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]  ***
     ********************************************/
    indx0 = 0;
    indx2 = 0;
    nn = 0;
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        if (!zero_row3[nb][indx2]) {
          kk = 0;
          indx = 2 * i + indx0;
          shift = 2 * nz * nn;
          for (k = 0; k < nz; ++k) {
            tmp2[kk + shift] = a[indx];
            tmp2[kk + 1 + shift] = a[indx + 1];
            kk += 2;
            indx += nxh2;
          }
          ++nn;
        }
        ++indx2;
      }
      indx0 += nxhz2;
    }

    gdevice_batch_cfftz_tmpz(true, nz, nn, n2ft3d, tmp2, d3db::tmpz);

    indx0 = 0;
    indx2 = 0;
    nn = 0;
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        if (!zero_row3[nb][indx2]) {
          kk = 0;
          indx = 2 * i + indx0;
          shift = 2 * nz * nn;
          for (k = 0; k < nz; ++k) {
            a[indx] = tmp2[kk + shift];
            a[indx + 1] = tmp2[kk + 1 + shift];
            kk += 2;
            indx += nxh2;
          }
          ++nn;
        }
        ++indx2;
      }
      indx0 += nxhz2;
    }

    /*
    indx0=0;
    indx2=0;
    for (q=0; q<nq; ++q)
    {
       for (i=0; i<nxh; ++i)
       {
          if (!zero_row3[nb][indx2])
          {
             kk   = 0;
             indx = 2*i+indx0;
             for (k=0; k<nz; ++k)
             {
                tmp2[kk]   = a[indx];
                tmp2[kk+1] = a[indx+1];
                kk   += 2;
                indx += nxh2;
             }

             dcfftf_(&nz,tmp2,tmpz);

             kk   = 0;
             indx = 2*i+indx0;
             for (k=0; k<nz; ++k)
             {
                a[indx]   = tmp2[kk];
                a[indx+1] = tmp2[kk+1];
                kk   += 2;
                indx += nxh2;
             }
          }
          ++indx2;
       }
       indx0 += nxhz2;
    }
    */

  }

  /*************************
   **** hilbert mapping ****
   *************************/
  else {
    /********************************************
     ***     do fft along nx dimension        ***
     ***   A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]  ***
     ********************************************/
    gdevice_batch_cfftx_tmpx(true, nx, nq1, n2ft3d, a, d3db::tmpx);
    /*
#if (defined NWPW_SYCL) || (defined NWPW_CUDA)
    gdevice_batch_cfftx(true,nx,nq1,n2ft3d,a);
#else
    indx = 0;
    for (q=0; q<nq1; ++q)
    {
       drfftf_(&nx,a+indx,d3db::tmpx);
       indx += nxh2;
    }
    d3db::cshift_fftf(nx,nq1,1,1,a);
#endif
    */
    d3db::c_ptranspose_ijk(nb, 0, a, tmp2, tmp3);

    /********************************************
     ***     do fft along ny dimension        ***
     ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
     ********************************************/
    gdevice_batch_cffty_tmpy_zero(true, ny, nq2, n2ft3d, a, d3db::tmpy,
                                  zero_row2[nb]);
    /*
#if (defined NWPW_SYCL) || (defined NWPW_CUDA)
    gdevice_batch_cffty(true,ny,nq2,n2ft3d,a);
#else
    indx = 0;
    for (q=0; q<nq2; ++q)
    {
       if (!zero_row2[nb][q])
          dcfftf_(&ny,a+indx,d3db::tmpy);
       indx += (2*ny);
    }
#endif
    */
    d3db::c_ptranspose_ijk(nb, 1, a, tmp2, tmp3);

    /********************************************
     ***     do fft along nz dimension        ***
     ***   A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]  ***
     ********************************************/
    gdevice_batch_cfftz_tmpz_zero(true, nz, nq3, n2ft3d, a, d3db::tmpz,
                                  zero_row3[nb]);
    /*
#if (defined NWPW_SYCL) || (defined NWPW_CUDA)
    gdevice_batch_cfftz(true,nz,nq3,n2ft3d,a);
#else
    indx = 0;
    for (q=0; q<nq3; ++q)
    {
       if (!zero_row3[nb][q])
          dcfftf_(&nz,a+indx,d3db::tmpz);
       indx += 2*nz;
    }
#endif
    */
  }

  delete[] tmp3;
  delete[] tmp2;
}

/********************************
 *                              *
 *     PGrid::c_unpack_start    *
 *                              *
 ********************************/
void PGrid::c_unpack_start(const int nb, double *tmp1, double *tmp2,
                           const int request_indx, const int msgtype) {
  if (balanced)
    mybalance->c_unbalance_start(nb, tmp1, request_indx, msgtype);
}

/********************************
 *                              *
 *     PGrid::c_unpack_mid      *
 *                              *
 ********************************/
void PGrid::c_unpack_mid(const int nb, double *tmp1, double *tmp2,
                         const int request_indx, const int msgtype) {
  if (balanced)
    mybalance->c_unbalance_end(nb, tmp1, request_indx);

  std::memcpy(tmp2, tmp1, 2 * (nida[nb] + nidb2[nb]) * sizeof(double));
  std::memset(tmp1, 0, n2ft3d * sizeof(double));

  c_bindexcopy((nida[nb] + nidb2[nb]), packarray[nb], tmp2, tmp1);
  // c_bindexcopy(nida[nb]+nidb[nb],packarray[nb],tmp2,tmp1);

  d3db::c_timereverse_start(tmp1, zplane_tmp1, zplane_tmp2, request_indx,
                            msgtype);
}

/********************************
 *                              *
 *     PGrid::c_unpack_end      *
 *                              *
 ********************************/
void PGrid::c_unpack_end(const int nb, double *tmp1, double *tmp2,
                         const int request_indx) {
  d3db::c_timereverse_end(tmp1, zplane_tmp1, zplane_tmp2, request_indx);
}

/********************************
 *                              *
 *        PGrid::pfftbz         *
 *                              *
 ********************************/
void PGrid::pfftbz(const int nb, double *tmp1, double *tmp2, int request_indx) {

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) {
    int kk, indx, shift;
    int nxh = nx / 2 + 1;
    int nxh2 = nx + 2;
    int nxhz2 = nxh2 * nz;

    /***************************************************
     ***     do fft along kz dimension               ***
     ***     A(kx,nz,ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
     ***************************************************/
    int indx0 = 0;
    int indx2 = 0;
    int nn = 0;
    for (auto q = 0; q < nq; ++q) {
      for (auto i = 0; i < nxh; ++i) {
        if (!zero_row3[nb][indx2]) {
          kk = 0;
          indx = 2 * i + indx0;
          shift = 2 * nz * nn;
          for (auto k = 0; k < nz; ++k) {
            tmp2[kk + shift] = tmp1[indx];
            tmp2[kk + 1 + shift] = tmp1[indx + 1];
            kk += 2;
            indx += nxh2;
          }
          nn += 1;
        }
        ++indx2;
      }
      indx0 += nxhz2;
    }

    gdevice_batch_cfftz_tmpz(false, nz, nn, n2ft3d, tmp2, d3db::tmpz);
    // for (auto i=0; i<nn; ++i)
    //    dcfftb_(&nz,tmp2+2*nz*i,d3db::tmpz);

    indx0 = 0;
    indx2 = 0;
    nn = 0;
    for (auto q = 0; q < nq; ++q) {
      for (auto i = 0; i < nxh; ++i) {
        if (!zero_row3[nb][indx2]) {
          kk = 0;
          indx = 2 * i + indx0;
          shift = 2 * nz * nn;
          for (auto k = 0; k < nz; ++k) {
            tmp1[indx] = tmp2[kk + shift];
            tmp1[indx + 1] = tmp2[kk + 1 + shift];
            kk += 2;
            indx += nxh2;
          }
          nn += 1;
        }
        ++indx2;
      }
      indx0 += nxhz2;
    }

    /*
     int indx0=0;
     int indx2=0;
     for (auto q=0; q<nq; ++q)
     {
        for (auto i=0; i<nxh; ++i)
        {
           if (!zero_row3[nb][indx2])
           {
              kk   = 0;
              indx = 2*i+indx0;
              for (auto k=0; k<nz; ++k)
              {
                 tmp2[kk]   = tmp1[indx];
                 tmp2[kk+1] = tmp1[indx+1];
                 kk   += 2;
                 indx += nxh2;
              }

              dcfftb_(&nz,tmp2,d3db::tmpz);

              kk   = 0;
              indx = 2*i+indx0;
              for (auto k=0; k<nz; ++k)
              {
                 tmp1[indx]   = tmp2[kk];
                 tmp1[indx+1] = tmp2[kk+1];
                 kk   += 2;
                 indx += nxh2;
              }
           }
           ++indx2;
        }
        indx0 += nxhz2;
     }
     */

    /***********************************************
     ***         Do a ptranspose of A            ***
     ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
     ************************************************/
    d3db::c_ptranspose1_jk_start(nb, tmp1, tmp2, tmp1, request_indx, 44);
  }
  /*************************
   **** hilbert mapping ****
   *************************/
  else {

    /************************************************
     ***     do fft along kz dimension            ***
     ***   A(nz,kx,ky) <- fft1d^(-1)[A(kz,kx,ky)] ***
     ************************************************/
    gdevice_batch_cfftz_tmpz_zero(false, nz, nq3, n2ft3d, tmp1, d3db::tmpz,
                                  zero_row3[nb]);

    d3db::c_ptranspose_ijk_start(nb, 2, tmp1, tmp2, tmp1, request_indx, 45);
    // d3db::c_ptranspose_ijk(nb,2,tmp1,tmp2,tmp1);
  }
}

/********************************
 *                              *
 *        PGrid::pfftby         *
 *                              *
 ********************************/
void PGrid::pfftby(const int nb, double *tmp1, double *tmp2, int request_indx) {

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) {
    int jj, indx, shift;
    int nxh = nx / 2 + 1;
    int nxh2 = nx + 2;
    int nxhy2 = nxh2 * ny;

    /***********************************************
     ***         Do a ptranspose of A            ***
     ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
     ************************************************/
    d3db::c_ptranspose1_jk_end(nb, tmp2, tmp1, request_indx);

    /********************************************
     ***     do fft along ny dimension        ***
     ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
     ********************************************/
    int indx0 = 0;
    int indx2 = 0;
    int nn = 0;
    for (auto q = 0; q < nq; ++q) {
      for (auto i = 0; i < nxh; ++i) {
        if (!zero_row2[nb][indx2]) {
          jj = 0;
          indx = 2 * i + indx0;
          shift = 2 * ny * nn;
          for (auto j = 0; j < ny; ++j) {
            tmp1[jj + shift] = tmp2[indx];
            tmp1[jj + 1 + shift] = tmp2[indx + 1];
            jj += 2;
            indx += nxh2;
          }
          nn += 1;
        }
        ++indx2;
      }
      indx0 += nxhy2;
    }

    gdevice_batch_cffty_tmpy(false, ny, nn, n2ft3d, tmp1, d3db::tmpy);
    // for (auto i=0; i<nn; ++i)
    //    dcfftb_(&ny,tmp1+2*ny*i,d3db::tmpy);

    indx0 = 0;
    indx2 = 0;
    nn = 0;
    for (auto q = 0; q < nq; ++q) {
      for (auto i = 0; i < nxh; ++i) {
        if (!zero_row2[nb][indx2]) {
          jj = 0;
          indx = 2 * i + indx0;
          shift = 2 * ny * nn;
          for (auto j = 0; j < ny; ++j) {
            tmp2[indx] = tmp1[jj + shift];
            tmp2[indx + 1] = tmp1[jj + 1 + shift];
            jj += 2;
            indx += nxh2;
          }
          nn += 1;
        }
        ++indx2;
      }
      indx0 += nxhy2;
    }

    /*
    int indx0=0;
    int indx2=0;
    for (auto q=0; q<nq; ++q)
    {
        for (auto i=0; i<nxh; ++i)
        {
           if (!zero_row2[nb][indx2])
           {
              jj   = 0;
              indx = 2*i+indx0;
              for (auto j=0; j<ny; ++j)
              {
                 tmp1[jj]   = tmp2[indx];
                 tmp1[jj+1] = tmp2[indx+1];
                 jj   += 2;
                 indx += nxh2;
              }

              dcfftb_(&ny,tmp1,tmpy);

              jj   = 0;
              indx = 2*i+indx0;
              for (auto j=0; j<ny; ++j)
              {
                 tmp2[indx]   = tmp1[jj];
                 tmp2[indx+1] = tmp1[jj+1];
                 jj   += 2;
                 indx += nxh2;
              }
           }
           ++indx2;
        }
        indx0 += nxhy2;
    }
    */
  }
  /*************************
   **** hilbert mapping ****
   *************************/
  else {
    d3db::c_ptranspose_ijk_end(nb, 2, tmp2, tmp1, request_indx);

    /********************************************
     ***     do fft along ny dimension        ***
     ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
     ********************************************/
    gdevice_batch_cffty_tmpy_zero(false, ny, nq2, n2ft3d, tmp2, d3db::tmpy,
                                  zero_row2[nb]);

    d3db::c_ptranspose_ijk_start(nb, 3, tmp2, tmp1, tmp2, request_indx, 46);
    // d3db::c_ptranspose_ijk(nb,3,tmp2,tmp1,tmp2);
  }
}

/********************************
 *                              *
 *        PGrid::pfftbx         *
 *                              *
 ********************************/
void PGrid::pfftbx(const int nb, double *tmp1, double *tmp2, int request_indx) {
  int i, j, k, jj, kk, q, indx, indx0, indx2, nxh, nxh2, nxhy, nxhy2, nxhz,
      nxhz2;

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) {
    /************************************************
     ***     do fft along kx dimension            ***
     ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
     ************************************************/
    gdevice_batch_cfftx_tmpx(false, nx, ny * nq, n2ft3d, tmp2, d3db::tmpx);
    d3db::zeroend_fftb(nx, ny, nq, 1, tmp2);
    std::memcpy(tmp1, tmp2, n2ft3d * sizeof(double));
  }
  /*************************
   **** hilbert mapping ****
   *************************/
  else {
    d3db::c_ptranspose_ijk_end(nb, 3, tmp1, tmp2, request_indx);

    /************************************************
     ***     do fft along kx dimension            ***
     ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
     ************************************************/
    gdevice_batch_cfftx_tmpx(false, nx, nq1, n2ft3d, tmp1, d3db::tmpx);
    d3db::zeroend_fftb(nx, nq1, 1, 1, tmp1);
    if (n2ft3d_map < n2ft3d)
      std::memset(tmp1 + n2ft3d_map, 0, (n2ft3d - n2ft3d_map) * sizeof(double));
  }
}

/********************************
 *                              *
 *       PGrid::pfftb_step      *
 *                              *
 ********************************/
void PGrid::pfftb_step(const int step, const int nb, double *a, double *tmp1,
                       double *tmp2, const int request_indx) {
  if (step == 0) {
    // parall->astart(request_indx,parall->np_i());

    // unpack start, tmp1-->tmp1
    std::memcpy(tmp1, a, 2 * (nida[nb] + nidb[nb]) * sizeof(double));
    this->c_unpack_start(nb, tmp1, tmp2, request_indx, 47);
  } else if (step == 1) {
    // unpack mid
    this->c_unpack_mid(nb, tmp1, tmp2, request_indx, 48);
  } else if (step == 2) {
    // unpack end; mem-->dev,  out=tmp1
    this->c_unpack_end(nb, tmp1, tmp2, request_indx);
  } else if (step == 3) {
    // pfftbz dev-->dev->mem,  tmp1->tmp1
    this->pfftbz(nb, tmp1, tmp2, request_indx);
  } else if (step == 4) {
    // pfftby mem->dev-->dev->mem
    // in=tmp1, tmp2->tmp1, tmp1=in , tmp2=tmp
    pfftby(nb, tmp1, tmp2, request_indx);
  } else if (step == 5) {
    // pfftbx mem->dev->dev->mem
    pfftbx(nb, tmp1, tmp2, request_indx);
    // parall->aend(request_indx);
  }
}

/********************************
 *                              *
 *    PGrid:cr_pfft3b_queuein   *
 *                              *
 ********************************/
void PGrid::cr_pfft3b_queuein(const int nb, double *a) {
  int shift1, shift2;
  int np = parall->np_i();

  for (auto q = 0; q < aqsize; ++q) {
    int indx = aqindx[q];
    int status = aqstatus[indx] + 1;
    shift1 = n2ft3d * (2 * indx);
    shift2 = n2ft3d * (2 * indx + 1);
    pfftb_step(status, nb, a, atmp + shift1, atmp + shift2, indx + 3);
    ++aqstatus[indx];
  }

  ++alast_index;
  if (alast_index >= aqmax)
    alast_index = 0;
  ++aqsize;
  aqindx[aqsize - 1] = alast_index;
  aqstatus[alast_index] = 0;

  // status = 0;
  shift1 = n2ft3d * (2 * alast_index);
  shift2 = n2ft3d * (2 * alast_index + 1);

  pfftb_step(0, nb, a, atmp + shift1, atmp + shift2, alast_index + 3);
}

/********************************
 *                              *
 *    PGrid:cr_pfft3b_queueout  *
 *                              *
 ********************************/
void PGrid::cr_pfft3b_queueout(const int nb, double *a) {
  int shift1, shift2;
  int indx1 = aqindx[0];

  while (aqstatus[indx1] < 5) {

    for (auto q = 0; q < aqsize; ++q) {
      int indx = aqindx[q];
      int status = aqstatus[indx] + 1;
      shift1 = n2ft3d * (2 * indx);
      shift2 = n2ft3d * (2 * indx + 1);
      pfftb_step(status, nb, a, atmp + shift1, atmp + shift2, indx + 3);
      ++aqstatus[indx];
    }
  }
  double scal1 = 1.0 / ((double)((nx) * (ny) * (nz)));
  double enrr0 = scal1 * d3db::rr_dot(atmp, atmp);

  shift1 = n2ft3d * (2 * indx1);
  std::memcpy(a, atmp + shift1, n2ft3d * sizeof(double));
  --aqsize;
  for (auto q = 0; q < aqsize; ++q)
    aqindx[q] = aqindx[q + 1];
}

/********************************
 *                              *
 * PGrid:cr_pfft3b_queuefilled  *
 *                              *
 ********************************/
int PGrid::cr_pfft3b_queuefilled() { return (aqsize >= aqmax); }

/********************************
 *                              *
 *        PGrid::pfftfx         *
 *                              *
 ********************************/
void PGrid::pfftfx(const int nb, double *a, double *tmp1, double *tmp2,
                   int request_indx) {

  /**** slab mapping ****/
  if (maptype == 1) {
    int nxh2 = nx + 2;

    // do fft along nx dimension
    // A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]
    /*
    int indx = 0;
    for (auto q=0; q<nq; ++q)
    for (auto j=0; j<ny; ++j)
    {
       drfftf_(&nx,tmp1+indx,tmpx);
       indx += nxh2;
    }
    d3db::cshift_fftf(nx,ny,nq,1,tmp1);
    */
    gdevice_batch_cfftx_tmpx(true, nx, ny * nq, n2ft3d, a, d3db::tmpx);
    std::memcpy(tmp1, a, n2ft3d * sizeof(double));
  }

  /**** hilbert mapping ****/
  else {
    // do fft along nx dimension
    // A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]
    gdevice_batch_cfftx_tmpx(true, nx, nq1, n2ft3d, a, d3db::tmpx);

    d3db::c_ptranspose_ijk_start(nb, 0, a, tmp1, tmp2, request_indx, 40);
  }
}

/********************************
 *                              *
 *        PGrid::pfftfy         *
 *                              *
 ********************************/
void PGrid::pfftfy(const int nb, double *tmp1, double *tmp2, int request_indx) {
  int i, j, k, jj, kk, q, indx, indx0, indx2, nxh, nxh2, nxhy, nxhy2, nxhz,
      nxhz2;

  /**** slab mapping ****/
  if (maptype == 1) {
    int jj, indx, nn, shift;
    int nxh = nx / 2 + 1;
    int nxh2 = nx + 2;
    int nxhy2 = nxh2 * ny;

    // do fft along ny dimension
    // A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]
    int indx0 = 0;
    int indx2 = 0;
    nn = 0;
    for (auto q = 0; q < nq; ++q) {
      for (auto i = 0; i < nxh; ++i) {
        if (!zero_row2[nb][indx2]) {
          jj = 0;
          indx = 2 * i + indx0;
          shift = 2 * ny * nn;
          for (auto j = 0; j < ny; ++j) {
            tmp2[jj + shift] = tmp1[indx];
            tmp2[jj + 1 + shift] = tmp1[indx + 1];
            jj += 2;
            indx += nxh2;
          }
          ++nn;
        }
        ++indx2;
      }
      indx0 += nxhy2;
    }

    gdevice_batch_cffty_tmpy(true, ny, nn, n2ft3d, tmp2, d3db::tmpy);

    indx0 = 0;
    indx2 = 0;
    nn = 0;
    for (auto q = 0; q < nq; ++q) {
      for (auto i = 0; i < nxh; ++i) {
        if (!zero_row2[nb][indx2]) {
          jj = 0;
          indx = 2 * i + indx0;
          shift = 2 * ny * nn;
          for (auto j = 0; j < ny; ++j) {
            tmp1[indx] = tmp2[jj + shift];
            tmp1[indx + 1] = tmp2[jj + 1 + shift];
            jj += 2;
            indx += nxh2;
          }
          ++nn;
        }
        ++indx2;
      }
      indx0 += nxhy2;
    }

    /*
    int indx0=0;
    int indx2=0;
    for (auto q=0; q<nq; ++q)
    {
       for (auto i=0; i<nxh; ++i)
       {
          if (!zero_row2[nb][indx2])
          {
             jj   = 0;
             indx = 2*i+indx0;
             for (auto j=0; j<ny; ++j)
             {
                tmp2[jj]   = tmp1[indx];
                tmp2[jj+1] = tmp1[indx+1];
                jj   += 2;
                indx += nxh2;
             }

             dcfftf_(&ny,tmp2,tmpy);

             jj   = 0;
             indx = 2*i+indx0;
             for (j=0; j<ny; ++j)
             {
                tmp1[indx]   = tmp2[jj];
                tmp1[indx+1] = tmp2[jj+1];
                jj   += 2;
                indx += nxh2;
             }
          }
          ++indx2;
       }
       indx0 += nxhy2;
     }
     */

    // Do a transpose of A
    // A(ky,nz,ky) <- A(kx,ky,nz)
    d3db::c_ptranspose2_jk_start(nb, tmp1, tmp2, tmp1, request_indx, 41);

  }

  /**** hilbert mapping ****/
  else {
    // in=tmp1, out=tmp2
    d3db::c_ptranspose_ijk_end(nb, 0, tmp1, tmp2, request_indx);

    // do fft along ny dimension
    // A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]
    gdevice_batch_cffty_tmpy_zero(true, ny, nq2, n2ft3d, tmp1, d3db::tmpy,
                                  zero_row2[nb]);

    // in=tmp2, out=tmp2
    d3db::c_ptranspose_ijk_start(nb, 1, tmp1, tmp2, tmp1, request_indx, 42);
  }
}

/********************************
 *                              *
 *        PGrid::pfftfz         *
 *                              *
 ********************************/
void PGrid::pfftfz(const int nb, double *tmp1, double *tmp2, int request_indx) {
  /**** slab mapping ****/
  if (maptype == 1) {

    d3db::c_ptranspose2_jk_end(nb, tmp2, tmp1, request_indx);

    int kk, indx, nn, shift;
    int nxh = nx / 2 + 1;
    int nxh2 = nx + 2;
    int nxhz2 = nxh2 * nz;

    // do fft along nz dimension
    // A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]
    int indx0 = 0;
    int indx2 = 0;
    nn = 0;
    for (auto q = 0; q < nq; ++q) {
      for (auto i = 0; i < nxh; ++i) {
        if (!zero_row3[nb][indx2]) {
          kk = 0;
          indx = 2 * i + indx0;
          shift = 2 * nz * nn;
          for (auto k = 0; k < nz; ++k) {
            tmp1[kk + shift] = tmp2[indx];
            tmp1[kk + 1 + shift] = tmp2[indx + 1];
            kk += 2;
            indx += nxh2;
          }
          ++nn;
        }
        ++indx2;
      }
      indx0 += nxhz2;
    }

    gdevice_batch_cfftz_tmpz(true, nz, nn, n2ft3d, tmp1, d3db::tmpz);

    indx0 = 0;
    indx2 = 0;
    nn = 0;
    for (auto q = 0; q < nq; ++q) {
      for (auto i = 0; i < nxh; ++i) {
        if (!zero_row3[nb][indx2]) {
          kk = 0;
          indx = 2 * i + indx0;
          shift = 2 * nz * nn;
          for (auto k = 0; k < nz; ++k) {
            tmp2[indx] = tmp1[kk + shift];
            tmp2[indx + 1] = tmp1[kk + 1 + shift];
            kk += 2;
            indx += nxh2;
          }
          ++nn;
        }
        ++indx2;
      }
      indx0 += nxhz2;
    }

    /*
    int indx0=0;
    int indx2=0;
    for (auto q=0; q<nq; ++q)
    {
       for (auto i=0; i<nxh; ++i)
       {
          if (!zero_row3[nb][indx2])
          {
             kk   = 0;
             indx = 2*i+indx0;
             for (auto k=0; k<nz; ++k)
             {
                tmp2[kk]   = tmp1[indx];
                tmp2[kk+1] = tmp1[indx+1];
                kk   += 2;
                indx += nxh2;
             }

             dcfftf_(&nz,tmp2,tmpz);

             kk   = 0;
             indx = 2*i+indx0;
             for (auto k=0; k<nz; ++k)
             {
                tmp1[indx]   = tmp2[kk];
                tmp1[indx+1] = tmp2[kk+1];
                kk   += 2;
                indx += nxh2;
             }
          }
          ++indx2;
       }
       indx0 += nxhz2;
    }
    */

  }

  /**** hilbert mapping ****/
  else {
    // in=tmp1, out=tmp2
    d3db::c_ptranspose_ijk_end(nb, 1, tmp2, tmp1, request_indx);

    // do fft along nz dimension
    // A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]
    gdevice_batch_cfftz_tmpz_zero(true, nz, nq3, n2ft3d, tmp2, d3db::tmpz,
                                  zero_row3[nb]);
  }
}

/********************************
 *                              *
 *       PGrid::pfftf_step      *
 *                              *
 ********************************/
void PGrid::pfftf_step(const int step, const int nb, double *a, double *tmp1,
                       double *tmp2, int request_indx) {
  if (step == 0) {
    // pfftfx mem-->device, in=a out=tmp2
    pfftfx(nb, a, tmp1, tmp2, request_indx);
  } else if (step == 1) {
    // pfftfy device, in=tmp1
    pfftfy(nb, tmp1, tmp2, request_indx);
  } else if (step == 2) {
    // pfftfz device-->mem
    pfftfz(nb, tmp1, tmp2, request_indx);
    this->c_pack_start(nb, tmp2, tmp1, request_indx, 47);
  } else if (step == 3) {
    // pfftf final
    this->c_pack_end(nb, tmp2, request_indx);
  }
}

/********************************
 *                              *
 *       PGrid:c_pack_start     *
 *                              *
 ********************************/
void PGrid::c_pack_start(const int nb, double *a, double *tmp1,
                         const int request_indx, const int msgtype) {
  // int one=1;

  // DCOPY_PWDFT(n2ft3d,a,one,tmp,one);
  std::memcpy(tmp1, a, n2ft3d * sizeof(double));
  std::memset(a, 0, n2ft3d * sizeof(double));

  c_aindexcopy(nida[nb] + nidb2[nb], packarray[nb], tmp1, a);

  if (balanced)
    mybalance->c_balance_start(nb, a, request_indx, msgtype);

  return;
}

/********************************
 *                              *
 *       PGrid:c_pack_end       *
 *                              *
 ********************************/
void PGrid::c_pack_end(const int nb, double *tmp1, const int request_indx) {

  if (balanced)
    mybalance->c_balance_end(nb, tmp1, request_indx);

  return;
}

/********************************
 *                              *
 *    PGrid:rc_pfft3f_queuein   *
 *                              *
 ********************************/
void PGrid::rc_pfft3f_queuein(const int nb, double *b) {
  int shift1, shift2;
  int np = parall->np_i();

  for (auto q = 0; q < bqsize; ++q) {
    int indx = bqindx[q];
    int status = bqstatus[indx] + 1;
    shift1 = n2ft3d * (2 * indx);
    shift2 = n2ft3d * (2 * indx + 1);
    pfftf_step(status, nb, b, btmp + shift1, btmp + shift2, indx + 3);
    ++bqstatus[indx];
  }

  ++blast_index;
  if (blast_index >= bqmax)
    blast_index = 0;
  ++bqsize;
  bqindx[bqsize - 1] = blast_index;
  bqstatus[blast_index] = 0;

  // status = 0;
  shift1 = n2ft3d * (2 * blast_index);
  shift2 = n2ft3d * (2 * blast_index + 1);

  pfftf_step(0, nb, b, btmp + shift1, btmp + shift2, blast_index + 3);
}

/********************************
 *                              *
 *    PGrid:rc_pfft3f_queueout  *
 *                              *
 ********************************/
void PGrid::rc_pfft3f_queueout(const int nb, double *b) {
  int shift1, shift2;
  int indx1 = bqindx[0];

  while (bqstatus[indx1] < 5) {

    for (auto q = 0; q < bqsize; ++q) {
      int indx = bqindx[q];
      int status = bqstatus[indx] + 1;
      shift1 = n2ft3d * (2 * indx);
      shift2 = n2ft3d * (2 * indx + 1);
      pfftf_step(status, nb, b, btmp + shift1, btmp + shift2, indx + 3);
      ++bqstatus[indx];
    }
  }
  double scal1 = 1.0 / ((double)((nx) * (ny) * (nz)));
  double enrr0 = scal1 * d3db::rr_dot(btmp, btmp);

  shift2 = n2ft3d * (2 * indx1 + 1);
  std::memcpy(b, btmp + shift2, n2ft3d * sizeof(double));
  --bqsize;
  for (auto q = 0; q < bqsize; ++q)
    bqindx[q] = bqindx[q + 1];
}

/********************************
 *                              *
 * PGrid:rc_pfft3f_queuefilled  *
 *                              *
 ********************************/
int PGrid::rc_pfft3f_queuefilled() { return (bqsize >= bqmax); }

/********************************
 *                              *
 *     PGrid:tc_pack_copy       *
 *                              *
 ********************************/
void PGrid::tc_pack_copy(const int nb, double *a, double *b) {
  int i, ii;
  int ng = nida[nb] + nidb[nb];

  ii = 0;
  for (i = 0; i < ng; ++i) {
    b[ii] = a[i];
    b[ii + 1] = 0.0;
    ii += 2;
  }
}

/********************************
 *                              *
 *      PGrid:tcc_pack_Mul      *
 *                              *
 ********************************/
void PGrid::tcc_pack_Mul(const int nb, const double *a, const double *b, double *c)
{
   int i, ii;
   int ng = nida[nb]+nidb[nb];
 
   ii = 0;
   for (i=0; i<ng; ++i)
   {
      c[ii]   = b[ii]*  a[i];
      c[ii+1] = b[ii+1]*a[i];
      ii += 2;
   }
}

/********************************
 *                              *
 *      PGrid:tcc_pack_aMul     *
 *                              *
 ********************************/
void PGrid::tcc_pack_aMul(const int nb, const double alpha, const double *a, const double *b, double *c)
{
   int i, ii;
   int ng = nida[nb]+nidb[nb];

   ii = 0;
   for (i=0; i<ng; ++i)
   {
      c[ii]   = alpha*b[ii]*  a[i];
      c[ii+1] = alpha*b[ii+1]*a[i];
      ii += 2;
   }
}

/********************************
 *                              *
 *      PGrid:tc_pack_Mul       *
 *                              *
 ********************************/
void PGrid::tc_pack_Mul(const int nb, const double *a, double *c) {
  int i, ii;
  int ng = nida[nb] + nidb[nb];

  ii = 0;
  for (i = 0; i < ng; ++i) {
    c[ii] = c[ii] * a[i];
    c[ii + 1] = c[ii + 1] * a[i];
    ii += 2;
  }
}

/********************************
 *                              *
 *    PGrid:tcc_pack_aMulAdd    *
 *                              *
 ********************************/
void PGrid::tcc_pack_aMulAdd(const int nb, const double alpha, const double *a, const double *b, double *c) 
{
   int i, ii;
   int ng = nida[nb] + nidb[nb];
 
   ii=0;
   for (i=0; i<ng; ++i) 
   {
      c[ii]   += alpha*b[ii] * a[i];
      c[ii+1] += alpha*b[ii+1]*a[i];
      ii += 2;
   }
}

/********************************
 *                              *
 *      PGrid:tcc_pack_iMul     *
 *                              *
 ********************************/
void PGrid::tcc_pack_iMul(const int nb, const double *a, const double *b,
                          double *c) {
  int i, ii;
  int ng = nida[nb] + nidb[nb];

  ii = 0;
  for (i = 0; i < ng; ++i) {
    c[ii] = -b[ii + 1] * a[i];
    c[ii + 1] = b[ii] * a[i];
    ii += 2;
  }
}

/*******************************************
 *                                         *
 *     PGrid:tcr_pack_iMul_unpack_fft      *
 *                                         *
 *******************************************/
void PGrid::tcr_pack_iMul_unpack_fft(const int nb, const double *a, const double *b, double *c) 
{
   int i, ii;
   int ng = nida[nb] + nidb[nb];

   ii = 0;
   for (i=0; i<ng; ++i)
   {
      c[ii]   = -b[ii+1]* a[i];
      c[ii+1] = b[ii]   * a[i];
      ii += 2;
   }
   this->c_unpack(nb,c);
   this->cr_pfft3b(nb,c);
   this->d3db::r_zero_ends(c);
}

/********************************
 *                              *
 *     PGrid:tc_pack_iMul       *
 *                              *
 ********************************/
void PGrid::tc_pack_iMul(const int nb, const double *a, double *c) {
  int i, ii;
  int ng = nida[nb] + nidb[nb];
  double x, y;

  ii = 0;
  for (i = 0; i < ng; ++i) {
    x = c[ii];
    y = c[ii + 1];

    c[ii] = -y * a[i];
    c[ii + 1] = x * a[i];
    ii += 2;
  }
}

/********************************
 *                              *
 *    PGrid:ttc_pack_MulSum2    *
 *                              *
 ********************************/
void PGrid::tcc_pack_MulSum2(const int nb, const double *a, const double *b,
                             double *c) {
  int i, ii;
  int ng = nida[nb] + nidb[nb];

  ii = 0;
  for (i = 0; i < ng; ++i) {
    c[ii] += b[ii] * a[i];
    c[ii + 1] += b[ii + 1] * a[i];
    ii += 2;
  }
}

/********************************
 *                              *
 *      PGrid:cc_pack_Sum2      *
 *                              *
 ********************************/
void PGrid::cc_pack_Sum2(const int nb, const double *a, double *b) {
  int i;
  int ng = 2 * (nida[nb] + nidb[nb]);

  for (i = 0; i < ng; ++i)
    b[i] += a[i];
}

/********************************
 *                              *
 *     PGrid:cccc_pack_Sum      *
 *                              *
 ********************************/
void PGrid::cccc_pack_Sum(const int nb, const double *a, const double *b,
                          const double *c, double *d) {
  int i;
  int ng = 2 * (nida[nb] + nidb[nb]);

  for (i = 0; i < ng; ++i)
    d[i] = (a[i] + b[i] + c[i]);
}

/********************************
 *                              *
 *     PGrid:c_pack_addzero     *
 *                              *
 ********************************/
void PGrid::c_pack_addzero(const int nb, const double vzero, double *a) {
  int pzero = ijktop(0, 0, 0);
  if (pzero == parall->taskid_i())
    a[0] += vzero;
}

/********************************
 *                              *
 *   PGrid:c_pack_noimagzero    *
 *                              *
 ********************************/

void PGrid::c_pack_noimagzero(const int nb, double *a) {
  int pzero = ijktop(0, 0, 0);
  if (pzero == parall->taskid_i())
    a[1] = 0.0;
}

/********************************
 *                              *
 *     PGrid:c_pack_zero        *
 *                              *
 ********************************/
void PGrid::c_pack_zero(const int nb, double *b) {
  int i;
  int ng = 2 * (nida[nb] + nidb[nb]);

  for (i = 0; i < ng; ++i)
    b[i] = 0.0;
}

/********************************
 *                              *
 *       PGrid:c_pack_SMul      *
 *                              *
 ********************************/
void PGrid::c_pack_SMul(const int nb, const double alpha, double *b) {
  int i;
  int ng = 2 * (nida[nb] + nidb[nb]);

  for (i = 0; i < ng; ++i)
    b[i] *= alpha;
}

/********************************
 *                              *
 *     PGrid:cc_pack_SMul       *
 *                              *
 ********************************/
void PGrid::cc_pack_SMul(const int nb, const double alpha, const double *a,
                         double *b) {
  int i;
  int ng = 2 * (nida[nb] + nidb[nb]);

  for (i = 0; i < ng; ++i)
    b[i] = alpha * a[i];
}

/********************************
 *                              *
 *      PGrid:cc_pack_daxpy     *
 *                              *
 ********************************/
void PGrid::cc_pack_daxpy(const int nb, const double alpha, const double *a,
                          double *b) {
  int i;
  int ng = 2 * (nida[nb] + nidb[nb]);

  for (i = 0; i < ng; ++i)
    b[i] += alpha * a[i];
}

/********************************
 *                              *
 *   PGrid:cct_pack_iconjgMul   *
 *                              *
 ********************************/
void PGrid::cct_pack_iconjgMul(const int nb, const double *a, const double *b, double *c) 
{
   for (int i=0; i<(nida[nb]+nidb[nb]); ++i)
      c[i] = a[2*i]*b[2*i+1] - a[2*i+1]*b[2*i];
}

/********************************
 *                              *
 *  PGrid:cct_pack_iconjgMulb   *
 *                              *
 ********************************/
void PGrid::cct_pack_iconjgMulb(const int nb, const double *a, const double *b, double *c) 
{
   for (int i=0; i<(nida[nb]+nidb[nb]); ++i)
      c[i] = a[2*i+1]*b[2*i] - a[2*i]*b[2*i+1];
}

/**********************************
 *                                *
 *    PGrid::regenerate_r_grid    *
 *                                *
 **********************************/
void PGrid::regenerate_r_grid() {
  int nxh = nx / 2;
  int nyh = ny / 2;
  int nzh = nz / 2;
  double a[9];
  for (auto i = 0; i < 3; ++i) {
    a[i] = lattice->unita1d(0 + i) / ((double)nx);
    a[3 + i] = lattice->unita1d(3 + i) / ((double)ny);
    a[6 + i] = lattice->unita1d(6 + i) / ((double)nz);
  }

  r_nzero(3, r_grid);

  /* grid points in coordination space */
  for (auto k3 = (-nzh); k3 < nzh; ++k3)
    for (auto k2 = (-nyh); k2 < nyh; ++k2)
      for (auto k1 = (-nxh); k1 < nxh; ++k1) {
        int i = k1 + nxh;
        int j = k2 + nyh;
        int k = k3 + nzh;
        int indx = ijktoindex2(i, j, k);
        int p = ijktop2(i, j, k);

        if (p == parall->taskid_i()) {
          r_grid[3 * indx] = a[0] * k1 + a[3] * k2 + a[6] * k3;
          r_grid[3 * indx + 1] = a[1] * k1 + a[4] * k2 + a[7] * k3;
          r_grid[3 * indx + 2] = a[2] * k1 + a[5] * k2 + a[8] * k3;
        }
      }
}

/************************************
 *                                  *
 *    PGrid::generate_r_sym_grid    *
 *                                  *
 ************************************/
void PGrid::generate_r_sym_grid(double *r_sym_grid) {
  int nxh = nx / 2;
  int nyh = ny / 2;
  int nzh = nz / 2;
  double a[9];
  for (auto i = 0; i < 3; ++i) {
    a[i] = lattice->unita1d(0 + i) / ((double)nx);
    a[3 + i] = lattice->unita1d(3 + i) / ((double)ny);
    a[6 + i] = lattice->unita1d(6 + i) / ((double)nz);
  }

  r_nzero(3, r_sym_grid);

  /* grid points in coordination space */
  for (auto k3 = (-nzh + 1); k3 < nzh; ++k3)
    for (auto k2 = (-nyh + 1); k2 < nyh; ++k2)
      for (auto k1 = (-nxh + 1); k1 < nxh; ++k1) {
        int i = k1 + nxh;
        int j = k2 + nyh;
        int k = k3 + nzh;
        int indx = ijktoindex2(i, j, k);
        int p = ijktop2(i, j, k);

        if (p == parall->taskid_i()) {
          r_sym_grid[3 * indx] = a[0] * k1 + a[3] * k2 + a[6] * k3;
          r_sym_grid[3 * indx + 1] = a[1] * k1 + a[4] * k2 + a[7] * k3;
          r_sym_grid[3 * indx + 2] = a[2] * k1 + a[5] * k2 + a[8] * k3;
        }
      }
}

/************************************
 *                                  *
 *    PGrid::generate_r_sym_mask    *
 *                                  *
 ************************************/
void PGrid::generate_r_sym_mask(double *rmask) {
  int nxh = nx / 2;
  int nyh = ny / 2;
  int nzh = nz / 2;
  r_zero(rmask);

  /* grid points in coordination space */
  for (auto k3 = (-nzh); k3 < nzh; ++k3)
    for (auto k2 = (-nyh); k2 < nyh; ++k2)
      for (auto k1 = (-nxh); k1 < nxh; ++k1) {
        int i = k1 + nxh;
        int j = k2 + nyh;
        int k = k3 + nzh;
        int indx = ijktoindex2(i, j, k);
        int p = ijktop2(i, j, k);

        if (p == parall->taskid_i())
          rmask[indx] = 1.0;
      }
}

/************************************
 *                                  *
 *       PGrid::c_Laplacian         *
 *                                  *
 ************************************/
void PGrid::c_Laplacian(const int nb, double *w) {
  int npack0 = this->npack(nb);
  double *Gx = this->Gpackxyz(nb, 0);
  double *Gy = this->Gpackxyz(nb, 1);
  double *Gz = this->Gpackxyz(nb, 2);
  int kk = 0;
  for (auto k = 0; k < npack0; ++k) {
    auto gg = Gx[k] * Gx[k] + Gy[k] * Gy[k] + Gz[k] * Gz[k];
    w[kk] *= (-gg);
    w[kk + 1] *= (-gg);
    kk += 2;
  }
}

/************************************
 *                                  *
 *       PGrid::cc_Laplacian        *
 *                                  *
 ************************************/
void PGrid::cc_Laplacian(const int nb, const double *w0, double *w) {
  int npack0 = this->npack(nb);
  double *Gx = this->Gpackxyz(nb, 0);
  double *Gy = this->Gpackxyz(nb, 1);
  double *Gz = this->Gpackxyz(nb, 2);
  int kk = 0;
  for (auto k = 0; k < npack0; ++k) {
    auto gg = Gx[k] * Gx[k] + Gy[k] * Gy[k] + Gz[k] * Gz[k];
    w[kk] = -w0[kk] * gg;
    w[kk + 1] = -w0[kk + 1] * gg;
    kk += 2;
  }
}

/************************************
 *                                  *
 *      PGrid::rr_Laplacian         *
 *                                  *
 ************************************/
void PGrid::rr_Laplacian(const int nb, const double *w0, double *w) {
  double scal1 = 1.0 / ((double)((nx) * (ny) * (nz)));

  rr_SMul(scal1, w0, w);
  this->rc_pfft3f(nb, w);
  this->c_pack(nb, w);

  this->c_Laplacian(nb, w);

  this->c_unpack(nb, w);
  this->cr_pfft3b(nb, w);
}

/************************************
 *                                  *
 *       PGrid::rr_Helmholtz        *
 *                                  *
 ************************************/
void PGrid::rr_Helmholtz(const int nb, const double *k2, const double *w,
                         double *Hw) {
  this->rr_Laplacian(nb, w, Hw);
  this->rrr_Mul2Add(k2, w, Hw);
}

/************************************
 *                                  *
 *   PGrid::rrr_solve_Helmholtz     *
 *                                  *
 ************************************/
/* The routine solves the inhomogeneous Helmholtz equation

 ∇^2 w(r) + k2(r)w(r) = b(r)  using CG.

 Entry - nb: 0-density grid, 1-wvfnc grid
         k2(r): wavenumber function
         b(r): source term

 Entry/Exit - w0(r) : initial guess / output solution

*/
void PGrid::rrr_solve_Helmholtz(const int nb, const double *k2, const double *b,
                                double *w) {
  double alpha, alpha0;
  double delta = 1.0;
  double *Aw = r_alloc();
  double *R = r_alloc();
  double *HR = r_alloc();

  double omega = this->lattice->omega();
  double scal1 = 1.0 / ((double)((nx) * (ny) * (nz)));
  double dv = omega * scal1;

  this->r_zero(w);

  int it = 0;
  while ((it < 10000) && (delta > 0.01)) {
    // compute the Aw and the residual
    this->rr_Helmholtz(nb, k2, w, Aw);
    this->rrr_Minus(b, Aw, R);

    this->rr_Helmholtz(nb, k2, R, HR);
    delta = rr_dot(R, R) * dv;
    // alpha = -1.0e-6;
    alpha0 = delta / rr_dot(R, HR);
    alpha = 1.0e-5 * alpha0;

    std::cout << "it=" << it << " delta=" << delta << " alpha0=" << alpha0
              << " alpha=" << alpha << std::endl;

    rr_daxpy(alpha, R, w);
    ++it;
  }

  r_dealloc(HR);
  r_dealloc(R);
  r_dealloc(Aw);
}

/************************************
 *                                  *
 *   PGrid::rrrr_FD_gradient        *
 *                                  *
 ************************************/
void PGrid::rrrr_FD_gradient(const double *rho, double *rhox, double *rhoy,
                             double *rhoz) {
  double ua[9];
  for (auto i = 0; i < 3; ++i) {
    ua[i] = lattice->unita1d(0 + i) / ((double)nx);
    ua[3 + i] = lattice->unita1d(3 + i) / ((double)ny);
    ua[6 + i] = lattice->unita1d(6 + i) / ((double)nz);
  }
  double dx = std::sqrt(ua[0] * ua[0] + ua[1] * ua[1] + ua[2] * ua[2]);
  double dy = std::sqrt(ua[3] * ua[3] + ua[4] * ua[4] + ua[5] * ua[5]);
  double dz = std::sqrt(ua[6] * ua[6] + ua[7] * ua[7] + ua[8] * ua[8]);

  this->rrrr_periodic_gradient(rho, rhox, rhoy, rhoz);
  for (auto i = 0; i < n2ft3d; ++i) {
    rhox[i] /= dx;
    rhoy[i] /= dy;
    rhoz[i] /= dz;
  }
}

/************************************
 *                                  *
 *   PGrid::rrrr_FD_laplacian       *
 *                                  *
 ************************************/
void PGrid::rrrr_FD_laplacian(const double *rho, double *rhoxx, double *rhoyy,
                              double *rhozz) {
  double ua[9];
  for (auto i = 0; i < 3; ++i) {
    ua[i] = lattice->unita1d(0 + i) / ((double)nx);
    ua[3 + i] = lattice->unita1d(3 + i) / ((double)ny);
    ua[6 + i] = lattice->unita1d(6 + i) / ((double)nz);
  }
  double dxx = (ua[0] * ua[0] + ua[1] * ua[1] + ua[2] * ua[2]);
  double dyy = (ua[3] * ua[3] + ua[4] * ua[4] + ua[5] * ua[5]);
  double dzz = (ua[6] * ua[6] + ua[7] * ua[7] + ua[8] * ua[8]);

  this->rrrr_periodic_laplacian(rho, rhoxx, rhoyy, rhozz);
  for (auto i = 0; i < n2ft3d; ++i) {
    rhoxx[i] /= (dxx);
    rhoyy[i] /= (dyy);
    rhozz[i] /= (dzz);
  }
}

} // namespace pwdft

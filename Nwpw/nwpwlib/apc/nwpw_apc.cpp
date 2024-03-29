/* nwpw_apc.cpp -
   Author - Eric Bylaska
*/

#include <cmath>
#include <cstring>
#include <iostream>

//#include "gdevice.hpp"
#include "nwpw_timing.hpp"

#include "blas.h"

#include "nwpw_apc.hpp"

namespace pwdft {

/* Constructors */

/*******************************************
 *                                         *
 *            nwpw_apc::nwpw_apc           *
 *                                         *
 *******************************************/
nwpw_apc::nwpw_apc(Ion *myionin, Pneb *mypnebin, Strfac *mystrfacin,
                   Control2 &control, std::ostream &coutput) {

  myion = myionin;
  mypneb = mypnebin;
  mystrfac = mystrfacin;
  apc_on = control.APC_on();
  v_apc_on = false;
  born_on = false;

  bool oprint =
      ((mypneb->PGrid::parall->is_master()) && (control.print_level("medium")));

  if (apc_on) {
    Gc = control.APC_Gc();
    nga = control.APC_nga();
    if (nga <= 0) {
      nga = 3;
      gamma = new double[nga];
      gamma[0] = 0.6;
      gamma[1] = 0.9;
      gamma[2] = 1.35;
    } else {
      gamma = new double[nga];
      for (auto i = 0; i < nga; ++i)
        gamma[i] = control.APC_gamma(i);
    }
    ngs = nga * (myion->nion);

    /* allocate APC memory from heap  */
    A = new double[4 * ngs * ngs];
    Am = new double[ngs * ngs];
    b = new double[4 * ngs];
    q = new double[ngs];
    u = new double[ngs];
    // utmp = new double [ngs];

    qion = new double[myion->nion];
    uion = new double[myion->nion];

    int npack0 = mypneb->npack(0);
    w = new double[npack0];
    gaus = new double[nga * npack0];
    vtmp = new double[2 * npack0];

    Eapc = 0.0;
    Papc = 0.0;

    /* define weight function */
    double gg, xx;
    double fourpi = 16.0 * atan(1.0);
    double *Gx = mypneb->Gpackxyz(0, 0);
    double *Gy = mypneb->Gpackxyz(0, 1);
    double *Gz = mypneb->Gpackxyz(0, 2);
    for (auto i = 0; i < npack0; ++i) {
      gg = (Gx[i] * Gx[i] + Gy[i] * Gy[i] + Gz[i] * Gz[i]);
      w[i] = 0.0;
      if ((gg > 1.0e-6) && (gg < (Gc * Gc))) {
        xx = (gg - Gc * Gc);
        w[i] = fourpi * xx * xx / (gg * Gc * Gc);
      }
    }

    /* define Gaussians in G-space */
    double coef = 1.0 / mypneb->lattice->omega();
    for (auto n = 0; n < nga; ++n) {
      xx = gamma[n] * gamma[n] / 4.0;
      for (auto i = 0; i < npack0; ++i) {
        gg = (Gx[i] * Gx[i] + Gy[i] * Gy[i] + Gz[i] * Gz[i]);
        gaus[n * npack0 + i] = coef * exp(-xx * gg);
      }
    }

    /* turn on self-consistency */
    v_apc_on = control.born_relax();
    for (auto ii = 0; ii < myion->nion; ++ii) {
      qion[ii] = control.APC_q(ii);
      uion[ii] = control.APC_u(ii);
      if (std::abs(uion[ii]) > 1.0e-9)
        v_apc_on = true;
    }

    /* write out APC header */
    if (oprint) {
      coutput << std::endl;
      coutput << " initializing nwpw_APC object" << std::endl;
      coutput << " ----------------------------" << std::endl;
      coutput << " nga = " << std::setw(3) << nga << " ngs = " << std::setw(5)
              << ngs << std::endl
              << " Gc  = " << std::fixed << std::setprecision(5) << std::setw(9)
              << Gc << std::endl;
      for (auto i = 0; i < nga; ++i)
        coutput << " APC gamma: " << i << " " << gamma[i] << std::endl;
      for (auto ii = 0; ii < myion->nion; ++ii)
        if (std::abs(uion[ii]) > 1.0e-9)
          coutput << " APC u: " << std::setw(8) << ii + 1 << "    "
                  << std::setw(12) << std::fixed << std::setprecision(5)
                  << uion[ii] << std::endl;
      if (v_apc_on)
        coutput << " - self-consistent" << std::endl;
      else
        coutput << " - not self-consistent" << std::endl;
    }

    born_on = control.born_on();
    if (born_on)
      myborn = new nwpw_born(myion, mypneb->PGrid::parall, control, coutput);
  }
}

/*******************************************
 *                                         *
 *            nwpw_apc::gen_APC            *
 *                                         *
 *******************************************/
void nwpw_apc::gen_APC(double *dng, bool move) {
  if (apc_on) {
    int npack0 = mypneb->npack(0);
    int ispin = mypneb->ispin;

    int i, j, indx, indxt;
    double fourpi = 16.0 * atan(1.0);
    double *Gx = mypneb->Gpackxyz(0, 0);
    double *Gy = mypneb->Gpackxyz(0, 1);
    double *Gz = mypneb->Gpackxyz(0, 2);
    double omega = mypneb->lattice->omega();

    /* allocate temporary memory from heap */
    double *exi = new double[2 * npack0];
    double *exj = new double[2 * npack0];
    double *gaus_i = new double[2 * npack0];
    double *gaus_j = new double[2 * npack0];
    double *xtmp = new double[npack0];

    std::memset(exi, 0, 2 * npack0 * sizeof(double));

    /* calculate N = dng(G=0)*omega */
    double N = ((double)(mypneb->ne[0] + mypneb->ne[ispin - 1]));

    /* calculate c_i = omega*gaus_i(G=0) = 1.0d0 */

    /* calculate b_i = omega*Sum(G) w(G)*Re(dcongj(dng(G))*gaus_i(G)) */
    for (auto ii = 0; ii < myion->nion; ++ii) {
      mystrfac->strfac_pack(0, ii, exi);

      for (auto iii = 0; iii < nga; ++iii) {
        i = iii + ii * nga;

        /* gaus_i(G)) */
        /* w(G)*gaus_i(G)) */
        mypneb->tcc_pack_Mul(0, &gaus[npack0 * iii], exi, gaus_i);
        mypneb->tc_pack_Mul(0, w, gaus_i);

        /* bi = omega*Sum(G) w(G)*Re(dcongj(dng(G))*gaus_i(G)) */
        b[i] = omega * mypneb->cc_pack_dot(0, dng, gaus_i);
        // std::cout << i << " " << b[i] << std::endl;

        if (move) {
          mypneb->cct_pack_iconjgMulb(0, dng, gaus_i, xtmp);
          b[i + ngs] = omega * mypneb->tt_pack_dot(0, Gx, xtmp);
          b[i + 2 * ngs] = omega * mypneb->tt_pack_dot(0, Gy, xtmp);
          b[i + 3 * ngs] = omega * mypneb->tt_pack_dot(0, Gz, xtmp);
          // std::cout << "i,b,db= " << std::right << std::setw(4) << i << " "
          // << std::fixed << std::setw(15) << std::setprecision(11) <<  b[i]
          //                << " " << b[i+ngs]
          //                << " " << b[i+2*ngs]
          //                << " " << b[i+3*ngs] << std::endl;
        }
      }
    }

    /* calculate A_ij = omega*Sum(G) w(G)*dcongj(gaus_i(G))*gaus_j(G)) */
    for (auto ii = 0; ii < myion->nion; ++ii) {
      mystrfac->strfac_pack(0, ii, exi);

      for (auto jj = ii; jj < myion->nion; ++jj) {
        mystrfac->strfac_pack(0, jj, exj);

        for (auto iii = 0; iii < nga; ++iii) {
          /* gaus_i(G)) */
          /* w(G)*gaus_i(G)) */
          mypneb->tcc_pack_Mul(0, &gaus[npack0 * iii], exi, gaus_i);
          mypneb->tc_pack_Mul(0, w, gaus_i);

          for (auto jjj = 0; jjj < nga; ++jjj) {
            /* gaus_j(G)) */
            mypneb->tcc_pack_Mul(0, &gaus[npack0 * jjj], exj, gaus_j);

            /* omega*Sum(G) w(G)*gaus_i(G)*gaus_j(G)) */
            double e1 = omega * mypneb->cc_pack_dot(0, gaus_i, gaus_j);
            i = iii + ii * nga;
            j = jjj + jj * nga;

            indx = i + j * ngs;
            indxt = j + i * ngs;

            A[indx] = e1;
            A[indxt] = e1;

            if (move) {
              mypneb->cct_pack_iconjgMulb(0, gaus_i, gaus_j, xtmp);
              A[indx + ngs * ngs] = omega * mypneb->tt_pack_dot(0, Gx, xtmp);
              A[indx + 2 * ngs * ngs] =
                  omega * mypneb->tt_pack_dot(0, Gy, xtmp);
              A[indx + 3 * ngs * ngs] =
                  omega * mypneb->tt_pack_dot(0, Gz, xtmp);
              if (indx != indxt) {
                A[indxt + ngs * ngs] = -A[indx + ngs * ngs];
                A[indxt + 2 * ngs * ngs] = -A[indx + 2 * ngs * ngs];
                A[indxt + 3 * ngs * ngs] = -A[indx + 3 * ngs * ngs];
              }
              // std::cout << "i,j,dA=" << std::setw(4) << indx << " " <<
              // std::setw(4) << indxt
              //           << std::fixed << std::setw(15) <<
              //           std::setprecision(11)
              //           << A[indx] << " " << A[indx+ngs*ngs] << " "
              //           << A[indx+2*ngs*ngs] << " " << A[indx+3*ngs*ngs] <<
              //           std::endl;
            }
          }
        }
      }
    }

    /* perform matrix operations in serial */
    std::memset(Am, 0, ngs * ngs * sizeof(double));
    std::memset(q, 0, ngs * sizeof(double));

    if (mypneb->PGrid::parall->is_master()) {
      int ierr, rank;

      int lwork = 5 * ngs * ngs;
      double work[lwork];
      double AAA[ngs * ngs];
      double rcond = 1.0e-9;

      std::memcpy(AAA, A, ngs * ngs * sizeof(double));
      // for (i=0; i<ngs*ngs; ++i) AAA[i]        = A[i];
      for (i = 0; i < ngs; ++i)
        Am[i + i * ngs] = 1.0;

      DGELSS_PWDFT(ngs, ngs, ngs, AAA, ngs, Am, ngs, q, rcond, rank, work,
                   lwork, ierr);

      /* calculate q_i */
      double sum = 0.0;
      double sum1 = 0.0;
      for (i = 0; i < ngs; ++i)
        for (j = 0; j < ngs; ++j) {
          indx = i + j * ngs;
          sum += Am[indx] * b[j];
          sum1 += Am[indx];
        }
      sum = (sum - N) / sum1;

      std::memset(q, 0, ngs * sizeof(double));
      for (i = 0; i < ngs; ++i) {
        sum1 = 0.0;
        for (j = 0; j < ngs; ++j) {
          indx = i + j * ngs;
          sum1 += Am[indx] * (b[j] - sum);
        }
        q[i] = sum1;
      }
    }

    /* synchronization */
    mypneb->PGrid::parall->Vector_SumAll(0, ngs, q);
    mypneb->PGrid::parall->Vector_SumAll(0, ngs * ngs, Am);

    /* deallocate temporary memory from heap */
    delete[] xtmp;
    delete[] gaus_j;
    delete[] gaus_i;
    delete[] exj;
    delete[] exi;

  } /*apc_on*/
}

/*******************************************
 *                                         *
 *            nwpw_apc::dngen_APC          *
 *                                         *
 *******************************************/
void nwpw_apc::dngen_APC(double *dn, bool move) {
  // define parameters
  int ispin = mypneb->ispin;
  int n2ft3d = mypneb->n2ft3d;
  double scal1 = 1.0 / ((double)((mypneb->nx) * (mypneb->ny) * (mypneb->nz)));

  // allocate temporary memory from heap
  double *tmp = mypneb->r_alloc();
  double *dng = mypneb->c_pack_allocate(0);

  // convert dn to dng
  mypneb->rrr_Sum(dn, &dn[(ispin - 1) * n2ft3d], tmp);
  mypneb->r_SMul(scal1, tmp);
  mypneb->rc_fft3d(tmp);
  mypneb->c_pack(0, tmp);
  mypneb->cc_pack_copy(0, tmp, dng);

  // generate APC
  gen_APC(dng, move);

  // deallocate temporary memory from heap
  mypneb->c_pack_deallocate(dng);
  mypneb->r_dealloc(tmp);
}

/*******************************************
 *                                         *
 *          private sumAm_APC              *
 *                                         *
 *******************************************/
static double sumAm_APC(const int ngs, const double *Am) {
  double sum1 = 0.0;
  for (auto i = 0; i < ngs * ngs; ++i)
    sum1 += Am[i];
  return sum1;
}

/*******************************************
 *                                         *
 *          private Amtimesu_APC           *
 *                                         *
 *******************************************/
static double Amtimesu_APC(const int ngs, const double *Am, const double *uq) {
  double sum0 = 0.0;
  for (auto i = 0; i < ngs; ++i)
    for (auto j = 0; j < ngs; ++j)
      sum0 += Am[i + j * ngs] * uq[j];
  return sum0;
}

/*******************************************
 *                                         *
 *          private Vfac_APC               *
 *                                         *
 *******************************************/
static double Vfac_APC(const int ngs, double *Am, const double *uq,
                       const int i) {
  double vsum0 = 0.0;
  double *AA = &Am[i];
  for (auto j = 0; j < ngs; ++j) {
    vsum0 += AA[0] * uq[j];
    AA += ngs;
  }
  // for (auto j=0; j<ngs; ++j)
  //    sum0 = sum0 +  Am[i+j*ngs]*uq[j];

  return vsum0;
}

/*******************************************
 *                                         *
 *            nwpw_apc::VQ_APC             *
 *                                         *
 *******************************************/
/*
    This routine calculates dE/drho(G) where E is a function of model q. In
  order to use this routine Am = inverse A must be calculated for the current
  geometry.

    Entry - nion: number of qm atoms
          - uq: dE/dq(ii) - derivative of E wrt to model charges q(ii)
    Exit - VQ(G) = dE/dq(ii)*dq(ii)/drho(G)

    Note - pspw_gen_APC needs to be called to generate Am=inv(A) before
           this routine is called.
*/

void nwpw_apc::VQ_APC(double *uq, double *VQ) {
  int npack0 = mypneb->npack(0);

  double omega = mypneb->lattice->omega();
  double sumAm = sumAm_APC(ngs, Am);
  double sumAmU = Amtimesu_APC(ngs, Am, uq);
  for (auto i = 0; i < ngs; ++i)
    uq[i] -= (sumAmU / sumAm);

  /* allocate temporary memory from heap */
  double *exi = new double[2 * npack0];
  double *gaus_i = new double[2 * npack0];

  for (auto ii = 0; ii < myion->nion; ++ii) {
    mystrfac->strfac_pack(0, ii, exi);

    for (auto iii = 0; iii < nga; ++iii) {
      int i = iii + ii * nga;
      double afac = omega * Vfac_APC(ngs, Am, uq, i);

      /* gaus_i(G)) */
      mypneb->PGrid::tcc_pack_Mul(0, &gaus[npack0 * iii], exi, gaus_i);

      /* VQ(G) += afac*w(G)*gaus_i(G)) */
      mypneb->PGrid::tcc_pack_aMulAdd(0, afac, w, gaus_i, VQ);
    }
  }
  mypneb->PGrid::c_pack_addzero(0, sumAmU / sumAm, VQ);

  /* deallocate temporary memory from heap */
  delete[] gaus_i;
  delete[] exi;
}

/*     dq(ia,ib)/dR(lb) = Sum(ja,jb) Am(ia,ib;ja,jb) * ( db(ja,jb)/dR(lb) -
 * Sum(ka,kb) (dA(ja,jb;ka,kb)/dR(lb) * q(ka,kb)) ) = Sum(ja,jb) Am(ia,ib;ja,jb)
 * * db(ja,lb)*delta(jb,lb)
 *                       + Sum(ja,jb) Am(ia,ib;ja,jb) * Sum(ka,kb)
 * dA(ja,lb;ka,kb)*q(ka,kb)*delta(jb,lb)
 *                       - Sum(ja,jb) Am(ia,ib;ja,jb) * Sum(ka,kb)
 * dA(ja,jb;ka,lb)*q(ka,kb)*delta(kb,lb)
 *
 *                       = Sum(ja) Am(ia,ib;ja,lb) * db(ja,lb)
 *                       + Sum(ja) Am(ia,ib;ja,lb) * Sum(ka,kb)
 * dA(ja,lb;ka,kb)*q(ka,kb)
 *                       - Sum(ja,jb) Am(ia,ib;ja,jb) * Sum(ka)
 * dA(ja,jb;ka,lb)*q(ka,lb)
 */

/*******************************************
 *                                         *
 *          private generate_dQdR          *
 *                                         *
 *******************************************/
static double generate_dQdR(const int lb, const int nga, const int nion,
                            const double *db, const double *q, const double *uq,
                            const double *dA, const double *Am,
                            const double sumAm, double *dbtmp, double *dAtmp) {
  int ngs = nga * nion;

  /* calculate dbtmp */
  for (auto ja = 0; ja < nga; ++ja) {
    double tmp = 0.0;
    for (auto kb = 0; kb < nion; ++kb)
      for (auto ka = 0; ka < nga; ++ka)
        tmp += dA[ja + lb * nga + (ka + kb * nga) * ngs] * q[ka + kb * nga];
    dbtmp[ja] = tmp;
  }

  /* calculate dAtmp */
  for (auto jb = 0; jb < nion; ++jb)
    for (auto ja = 0; ja < nga; ++ja) {
      double tmp = 0.0;
      for (auto ka = 0; ka < nga; ++ka)
        tmp += dA[ja + jb * nga + (ka + lb * nga) * ngs] * q[ka + lb * nga];
      dAtmp[ja + jb * nga] = tmp;
    }

  double ff = 0.0;
  double sumdQdR = 0.0;
  for (auto ib = 0; ib < nion; ++ib)
    for (auto ia = 0; ia < nga; ++ia) {
      double tmp = 0.0;
      for (auto ja = 0; ja < nga; ++ja)
        tmp += Am[ia + ib * nga + (ja + lb * nga) * ngs] *
               (db[ja + lb * nga] + dbtmp[ja]);

      double tmp2 = 0.0;
      for (auto jb = 0; jb < nion; ++jb)
        for (auto ja = 0; ja < nga; ++ja)
          tmp2 +=
              Am[ia + ib * nga + (ja + jb * nga) * ngs] * dAtmp[ja + jb * nga];
      ff = ff + uq[ia + ib * nga] * (tmp - tmp2);
      sumdQdR = sumdQdR + (tmp - tmp2);
    }
  double fac = sumdQdR / sumAm;

  for (auto ib = 0; ib < nion; ++ib)
    for (auto ia = 0; ia < nga; ++ia) {
      double tmp = 0.0;
      for (auto jb = 0; jb < nion; ++jb)
        for (auto ja = 0; ja < nga; ++ja)
          tmp += Am[ia + ib * nga + (ja + jb * nga) * ngs] * fac;
      ff -= tmp * uq[ia + ib * nga];
    }

  return ff;
}

/*  B  dq(ia,ib)/dR(lb) = Sum(ja,jb) Am(ia,ib;ja,jb) * ( db(ja,jb)/dR(lb) -
 * Sum(ka,kb) (dA(ja,jb;ka,kb)/dR(lb) * q(ka,kb)) ) = Sum(ja,jb) Am(ia,ib;ja,jb)
 * * db(ja,lb)*delta(jb,lb)
 *                       + Sum(ja,jb) Am(ia,ib;ja,jb) * Sum(ka,kb)
 * dA(ja,lb;ka,kb)*q(ka,kb)*delta(jb,lb)
 *                       - Sum(ja,jb) Am(ia,ib;ja,jb) * Sum(ka,kb)
 * dA(ja,jb;ka,lb)*q(ka,kb)*delta(kb,lb)
 *
 *                       = Sum(ja) Am(ia,ib;ja,lb) * db(ja,lb)
 *                       + Sum(ja) Am(ia,ib;ja,lb) * Sum(ka,kb)
 * dA(ja,lb;ka,kb)*q(ka,kb)
 *                       - Sum(ja,jb) Am(ia,ib;ja,jb) * Sum(ka)
 * dA(ja,jb;ka,lb)*q(ka,lb)
 */

/*******************************************
 *                                         *
 *          private generate_dQdR0         *
 *                                         *
 *******************************************/
static double generate_dQdR0(const int lb, const int nga, const int nion,
                             const double *db, const double *q,
                             const double *uq, const double *dA,
                             const double *Am, const double sumAm,
                             double *dbtmp, double *dAtmp) {
  int ngs = nga * nion;

  /* calculate dbtmp */
  for (auto ja = 0; ja < nga; ++ja) {
    double tmp = 0.0;
    for (auto kb = 0; kb < nion; ++kb)
      for (auto ka = 0; ka < nga; ++ka)
        tmp += dA[ja + lb * nga + (ka + kb * nga) * ngs] * q[ka + kb * nga];
    dbtmp[ja] = tmp;
  }

  /* calculate dAtmp */
  for (auto jb = 0; jb < nion; ++jb)
    for (auto ja = 0; ja < nga; ++ja) {
      double tmp = 0.0;
      for (auto ka = 0; ka < nga; ++ka)
        tmp += dA[ja + jb * nga + (ka + lb * nga) * ngs] * q[ka + lb * nga];
      dAtmp[ja + jb * nga] = tmp;
    }

  double ff = 0.0;
  double sumdQdR = 0.0;
  for (auto ib = 0; ib < nion; ++ib)
    for (auto ia = 0; ia < nga; ++ia) {
      double tmp = 0.0;
      for (auto ja = 0; ja < nga; ++ja)
        tmp += Am[ia + ib * nga + (ja + lb * nga) * ngs] *
               (db[ja + lb * nga] + dbtmp[ja]);
      // tmp +=  Am[ia+ib*nga+(ja+lb*nga)*ngs]*(db[ja+lb*nga]);
      // tmp +=  Am[ia+ib*nga+(ja+lb*nga)*ngs]*(dbtmp[ja]);

      double tmp2 = 0.0;
      for (auto jb = 0; jb < nion; ++jb)
        for (auto ja = 0; ja < nga; ++ja)
          tmp2 +=
              Am[ia + ib * nga + (ja + jb * nga) * ngs] * dAtmp[ja + jb * nga];
      ff = ff + uq[ia + ib * nga] * (tmp - tmp2);
      sumdQdR = sumdQdR + (tmp - tmp2);
    }
  double fac = sumdQdR / sumAm;

  for (auto ib = 0; ib < nion; ++ib)
    for (auto ia = 0; ia < nga; ++ia) {
      double tmp = 0.0;
      for (auto jb = 0; jb < nion; ++jb)
        for (auto ja = 0; ja < nga; ++ja)
          tmp += Am[ia + ib * nga + (ja + jb * nga) * ngs] * fac;
      ff -= tmp * uq[ia + ib * nga];
    }

  return ff;
}

/*******************************************
 *                                         *
 *            nwpw_apc::dQdR_APC           *
 *                                         *
 *******************************************/
/*
    This routine calculates dq/Rion where E is a function of model q. In order
  to use this routine Am = inverse A must be calculated
  for the current geometry.

    Entry - nion: number of qm atoms
          - uq: dE/dq(ii) - derivative of E wrt to model charges q(ii)
    Exit - fion = dE/dq(ii)*dq(ii)/dR = uq(ii) * dq(ii)/dR

    Note - pspw_gen_APC needs to be called to Am=inv(A) before
           this routine is called.
*/
void nwpw_apc::dQdR_APC(double *uq, double *fion) {
  int taskid = mypneb->PGrid::parall->taskid();
  int np = mypneb->PGrid::parall->np();
  int nion = myion->nion;
  int nion3 = 3 * nion;
  int ione = 1;
  double rone = 1.0;

  /* allocate memory from stack */
  double dbtmp[nga];
  double dAtmp[ngs];

  double ftmp[nion3];
  std::memset(ftmp, 0, nion3 * sizeof(double));

  double sumAm = sumAm_APC(ngs, Am);

  for (auto ii = taskid; ii < nion; ii += np) {
    ftmp[3 * ii] = generate_dQdR(ii, nga, nion, &b[ngs], q, uq, &A[ngs * ngs],
                                 Am, sumAm, dbtmp, dAtmp);
    ftmp[3 * ii + 1] =
        generate_dQdR(ii, nga, nion, &b[2 * ngs], q, uq, &A[2 * ngs * ngs], Am,
                      sumAm, dbtmp, dAtmp);
    ftmp[3 * ii + 2] =
        generate_dQdR(ii, nga, nion, &b[3 * ngs], q, uq, &A[3 * ngs * ngs], Am,
                      sumAm, dbtmp, dAtmp);
  }
  mypneb->PGrid::parall->Vector_SumAll(0, nion3, ftmp);
  DAXPY_PWDFT(nion3, rone, ftmp, ione, fion, ione);
}

/*******************************************
 *                                         *
 *           nwpw_APC::V_APC_cdft          *
 *                                         *
 *******************************************/
void nwpw_apc::V_APC_cdft(double *dng, double *zv, double *vapc, bool move,
                          double *fion) {
  // generate APC charges qion
  gen_APC(dng, move);
  for (auto ii = 0; ii < myion->nion; ++ii)
    qion[ii] = -Qtot_APC(ii) + zv[myion->katm[ii]];

  // set u from APC uion
  for (auto ii = 0; ii < (myion->nion); ++ii)
    for (auto k = 0; k < nga; ++k)
      u[ii * nga + k] = -uion[ii];

  if (move)
    dQdR_APC(u, fion);

  // Calculate Eapc - cdft,qmmm,cosmo,born??
  Eapc = 0.0;
  for (auto ii = 0; ii < (myion->nion); ++ii)
    Eapc += qion[ii] * uion[ii];

  // Calculate Vapc
  int npack0 = mypneb->npack(0);
  std::memset(vtmp, 0, 2 * npack0 * sizeof(double));
  VQ_APC(u, vtmp);
  mypneb->cc_pack_daxpy(0, 1.0, vtmp, vapc);

  // Calculate Papc
  Papc = mypneb->cc_pack_dot(0, dng, vtmp);

  // F =  -sum(i,j) 0.5*q(i)*(dM/dR)*q(j) - sum u(i)*dq(i)/dR
  // fion = -sum u(i)*dq(i)/dR
}

/*******************************************
 *                                         *
 *           nwpw_apc::V_APC_born          *
 *                                         *
 *******************************************/
void nwpw_apc::V_APC_born(double *dng, double *zv, double *vapc, bool move,
                          double *fion) {
  // generate APC charges qion
  gen_APC(dng, move);
  for (auto ii = 0; ii < myion->nion; ++ii)
    qion[ii] = -Qtot_APC(ii) + zv[myion->katm[ii]];

  myborn->dVdq(qion, uion);

  // Eapc = 0.0;
  // for (auto ii=0; ii<myion->nion; ++ii)
  //    Eapc += qion[ii]*uion[ii];

  // set u from APC uion
  for (auto ii = 0; ii < myion->nion; ++ii)
    for (auto k = 0; k < nga; ++k)
      u[ii * nga + k] = uion[ii];

  // double elocal = mypneb->cc_pack_dot(0,dng,vapc);

  // if (move) std::memcpy(utmp,u,ngs*sizeof(double));
  if (move) {
    myborn->fion(qion, fion);
    dQdR_APC(u, fion);
  }

  // Calculate Vapc
  int npack0 = mypneb->npack(0);
  std::memset(vtmp, 0, 2 * npack0 * sizeof(double));
  VQ_APC(u, vtmp);
  mypneb->cc_pack_daxpy(0, 1.0, vtmp, vapc);

  // Calculate Eborn
  Eapc = myborn->energy(qion);

  // Calculate Pborn
  Papc = mypneb->cc_pack_dot(0, dng, vtmp);
  // double Papc1 = mypneb->cc_pack_dot(0,dng,vapc);

  // std::cout << "Elocal=" << std::fixed << std::setw(15) <<
  // std::setprecision(11) << elocal << std::endl; std::cout << std::endl;
  // std::cout << "PAPC=" << std::fixed << std::setw(15) <<
  // std::setprecision(11) << Papc << std::endl; std::cout << "PAPC1=" <<
  // std::fixed << std::setw(15) << std::setprecision(11) << Papc1 << std::endl;
  // std::cout << "vapc=" << std::fixed << std::setw(15) <<
  // std::setprecision(11) << vapc[0] << " " << vapc[1] << std::endl;

  // if (move)
  //{
  //    myborn->fion(qion,fion);
  //    dQdR_APC(utmp,fion);
  // }
}

/*******************************************
 *                                         *
 *            nwpw_apc::V_APC              *
 *                                         *
 *******************************************/
void nwpw_apc::V_APC(double *dng, double *zv, double *vapc, bool move,
                     double *fion) {
  if ((apc_on) && (v_apc_on)) {
    if (born_on)
      this->V_APC_born(dng, zv, vapc, move, fion);
    else
      this->V_APC_cdft(dng, zv, vapc, move, fion);
  }
}

/*******************************************
 *                                         *
 *            nwpw_apc::f_APC_cdft         *
 *                                         *
 *******************************************/
void nwpw_apc::f_APC_cdft(double *dng, double *zv, double *fion) {
  // generate APC charges qion
  gen_APC(dng, true);
  for (auto ii = 0; ii < (myion->nion); ++ii)
    qion[ii] = -Qtot_APC(ii) + zv[myion->katm[ii]];

  // double uintern[myion->nion];
  // memset(uintern,0,myion->nion*sizeof(double));
  // for (auto ii=0; ii<(myion->nion); ++ii)
  // for (auto jj=0; jj<(myion->nion); ++jj)
  //{
  //    double dx = myion->rion1[3*ii]  -myion->rion1[3*jj];
  //    double dy = myion->rion1[3*ii+1]-myion->rion1[3*jj+1];
  //    double dz = myion->rion1[3*ii+2]-myion->rion1[3*jj+2];
  //    double dist2 = dx*dx + dy*dy + dz*dz;
  //    if (dist2>1.0e-4)
  //    {
  //       double f = std::sqrt(dist2);
  //
  //      uintern[ii] += 0.5*q[jj]/f;
  //       uintern[jj] += 0.5*q[ii]/f;
  //    }
  // }

  // set u from APC uion
  for (auto ii = 0; ii < (myion->nion); ++ii)
    for (auto k = 0; k < nga; ++k)
      u[ii * nga + k] = -uion[ii];

  // double sumAm  = sumAm_APC(ngs,Am);
  // double sumAmU = Amtimesu_APC(ngs,Am,u);
  // for (auto i=0; i<ngs; ++i)
  //    u[i] -= (sumAmU/sumAm);

  // fion = -sum u(i)*dq(i)/dR
  // std::cout << "fion before: " << fion[0] << " "
  //                             << fion[1] << " "
  //                             << fion[2] << " " << std::endl
  //          << "           : " << fion[3] << " "
  //                             << fion[4] << " "
  //                             << fion[5] << " " << std::endl
  //          << "           : " << fion[6] << " "
  //                             << fion[7] << " "
  //                             << fion[8] << " " << std::endl;
  dQdR_APC(u, fion);
  // std::cout << "fion +apc  : " << fion[0] << " "
  //                              << fion[1] << " "
  //                              << fion[2] << " " << std::endl
  //           << "           : " << fion[3] << " "
  //                              << fion[4] << " "
  //                              << fion[5] << " " << std::endl
  //           << "           : " << fion[6] << " "
  //                              << fion[7] << " "
  //                              << fion[8] << " " << std::endl;
}

/*******************************************
 *                                         *
 *            nwpw_apc::f_APC_born         *
 *                                         *
 *******************************************/
void nwpw_apc::f_APC_born(double *dng, double *zv, double *fion) {
  // generate APC charges qion
  gen_APC(dng, true);
  for (auto ii = 0; ii < myion->nion; ++ii)
    qion[ii] = -Qtot_APC(ii) + zv[myion->katm[ii]];

  myborn->dVdq(qion, uion);

  // set u from APC uion
  for (auto ii = 0; ii < (myion->nion); ++ii)
    for (auto k = 0; k < nga; ++k)
      u[ii * nga + k] = uion[ii];

  myborn->fion(qion, fion);
  dQdR_APC(u, fion);
}

/*******************************************
 *                                         *
 *            nwpw_apc::f_APC              *
 *                                         *
 *******************************************/

void nwpw_apc::f_APC(double *dng, double *zv, double *fion) {
  if ((apc_on) && (v_apc_on)) {
    if (born_on)
      this->f_APC_born(dng, zv, fion);
    else
      this->f_APC_cdft(dng, zv, fion);
  }
}

/*******************************************
 *                                         *
 *            nwpw_apc::Qtot_APC           *
 *                                         *
 *******************************************/
double nwpw_apc::Qtot_APC(const int ii) {
  double qq = 0.0;
  if (nga > 0)
    for (auto n = 0; n < nga; ++n)
      qq += q[n + ii * nga];
  return qq;
}

/*******************************************
 *                                         *
 *         nwpw_apc::shortprint_APC        *
 *                                         *
 *******************************************/
std::string nwpw_apc::shortprint_APC() {
  std::stringstream stream;

  stream << std::endl;
  stream << " APC Potential:" << std::endl;
  for (auto ii = 0; ii < myion->nion; ++ii) {
    stream << std::setw(14) << std::fixed << std::setprecision(9) << uion[ii];
    if ((((ii + 1) % 10) == 0) && (ii != (myion->nion - 1)))
      stream << std::endl;
  }
  stream << std::endl << std::endl;
  stream << " APC Point Charges:" << std::endl;
  for (auto ii = 0; ii < myion->nion; ++ii) {
    stream << std::setw(14) << std::fixed << std::setprecision(9) << qion[ii];
    if ((((ii + 1) % 10) == 0) && (ii != (myion->nion - 1)))
      stream << std::endl;
  }
  stream << std::endl << std::endl;

  return stream.str();
}

/*******************************************
 *                                         *
 *            nwpw_apc::print_APC          *
 *                                         *
 *******************************************/
std::string nwpw_apc::print_APC(const double *zv) {
  std::stringstream stream;

  stream << std::endl
         << "*************************************************************"
         << std::endl
         << "**                                                         **"
         << std::endl
         << "**          PSPW Atomic Point Charge (APC) Analysis        **"
         << std::endl
         << "**                                                         **"
         << std::endl
         << "**   Point charge analysis based on paper by P.E. Blochl   **"
         << std::endl
         << "**         (J. Chem. Phys. vol 103, page 7422, 1995)       **"
         << std::endl
         << "**                                                         **"
         << std::endl
         << "*************************************************************"
         << std::endl
         << std::endl
         << " nwpw_APC object" << std::endl
         << " ---------------" << std::endl
         << " nga = " << std::setw(3) << nga << " ngs = " << std::setw(5) << ngs
         << std::endl
         << " Gc  = " << std::fixed << std::setprecision(5) << std::setw(9)
         << Gc << std::endl;
  for (auto i = 0; i < nga; ++i)
    stream << " APC gamma: " << i << " " << gamma[i] << std::endl;
  if (v_apc_on)
    stream << " - self-consistent" << std::endl;
  else
    stream << " - not self-consistent" << std::endl;

  if (v_apc_on) {
    stream << std::endl;
    stream << " APC Potential:" << std::endl;
    for (auto ii = 0; ii < myion->nion; ++ii) {
      stream << std::setw(14) << std::fixed << std::setprecision(9) << uion[ii];
      if ((((ii + 1) % 10) == 0) && (ii != (myion->nion - 1)))
        stream << std::endl;
    }
    stream << std::endl << std::endl;
    stream << " APC Point Charges:" << std::endl;
    for (auto ii = 0; ii < myion->nion; ++ii) {
      stream << std::setw(14) << std::fixed << std::setprecision(9) << qion[ii];
      if ((((ii + 1) % 10) == 0) && (ii != (myion->nion - 1)))
        stream << std::endl;
    }
    stream << std::endl;
  }

  stream << std::endl << std::endl;
  stream << " charge analysis on each atom" << std::endl
         << " ----------------------------" << std::endl
         << std::endl
         << "      no  atom        Qelc        Qion      Qtotal" << std::endl
         << "   -----  ----     -------     -------     -------" << std::endl;
  double sume = 0.0;
  double sumi = 0.0;
  for (auto ii = 0; ii < myion->nion; ++ii) {
    int ia = myion->katm[ii];
    double sum = 0.0;
    for (auto i = 0; i < nga; ++i)
      sum += q[i + ii * nga];
    sume -= sum;
    sumi += zv[ia];
    stream << std::setw(8) << ii + 1 << "    " << std::left << std::setw(2)
           << myion->symbol(ii) << std::right << std::setw(12) << std::fixed
           << std::setprecision(3) << -sum << std::setw(12) << std::fixed
           << std::setprecision(3) << zv[ia] << std::setw(12) << std::fixed
           << std::setprecision(3) << zv[ia] - sum << std::endl;
  }
  stream << std::setw(14) << "Total Q" << std::setw(12) << std::fixed
         << std::setprecision(3) << sume << std::setw(12) << std::fixed
         << std::setprecision(3) << sumi << std::setw(12) << std::fixed
         << std::setprecision(3) << sume + sumi << std::endl;

  stream << std::endl
         << std::endl
         << " Gaussian coefficients of model density" << std::endl
         << " --------------------------------------" << std::endl
         << std::endl;
  stream << "      no  atom";
  stream << std::setw(7) << "g=" << std::fixed << std::setprecision(3) << 0.0;
  for (auto i = 0; i < nga; ++i)
    stream << std::setw(7) << "g=" << std::fixed << std::setprecision(3)
           << gamma[i];
  stream << std::endl;
  stream << "   -----  ----";
  for (auto i = 0; i < (nga + 1); ++i)
    stream << std::setw(12) << "-------";
  stream << std::endl;
  for (auto ii = 0; ii < myion->nion; ++ii) {
    int ia = myion->katm[ii];
    stream << std::setw(8) << ii + 1 << std::setw(6) << myion->symbol(ii)
           << std::setw(12) << std::fixed << std::setprecision(3) << zv[ia];
    for (auto i = 0; i < nga; ++i)
      stream << std::setw(12) << std::fixed << std::setprecision(3)
             << -q[i + ii * nga];
    stream << std::endl;
  }

  /* include born output */
  if (born_on) {
    double qion[myion->nion];
    for (auto ii = 0; ii < myion->nion; ++ii) {
      int ia = myion->katm[ii];
      double sum = 0.0;
      for (auto i = 0; i < nga; ++i)
        sum += q[i + ii * nga];
      qion[ii] = zv[ia] - sum;
    }
    stream << myborn->Qprint(qion);
  }

  stream << std::endl << std::endl;

  return stream.str();
}

} // namespace pwdft

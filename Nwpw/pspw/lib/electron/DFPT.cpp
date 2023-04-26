
#include "Coulomb12.hpp"
#include "Ion.hpp"
#include "Kinetic.hpp"
#include "Pneb.hpp"
#include "Pseudopotential.hpp"
#include "Strfac.hpp"
#include "blas.h"
#include "exchange_correlation.hpp"
//#include        "v_exc.hpp"

#include "psi_H.hpp"

#include "DFPT.hpp"

namespace pwdft {

/********************************************
 *                                          *
 *     DFPT_Operators::DFPT_Operators       *
 *                                          *
 ********************************************/
DFPT_Operators::DFPT_Operators(Pneb *mygrid0, Kinetic_Operator *myke0,
                               Coulomb12_Operator *mycoul0, XC_Operator *myxc0,
                               Pseudopotential *mypsp0, const double eps0) {
  mygrid = mygrid0;
  myke = myke0;
  mypsp = mypsp0;
  mycoulomb12 = mycoul0;
  myxc = myxc0;
  eps = eps0;

  ispin = mygrid->ispin;
  neall = mygrid->neq[0] + mygrid->neq[1];

  aperiodic = mycoulomb12->has_coulomb2;
  periodic = mycoulomb12->has_coulomb1;

  /* allocate memory */
  Hpsi1 = mygrid->g_allocate(1);
  psi1 = mygrid->g_allocate(1);
  b1 = mygrid->g_allocate(1);
  psi1_r = mygrid->h_allocate();

  dn1 = mygrid->r_nalloc(ispin);
  dnallp = mygrid->r_nalloc(ispin);
  dnallm = mygrid->r_nalloc(ispin);
  rho1 = mygrid->r_alloc();
  dng1 = mygrid->c_pack_allocate(0);

  vxc1 = mygrid->r_nalloc(ispin);
  xcpp = mygrid->r_nalloc(ispin);
  xcpm = mygrid->r_nalloc(ispin);
  if (periodic)
    vc1 = mygrid->c_pack_allocate(0);
  if (aperiodic)
    vc1 = mygrid->r_alloc();

  tmp = mygrid->r_alloc();

  hmltmp = mygrid->m_allocate(-1, 1);

  omega = mygrid->lattice->omega();
  scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
  scal2 = 1.0 / omega;
  dv = omega * scal1;

  n2ft3d = (mygrid->n2ft3d);
  shift1 = 2 * (mygrid->npack(1));
  npack1 = shift1;
  shift2 = (mygrid->n2ft3d);
}

/********************************************
 *                                          *
 *           DFPT_Operators::start          *
 *                                          *
 ********************************************/
void DFPT_Operators::start(double *psi, double *psi_r, double *dnall,
                           double *dEpertdpsi_r) {
  psi0 = psi;
  psi0_r = psi_r;
  dnall0 = dnall;
  dEpertdpsi0_r = dEpertdpsi_r;
}

/********************************************
 *                                          *
 *           DFPT_Operators::run            *
 *                                          *
 ********************************************/
void DFPT_Operators::run() {
  ++counter;

  /* convert psi1(G) to psi1(r) */
  int indx1 = 0;
  int indx2 = 0;
  for (int i = 0; i < neall; ++i) {
    mygrid->cc_pack_copy(1, psi1 + indx1, psi1_r + indx2);
    mygrid->c_unpack(1, psi1_r + indx2);
    mygrid->cr_fft3d(psi1_r + indx2);
    indx1 += shift1;
    indx2 += shift2;
  }

  /* generate dn1 */
  mygrid->hhr_aSumMul(scal2, psi0_r, psi1_r, dn1);

  /* generate rho1 and dng */
  mygrid->rrr_Sum(dn1, dn1 + (ispin - 1) * n2ft3d, rho1);
  mygrid->rr_SMul(scal1, rho1, tmp);
  mygrid->rc_fft3d(tmp);
  mygrid->c_pack(0, tmp);
  mygrid->cc_pack_copy(0, tmp, dng1);

  /* generate coulomb potential */
  if (periodic)
    mycoulomb12->mycoulomb1->vcoulomb(dng1, vc1);
  if (aperiodic) {
    mygrid->rrr_Sum(dn1, dn1 + (ispin - 1) * n2ft3d, rho1);
    mycoulomb12->mycoulomb2->vcoulomb(rho1, vc1);
  }

  /*  generate dnallp and dnallm - used for generate xcp1 */
  for (int ms = 0; ms < ispin; ++ms) {
    mygrid->rrr_SMulAdd(eps, dn1 + ms * n2ft3d, dnall0 + ms * n2ft3d,
                        dnallp + ms * n2ft3d);
    mygrid->rrr_SMulAdd(-eps, dn1 + ms * n2ft3d, dnall0 + ms * n2ft3d,
                        dnallm + ms * n2ft3d);
  }

  /* generate exchange-correlation potential - vxc1 = 0.5*eps*(xcpp - xcpm) */
  myxc->v_exc_all(ispin, dnallp, xcpp, tmp);
  myxc->v_exc_all(ispin, dnallm, xcpm, tmp);
  for (int ms = 0; ms < ispin; ++ms) {
    mygrid->arrr_Minus((0.5 / eps), xcpp + ms * n2ft3d, xcpm + ms * n2ft3d,
                       vxc1 + ms * n2ft3d);
  }

  /* generate b1 -  b1(ms,i) = rc_fft[(vxc1(ms)+vc1)*psi0(ms,i) +
   * dEpertdpsi(ms,i)] */
  indx1 = 0;
  for (int ms = 0; ms < ispin; ++ms)
    for (int i = 0; i < (mygrid->neq[ms]); ++i) {
      mygrid->rrrrr_SumMulAdd(vc1, vxc1 + ms * n2ft3d, psi0_r + indx2,
                              dEpertdpsi0_r + indx2, tmp);
      mygrid->r_SMul(scal1, tmp);
      mygrid->rc_fft3d(tmp);
      mygrid->c_pack(1, tmp);
      mygrid->cc_pack_copy(1, tmp, b1 + indx1);
      indx1 += shift1;
    }

  /* project out psi0 from b1 - b1 = Pe*b1 = (1-|psi0><psi0|)*b1 */
  mygrid->ffm_Multiply(-1, psi0, b1, hmltmp);
  mygrid->fmf_Multiply(-1, psi0, hmltmp, -1.0, b1, 1.0);
}

} // namespace pwdft

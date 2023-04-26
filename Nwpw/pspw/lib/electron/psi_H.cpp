

#include "Kinetic.hpp"
#include "PGrid.hpp"
#include "Pseudopotential.hpp"

namespace pwdft {

/*************************************
 *                                   *
 *             psi_H                 *
 *                                   *
 *************************************

   This routine calculates

        Hpsi_k = KE*psi_k + Vnl*psi_k +VSic*psi_k + FFT[(vall+xcp)*psi_r]

   where vall = iFFT[Vl+Vc] + Vfield

    Entry - ispin,neq           - spin and number of electron in each spin
            psi_k,psi_r         - orbitals in k-space and r-space
            vl                  - local psp in k-space
            vc                  - coulomb potential in k-space
            xcp                 - xc potential in r-space
            move                - flag to compute ionic forces
    Exit - Hpsi - gradient in k-space
           fion   - ionic forces
*/

void psi_H(Pneb *mygrid, Kinetic_Operator *myke, Pseudopotential *mypsp,
           double *psi, double *psi_r, double *vl, double *vc, double *xcp,
           double *Hpsi, bool move, double *fion)

{
  int indx1 = 0;
  int indx2 = 0;
  int indx1n = 0;
  int indx2n = 0;
  int ispin = mygrid->ispin;
  int shift1 = 2 * (mygrid->npack(1));
  int shift2 = (mygrid->n2ft3d);
  int n2ft3d = (mygrid->n2ft3d);
  int ms = 0;
  int n1 = mygrid->neq[0];
  int n2 = mygrid->neq[0] + mygrid->neq[1];
  ;

  bool done = false;

  double omega = mygrid->lattice->omega();
  double scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
  double scal2 = 1.0 / omega;

  /* allocate temporary memory */
  double *vall = mygrid->r_alloc();
  double *vpsi = mygrid->r_alloc();
  double *tmp = mygrid->r_alloc();

  /* apply k-space operators */
  myke->ke(psi, Hpsi);

  /* apply non-local PSP  - Expensive */
  mypsp->v_nonlocal_fion(psi, Hpsi, move, fion);

  /* apply r-space operators  - Expensive*/
  mygrid->cc_pack_SMul(0, scal2, vl, vall);
  mygrid->cc_pack_Sum2(0, vc, vall);
  mygrid->c_unpack(0, vall);
  mygrid->cr_fft3d(vall);

  /* add v_field to vall */
  if (mypsp->myefield->efield_on)
    mygrid->rr_Sum(mypsp->myefield->v_field, vall);

  /*
     for (ms=0; ms<ispin; ++ms)
     {
        mygrid->rrr_Sum(vall,xcp+ms*n2ft3d,tmp);
        for (int i=0; i<(mygrid->neq[ms]); ++i)
        {
           mygrid->rrr_Mul(tmp,psi_r+indx2,vpsi);
           mygrid->rc_fft3d(vpsi);
           mygrid->c_pack(1,vpsi);
           mygrid->cc_pack_daxpy(1,(-scal1),vpsi,Hpsi+indx1);

           indx1 += shift1;
           indx2 += shift2;
        }
     }
     */

  {
    nwpw_timing_function ftimer(1);

    mygrid->rrr_Sum(vall, xcp, tmp);
    while (!done) {
      if (indx1 < n2) {
        if (indx1 >= n1) {
          ms = 1;
          mygrid->rrr_Sum(vall, xcp + ms * n2ft3d, tmp);
        }

        mygrid->rrr_Mul(tmp, psi_r + indx1n, vpsi);

        mygrid->rc_pfft3f_queuein(1, vpsi);
        indx1n += shift2;
        ++indx1;
      }

      if ((mygrid->rc_pfft3f_queuefilled()) || (indx1 >= n2)) {
        mygrid->rc_pfft3f_queueout(1, vpsi);
        mygrid->cc_pack_daxpy(1, (-scal1), vpsi, Hpsi + indx2n);
        indx2n += shift1;
        ++indx2;
      }
      done = ((indx1 >= n2) && (indx2 >= n2));
    }
  }

  /* deallocate temporary memory */
  mygrid->r_dealloc(tmp);
  mygrid->r_dealloc(vpsi);
  mygrid->r_dealloc(vall);
}

/*************************************
 *                                   *
 *             psi_Hv4               *
 *                                   *
 *************************************
   This routine calculates

        Hpsi_k = KE*psi_k + Vnl*psi_k +VSic*psi_k + FFT[(vall+xcp)*psi_r]

   where vall = iFFT[Vsr_l]+ vc + vlr_l + Vfield

   Entry - ispin,neq           - spin and number of electron in each spin
           psi_k,psi_r         - orbitals in k-space and r-space
           vsr_l               - short-range local psp in k-space
           vlr_l               - long-range local psp in r-space
           vc                  - coulomb potential in r-space
           xcp                 - xc potential in r-space
           move                - flag to compute ionic forces
   Exit  - Hpsi - gradient in k-space
           fion   - ionic forces
*/

void psi_Hv4(Pneb *mygrid, Kinetic_Operator *myke, Pseudopotential *mypsp,
             double *psi, double *psi_r, double *vsr_l, double *vlr_l,
             double *vc, double *xcp, double *Hpsi, bool move, double *fion) {
  int indx1 = 0;
  int indx2 = 0;
  int ispin = mygrid->ispin;
  int shift1 = 2 * (mygrid->npack(1));
  int shift2 = (mygrid->n2ft3d);
  int n2ft3d = (mygrid->n2ft3d);

  double omega = mygrid->lattice->omega();
  double scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
  double scal2 = 1.0 / omega;

  /* allocate temporary memory */
  double *vall = mygrid->r_alloc();
  double *vpsi = mygrid->r_alloc();
  double *tmp = mygrid->r_alloc();

  /* apply k-space operators */
  myke->ke(psi, Hpsi);

  /* apply non-local PSP  - Expensive */
  mypsp->v_nonlocal_fion(psi, Hpsi, move, fion);

  /* add up k-space potentials, vall = scal2*vsr_l */
  mygrid->cc_pack_SMul(0, scal2, vsr_l, vall);
  mygrid->c_unpack(0, vall);
  mygrid->cr_fft3d(vall);

  /* add vall += vlr_l + vc */
  mygrid->rrr_Sum2Add(vlr_l, vc, vall);

  /* add v_field to vall */
  if (mypsp->myefield->efield_on)
    mygrid->rr_Sum(mypsp->myefield->v_field, vall);

  /* apply r-space operators  - Expensive*/
  {
    nwpw_timing_function ftimer(1);

    for (int ms = 0; ms < ispin; ++ms) {
      mygrid->rrr_Sum(vall, xcp + ms * n2ft3d, tmp);
      for (int i = 0; i < (mygrid->neq[ms]); ++i) {
        mygrid->rrr_Mul(tmp, psi_r + indx2, vpsi);
        mygrid->rc_fft3d(vpsi);
        mygrid->c_pack(1, vpsi);
        mygrid->cc_pack_daxpy(1, (-scal1), vpsi, Hpsi + indx1);

        indx1 += shift1;
        indx2 += shift2;
      }
    }
  }

  /* deallocate temporary memory */
  mygrid->r_dealloc(tmp);
  mygrid->r_dealloc(vpsi);
  mygrid->r_dealloc(vall);
}

/*************************************
 *                                   *
 *             psi1_H_DFPT           *
 *                                   *
 *************************************

   This routine is used for DFPT and calculates

        Hks*psi1_k = KE*psi1_k + Vnl*psi1_k + rc_FFT[vall_r*psi1_r]

    Entry - psi1_k,psi1_r - orbitals in k-space and r-space
            vall_r        - local psp + coulomb potential and xc and e-field in
 r-space hml           - <psi|H|psi> matrix Exit - Hpsi1 - gradient in k-space
*/

void psi1_H_DFPT(Pneb *mygrid, Kinetic_Operator *myke, Pseudopotential *mypsp,
                 double *psi1, double *psi1_r, double *vall_r, double *hml,
                 double *Hpsi1)

{
  int indx1 = 0;
  int indx2 = 0;
  int ispin = mygrid->ispin;
  int shift1 = 2 * (mygrid->npack(1));
  int shift2 = (mygrid->n2ft3d);
  int n2ft3d = (mygrid->n2ft3d);

  double scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));

  /* allocate temporary memory */
  double *vpsi1 = mygrid->r_alloc();

  /* apply k-space operators */
  myke->ke(psi1, Hpsi1);

  /* apply non-local PSP  - Expensive */
  mypsp->v_nonlocal(psi1, Hpsi1);

  /* apply r-space operators  - Expensive*/
  for (int ms = 0; ms < ispin; ++ms) {
    for (int i = 0; i < (mygrid->neq[ms]); ++i) {
      mygrid->rrr_Mul(vall_r + ms * n2ft3d, psi1_r + indx2, vpsi1);
      mygrid->rc_fft3d(vpsi1);
      mygrid->c_pack(1, vpsi1);
      mygrid->cc_pack_daxpy(1, (scal1), vpsi1,
                            Hpsi1 + indx1); // generates  -Hpsi1

      indx1 += shift1;
      indx2 += shift2;
    }
  }

  /* generate -(H - hml0)*psi1 = -H*psi1 + hml0*psi1) */
  mygrid->fmf_Multiply(1, psi1, hml, 1.0, Hpsi1, 1.0);

  /* deallocate temporary memory */
  mygrid->r_dealloc(vpsi1);
}

/*************************************
 *                                   *
 *             gen_vall_DFPT         *
 *                                   *
 *************************************

   This routine calculates vall_r = iFFT[vl+vc] + xcp + Vfield

    Entry - vl    - local psp in k-space
            vc    - coulomb potential in k-space
            xcp   - xc potential in r-space
    Exit - vall_r - spin potential in r-space
*/

void gen_vall_DFPT(Pneb *mygrid, Pseudopotential *mypsp, double *vl, double *vc,
                   double *xcp, double *vall_r) {
  int ispin = mygrid->ispin;
  int n2ft3d = (mygrid->n2ft3d);

  double omega = mygrid->lattice->omega();
  double scal2 = 1.0 / omega;

  /* allocate temporary memory */
  double *tmp = mygrid->r_alloc();

  /* generate vall_r */
  mygrid->cc_pack_SMul(0, scal2, vl, tmp);
  mygrid->cc_pack_Sum2(0, vc, tmp);
  mygrid->c_unpack(0, tmp);
  mygrid->cr_fft3d(tmp);

  /* add v_field to tmp */
  if (mypsp->myefield->efield_on)
    mygrid->rr_Sum(mypsp->myefield->v_field, tmp);

  /* generate vall_r */
  for (int ms = 0; ms < ispin; ++ms)
    mygrid->rrr_Sum(tmp, xcp + ms * n2ft3d, vall_r + ms * n2ft3d);

  /* deallocate temporary memory */
  mygrid->r_dealloc(tmp);
}

/*************************************
 *                                   *
 *          gen_vall_DFPT_v4         *
 *                                   *
 *************************************

   This routine calculates vall_r = iFFT[vsr_l] + vlr_l + vc + vlr + xcp +
 Vfield

    Entry - vsr_l - short-range local psp in k-space
            vlr_l - long-range local psp in r-space
            vc    - coulomb potential in r-space
            xcp   - xc potential in r-space
    Exit - vall_r - spin potential in r-space
*/

void gen_vall_v4_DFPT(Pneb *mygrid, Pseudopotential *mypsp, double *vsr_l,
                      double *vlr_l, double *vc, double *xcp, double *vall_r) {
  int ispin = mygrid->ispin;
  int n2ft3d = (mygrid->n2ft3d);

  double omega = mygrid->lattice->omega();
  double scal2 = 1.0 / omega;

  /* allocate temporary memory */
  double *tmp = mygrid->r_alloc();

  /* add up k-space potentials, tmp = scal2*vsr_l */
  mygrid->cc_pack_SMul(0, scal2, vsr_l, tmp);
  mygrid->c_unpack(0, tmp);
  mygrid->cr_fft3d(tmp);

  /* add tmp += vlr_l + vc */
  mygrid->rrr_Sum2Add(vlr_l, vc, tmp);

  /* add v_field to vall */
  if (mypsp->myefield->efield_on)
    mygrid->cc_pack_Sum2(0, mypsp->myefield->v_field, tmp);

  /* generate vall_r */
  for (int ms = 0; ms < ispin; ++ms)
    mygrid->rrr_Sum(tmp, xcp + ms * n2ft3d, vall_r + ms * n2ft3d);

  /* deallocate temporary memory */
  mygrid->r_dealloc(tmp);
}

} // namespace pwdft

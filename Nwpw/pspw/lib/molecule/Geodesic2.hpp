#ifndef _GEODESIC2_HPP_
#define _GEODESIC2_HPP_

#pragma once

#include "Electron.hpp"
#include "Molecule.hpp"
#include "Pneb.hpp"
#include <cmath>

//#include	"util.hpp"

namespace pwdft {

class Geodesic2 {

  int minimizer;
  Molecule *mymolecule;
  Electron_Operators *myelectron;

  // tmp space of m matrices
  double *R, *A, *MM, *NN, *TT;

  // tmp space of m4 matrices
  double *T4, *V4, *W4, *A4, *B4, *R4;

  // tmp space for ffm multiplies
  double *Q, *Hold, *S;
  int nemax, npack1;

public:
  Pneb *mygrid;

  /* Constructors */
  Geodesic2(int minimizer0, Molecule *mymolecule0) {
    mymolecule = mymolecule0;
    minimizer = minimizer0;
    myelectron = mymolecule->myelectron;
    mygrid = mymolecule->mygrid;

    nemax = mygrid->neq[0] + mygrid->neq[1];
    npack1 = mygrid->PGrid::npack(1);

    Hold = mygrid->g_allocate(1);
    Q = mygrid->g_allocate(1);

    R = mygrid->m_allocate(-1, 1);
    A = mygrid->m_allocate(-1, 1);
    MM = mygrid->m_allocate(-1, 1);
    NN = mygrid->m_allocate(-1, 1);
    TT = mygrid->m_allocate(-1, 1);

    T4 = mygrid->m4_allocate(-1, 1);
    V4 = mygrid->m4_allocate(-1, 1);
    W4 = mygrid->m4_allocate(-1, 1);
    A4 = mygrid->m4_allocate(-1, 1);
    B4 = mygrid->m4_allocate(-1, 1);
    R4 = mygrid->m4_allocate(-1, 1);

    S = new double[2 * (mygrid->ne[0] + mygrid->ne[1])];
    // tmp space
  }

  /* destructor */
  ~Geodesic2() {
    delete[] Hold;
    delete[] Q;

    delete[] R;
    delete[] A;
    delete[] MM;
    delete[] NN;

    delete[] T4;
    delete[] V4;
    delete[] W4;
    delete[] A4;
    delete[] B4;
    delete[] R4;

    delete[] S;
  }

  /*****************************************
   *                                       *
   *               start                   *
   *                                       *
   *****************************************

   This function determines the pxp matrices R and YA, and
   the orthogonal nxp matrix Q.   Q and R are determined from
   the QR decomposition of the projected direction (I-YY^t)H, and
   YH is defined as the Lagrange Multiplier pxp matrix Y^tH.
   */

  double start(double *Y, double *H, double *max_sigma, double *min_sigma) {

    int nn = 2 * npack1 * nemax;
    int one = 1;
    double mrone = 1.0;

    //**** Hold <-- H ****
    std::memcpy(Hold, H, 2 * nemax * npack1 * sizeof(double));

    //**** calculate A=<Y|H> ****
    mygrid->ffm_Multiply(-1, Y, H, A);

    //**** calculate Q=(I-YYt)H - should not be necessary but just in case ****
    mygrid->fmf_Multiply(-1, Y, A, 1.0, Q, 0.0);
    DAXPY_PWDFT(nn, mrone, H, one, Q, one);
    DSCAL_PWDFT(nn, mrone, Q, one);

    //**** calculate QR using Modified Gram-Schmidt ****
    mygrid->fm_QR(-1, Q, R);

    //**** generate T from A and R ****
    mygrid->mmm4_AR_to_T4(-1, A, R, T4);

    //**** Factor T--> V,W,and S ****
    mygrid->m4_FactorSkew(-1, T4, V4, W4, S);

    /* calculate  and return 2*<A|H|psi> */
    return (2.0 * myelectron->eorbit(H));
  }

  /*****************************************
   *                                       *
   *               get                     *
   *                                       *
   *****************************************
   This routine calculates

   Ynew = Yold*M(t) + Q*N(t)

   where
       -    -               - -
      | M(t) | = Exp(t*T)* | I |
      | N(t) |             | 0 |
       -    -               - -
   */
   void get(double t, double *Yold, double *Ynew) 
   {
      this->get_MandN(t, MM, NN);
      mygrid->fmf_Multiply(-1, Yold, MM, 1.0, Ynew, 0.0);
      mygrid->fmf_Multiply(-1, Q, NN, 1.0, Ynew, 1.0);
   }

  /*****************************************
   *                                       *
   *               transport               *
   *                                       *
   *****************************************
   This routine calculates

   Hnew = Hold*M(t)    + Yold*R^t*N(t)

   where
      -    -               - -
     | M(t) | = Exp(t*T)* | I |
     | N(t) |             | 0 |
      -    -               - -
   */
   void transport(double t, double *Yold, double *Hnew) 
   {

      this->get_MandN(t, MM, NN);

      //*** TT(t) = -R^t*NN(t) ***
      mygrid->mmm_Multiply2(-1, R, NN, -1.0, TT, 0.0);

      //*** Hnew <-- Hold*M(t) + Yold*TT(t) ***
      mygrid->fmf_Multiply(-1, Hold, MM, 1.0, Hnew, 0.0);
      mygrid->fmf_Multiply(-1, Yold, TT, 1.0, Hnew, 1.0);
   }

  /*****************************************
   *                                       *
   *             get_MandN                 *
   *                                       *
   *****************************************
   This routine returns
       -    -               - -
      | M(t) | = Exp(t*T)* | I |
      | N(t) |             | 0 |
       -    -               - -
  where

     T =  U*Sigma*U^H, with U=(V+iW)

     is a skew matrix that is decomposed into V,W,and Sigma

   */
   void get_MandN(const double t, double *M, double *N) 
   {
      mygrid->m4_RotationSkew(-1, t, V4, W4, S, A4, B4, R4);
      mygrid->m4_R4_to_MN(-1, R4, M, N);
   }

   void psi_1transport(double t, double *H0) 
   {
      this->transport(t, mymolecule->psi1, H0);
   }

   double energy(double t) 
   {
      this->get(t, mymolecule->psi1, mymolecule->psi2);
      return (mymolecule->psi2_energy());
   }

   double denergy(double t) 
   {
      this->transport(t, mymolecule->psi1, mymolecule->psi2);
      return (2.0 * mymolecule->psi2_eorbit());
   }

   void psi_final(double t) 
   { 
      this->get(t, mymolecule->psi1, mymolecule->psi2); 
   }

};

} // namespace pwdft

#endif

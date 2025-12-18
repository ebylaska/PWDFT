#ifndef _BAND_GEODESIC_HPP_
#define _BAND_GEODESIC_HPP_
/*
 * band_Geodesic.hpp
 * 
 * This header file defines the band_Geodesic class within the pwdft namespace, representing a Grassmann manifold used in quantum chemical computations, particularly in plane wave density functional theory (PWDFT) calculations.
 * 
 * Key Features:
 * - The class handles computations related to Grassmann manifold projections, which are used in optimizing electronic structure calculations.
 * - Integrates with the Solid and cElectron_Operators classes to perform operations related to electronic structure and energy calculations.
 * - Provides memory management functionalities for temporary storage during computational operations, including matrix transforms and singular value decomposition (SVD).
 * - Implements methods for starting, updating, transporting, and finalizing state data in quantum chemical simulations.
 * - Utilizes specialized grid operations from the Cneb class to efficiently manage calculations involving matrices and other high-dimensional data structures.
 * 
 * Dependencies:
 * - Solid.hpp: Interface with material properties and electron operations.
 * - cElectron.hpp: Provides electron-related functionalities used in quantum chemical calculations.
 * - Cneb.hpp: Supplies grid operations essential for matrix and vector calculations within the band_Geodesic workflow.
 * 
 * The band_Geodesic class facilitates efficient and modular handling of Grassmann manifold operations necessary for PWDFT calculations, supporting both memory management and numerical operations crucial to high-fidelity simulation and analysis.
 * 
 * NOTE: Ensure all allocated resources are carefully managed to prevent memory leaks and optimize computational efficiency.
 */

#pragma once

#include "cElectron.hpp"
#include "Solid.hpp"
#include "Cneb.hpp"
#include <cmath>

//#include	"util.hpp"

namespace pwdft {

//This is a Grassmann manifold
class band_Geodesic {

  int minimizer;
  Solid *mysolid;
  cElectron_Operators *myelectron;

  double *U, *Vt, *S;

  // tmp space for ffm multiplies
  double *tmp1, *tmp2, *tmp3, *tmpC, *tmpS;

public:
  Cneb *mygrid;

  /* Constructors */
  band_Geodesic(int minimizer0, Solid *mysolid0) 
  {
     mysolid = mysolid0;
     minimizer = minimizer0;
     myelectron = mysolid->myelectron;
     mygrid = mysolid->mygrid;
     U = mygrid->g_allocate_nbrillq_all();
     //Vt = mygrid->w_allocate(-1, 1);
     Vt = mygrid->w_allocate_nbrillq_all();
     S = new double[mygrid->nbrillq*(mygrid->ne[0] + mygrid->ne[1])];
    
     // tmp space
     tmp1 = mygrid->w_allocate(0,1);
     tmp2 = mygrid->w_allocate(0,1);
     tmp3 = mygrid->w_allocate(0,1);
     tmpC = new double[mygrid->ne[0] + mygrid->ne[1]];
     tmpS = new double[mygrid->ne[0] + mygrid->ne[1]];
  }

  /* destructor */
  ~band_Geodesic() {
    delete[] tmpS;
    delete[] tmpC;
    delete[] tmp3;
    delete[] tmp2;
    delete[] tmp1;
    delete[] S;
    delete[] Vt;
    delete[] U;
  }

   double start(double *A, double *max_sigma, double *min_sigma) 
   {
      double *V = mygrid->w_allocate_nbrillq_all();
      mygrid->ggw_SVD(A, U, S, V);

      int npack1 = mygrid->CGrid::npack1_max(); 
      int npack2 = 2*mygrid->CGrid::npack1_max();
      int shift2 = (mygrid->neq[0]+mygrid->neq[1])*npack2;
      int shift1 = 2*(mygrid->ne[0]*mygrid->ne[0] + mygrid->ne[1]*mygrid->ne[1]);
      int neall = mygrid->ne[0] + mygrid->ne[1];

    
      for (auto nbq=0; nbq<mygrid->nbrillq; ++nbq)
      {
         double mmsig = 9.99e9;
         double msig = 0.0;

         double *Vk  = V  + nbq*shift1;
         double *Vtk = Vt + nbq*shift1;

         for (int i=0; i<neall; ++i) 
         {
            if (std::fabs(S[i]) > msig)
               msig = fabs(S[i]);
            if (std::fabs(S[i]) < mmsig)
               mmsig = fabs(S[i]);
         }
         *max_sigma = msig;
         *min_sigma = mmsig;
       
         /* calculate Vt */
         mygrid->Cneb::ww_hermit_transpose(-1, Vk, Vtk);
       
         // double *tmp1 = mygrid->m_allocate(-1,1);
         // mygrid->mmm_Multiply(-1,Vt,V,1.0,tmp1,0.0);
         // util_matprint("Vt*V",4,tmp1);
       
         // mygrid->ggm_sym_Multiply(U,U,tmp1);
         // util_matprint("Ut*U",4,tmp1);
         // delete [] tmp1;
      }
    
      delete[] V;
    
      /* calculate  and return 2*<A|H|psi> */
      return (2.0 * myelectron->eorbit(A));
   }


   void get(double t, double *Yold, double *Ynew) 
   {
      int npack1 =   mygrid->CGrid::npack1_max(); 
      int npack2 = 2*mygrid->CGrid::npack1_max();
      int shift2 = (mygrid->neq[0]+mygrid->neq[1])*npack2;
      int shift1 = 2*(mygrid->ne[0]*mygrid->ne[0]+mygrid->ne[1]*mygrid->ne[1]);
      int neall = mygrid->ne[0] + mygrid->ne[1];

      double rone[2]  = {1.0,0.0};
      double rmone[2]  = {-1.0,0.0};
      double rzero[2] = {0.0,0.0};

      for (auto nbq=0; nbq<(mygrid->nbrillq); ++nbq)
      {
         double *Sk = S+nbq*neall;
         double *Vtk  = Vt + nbq*shift1;
         double *Uk   = U  + nbq*shift2;
         double *Yoldk = Yold  + nbq*shift2;
         double *Ynewk = Ynew  + nbq*shift2;

         mygrid->ww_SCtimesVtrans(-1, t, Sk, Vtk, tmp1, tmp3, tmpC, tmpS);

       
         /* Ynew = Yold*V*cos(Sigma*t)*Vt + U*sin(Sigma*t)*Vt */
         mygrid->www_Multiply2(-1, Vtk, tmp1, rone, tmp2, rzero);

         mygrid->fwf_Multiply(-1, Yoldk, tmp2, rone, Ynewk, rzero);
         mygrid->fwf_Multiply(-1, Uk, tmp3, rone, Ynewk, rone);
       
         /* ortho check  - need to figure out what causes this to happen */
         double sum2 = mygrid->gg_traceall(Ynew, Ynew);
         double sum1 = mygrid->ne[0] + mygrid->ne[1];
         if ((mygrid->ispin) == 1)
            sum1 *= 2;
         if (std::fabs(sum2 - sum1) > 1.0e-10) 
         {
           // if (myparall->is_master()) std::cout << " Warning - Gram-Schmidt being
           // performed on psi2" << std::endl; std::cout << " Warning - Gram-Schmidt
           // being performed on psi2, t="
           //           << t << " sum1=" << sum1 << " sum2=" << sum2 <<  " error=" <<
           //           fabs(sum2-sum1) << std::endl;
            mygrid->g_ortho(Ynew);
         }
      }
   }

   void transport(double t, double *Yold, double *Ynew) 
   {
      double rone[2]  = {1.0,0.0};
      double rmone[2]  = {-1.0,0.0};
      double rzero[2] = {0.0,0.0};

      int npack1 =   mygrid->CGrid::npack1_max(); 
      int npack2 = 2*mygrid->CGrid::npack1_max();
      int shift2 = (mygrid->neq[0]+mygrid->neq[1])*npack2;
      int shift1 = 2*(mygrid->ne[0]*mygrid->ne[0]+mygrid->ne[1]*mygrid->ne[1]);
      int neall = mygrid->ne[0] + mygrid->ne[1];

      for (auto nbq=0; nbq<mygrid->nbrillq; ++nbq)
      {
          //std::cout << "transport nbq=" << nbq << std::endl;

         double *Sk = S+nbq*neall;
         double *Vtk  = Vt  + nbq*shift1;
         double *Uk   = U  + nbq*shift2;
         double *Yoldk = Yold  + nbq*shift2;
         double *Ynewk = Ynew  + nbq*shift2;

         //std::cout << "transport2 nbq=" << nbq << std::endl;
         mygrid->ww_SCtimesVtrans2(-1, t, Sk, Vtk, tmp1, tmp3, tmpC, tmpS);
     
         /* tHnew = (-Yold*V*sin(Sigma*t) + U*cos(Sigma*t))*Sigma*Vt */
         //std::cout << "transport3 nbq=" << nbq << std::endl;
         mygrid->www_Multiply2(-1, Vtk, tmp1, rone, tmp2, rzero);
         //std::cout << "transport4 nbq=" << nbq << std::endl;
         //std::cout << "tmp2=" ;
         //for (auto k=0; k<(2*mygrid->ne[0]*2*mygrid->ne[0]); ++k)
         //    std::cout << tmp2[k] << " " ;
         //std::cout << std::endl;

         mygrid->fwf_Multiply(-1, Yoldk, tmp2, rmone, Ynewk, rzero);
         //std::cout << "transport5 nbq=" << nbq << std::endl;
         mygrid->fwf_Multiply(-1, Uk, tmp3, rone, Ynewk, rone);
         //std::cout << "transport6 nbq=" << nbq << std::endl;
      }
   }

   void psi_1transport(double t, double *H0) 
   {
      this->transport(t, mysolid->psi1, H0);
   }

   void Gtransport(double t, double *Yold, double *tG) 
   {
      double rone[2]  = {1.0,0.0};
      double rmone[2]  = {-1.0,0.0};
      double rzero[2] = {0.0,0.0};

      int npack1 =   mygrid->CGrid::npack1_max(); 
      int npack2 = 2*mygrid->CGrid::npack1_max();
      int shift2 = (mygrid->neq[0]+mygrid->neq[1])*npack2;
      int shift1 = 2*(mygrid->ne[0]*mygrid->ne[0]+mygrid->ne[1]*mygrid->ne[1]);
      int neall = mygrid->ne[0] + mygrid->ne[1];

      for (auto nbq=0; nbq<mygrid->nbrillq; ++nbq)
      {
         double *Sk = S+nbq*neall;
         double *Vtk  = Vt  + nbq*shift1;
         double *Uk  = U  + nbq*shift2;
         double *tGk = tG + nbq*shift2;
         double *Yoldk = Yold  + nbq*shift2;

         mygrid->ffw_sym_Multiply(-1,Uk,tGk,tmp2);
         mygrid->ffw_Multiply(-1, Uk, tGk, tmp2);
         mygrid->ww_SCtimesVtrans3(-1, t, Sk, tmp2, tmp1, tmp3, tmpC, tmpS);
         mygrid->www_Multiply2(-1, Vtk, tmp1, rone, tmp2, rzero);
        
         mygrid->fwf_Multiply(-1, Yoldk, tmp2, rmone, tGk, rone);
         mygrid->fwf_Multiply(-1, Uk, tmp3, rmone, tGk, rone);
      }
   }

   void psi_1Gtransport(double t, double *H0) 
   {
      this->Gtransport(t, mysolid->psi1, H0);
   }

   double energy(double t) 
   {
     this->get(t, mysolid->psi1, mysolid->psi2);
     return (mysolid->psi2_energy());
   }

   double energy0(double t) 
   {
     this->get(t, mysolid->psi1, mysolid->psi2);
     return (mysolid->psi2_energy0());
   }

   double denergy(double t) 
   {
      this->transport(t, mysolid->psi1, mysolid->psi2);
      return (2.0 * mysolid->psi2_eorbit());
   }

  void psi_final(double t) { this->get(t, mysolid->psi1, mysolid->psi2); }
};

} // namespace pwdft

#endif

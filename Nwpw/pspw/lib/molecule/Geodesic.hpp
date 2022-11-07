#ifndef	_GEODESIC_HPP_
#define _GEODESIC_HPP_

#include	<cmath>
#include	"Pneb.hpp"
#include	"Electron.hpp"
#include	"Molecule.hpp"

//#include	"util.hpp"

namespace pwdft {


class	Geodesic {

   int minimizer;
   Molecule            *mymolecule;
   Electron_Operators  *myelectron;

   double *U, *Vt, *S;
   
   // tmp space for ffm multiplies
   double *tmp1,*tmp2,*tmp3,*tmpC,*tmpS;

public:
   Pneb                *mygrid;

   /* Constructors */
   Geodesic(int minimizer0, Molecule *mymolecule0) {
      mymolecule = mymolecule0;
      minimizer  = minimizer0;
      myelectron = mymolecule->myelectron;
      mygrid     = mymolecule->mygrid;
      U  = mygrid->g_allocate(1);
      Vt = mygrid->m_allocate(-1,1);
      S = new double[mygrid->ne[0]+mygrid->ne[1]];

      // tmp space 
      tmp1 = mygrid->m_allocate(-1,1);
      tmp2 = mygrid->m_allocate(-1,1);
      tmp3 = mygrid->m_allocate(-1,1);
      tmpC = new double[mygrid->ne[0]+mygrid->ne[1]];
      tmpS = new double[mygrid->ne[0]+mygrid->ne[1]];
   }


   /* destructor */
   ~Geodesic() {
         delete [] tmpS;
         delete [] tmpC;
         delete [] tmp3;
         delete [] tmp2;
         delete [] tmp1;
         delete [] S;
         delete [] Vt;
         delete [] U;
    }


   double start(double *A, double *max_sigma, double *min_sigma) {
      double *V = mygrid->m_allocate(-1,1);
      mygrid->ggm_SVD(A,U,S,V);

      int neall = mygrid->ne[0] + mygrid->ne[1];
      double mmsig=9.99e9;
      double msig =0.0;
      for (int i=0; i<neall; ++i) {
         if (std::fabs(S[i])>msig)   msig = fabs(S[i]);
         if (std::fabs(S[i])<mmsig) mmsig = fabs(S[i]);
      }
      *max_sigma = msig;
      *min_sigma = mmsig;

      /* calculate Vt */
      mygrid->mm_transpose(-1,V,Vt);

      //double *tmp1 = mygrid->m_allocate(-1,1);
      //mygrid->mmm_Multiply(-1,Vt,V,1.0,tmp1,0.0);
      //util_matprint("Vt*V",4,tmp1);

      //mygrid->ggm_sym_Multiply(U,U,tmp1);
      //util_matprint("Ut*U",4,tmp1);
      //delete [] tmp1;

      delete [] V;

      /* calculate  and return 2*<A|H|psi> */
      return(2.0*myelectron->eorbit(A));

    }
 
    void get(double t, double *Yold, double *Ynew) {
       mygrid->mm_SCtimesVtrans(-1,t,S,Vt,tmp1,tmp3,tmpC,tmpS);

       /* Ynew = Yold*V*cos(Sigma*t)*Vt + U*sin(Sigma*t)*Vt */
       mygrid->mmm_Multiply2(-1,Vt,tmp1,1.0,tmp2,0.0);
       mygrid->fmf_Multiply(-1,Yold,tmp2,1.0,Ynew,0.0);
       mygrid->fmf_Multiply(-1,U,tmp3,1.0,Ynew,1.0);


       /* ortho check  - need to figure out what causes this to happen */
       double sum2  = mygrid->gg_traceall(Ynew,Ynew);
       double sum1  = mygrid->ne[0] + mygrid->ne[1];
       if ((mygrid->ispin)==1) sum1 *= 2;
       if (std::fabs(sum2-sum1)>1.0e-10)
       {
          //if (myparall->is_master()) std::cout << " Warning - Gram-Schmidt being performed on psi2" << std::endl;
          //std::cout << " Warning - Gram-Schmidt being performed on psi2, t=" 
          //          << t << " sum1=" << sum1 << " sum2=" << sum2 <<  " error=" <<  fabs(sum2-sum1) << std::endl;
          mygrid->g_ortho(Ynew);
       }
    }


    void transport(double t, double *Yold, double *Ynew) {
       mygrid->mm_SCtimesVtrans2(-1,t,S,Vt,tmp1,tmp3,tmpC,tmpS);

       /* tHnew = (-Yold*V*sin(Sigma*t) + U*cos(Sigma*t))*Sigma*Vt */
       mygrid->mmm_Multiply2(-1,Vt,tmp1,1.0,tmp2,0.0);
       mygrid->fmf_Multiply(-1,Yold,tmp2,-1.0,Ynew,0.0);
       mygrid->fmf_Multiply(-1,U,tmp3,1.0,Ynew,1.0);
    }

    void psi_1transport(double t, double *H0) {
        this->transport(t, mymolecule->psi1,H0);
    }

    void Gtransport(double t, double *Yold, double *tG) {
       //mygrid->ffm_sym_Multiply(-1,U,tG,tmp2);
       mygrid->ffm_Multiply(-1,U,tG,tmp2);
       mygrid->mm_SCtimesVtrans3(-1,t,S,tmp2,tmp1,tmp3,tmpC,tmpS);
       mygrid->mmm_Multiply2(-1,Vt,tmp1,1.0,tmp2,0.0);

       mygrid->fmf_Multiply(-1,Yold,tmp2,-1.0,tG,1.0);
       mygrid->fmf_Multiply(-1,U,tmp3,-1.0,tG,1.0);
    }

    void psi_1Gtransport(double t, double *H0) {
        this->Gtransport(t, mymolecule->psi1,H0);
    }

    double energy(double t) {
       this->get(t,mymolecule->psi1,mymolecule->psi2);
       return(mymolecule->psi2_energy());
    }

    double denergy(double t) {
       this->transport(t,mymolecule->psi1,mymolecule->psi2);
       return(2.0*mymolecule->psi2_eorbit());
    }

    void psi_final(double t) {
       this->get(t,mymolecule->psi1,mymolecule->psi2);
    }
};

}

#endif

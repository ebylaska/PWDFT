

#include "cKinetic.hpp"
#include "CGrid.hpp"
#include "CPseudopotential.hpp"

namespace pwdft {

/*************************************
 *                                   *
 *             cpsi_H                *
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

void cpsi_H(Cneb *mygrid, cKinetic_Operator *myke, CPseudopotential *mypsp,
            double *psi, double *psi_r, double *vl, double *vc, double *xcp,
            double *Hpsi, bool move, double *fion)

{
   int indx1 = 0;
   int indx2 = 0;
   int indx1n = 0;
   int indx2n = 0;
   int ispin = mygrid->ispin;
   int shift1 = 2*(mygrid->npack1_max());
   int shift2 = (mygrid->n2ft3d);
   int n2ft3d = (mygrid->n2ft3d);
   int nfft3d = (mygrid->nfft3d);
   int n1 = mygrid->neq[0];
   int n2 = mygrid->neq[0] + mygrid->neq[1];
   int nn = n2*mygrid->nbrillq;
 
   bool done = false;
 
   double omega = mygrid->lattice->omega();
   double scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));
   double scal2 = 1.0 / omega;
 
   /* allocate temporary memory */
   double *vall = mygrid->c_alloc();
   double *vpsi = mygrid->c_alloc();
   double *tmp  = mygrid->c_alloc();

 
   /* apply k-space operators */
   myke->ke(psi, Hpsi);

   /* apply non-local PSP  - Expensive */
   mypsp->v_nonlocal_fion(psi, Hpsi, move, fion);

   /* apply r-space operators  - Expensive*/
   mygrid->cc_pack_SMul(0,scal2,vl,vall);
   mygrid->cc_pack_Sum2(0,vc,vall);

   mygrid->c_unpack(0,vall);
   mygrid->cr_pfft3b(0,vall);


   for (auto nbq=0; nbq<mygrid->nbrillq; ++nbq)
   {
     for (auto ms=0; ms<ispin; ++ms)
     {
        mygrid->rcc_Sum(xcp+ms*nfft3d,vall,tmp);
        for (int i=0; i<(mygrid->neq[ms]); ++i)
        {
           std::memcpy(vpsi,tmp,n2ft3d*sizeof(double));
           mygrid->bb_Mul(psi_r+indx1n,vpsi);
           mygrid->rc_fft3d(vpsi);
           mygrid->c_pack(1+nbq,vpsi);
           mygrid->cc_pack_daxpy(1+nbq,(-scal1),vpsi,Hpsi+indx1);
 
           indx1n += shift2;
           indx1  += shift1;
        }
     }
   }

/*
   { nwpw_timing_function ftimer(1);
 
     int ms = 0;

     while (!done) 
     {
        if (indx1<nn) 
        {
           int nbq1 = (indx1/n2) + 1;

           if ((indx1 - (nbq1-1)*n2) < n1)
              ms = 0;
           else 
              ms = 1;

           mygrid->rcc_Sum(xcp+ms*nfft3d,vall,tmp);
           std::memcpy(vpsi,tmp,n2ft3d*sizeof(double));

           mygrid->bb_Mul(psi_r+indx1n,vpsi);

           mygrid->rc_pfft3f_queuein(nbq1,1,vpsi);
           indx1n += shift2;
           ++indx1;
        }
        
        if ((mygrid->rc_pfft3f_queuefilled()) || (indx1 >= nn)) 
        {
           int nbq2 = (indx2/n2) + 1;
           mygrid->rc_pfft3f_queueout(nbq2,1,vpsi);

           mygrid->cc_pack_daxpy(nbq2,(-scal1),vpsi,Hpsi+indx2n);

           indx2n += shift1;
           ++indx2;
        }
        done = ((indx1 >= nn) && (indx2 >= nn));
     }
   }
   */
   /*{
      double hml[4*22*22*4];
      mygrid->ggw_sym_Multiply(psi, Hpsi, hml);
      std::cout << "C hml=" << hml[0]  << " " << hml[1] << " " 
                            << hml[46] << " " << hml[47] << " " 
                            << hml[92] << " " << hml[93] <<  std::endl;
   }
   */
   
   /* deallocate temporary memory */
   mygrid->r_dealloc(tmp);
   mygrid->r_dealloc(vpsi);
   mygrid->r_dealloc(vall);
}


} // namespace pwdft

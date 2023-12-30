

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
   int shift1 = 2 * (mygrid->npack(1));
   int shift2 = (mygrid->n2ft3d);
   int n2ft3d = (mygrid->n2ft3d);
   int ms = 0;
   int n1 = mygrid->neq[0];
   int n2 = mygrid->neq[0] + mygrid->neq[1];
 
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

   std::cout << "ke_psi = " ;
   for (auto i=0; i<20; ++i)
      std::cout << Hpsi[i] << " ";
   std::cout << std::endl << std::endl;
 
   /* apply non-local PSP  - Expensive */
   mypsp->v_nonlocal_fion(psi, Hpsi, move, fion);

   std::cout << "vnl_psi = " ;
   for (auto i=0; i<20; ++i)
      std::cout << Hpsi[i] << " ";
   std::cout << std::endl << std::endl;
 
   /* apply r-space operators  - Expensive*/
   mygrid->cc_pack_SMul(0, scal2, vl, vall);
   mygrid->cc_pack_Sum2(0,vc,vall);
   mygrid->c_unpack(0,vall);
   mygrid->cr_fft3d(vall);
 
   /* add v_field to vall */
   //if (mypsp->myefield->efield_on)
   //   mygrid->rr_Sum(mypsp->myefield->v_field,vall);
 
   
   /*
    std::cout << "INTO psiH rc_fftd  -------------------" << std::endl;
 
     //OLD - stable?
     for (ms=0; ms<ispin; ++ms)
     {
        mygrid->rrr_Sum(vall,xcp+ms*n2ft3d,tmp);
        for (int i=0; i<(mygrid->neq[ms]); ++i)
        {
            std::cout << "psiH i=" << i << " -------------------" << std::endl;
           mygrid->rrr_Mul(tmp,psi_r+indx2,vpsi);
           mygrid->rc_fft3d(vpsi);
           mygrid->c_pack(1,vpsi);
           mygrid->cc_pack_daxpy(1,(-scal1),vpsi,Hpsi+indx1);
 
           indx1 += shift1;
           indx2 += shift2;
        }
     }
    std::cout << "OUT psiH rc_fftd  -------------------" << std::endl;
    */
  
  
   { nwpw_timing_function ftimer(1);
 
     mygrid->rrr_Sum(vall,xcp,tmp);
     while (!done) 
     {
        if (indx1<n2) 
        {
           if (indx1>=n1) 
           {
              ms = 1;
              mygrid->rrr_Sum(vall,xcp+ms*n2ft3d,tmp);
           }
           
           mygrid->rrr_Mul(tmp,psi_r+indx1n,vpsi);
           
           mygrid->rc_pfft3f_queuein(1,vpsi);
           indx1n += shift2;
           ++indx1;
        }
        
        if ((mygrid->rc_pfft3f_queuefilled()) || (indx1 >= n2)) 
        {
           mygrid->rc_pfft3f_queueout(1,vpsi);
           mygrid->cc_pack_daxpy(1,(-scal1),vpsi,Hpsi+indx2n);
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


} // namespace pwdft

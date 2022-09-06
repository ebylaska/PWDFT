

#include        "PGrid.hpp"
#include        "Kinetic.hpp"
#include        "Pseudopotential.hpp"

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

void psi_H(Pneb *mygrid, 
           Kinetic_Operator *myke,
           Pseudopotential  *mypsp,
           double *psi, double *psi_r,
           double *vl, double *vc, double *xcp,
           double *Hpsi,
           bool move, double *fion)

{
   int indx1 = 0;
   int indx2 = 0;
   int ispin = mygrid->ispin;
   int shift1 = 2*(mygrid->npack(1));
   int shift2 = (mygrid->n2ft3d);
   int n2ft3d = (mygrid->n2ft3d);

   double omega = mygrid->lattice->omega();
   double scal1 = 1.0/((double) ((mygrid->nx)*(mygrid->ny)*(mygrid->nz)));
   double scal2 = 1.0/omega;

   /* allocate temporary memory */
   double *vall= mygrid->r_alloc();
   double *vpsi= mygrid->r_alloc();
   double *tmp = mygrid->r_alloc();

   /* apply k-space operators */
   myke->ke(psi,Hpsi);

   /* apply non-local PSP  - Expensive */
   mypsp->v_nonlocal_fion(psi,Hpsi,move,fion);

   /* apply r-space operators  - Expensive*/
   mygrid->cc_pack_SMul(0,scal2,vl,vall);
   mygrid->cc_pack_Sum2(0,vc,vall);
   mygrid->c_unpack(0,vall);
   mygrid->cr_fft3d(vall);
   for (int ms=0; ms<ispin; ++ms)
   {
      mygrid->rrr_Sum(vall,&xcp[ms*n2ft3d],tmp);
      for (int i=0; i<(mygrid->neq[ms]); ++i)
      {
         mygrid->rrr_Mul(tmp,&psi_r[indx2],vpsi);
         mygrid->rc_fft3d(vpsi);
         mygrid->c_pack(1,vpsi);
         mygrid->cc_pack_daxpy(1,(-scal1),vpsi,&Hpsi[indx1]);

         indx1 += shift1;
         indx2 += shift2;
      }
   }

   /* deallocate temporary memory */
   mygrid->r_dealloc(tmp);
   mygrid->r_dealloc(vpsi);
   mygrid->r_dealloc(vall);
}

}

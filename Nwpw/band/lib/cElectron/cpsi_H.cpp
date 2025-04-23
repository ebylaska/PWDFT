

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
            double *Hpsi, bool move, double *fion, double *occ)

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
   //mypsp->v_nonlocal_fion(psi, Hpsi, move, fion);
   mypsp->v_nonlocal_fion(psi, Hpsi, move, fion, occ);

   /* apply r-space operators  - Expensive*/
   mygrid->cc_pack_SMul(0,scal2,vl,vall);
   mygrid->cc_pack_Sum2(0,vc,vall);

   mygrid->c_unpack(0,vall);
   mygrid->cr_pfft3b(0,vall);


/*
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
*/

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
   
   /*{
      double hml[4*22*22*4];
      mygrid->ggw_sym_Multiply(psi, Hpsi, hml);
      std::cout << "C hml=" << hml[0]  << " " << hml[1] << " " 
                            << hml[46] << " " << hml[47] << " " 
                            << hml[92] << " " << hml[93] <<  std::endl;
   }
   */

   /* deallocate temporary memory */
   mygrid->c_dealloc(tmp);
   mygrid->c_dealloc(vpsi);
   mygrid->c_dealloc(vall);
}

/*************************************
 *                                   *
 *             cpsi_H_orb            *
 *                                   *
 *************************************/
/**
 * @brief Applies the Hamiltonian operator H to a wavefunction orbital.
 *
 * Computes the action of the Hamiltonian on a given orbital by combining:
 *   - Kinetic energy in reciprocal space
 *   - Non-local pseudopotential (NLPP) in reciprocal space
 *   - Local potential in real space (FFT-transformed and accumulated)
 *
 * The orbital (`orb`) and result (`Horb`) are in packed reciprocal-space
 * format and should be allocated using `mygrid->c_pack_allocate(nbq1)`.
 * The real-space buffers (`orb_r`, `vall_r`, `vpsi`) are full complex
 * grids allocated with `mygrid->c_alloc()` and must be of size `2 * nfft3d`.
 *
 * @param[in]  nbq1     Brillouin zone point index
 * @param[in]  mygrid   Grid manager (FFT, memory, packing)
 * @param[in]  myke     Kinetic energy operator
 * @param[in]  mypsp    Pseudopotential operator (NLPP handler)
 * @param[in]  orb      Input orbital (complex, packed reciprocal space)
 * @param[in]  orb_r    Orbital transformed to real space (complex grid)
 * @param[in]  vall_r   Local potential in real space (complex grid)
 * @param[in,out] Horb  Output: H * orb (packed format, overwritten)
 *
 * @note `orb` and `Horb` must be allocated via `c_pack_allocate(nbq1)`.
 * @note `orb_r`, `vall_r`, and `vpsi` must be allocated via `c_alloc()`
 *       and have size `2 * nfft3d`.
 */
void cpsi_H_orb(const int nbq1, 
                Cneb *mygrid, cKinetic_Operator *myke, CPseudopotential *mypsp,
                double *orb, double *orb_r, double *vall_r, double *Horb)
{
   int n2ft3d = (mygrid->n2ft3d);
   double scal1 = 1.0 / ((double)((mygrid->nx) * (mygrid->ny) * (mygrid->nz)));

   /* allocate temporary memory */
   double *vpsi = mygrid->c_alloc();
 
   /* apply k-space operators */
   myke->ke_orb(nbq1,orb,Horb);
 
   /* apply non-local PSP  - Expensive */
   mypsp->v_nonlocal_orb(nbq1,orb, Horb);

   /* apply r-space operators  - Expensive*/
   // mygrid->ccc_Mul(vall_r,orb_r,vpsi);
   //std::memcpy(vpsi,vall_r,n2ft3d*sizeof(double));
   std::memcpy(vpsi,vall_r,(n2ft3d)*sizeof(double));
   mygrid->bb_Mul(orb_r,vpsi);

   mygrid->rc_fft3d(vpsi);
   mygrid->c_pack(nbq1,vpsi);
   mygrid->cc_pack_daxpy(nbq1,(-scal1),vpsi,Horb);


   /* deallocate temporary memory */
   mygrid->c_dealloc(vpsi);

}




} // namespace pwdft

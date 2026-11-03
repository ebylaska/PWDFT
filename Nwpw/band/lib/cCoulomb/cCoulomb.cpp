/* cCoulomb.cpp -
   Author - Eric Bylaska
*/

/**
 * class: cCoulomb_Operator
 *
 * brief: Class for computing the Coulomb operator in electronic structure calculations.
 *
 * The Coulomb operator is an integral operator that appears in many electronic
 * structure methods. This class provides methods to compute the Coulomb operator
 * and its associated energy, given the density in reciprocal space. The Coulomb
 * operator is defined as:
 *
 *   V(r) = ∫ d³r' ρ(r')/|r - r'|,
 *
 * where ρ(r') is the electron density at position r'.
 *
 * This class requires a Pneb object to define the grid on which the density is
 * represented, and a Control2 object to provide some options for the computation.
 *
 * author: Eric J. Bylaska
 * date: 5/9/2023
 */

/*
#include        <cstdio>
#include        <cstdlib>
#include        <iostream>
#include        <stdio.h>
#include	<string>
*/
#include "units.hpp"
#include "cCoulomb.hpp"
#include <cmath>

namespace pwdft {

/*******************************************
 *                                         *
 *   cCoulomb_Operator::cCoulomb_Operator  *
 *                                         *
 *******************************************/
cCoulomb_Operator::cCoulomb_Operator(Cneb *mygrid, Control2 &control) 
{
   int pzero, zero, taskid;
   double gg;
   double *Gx = mygrid->Gxyz(0);
   double *Gy = mygrid->Gxyz(1);
   double *Gz = mygrid->Gxyz(2);
   vg = new double[mygrid->npack(0)];
   double *tmp = new double[mygrid->nfft3d];
   double fourpi = 16.0 * atan(1.0);
   // double kk0 = 0.0; // where k0 = the uniform background potential
 
   mycneb = mygrid;
 
   taskid = mycneb->c3db::parall->taskid_i();
   pzero = mycneb->cijktop(0, 0, 0);
   zero = mycneb->cijktoindex(0, 0, 0);
 
   for (auto k=0; k<(mycneb->nfft3d); ++k) 
   {
      gg = Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k];
      if ((pzero == taskid) && (k == zero))
         tmp[k] = 0.0;
      else
         tmp[k] = fourpi/gg;
   }
   mycneb->t_pack(0, tmp);
   mycneb->tt_pack_copy(0, tmp, vg);

   delete[] tmp;
}

/*******************************************
 *                                         *
 *        cCoulomb_Operator::vcoulomb      *
 *                                         *
 *******************************************/
void cCoulomb_Operator::vcoulomb(const double *dng, double *vcout) 
{
   int npack0 = (mycneb->npack(0));
   int k1 = 0;
   for (auto k=0; k<npack0; ++k) 
   {
      vcout[k1]   = vg[k] * dng[k1];
      vcout[k1+1] = vg[k] * dng[k1+1];
      k1 += 2;
   }
}

/*******************************************
 *                                         *
 *        cCoulomb_Operator::ecoulomb      *
 *                                         *
 *******************************************/
double cCoulomb_Operator::ecoulomb(const double *dng) 
{
   int npack0 = (mycneb->npack(0));
   double ave = 0.0;

   int k1 = 0;
   for (auto k=0; k<npack0; ++k) 
   {
      ave += vg[k] * (dng[k1]*dng[k1] + dng[k1+1]*dng[k1+1]);
      k1  += 2;
   }
   ave = mycneb->c3db::parall->SumAll(1,ave);
   ave *= 0.5*(mycneb->lattice->omega());
 
   return ave;
}



/*******************************************
 *                                         *
 *        cCoulomb_Operator::euv           *
 *                                         *
 *******************************************/
/**
 * @brief Compute the Coulomb contribution to the stress tensor for
 *        the current density.
 *
 * Given the packed complex density @p dng and the reciprocal-space
 * metric of the lattice, this routine returns the 3×3 stress tensor
 * (row-major, flattened as @c stress[u+3*v]) arising from the
 * long-range Coulomb (Hartree) energy.  The result combines an
 * isotropic term proportional to the Coulomb energy @c ecoul with
 * an anisotropic reciprocal-space sum over the squared density.
 *
 * The computation proceeds in three stages:
 *
 * 1. The 3×3 reciprocal metric
 *    @f$h_{ij} = \frac{1}{2\pi}\,u_{ij}@f$ is formed from the
 *    lattice's unit reciprocal vectors, where @f$u_{ij}@f$ is
 *    @c lattice->unitg(i,j).  This @c hm array is stored at
 *    @c hm[i+3*j].
 *
 * 2. A per-G-vector weight
 *    @f$\tilde n_G^2 = (4\pi/G^2)^2\,|n_G|^2\,v_G^2@f$ is built
 *    by @c crr_pack_SqrMul2 into @c tmp2.
 *
 * 3. The tensor
 *    @f[
 *       B_{us} = \frac{\omega}{4\pi}
 *                \sum_G \tilde n_G^2 \, G_u \, G_s
 *    @f]
 *    is accumulated with @c rrr_pack_MulDot, exploiting
 *    @f$B_{us} = B_{su}@f$ so only the upper triangle is summed.
 *
 * The final stress is then
 * @f[
 *    \sigma_{uv} = -E_{\mathrm{coul}}\,h_{uv}
 *                  + \sum_s B_{us}\,h_{sv}
 * @f]
 * stored at @c stress[u+3*v].
 *
 * @param[in]  dng     Packed complex density amplitudes for the
 *                     current density, in the same layout expected
 *                     by @c crr_pack_SqrMul2 (interleaved real/imag).
 * @param[out] stress  Caller-provided array of nine doubles.  On
 *                     entry the contents are overwritten with zeros;
 *                     on exit it holds the 3×3 Coulomb stress tensor
 *                     in row-major order, @c stress[u+3*v].
 *
 * @note @b The @c hm array is the reciprocal metric, not the
 *       identity.  The isotropic term is @f$-E_{\mathrm{coul}}\,
 *       h_{uv}@f$, @e not @f$-E_{\mathrm{coul}}\,\delta_{uv}@f$.
 *       This is correct: every reciprocal-space @f$G@f$-vector
 *       component carries a factor of the reciprocal basis, so the
 *       strain derivative of the @f$1/G^2@f$ Coulomb kernel produces
 *       an @f$h_{uv}@f$ factor even in the term that looks
 *       isotropic.  The plain Hartree energy returned by
 *       @c ecoulomb() does @e not include this metric factor; the
 *       multiplication by @c hm here is what supplies it.  Readers
 *       comparing against textbook stress formulas should not
 *       "correct" this term to @f$\delta_{uv}@f$.
 * @note @p stress must point to at least nine writable doubles.
 * @note The reciprocal-space sums are over the full G-vector set of
 *       the plane-wave basis as stored in @c mycneb; the reduction
 *       is performed by @c rrr_pack_MulDot via @c SumAll, so the
 *       result is identical on every MPI rank.
 * @note The @c std::array<double,9> @c Bus is symmetric by
 *       construction; the lower triangle is filled by copying from
 *       the upper triangle after the accumulation loop.
 * @note The constant @c scal = 1/(2π) reflects the Fourier
 *       convention used throughout the Nwpw code.
 *
 * @see CGrid::crr_pack_SqrMul2
 * @see CGrid::rrr_pack_MulDot
 * @see cCoulomb_Operator::ecoulomb
 */
void cCoulomb_Operator::euv(const double *dng, double *stress)
{
   std::fill(stress, stress + 9, 0.0);

   constexpr double pi     = units::PI;
   constexpr double fourpi = 4.0*pi;
   constexpr double scal   = 1.0/(2.0*pi);

   int ksize    = mycneb->npack(0);
   double omega = mycneb->lattice->omega();
   double ss    = omega/fourpi;

   // tmp space
   std::vector<double> tmp2(ksize);

   // define hm
   double hm[9];
   for (size_t i=0; i<3; ++i)
   for (size_t j=0; j<3; ++j)
      hm[i+3*j] = scal*mycneb->lattice->unitg(i,j);

   // tmp2(G) = (n(G)**2) * (4*pi/G**2)**2
   mycneb->crr_pack_SqrMul2(0,dng,vg,tmp2.data());

   // Bus = Sum(G) (omega/4*pi)*tmp2(G)*Gu*Gs

   std::array<double, 9> Bus = {0.0}; // Fixed size symmetric tensor

   // An array of pointers to represent the 3 segments of G.
   // This is compatible with C++11/14/17 and avoids std::span errors.
   const double* g_segments[3];

   g_segments[0] = mycneb->Gpackxyz(0,0);                  // Start of segment 1 - Gx
   g_segments[1] = mycneb->Gpackxyz(0,1);                  // Start of segment 2 - Gy
   g_segments[2] = mycneb->Gpackxyz(0,2);                  // Start of segment 3 - Gz

   for (size_t u=0; u<3; ++u)
   for (size_t s=u; s<3; ++s)
   {
      // sum = Sum_k gu[k]*gs[k]*tmp2[k]
      double sum = mycneb->rrr_pack_MulDot(0,g_segments[u],g_segments[s],tmp2.data());
      Bus[u+3*s] = ss*sum;
   }
   for (size_t u=0;   u<3; ++u)
   for (size_t s=u+1; s<3; ++s)
      Bus[s+3*u] = Bus[u+3*s];

   double ecoul = this->ecoulomb(dng);
   for (size_t v=0; v<3; ++v)
   for (size_t u=0; u<3; ++u)
   {
      stress[u+3*v] = -ecoul*hm[u+3*v];
      for (size_t s=0; s<3; ++s)
         stress[u+3*v] += Bus[u+3*s]*hm[s+3*v];
   }

}






} // namespace pwdft

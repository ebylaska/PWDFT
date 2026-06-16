/* Kinetic.C -
   Author - Eric Bylaska
*/

/*


#include        <cmath>
#include        <cstdio>
#include        <cstdlib>
#include        <iostream>
#include        <stdio.h>
#include	<string>
*/

#include "units.hpp"
#include "Kinetic.hpp"
#include "PGrid.hpp"

namespace pwdft {

/*******************************************
 *                                         *
 *     Kinetic_Operator::Kinetic_Operator  *
 *                                         *
 *******************************************/
/**
 * @brief Constructor: Pre-calculates kinetic energy coefficients in reciprocal space.
 * 
 * This constructor initializes the operator by computing the kinetic energy 
 * term (-0.5 * |G|^2) for every point in the reciprocal lattice. It transforms 
 * the values from the 3D grid format into a packed 1D array (tg) for fast access.
 *
 * @param mygrid Pointer to the Pneb object containing the reciprocal lattice vectors.
 * 
 * @note This constructor allocates memory for 'tg' that must be managed by the class.
 */
Kinetic_Operator::Kinetic_Operator(Pneb *mygrid) : mypneb(mygrid)
{
   //mypneb = mygrid;
   if (!mypneb) 
   {
      // Stop the program before it tries to access a null pointer
      throw std::runtime_error("Kinetic_Operator: Passed Pneb pointer is null!");
   }

   size_t ksize = (mypneb->npack(1));

   // Accessing G-vectors from the grid
   double *Gx = mypneb->Gpackxyz(1,0);
   double *Gy = mypneb->Gpackxyz(1,1);
   double *Gz = mypneb->Gpackxyz(1,2);

   // Allocate the persistent kinetic coefficient array (tg)
   // This array lives as long as the Kinetic_Operator object lives.
   this->tg = new double[ksize];

   for (size_t k=0; k<ksize; ++k) 
   {
      double gg = Gx[k]*Gx[k] + Gy[k]*Gy[k] + Gz[k]*Gz[k];
      this->tg[k] = -0.5*gg;
   }
}


/*******************************************
 *                                         *
 *         Kinetic_Operator::ke            *
 *                                         *
 *******************************************/
/**
 * @brief Accumulates the kinetic energy gradient into the master gradient array.
 * 
 * This function performs an additive update to the `tpsi` array. It takes 
 * the kinetic-weighted components of the wavefunction (`psi`) and adds them 
 * to the existing values in `tpsi`. This allows multiple energy-component 
 * operators to contribute to a single global gradient buffer.
 *
 * @param[in]  psi    Pointer to the input wave function/basis coefficients.
 * @param[in,out] tpsi Pointer to the master gradient array (the accumulator).
 *                    This array is NOT cleared by this function; it expects 
 *                    pre-existing partial gradients which this function 
 *                    will add onto.
 * 
 * @note This function assumes that the `tpsi` buffer has been initialized 
 *       to zero at the start of the total energy/gradient calculation step.
 * @note Complexity: $O(N_{size} \times K_{size})$ where $N$ is number of shells 
 *       and $K$ is number of basis functions.
 */
void Kinetic_Operator::ke(double *psi, double *tpsi) 
{
   size_t nsize = (mypneb->neq[0] + mypneb->neq[1]);
   size_t ksize = (mypneb->npack(1));

   size_t k1 = 0;
   size_t k2 = 1;
   for (size_t n=0; n<nsize; ++n)
      for (size_t k=0; k<ksize; ++k) 
      {
         tpsi[k1] += tg[k]*psi[k1];
         tpsi[k2] += tg[k]*psi[k2];
         k1 += 2;
         k2 += 2;
      }
}

/*******************************************
 *                                         *
 *         Kinetic_Operator::ke_orb        *
 *                                         *
 *******************************************/
/**
 * @brief Applies the kinetic operator to a single orbital.
 * 
 * This function performs an in-place additive update to a single orbital 
 * coefficient array. It is useful for localized updates where only one 
 * specific orbital needs to be processed through the kinetic filter.
 *
 * @param[in]  orb  Pointer to the input orbital coefficients.
 * @param[in,out] torb Pointer to the output/accumulator array.
 * 
 * @note Complexity: $O(K_{size})$ where $K$ is the number of basis functions.
 */
void Kinetic_Operator::ke_orb(double *orb, double *torb) 
{
 
   size_t ksize = (mypneb->npack(1));

   size_t k1 = 0;
   size_t k2 = 1;
   for (size_t k=0; k<ksize; ++k) 
   {
      torb[k1] += tg[k] * orb[k1];
      torb[k2] += tg[k] * orb[k2];
      k1 += 2;
      k2 += 2;
   }
}


/*******************************************
 *                                         *
 *        Kinetic_Operator::ke_ave         *
 *                                         *
 *******************************************/
/**
 * @brief Calculates the expectation value of the kinetic energy.
 * 
 * This function computes the average kinetic energy contribution $\langle \hat{T} \rangle$ 
 * by integrating the kinetic operator over the wave function coefficients. 
 * It traverses the basis functions across all shells, applying occupation weights 
 * and accounting for both the "zero" (s-space) and "packed" (k-space) regions 
 * of the expansion.
 *
 * @param[in]  psi    Pointer to the array of wave function/basis coefficients.
 *                    Must be large enough to accommodate all indices $k$ up to $k_{size2}$.
 * @param[in]  occ    (Optional) Pointer to an array of occupation weights, 
 *                    one per index $n$. If nullptr, weight is assumed to be 1.0.
 * 
 * @return double      The calculated average kinetic energy (with a sign inversion 
 *                     applied as per the physical model convention).
 * 
 * @note The function handles both the "zero" (s-space) and "packed" (k-space) 
 *       regions of the summation, applying a factor of 2.0 for the packed region.
 * @note Complexity: $O(N_{size} \times K_{size})$ where $N$ is number of shells 
 *       and $K$ is number of basis functions.
 */
double Kinetic_Operator::ke_ave(double *psi, double *occ)
{
   size_t nsize  = (mypneb->neq[0] + mypneb->neq[1]);
   size_t ksize1 = (mypneb->nzero(1));
   size_t ksize2 = (mypneb->npack(1));

   double ave = 0.0;
   size_t k1 = 0;
   size_t k2 = 1;
   for (size_t n=0; n<nsize; ++n) 
   {
      double wght = (occ != nullptr) ? occ[n] : 1.0;

      for (size_t k=0; k<ksize1; ++k) 
      {
         ave += tg[k]*(psi[k1]*psi[k1] + psi[k2]*psi[k2])*wght;
         k1 += 2;
         k2 += 2;
      }
      for (size_t k=ksize1; k<ksize2; ++k) 
      {
         ave += 2.0*tg[k]*(psi[k1]*psi[k1] + psi[k2]*psi[k2])*wght;
         k1 += 2;
         k2 += 2;
      }
   }
   ave = mypneb->d3db::parall->SumAll(0, ave);
   if (mypneb->ispin == 1)
      ave *= 2.0;
   ave = -ave;

   return ave;
}


/***********************************************
 *                                             *
 *       Kinetic_Operator::ke_euv              *
 *                                             *
 ***********************************************/
/**
 * @brief Calculates the kinetic energy contribution to the lattice stress tensor.
 * 
 * This function computes the $A_{us}$ component of the kinetic energy via an Ewald-style 
 * summation over basis functions. It calculates the metric tensor ($H$) from the 
 * lattice parameters and performs a tensor contraction to derive the final 
 * stress tensor components in the laboratory frame.
 *
 * @param[in]  psi    Pointer to the array of wave function/basis coefficients.
 *                    Must be large enough to accommodate all indices $k$ up to $k_{size2}$.
 * @param[out] stress Pointer to a 9-element array representing the $3 \times 3$ 
 *                    stress tensor (initialized to zero by this function).
 * @param[in]  occ    (Optional) Pointer to an array of occupation weights, one per 
 *                    index $n$. If nullptr, weight is assumed to be 1.0.
 * 
 * @note The function handles both the "zero" (s-space) and "packed" (k-space) 
 *       regions of the summation, applying a factor of 2.0 for the packed region.
 * @note Complexity: $O(N_{size} \times K_{size})$ where $N$ is number of shells 
 *       and $K$ is number of basis functions.
 */
void Kinetic_Operator::ke_euv(double *psi, double *stress, double *occ) 
{
   //Initialize the master stress array to zero safely
   std::fill(stress, stress + 9, 0.0);

   double Aus[9] = {0.0};

   constexpr double pi     = units::PI;
   double scal = 1.0/(2.0*pi);
   double hm[9];
   for (size_t i=0; i<3; ++i)
   for (size_t j=0; j<3; ++j)
      hm[i+3*j] = scal*mypneb->lattice->unitg(i,j);

   size_t nsize  = (mypneb->neq[0] + mypneb->neq[1]);
   size_t ksize1 = (mypneb->nzero(1));
   size_t ksize2 = (mypneb->npack(1));

   for (size_t u=0; u<3; ++u)
   for (size_t s=u; s<3; ++s)
   {
      double *gu = mypneb->Gpackxyz(1,u);
      double *gs = mypneb->Gpackxyz(1,s);

      double ave = 0.0;
      size_t k1 = 0;
      size_t k2 = 1;
      for (size_t n=0; n<nsize; ++n) 
      {
         double wght = (occ != nullptr) ? occ[n] : 1.0;
 
         for (size_t k=0; k<ksize1; ++k) 
         {
            ave += wght*gu[k]*gs[k]*(psi[k1]*psi[k1] + psi[k2]*psi[k2]);
            k1 += 2;
            k2 += 2;
         }
         for (size_t k=ksize1; k<ksize2; ++k) 
         {
           ave += 2.0*wght*gu[k]*gs[k]*(psi[k1]*psi[k1] + psi[k2]*psi[k2]);
           k1 += 2;
           k2 += 2;
         }
      }

      ave = mypneb->d3db::parall->SumAll(0,ave);
      if (mypneb->ispin == 1) 
         ave *= 2.0;

      Aus[u+3*s] += ave;
   }
   
   for (size_t u=0; u<3; ++u)
   for (size_t s=u+1; s<3; ++s)
       Aus[s+3*u] = Aus[u+3*s];

   // Calculate stress = -Sum(s) hm(s,v)*Aus(u,s)
   for (size_t v=0; v<3; ++v)
   for (size_t u=0; u<3; ++u)
      for (size_t s=0; s<3; ++s)
         stress[u+3*v] -= Aus[u+3*s]*hm[s+3*v];

}




/***********************************************
 *                                             *
 *       Kinetic_Operator::ke_precondition     *
 *                                             *
 ***********************************************/
// **** My preconditioner ****
/**
 * @brief Performs frequency-dependent preconditioning for the kinetic energy operator.
 * 
 * This function implements a polynomial-based preconditioner designed to flatten 
 * the eigenvalue spectrum of the kinetic energy operator. It scales the components 
 * of the search direction (tpsi) based on their contribution to the kinetic energy 
 * expectation value <psi|T|psi>. By damping high-frequency G-vector components, 
 * it accelerates the convergence of iterative solvers (e.g., Davidson or CG).
 *
 * @param Ep      The current eigenvalue/energy (used for context/scaling).
 * @param neall   Total number of electronic bands to be processed.
 * @param psi     Pointer to the input wavefunction array (shape: neall * npack2).
 * @param tpsi    Pointer to the output preconditioned array (accumulated result).
 * 
 * @note Complexity: O(neall * npack1)
 * @warning Ensure that 'tpsi' is properly initialized before calling this function.
 */
void Kinetic_Operator::ke_precondition(const double Ep, const int neall, double *psi, double *tpsi)
{
   size_t npack1 = mypneb->npack(1);
   size_t npack2 = 2*npack1;

   // Use a local buffer on the stack if small, or better: 
   // pass in a pre-allocated 'tmp' from the caller to avoid 'new' overhead.
   std::vector<double> tmp_vec(npack2); 
   double *tmp = tmp_vec.data();

   for (size_t n=0; n<(size_t)neall; ++n)
   {
      double *worb = psi + n * npack2;
      
      // Compute T|psi>
      mypneb->tcc_pack_Mul(1, tg, worb, tmp);
      
      // Compute <psi|T|psi>
      double sum = mypneb->cc_pack_dot(1, tmp, worb);

      // Safety check to prevent NaN propagation
      if (std::abs(sum) < 1e-16) continue; 
      
      double inv_sum = 1.0 / sum;

      for (size_t k=0; k<npack1; ++k)    
      {
         // Use direct indexing: 2*k and 2*k + 1
         size_t idx1 = 2 * k;
         size_t idx2 = 2 * k + 1;

         double x = tg[k];
         // Magnitude of the wavefunction at this G-vector
         double density = (worb[idx1]*worb[idx1] + worb[idx2]*worb[idx2]);
         
         x = x * density * inv_sum;

         // Polynomial damping logic
         double cm = 27.00 + (18.00 + (12.00 + 8.00*x)*x)*x;
         double dm = (cm + 16.00 * x*x*x*x);
         
         cm = cm / dm;

         // Apply scaling to the output wavefunction
         tpsi[idx1] *= cm;
         tpsi[idx2] *= cm;
      }
   }
}


} // namespace pwdft

/* Coulomb.cpp -
   Author - Eric Bylaska
*/

/**
 * class: Coulomb_Operator
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
#include "Coulomb.hpp"
#include <cmath>

namespace pwdft {

/*******************************************
 *                                         *
 *     Coulomb_Operator::Coulomb_Operator  *
 *                                         *
 *******************************************/
Coulomb_Operator::Coulomb_Operator(Pneb *mygrid, Control2 &control) {
  int k, pzero, zero, taskid;
  double gg;
  double *Gx = mygrid->Gxyz(0);
  double *Gy = mygrid->Gxyz(1);
  double *Gz = mygrid->Gxyz(2);
  vg = new double[mygrid->npack(0)];
  double *tmp = new double[mygrid->nfft3d];
  double fourpi = 16.0 * atan(1.0);
  // double kk0 = 0.0; // where k0 = the uniform background potential

  mypneb = mygrid;

  taskid = mypneb->d3db::parall->taskid_i();
  pzero = mypneb->ijktop(0, 0, 0);
  zero = mypneb->ijktoindex(0, 0, 0);

  for (k = 0; k < (mypneb->nfft3d); ++k) {
    gg = Gx[k] * Gx[k] + Gy[k] * Gy[k] + Gz[k] * Gz[k];
    if ((pzero == taskid) && (k == zero))
      tmp[k] = 0.0;
    else
      tmp[k] = fourpi / gg;
  }
  mypneb->t_pack(0, tmp);
  mypneb->tt_pack_copy(0, tmp, vg);

  delete[] tmp;
}

/*******************************************
 *                                         *
 *        Coulomb_Operator::vcoulomb       *
 *                                         *
 *******************************************/
void Coulomb_Operator::vcoulomb(const double *dng, double *vcout) 
{
   int ksize = (mypneb->npack(0));
   int k1 = 0;
   for (auto k=0; k<ksize; ++k) 
   {
      vcout[k1]   = vg[k] * dng[k1];
      vcout[k1+1] = vg[k] * dng[k1+1];
      k1 += 2;
   }
}

/*******************************************
 *                                         *
 *        Coulomb_Operator::ecoulomb       *
 *                                         *
 *******************************************/
double Coulomb_Operator::ecoulomb(const double *dng) 
{
   int ksize1 = (mypneb->nzero(0));
   int ksize2 = (mypneb->npack(0));
   double ave = 0.0;

   int k1 = 0;
   for (auto k=0; k<ksize1; ++k) 
   {
      ave += vg[k] * (dng[k1]*dng[k1] + dng[k1+1]*dng[k1+1]);
      k1 += 2;
   }
   for (auto k=ksize1; k<ksize2; ++k) 
   {
      ave += 2.0 * vg[k] * (dng[k1]*dng[k1] + dng[k1+1]*dng[k1+1]);
      k1  += 2;
   }
   ave = mypneb->d3db::parall->SumAll(1,ave);
   ave *= 0.5*(mypneb->lattice->omega());
 
   return ave;
}

/*******************************************
 *                                         *
 *        Coulomb_Operator::euv            *
 *                                         *
 *******************************************/
/**
 * @brief Computates the electrostatic (Coulomb) contribution to the stress tensor.
 * 
 * This function calculates the derivative of the Coulomb energy with respect to the 
 * strain tensor. It utilizes the reciprocal lattice vectors and the Poisson kernel 
 * in G-space to construct the symmetric part of the stress tensor.
 *
 * Physics Implementation:
 *   The stress tensor contribution is derived from the derivative of the 
 *   Coulomb energy (E_c) with respect to the lattice strain. 
 *   Specifically, it implements the summation:
 *     stress_{uv} = -E_c * hm_{uv} + \sum_{G} [ (omega/4pi) * (4pi/G^2)^2 * n(G)^2 * (G_u G_s) * hm_{sv} ]
 *   where 'hm' represents the scaled reciprocal lattice vectors.
 *
 * @param[in]  dng     Pointer to the array of density coefficients (n(G)) in reciprocal space.
 * @param[out] stress  Pointer to a 9-element array representing the 3x3 stress tensor 
 *                     (stored in column-major order: index = u + 3*v).
 * 
 * @note The calculation involves computing the G-space kernel (4*pi/G^2) and 
 *       summing the contributions of the reciprocal lattice vectors.
 * @complexity O(N_G), where N_G is the number of G-vectors in the expansion.
 */
void Coulomb_Operator::euv(const double *dng, double *stress) 
{
   std::fill(stress, stress + 9, 0.0);

   constexpr double pi     = units::PI;
   constexpr double fourpi = 4.0*pi;
   constexpr double scal   = 1.0/(2.0*pi);

   int ksize    = mypneb->npack(0);
   double omega = mypneb->lattice->omega();
   double ss    = omega/fourpi;

   // tmp space
   std::vector<double> tmp2(ksize);

   // define hm
   double hm[9];
   for (size_t i=0; i<3; ++i)
   for (size_t j=0; j<3; ++j)
      hm[i+3*j] = scal*mypneb->lattice->unitg(i,j);


   // tmp2(G) = (n(G)**2) * (4*pi/G**2)**2  
   mypneb->ctt_pack_SqrMul2(0,dng,vg,tmp2.data());


   // Bus = Sum(G) (omega/4*pi)*tmp2(G)*Gu*Gs 

   std::array<double, 9> Bus = {0.0}; // Fixed size symmetric tensor

   // An array of pointers to represent the 3 segments of G.
   // This is compatible with C++11/14/17 and avoids std::span errors.
   const double* g_segments[3];

   g_segments[0] = mypneb->Gpackxyz(0,0);                  // Start of segment 1 - Gx
   g_segments[1] = mypneb->Gpackxyz(0,1);                  // Start of segment 2 - Gy
   g_segments[2] = mypneb->Gpackxyz(0,2);                  // Start of segment 3 - Gz


   for (size_t u=0; u<3; ++u)
   for (size_t s=u; s<3; ++s)
   {
      // sum = Sum_k gu[k]*gs[k]*tmp2[k]
      double sum = mypneb->ttt_pack_MulDot(0,g_segments[u],g_segments[s],tmp2.data());
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

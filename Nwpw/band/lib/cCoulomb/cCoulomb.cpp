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
   std::cout << "vg=" << vg[0] << " " 
                      << vg[1] << " " 
                      << vg[2] << " " 
                      << vg[3] << " " 
                      << vg[4] << " " 
                      << vg[5] << " " 
                      << vg[6] << " " 
                      << vg[7] << " " 
                      << vg[8] << " " 
                      << vg[9] << std::endl;
 
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

} // namespace pwdft

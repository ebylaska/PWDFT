/* ion_ion.cpp
   Author - Eric Bylaska
*/

#include        <cmath>

namespace pwdft {
using namespace pwdft;


/*******************************
 *                             *
 *        ion_ion_e            *
 *                             *
 *******************************/
/*  This function computes the (free-space) Coulomb energy between
   ion cores.

   Entry -
       nion  ---- number of ions
       Q[]   ---- charges of ions
       R[]   ---- coordinates of ions

*/


double ion_ion_e(const int nion, const double Q[], const double R[])
{
   double e1 = 0.0;
   for (auto jj=0; jj<nion; ++jj)
      for (auto ii=0; ii<(jj-1); ++ii)
      {
         double x = R[3*jj]-R[3*ii];
         double y = R[3*jj+1]-R[3*ii+1];
         double z = R[3*jj+2]-R[3*ii+2];
         double r = std::sqrt(x*x + y*y + z*z);
         if (r>1.0e-6)
            e1 += Q[ii]*Q[jj]/r;
      }

   return e1;
}

/*******************************
 *                             *
 *        ion_ion_f            *
 *                             *
 *******************************/
/*   This routine computes the (free-space) Coulomb forces between
    ion cores.
 
    Entry -
        nion  ---- number of ions
        Q[]   ---- charges on the ions
        R[]   ---- coordinates of ions
 
    Exit -
      F[]     ---- force vectors
*/

void ion_ion_f(const int nion, const double Q[], const double R[], double F[])
{
   if (nion>1)
   {
      for (auto jj=0; jj<nion; ++jj)
         for (auto ii=0; ii<(jj-1); ++ii)
         {
            double x = R[3*jj]-R[3*ii];
            double y = R[3*jj+1]-R[3*ii+1];
            double z = R[3*jj+2]-R[3*ii+2];
            double r = std::sqrt(x*x + y*y + z*z);
            if (r>1.0e-6)
            {
               double v = Q[ii]*Q[jj]/(r*r*r);
               F[3*ii]   += x*v;
               F[3*ii+1] += y*v;
               F[3*ii+2] += z*v;
               F[3*jj]   -= x*v;
               F[3*jj+1] -= y*v;
               F[3*jj+2] -= z*v;
            }
         }
   }
}

}

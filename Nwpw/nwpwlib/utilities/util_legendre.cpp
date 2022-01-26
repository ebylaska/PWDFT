
#include        <iostream>
#include        <cstdlib>
#include        <cmath>

#include        "util_legendre.hpp"

namespace pwdft {
using namespace pwdft;

/*********************************************
 *                                           *
 *             util_legendre_lm              *
 *                                           *
 *********************************************/
/*
* Compute the associated legendre polynomial w/ a Condon-Shortley phase?
*/
double util_legendre_lm(const int l, const int m, const double x)
{
   double lgndr;
   int mod_m = abs(m);      

   if (abs(x)>1.0) 
   {
      std::cout << "parameter out of range in legendre_lm" << std::endl;
      return 0.0;
   }
   if (mod_m > l) return 0.0;

   // *** find P(mod_m,mod_m) for mod_m=0 ***
   double p_mm = 1.0;

   // *** find P(mod_m,mod_m) for mod_m > 0 ***
   if (mod_m > 0)
   {
      double tmp  = sqrt((1.0-x)*(1.0+x));
      double fact = 1.0;
      for (auto i=0; i<mod_m; ++i)
      {
         p_mm *= -fact*tmp;
         fact += 2.0;
      }
   }

   // *** find P(l,mod_m) ***
   if (mod_m == l)
      lgndr = p_mm;

   // *** find P(mod_m+1,mod_m) ***
   else
   {
      double p_mp1m = x*(2*mod_m + 1.0)*p_mm;
      if (l==(mod_m+1))
         lgndr = p_mp1m;
      else
      {
         double tmp;
         for (auto i=(mod_m+2); i<=l; ++i)
         {
            tmp    = (x*(2*i-1)*p_mp1m-(i+mod_m-1)*p_mm)/((double) (i-mod_m));
            p_mm   = p_mp1m;
            p_mp1m = tmp;
         }
         lgndr = tmp;
      }
   }

   // *** negative m - this routine is only call with negative m from dtheta_lm and ddtheta_lm ***
   if (m<0)
   {
      double coeff = 1.0;
      for (auto i=1; i<=(2*mod_m); ++i)
         coeff /= ((double) (l - mod_m + i));
      if ((mod_m%2)==1)
         coeff *= -1.0;
      lgndr *= coeff;
   }

   return lgndr;
}

}

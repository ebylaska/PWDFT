
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

/*********************************************
 *                                           *
 *            util_rlegendre_lm              *
 *                                           *
 *********************************************/
/*
* Compute the associated legendre polynomial w/o a Condon-Shortley phase?
*/
double util_rlegendre_lm(const int l, const int m, const double x)
{
   double rlgndr;
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
      rlgndr = p_mm;

   // *** find P(mod_m+1,mod_m) ***
   else
   {
      double p_mp1m = x*(2*mod_m + 1.0)*p_mm;
      if (l==(mod_m+1))
         rlgndr = p_mp1m;
      else
      {
         double tmp;
         for (auto i=(mod_m+2); i<=l; ++i)
         {
            tmp    = (x*(2*i-1)*p_mp1m-(i+mod_m-1)*p_mm)/((double) (i-mod_m));
            p_mm   = p_mp1m;
            p_mp1m = tmp;
         }
         rlgndr = tmp;
      }
   }
   return rlgndr;
}


/***********************************
 *                                 *
 *       util_rlegendre_lm_div     *
 *                                 *
 ***********************************/
/*
  Compute the associated legendre polynomial divided by sin(theta) w/o
  a Condon-Shortley phase.
*/
double util_rlegendre_lm_div(const int l, const int m, const double x)
{
   double f;

   // *** check the arguments ***
   if ((m<0) || (m>l) || (abs(x)>1.0) || (m==0)) std::cout << "parameter ot of range in rlegendre_lm_div" << std::endl;

   // *** P(1,1,x)/dsqrt(1-x**2) ***
   double p_mm = 1.0;

   // *** P(m,m,x)/dsqrt(1-x**2)            ***
   // ***   = (2*m-1)*dsqrt(1-x**2)         ***
   // ***     *(P(m-1,m-1,x)/dsqrt(1-x**2)) ***
   double fact = 3.0;
   double tmp  = sqrt((1.0-x)*(1.0+x));

   for (auto i=2; i<=m; ++i)
   {
      p_mm *= fact*tmp;
      fact += 2.0;
   }

   // *** find P(l,m) ***
   if (m==l)
      f = p_mm;
   else
   {
      // *** find P(m+1,m) ***
      double p_mp1m = x*(2*m+1)*p_mm;

      if (l==(m+1)) 
         f = p_mp1m;
      else
      {
         for (auto i=(m+2); i<=l; ++i)
         {
            tmp = (x*(2*i-1)*p_mp1m - (i+m-1)*p_mm)/((double) (i-m));
            p_mm = p_mp1m;
            p_mp1m = tmp;
         }
         f = tmp;
      }
   }

   return f;
}



/***********************************
 *                                 *
 *       util_legendre_lm_div      *
 *                                 *
 ***********************************/
/*
   Purpose : calculates associated Legendre
             polynomial divided by sin(theta)
             for a scalar argument
*/
double util_legendre_lm_div(const int l, const int m, const double x)
{
   double f;

   // *** check the arguments ***
   if ((m<0) || (m>l) || (abs(x)>1.0) || (m==0)) std::cout << "parameter ot of range in legendre_lm_div" << std::endl;

   // *** P(1,1,x)/sqrt(1-x**2) ***
   double p_mm = -1.0;

   double fact = 3.0;
   double tmp = sqrt((1.0-x)*(1.0+x));

   for (auto i=2; i<=(m); ++i)
   {
      p_mm *= -fact*tmp;
      fact +=  2.0;
   }

   // *** find P(l,m) ***
   if (m==l)
      f = p_mm;
   else
   {
      // *** find P(m+1,m) ***
      double  p_mp1m = x*(2*m+1)*p_mm;

       if (l==(m+1)) 
            f = p_mp1m;
       else
       {
          for (auto i=(m+2); i<=l; ++i)
          {
             tmp = (x*(2*i-1)*p_mp1m - (i+m-1)*p_mm)/(i-m);
             p_mm = p_mp1m;
             p_mp1m = tmp;
          }
          f = tmp;
       }
   }
   return f;
}




/*********************************************
 *                                           *
 *              util_theta_lm                *
 *                                           *
 *********************************************/
/*
*    Name    : theta_lm
*
*
*    Purpose : calculates theta_lm for a scalar cos_theta
*              such that
*
*        Y_lm(cos_theta,phi)=theta_lm(cos_theta)*exp(i*m*phi)
*/
double util_theta_lm(const int l, const int m, const double cos_theta)
{
   double coeff;
   int mod_m = abs(m);
   double fourpi = 16.0*atan(1.0);

   if (m>l) std::cout << "parameter out of order in function theta_lm" << std::endl;

   // *** find coefficient ***
   if (mod_m==0)
      coeff= 1.0;
   else if (mod_m>0)
   {
      coeff= 1.0;
      for (auto i=1; i<=(2*mod_m); ++i)
          coeff /= ((double) (l-mod_m+i));
   }
   coeff *= (2*l+1)/fourpi;
   coeff = sqrt(coeff);
   double f = coeff*util_legendre_lm(l,mod_m,cos_theta);
   if ((m<0) && ((mod_m%2)==1))
      f = -f;

   return f;
}

/*********************************************
 *                                           *
 *              util_rtheta_lm               *
 *                                           *
 *********************************************/
/*
* 
*      Purpose : calculates rtheta_lm for a scalar cos_theta
*                such that
*                                                   {cos(|m|*phi)   m>0
*          T_lm(cos_theta,phi)=rtheta_lm(cos_theta)*{1              m==0
*                                                   {sin(|m|*phi)   m<0
*/

double util_rtheta_lm(const int l, const int m, const double cos_theta)
{
   double coeff;
   int mod_m    = abs(m);
   double twopi = 8.00*atan(1.0);

   if (mod_m>l) std::cout << "parameter out of order in function rtheta_lm" << std::endl;

   // *** find coefficient ***
   if (mod_m==0)
      coeff= 0.5;
   else
   {
      coeff= 1.0;
      for (auto i=1; i<=(2*mod_m); ++i)
         coeff /= ((double) (l-mod_m+i));
   }
   coeff *= (2*l+1)/twopi;
   coeff = sqrt(coeff);
   double f = coeff*util_rlegendre_lm(l,mod_m,cos_theta);

   return f;
}


/*********************************************
 *                                           *
 *              util_ytheta_lm               *
 *                                           *
 *********************************************/
/*
*      Purpose : calculates ytheta_lm for a scalar cos_theta
*                such that
* 
*          Y_lm(cos_theta,phi)=ytheta_lm(cos_theta)*exp(i*m*phi)
*/
double util_ytheta_lm(const int l, const int m, const double cos_theta)
{
   double coeff;
   int mod_m = abs(m);
   double fourpi = 16.0*atan(1.0);

   if (mod_m>l) std::cout << "parameter out of order in function ytheta_lm" << std::endl;

   // *** find coefficient ***
   if (mod_m==0) 
      coeff= 1.0;
   else
   {
      coeff= 1.0;
      for (auto i=1; i<=(2*mod_m); ++i)
         coeff /= ((double) (l-mod_m+i));
   }
   coeff *= (2*l+1)/fourpi;
   coeff = sqrt(coeff);
   double f = coeff*util_legendre_lm(l,mod_m,cos_theta);
   if ((m<0) && ((mod_m%2)==1)) f=-f;
   return f;
}


/*********************************************
 *                                           *
 *              util_rtheta2_lm              *
 *                                           *
 *********************************************/
/*
*      Purpose : calculates rtheta2_lm for a scalar cos_theta
*                such that
*                                                            {cos(|m|*phi)   m>0
*          dT_lm(cos_theta,phi)/dtheta=rtheta2_lm(cos_theta)*{1              m==0
*                                                            {sin(|m|*phi)   m<0
*/
double util_rtheta2_lm(const int l, const int m, const double cos_theta)
{
   double coeff;
   int mod_m = abs(m);
   double twopi = 8.0*atan(1.0);

   if ((mod_m+1)>l) return 0.0;

   // *** find coefficient ***
   if (mod_m==0)
      coeff= 0.5;
   else
   {
      coeff= 1.00;
      for (auto i=1; i<=(2*mod_m); ++i)
         coeff /= ((double) (l-mod_m+i));
   }
   coeff *= (2*l+1)/twopi;
   coeff = sqrt(coeff);

   double f = -coeff*util_rlegendre_lm(l,mod_m+1,cos_theta);
   return f;
}


/*********************************************
 *                                           *
 *              util_ytheta2_lm              *
 *                                           *
 *********************************************/
/*
*      Purpose : calculates ytheta_lm for a scalar cos_theta
*                such that
* 
*          dY_lm(cos_theta,phi)/dtheta=ytheta_lm(cos_theta)*exp(i*m*phi)
*/
double util_ytheta2_lm(const int l, const int m, const double cos_theta)
{
   double coeff;
   int mod_m = abs(m);
   double fourpi = 16.0*atan(1.0);

   if ((mod_m+1)>l) return 0.0;

   // *** find coefficient ***
   if (mod_m==0)
      coeff= 1.0;
   else
   {
      coeff= 1.00;
      for (auto i=1; i<=(2*mod_m); ++i)
         coeff /= ((double) (l-mod_m+i));
   }
   coeff *= (2*l+1)/fourpi;
   coeff = sqrt(coeff);

   double f = coeff*util_legendre_lm(l,mod_m+1,cos_theta);
   if ((m<0) && (((mod_m+1)%2)==1)) f=-f;

   return (-f);
}



/*********************************************
 *                                           *
 *            util_theta_lm_div              *
 *                                           *
 *********************************************/
/*
    Purpose : calculates theta_lm/sin(theta)
*/

double util_theta_lm_div(const int l, const int m, const double cos_theta)
{
   double coeff;
   int mod_m = abs(m);
   double fourpi = 16.0*atan(1.0);
      
   if ( m>l) std::cout << "parameter out of order in function theta_lm_div" << std::endl;
           
   // *** find coefficient ***
   if (mod_m==0)
      coeff = 1.0;
   else
   {
      coeff = 1.0;
      for (auto i=1; i<=(2*mod_m); ++i)
         coeff /= ((double) (l-mod_m+i));
   }
   coeff *= (2*l+1)/fourpi;
   coeff = sqrt(coeff);

   double f = coeff*util_legendre_lm_div(l,mod_m,cos_theta);
   if ((m<0) && ((mod_m%2)==1)) f = -f;
   return f;
}



/*********************************************
 *                                           *
 *            util_rtheta_lm_div             *
 *                                           *
 *********************************************/
/*
*      Purpose : calculates rtheta_lm_div for a scalar cos_theta
*                such that
*                                                             {cos(|m|*phi)   m>0
*      T_lm(cos_theta,phi)/sin_theta=rtheta_lm_div(cos_theta)*{1              m==0
*                                                             {sin(|m|*phi)   m<0
*/
double util_rtheta_lm_div(const int l, const int m, const double cos_theta)
{
   double coeff;
   int mod_m = abs(m);
   double twopi = 8.0*atan(1.0);

   if (mod_m>l) std::cout << "parameter out of order in function rtheta_lm_div" << std::endl;

   // *** find coefficient ***
   if (mod_m==0)
      coeff= 0.5;
   else
   {
      coeff= 1.0;
      for (auto i=1; i<(2*mod_m); ++i)
         coeff /= ((double) (l-mod_m+i));
         
   }
   coeff *= (2*l+1)/twopi;
   coeff = sqrt(coeff);

   double f = coeff*util_rlegendre_lm_div(l,mod_m,cos_theta);
   return f;
}

/*********************************************
 *                                           *
 *            util_ytheta_lm_div             *
 *                                           *
 *********************************************/
/*
   Purpose : calculates ytheta_lm/sin(theta)
*/
double util_ytheta_lm_div(const int l, const int m, const double cos_theta)
{
   double coeff;
   int mod_m = abs(m);
   double fourpi = 16.0*atan(1.0);

   if (m>l) std::cout << "parameter out of order in function theta_lm_div" << std::endl;

   // *** find coefficient ***
   if (mod_m==0) 
      coeff = 1.0;
   else 
   {
      coeff = 1.0;
      for (auto i=1; i<=(2*mod_m); ++i) 
         coeff /=  ((double) (l-mod_m+i));
   }
   coeff *= (2*l+1)/fourpi;
   coeff = sqrt(coeff);

   double f = coeff*util_legendre_lm_div(l,mod_m,cos_theta);
   if ((m<0) && ((mod_m%2)==1)) f = -f;
   return f;
}



}

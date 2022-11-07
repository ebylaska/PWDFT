
#include	<iostream>
#include	<cstdlib>
#include	<cmath>
#include        "Parallel.hpp"

#include	"util_gamma.hpp"

namespace pwdft {




/**************************************
 *                                    *
 *           util_ln_gamma            *
 *                                    *
 **************************************/

/* This is a slightly modified version of Log(Gamma) function program  */

double util_ln_gamma(const double x)
{
   double cof[6] = {76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
   double stp    = 2.5066282746310005;

   double ser = 1.000000000190015;
   double y   = x;
   for (auto j=0; j<6; ++j)
   {
      y   += 1.0;
      ser += cof[j]/y;
   }
   double tmp = x+5.50;

   tmp = (x+0.50)*log(tmp)-tmp;

   return (tmp + log(stp*ser/x));
}

/**************************************
 *                                    *
 *             util_gamma             *
 *                                    *
 **************************************/

double util_gamma(const double x)
{
   return (exp(util_ln_gamma(x)));
}

/**************************************
 *                                    *
 *             util_gser              *
 *                                    *
 **************************************/

// convert   subroutine gser(gamser,a,x,gln)
//
double util_gser(const double a, const double x)
{
   int ITMAX  = 1000;
   double EPS = 3.0e-12;

   if (x < 0.0)  { std::cout << "x < 0 in gser" << std::endl; return 0.0; }
   if (x == 0.0) { return 0.0; }   

   double gln = util_ln_gamma(a);
   double ap  = a;
   double sum = 1.0/a;
   double del = sum;
   bool done = false;

   int n=1;
   while ((n<=ITMAX) && (!done))
   {
      ap  += 1.0;
      del *= x/ap;
      sum += del;
      if (std::abs(del)<(std::abs(sum)*EPS)) done = true;
      ++n;
   }
   if (n>ITMAX) std::cout << "a too large, ITMAX too small in gser" << std::endl;

   return (sum*exp(-x+a*log(x)-gln));
}

/**************************************
 *                                    *
 *             util_gcf               *
 *                                    *
 **************************************/

// converting SUBROUTINE gcf(gammcf,a,x,gln)
//
double util_gcf(const double a, const double x)
{
   int ITMAX=1000;
   double EPS    = 3.0e-12;
   double FPMIN  = 1.0e-30;
   double undovl = -20.0*2.3025;

   double gln = util_ln_gamma(a);

   double b = x + 1.00 - a;
   double c = 1.00/FPMIN;
   double d = 1.00/b;
   double h = d;
   bool done = false;
   double gammcf;

   int i=1;
   while ((i<=ITMAX) && (!done))
   {
      double an = -i*(i-a);

      b += 2.00;
      d = an*d + b;

      if (std::abs(d) < FPMIN) d = FPMIN;
      c= b + an/c;
      if (std::abs(c) < FPMIN) c = FPMIN;
      d = 1.0/d;
      double  del = d*c;
      h *= del;
      if (std::abs(del-1.00)<EPS) done = true;
      ++i;
   }
   if (i>ITMAX) std::cout << "a too large, ITMAX too small in util_gcf" << std::endl;

   if ((-x+a*log(x)-gln)>undovl)
   {
      gammcf = exp(-x+a*log(x)-gln)*h;
   }
   else
   {
      gammcf = 0.0;
   }
   return gammcf;
}

/**************************************
 *                                    *
 *             util_gammap            *
 *                                    *
 **************************************/
double util_gammp(const double a, const double x)
{
   if ((x<0.0) || (a<=0.0)) { std::cout << "bad arguments in util_gammp" << std::endl; return 0.0;}

   double gammp;

   if (x < (a+1.0))
   {
      gammp = util_gser(a,x);
   }
   else
   {
      gammp = 1.0 - util_gcf(a,x);
   }
   return gammp;
}





}


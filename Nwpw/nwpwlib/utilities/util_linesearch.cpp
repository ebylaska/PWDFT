#include	<cstdlib>
#include	<cmath>
#include	"util_linesearch.hpp"

namespace pwdft {
using namespace pwdft;


/* util_linesearch common block */
static int  maxiter = 8;
static int  counter = 0;

/**************************************
 *                                    *
 *       util_linesearch_init         *
 *                                    *
 **************************************/
void util_linesearch_init()
{
    maxiter = 8;
    counter = 0;
}

/**************************************
 *                                    *
 *    util_linesearch_maxiter_set     *
 *                                    *
 **************************************/
void util_linesearch_maxiter_set(const int maxiter0) { maxiter = maxiter0; }

/**************************************
 *                                    *
 *      util_linesearch_maxiter       *
 *                                    *
 **************************************/
int util_linesearch_maxiter() { return maxiter; }

/**************************************
 *                                    *
 *      util_linesearch_counter       *
 *                                    *
 **************************************/
int util_linesearch_counter() { return counter; }

/**************************************
 *                                    *
 *           Order_Values             *
 *                                    *
 **************************************/
/* This subroutine makes f(indx(1)) < f(indx(2)) < f(indx(3)) < ....
   via a bubble sort
   Entry - n,f
   Exit - returns indx
*/

void Order_Values(int n, double f[], int indx[])
{
   int idum;

   for (int i=0; i<n; ++i) indx[i] = i;
   
   for (int i=0; i<(n-1); ++i)   
   for (int j=i+1; j<n; ++j)   
      if (f[indx[j]] < f[indx[i]])
      {
         idum    = indx[i];
         indx[i] = indx[j];
         indx[j] = idum;
      }
}

/**************************************
 *                                    *
 *           util_linesearch          *
 *                                    *
 **************************************/

double util_linesearch(const double t0, const double f0, const double df0, double deltat,
                       double (*FUNC1)(double), double (*DFUNC1)(double), 
                       const double tolerance, 
                       double *tmin_ret, double *deltaE_ret, 
                       const int stoptype)
{
   bool secant = true;
   bool notfinished;

   int    indx[3],iteration;

   double tmin;
   double t[3],f[3],df[3];
   double t_last, f_last;
   double t_first,f_first,df_first;
   double up,down,fmin,dfmin,deltaf;

   ++counter;

   t[0]  = t0; 
   f[0]  = f0; 
   df[0] = df0;
   t_last = t[0];
   f_last = f[0];
   t_first  = t[0];
   f_first  = f[0];
   df_first = df[0];

   f[1]  =  FUNC1(t[0]+deltat);
   df[1] = DFUNC1(t[0]+deltat);

   iteration = 1;

   /* make sure that f1 < f0 */
   while ((f[1] > f[0]) && (iteration <= maxiter))
   {
      deltat *= 0.5;
      f[1]  =  FUNC1(t[0]+deltat);
      df[1] = DFUNC1(t[0]+deltat);
      ++iteration;
   }
   t[1] = t[0] + deltat;
   t_last = t[1];
   f_last = f[1];

   /* use secant method to generate f[2] */
   deltat = -df[0]*(t[1]-t[0])/(df[1]-df[0]);
   t[2]   = t[0] + deltat;
   f[2]   =  FUNC1(t[2]);
   df[2]  = DFUNC1(t[2]);
   ++iteration;
   t_last = t[2];
   f_last = f[2];

   /* sort the function values */
   Order_Values(3,f,indx);

   deltaf = f[indx[1]] - f[indx[0]];

   if (stoptype==1)
      notfinished = (fabs(deltaf)>tolerance) &&  (iteration<=maxiter);
   else
      notfinished = (fabs(df[indx[0]]/df_first)>tolerance) && (iteration<=maxiter);

   while (notfinished) 
   {

      /* use secant interpolation to generate tmin */
      if (secant)
      {
         deltat = -df[indx[0]] * (t[indx[1]]-t[indx[0]]) / (df[indx[1]]-df[indx[0]]);
         tmin = t[indx[0]] + deltat;
         fmin  =  FUNC1(tmin);
         dfmin = DFUNC1(tmin);
         ++iteration;
         t_last = tmin;
         f_last = fmin;

         /* finish using secant method */
         if (fmin >= f[indx[0]])
         {
            secant = false;
            if (fmin < f[indx[2]])
            {
                t[indx[2]]  = tmin;
                f[indx[2]]  = fmin;
                df[indx[2]] = dfmin;
                Order_Values(3,f,indx);
            }
         }
      }


      /* use quadradic interpolation to generate tmin */
      if (!secant)
      {
         up  = (t[1]*t[1] - t[2]*t[2])*f[0]
             + (t[2]*t[2] - t[0]*t[0])*f[1]
             + (t[0]*t[0] - t[1]*t[1])*f[2];
         down = (t[1] - t[2])*f[0]
              + (t[2] - t[0])*f[1]
              + (t[0] - t[1])*f[2];

         /* check added by E.Apra */
         if (fabs(down) > (tolerance*tolerance))
         {
            tmin = 0.50*up/down;
            fmin  =  FUNC1(tmin);
            dfmin = DFUNC1(tmin);
            ++iteration;
            t_last = tmin;
            f_last = fmin;
         }
         /* parabola fit failed - exit loop */
         else
         {
            tmin=t[indx[2]];
            fmin=f[indx[2]]+tolerance;
            iteration = maxiter+1;
         }
      }

      /* tolerance check and replacement */
      if (fmin<f[indx[2]])
      {
         t[indx[2]]  = tmin;
         f[indx[2]]  = fmin;
         df[indx[2]] = dfmin;
         Order_Values(3,f,indx);
         deltaf = f[indx[1]] - f[indx[0]];
      }
      else
         deltaf=0.0;

      if (stoptype==1)
         notfinished = (fabs(deltaf)>tolerance) && (iteration<=maxiter);
      else
         notfinished = (fabs(df[indx[0]]/df_first)>tolerance) && (iteration<=maxiter);
   }


   /* make sure that tmin is last functions call */
   tmin = t[indx[0]];
   fmin = f[indx[0]];
   if (tmin!=t_last) fmin = FUNC1(tmin);

   *tmin_ret   = tmin;
   *deltaE_ret = (fmin-f_first);
  
   return fmin;
}

}

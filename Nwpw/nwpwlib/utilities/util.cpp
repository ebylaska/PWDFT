
//extern "C" {
//#include        "compressed_io.h"
//}
#include        "compressed_io.hpp"

#include	<cstdlib>
#include	<cmath>
#include        "Parallel.hpp"
#include	"util.hpp"

void c_aindexcopy(const int n, const int *indx, double *A, double *B)
{
   int ii,jj;
   ii = 0;
   for (int i=0; i<n; ++i)
   {
      jj      = 2*indx[i];
      B[ii]   = A[jj];
      B[ii+1] = A[jj+1];
      ii += 2;
   }
}

void c_bindexcopy(const int n, const int *indx, double *A, double *B)
{
   int ii,jj;
   ii = 0;
   for (int i=0; i<n; ++i)
   {
      jj      = 2*indx[i];
      B[jj]   = A[ii];
      B[jj+1] = A[ii+1];
      ii += 2;
   }
}
void c_bindexcopy_conjg(const int n, const int *indx, double *A, double *B)
{
   int ii,jj;
   ii = 0;
   for (int i=0; i<n; ++i)
   {
      jj      = 2*indx[i];
      B[jj]   = A[ii];
      B[jj+1] = -A[ii+1];
      ii += 2;
   }
}


void t_aindexcopy(const int n, const int *indx, double *A, double *B)
{
   for (int i=0; i<n; ++i)
      B[i] = A[indx[i]];
}

void t_bindexcopy(const int n, const int *indx, double *A, double *B)
{
   for (int i=0; i<n; ++i)
      B[indx[i]] = A[i];
}


void i_aindexcopy(const int n, const int *indx, int *A, int *B)
{
   for (int i=0; i<n; ++i)
      B[i] = A[indx[i]];
}




/**************************************
 *                                    *
 *            eigsrt                  *
 *                                    *
 **************************************/

void eigsrt(double *D, double *V, int n)
{
   int i,j,k;
   double p;

   for (i=0; i<(n-1); ++i)
   {
      k = i;
      p = D[i];
      for(j=i+1; j<n; ++j)
         if (D[j]>=p)
         {
            k = j;
            p = D[j];
         }

      if (k!=i)
      {
         D[k] = D[i];
         D[i] = p;
         for (j=0; j<n; ++j)
         {
            p = V[j+i*n];
            V[j+i*n] = V[j+k*n];
            V[j+k*n] = p;
         }
      }
   }
}



/**************************************
 *                                    *
 *         util_getfilling            *
 *                                    *
 **************************************/
void util_getfilling(int f, int nfft[], int *filling, double zvalue[])
{
   int h = (f%2);
   int f2 = (f-h)/2;
   int k=0;
   while ((k+1)*(k+2)*(k+3) <= (6*f2))
      ++k;

   f2 -= k*(k+1)*(k+2)/6;
   int j=0;
   while ((j+1)*(j+2) <= (2*f2))
      ++j;

   int i = f2 - j*(j+1)/2;

   filling[0] = i;
   if (i==0) 
   {
      filling[1] = j;
      filling[2] = k;
   }
   else
   {

      if (((j+1)%2)==0)
         filling[1] = (j+1)/2;
      else 
         filling[1] = -(j+1)/2;

      if (((k+1)%2)==0)
         filling[2] = (k+1)/2;
      else 
         filling[2] = -(k+1)/2;
   }

   if ((i==0) && (j==0) && (k==0))
   {
      filling[3] = 0;
      zvalue[0] = 1.0;
      zvalue[1] = 0.0;
   }
   else if (h==0)
   {
      filling[3] = 2;
      zvalue[0] = 1.0/sqrt(2.0);
      zvalue[1] = 0.0;
   }
   else
   {
      filling[3] = -2;
      zvalue[0] = 0.0;
      zvalue[1] = 1.0/sqrt(2.0);
   }

   /* modularize the filling */
   int inc2c =( nfft[0]/2 + 1);
   filling[0] = (filling[0]+inc2c)%inc2c;
   filling[1] = (filling[1]+nfft[1])%nfft[1];
   filling[2] = (filling[2]+nfft[2])%nfft[2];
}


/**************************************
 *                                    *
 *           util_random              *
 *                                    *
 **************************************/

/* returns a random number between 0 and 1

   Entry - seed - if zero set a seed with srand
   Exit - returns a random number between 0 and 1
   Uses - rand and srand stdlib functions
*/
double util_random(const int seed)
{
   if (seed>0) std::srand(((double) seed));
   return ( (double) std::rand()/RAND_MAX);
}


/**************************************
 *                                    *
 *           util_filefind            *
 *                                    *
 **************************************/
bool util_filefind(Parallel *myparall, char *fname)
{
   int ifound;

   if (myparall->is_master())
      ifound = cfileexists(fname);

   myparall->Brdcst_iValue(0,0,&ifound);

   return (ifound>0);
}


/**************************************
 *                                    *
 *           util_spline              *
 *                                    *
 **************************************/
void   util_spline(double *x, double *y, int n, double yp1, double ypn, double *y2, double *utmp)
{
   double sig,p,qn,un;

   if (yp1>0.99e30)
   {
      y2[0]   = 0.0;
      utmp[0] = 0.0;
   }
   else
   {
      y2[0]   = -0.5;
      utmp[0] = 3.0 / (x[1]-x[0]) * ( (y[1]-y[0]) / (x[1]-x[0]) - yp1);
   }
   for (auto i=1; i<(n-1); ++i)
   {
      sig = (x[i]-x[i-1]) / (x[i+1]-x[i-1]);
      p = sig*y2[i-1] + 2.00;
      y2[i] = (sig-1.00) / p;
      utmp[i] = (
                 6.00 *
                 (
                   (y[i+1]-y[i  ]) / (x[i+1]-x[i  ])
                 - (y[i  ]-y[i-1]) / (x[i  ]-x[i-1])
                 )
                 / (x[i+1]-x[i-1]) - sig*utmp[i-1]
               )
               / p;
   }
   
   if (ypn>0.99e30)
   {
      qn = 0.0;
      un = 0.0;
   }
   else
   {
      qn = 0.5;
      un = 3.00 / (x[n-1]-x[n-2])
              * ( ypn - (y[n-1]-y[n-2]) / (x[n-1]-x[n-2]) );
   
   }

   y2[n-1] = (un-qn*utmp[n-2]) / (qn*y2[n-2]+1.00);
   for (auto k=n-2; k>=0; --k)
      y2[k] = y2[k]*y2[k+1] + utmp[k];

}



/**************************************
 *                                    *
 *           util_splint              *
 *                                    *
 **************************************/
double util_splint(double *xa, double *ya, double *y2a, int n, int nx, double x)
{
   int khi = nx;
   int klo = nx-1;

   while ((xa[klo]>x) || (xa[khi]<x))
   {
      if (xa[klo]>x)
      {
        klo = klo - 1;
        khi = khi - 1;
      }
      if (xa[khi]<x)
      {
           klo = klo + 1;
           khi = khi + 1;
      }
   }

   double h = xa[khi]-xa[klo];
   double a = (xa[khi]-x)/h;
   double b = (x-xa[klo])/h;
   
   return  (a*ya[klo] + b*ya[khi] + ((a*a*a-a)*y2a[klo] + (b*b*b-b)*y2a[khi]) * h*h / 6.0);
}

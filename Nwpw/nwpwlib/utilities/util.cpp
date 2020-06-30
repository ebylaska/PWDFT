
extern "C" {
#include        "compressed_io.h"
}

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


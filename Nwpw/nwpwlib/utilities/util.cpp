
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


/* d3db.C
   Author - Eric Bylaska

	this class is used for defining 3d parallel maps
*/
/*
#include        <iostream>
#include        <cstdio>
#include        <stdio.h>
#include        <cstdlib>
using namespace std;
*/

extern "C" {
#include        "compressed_io.h"
}

#include	"Parallel.hpp"
#include	"util.hpp"
#include	"fft.h"
#include	"blas.h"
#include        <cmath>

#include	"d3db.hpp"

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

d3db::d3db(Parallel *inparall,const int inmaptype, const int nx, const int ny, const int nz)
   : Mapping3(inmaptype,inparall->np_i(),inparall->taskid_i(),nx,ny,nz)
{
   int i,j,k,index1,index2,it,proc_to,proc_from;
   int j1,j2,k1,k2,nyh,nzh;
   int phere,pto,pfrom;

   parall = inparall;

   if (maptype==1)
   {
      iq_to_i1 = new int*[1]; iq_to_i1[0] = new int[(nx/2+1)*ny*nq];
      iq_to_i2 = new int*[1]; iq_to_i2[0] = new int[(nx/2+1)*ny*nq];
      i1_start = new int*[1]; i1_start[0] = new int[nz+1];
      i2_start = new int*[1]; i2_start[0] = new int[nz+1];
      index1 = 0;
      index2 = 0;
      for (it=0; it<np; ++it)
      {
         proc_to   = (taskid+it)%np;
         proc_from = (taskid-it+np)%np;
         i1_start[0][it] = index1;
         i2_start[0][it] = index2;

         for (k=0; k<nz; ++k)
         for (j=0; j<ny; ++j)
         {
            /* packing scheme */
            phere = ijktop(0,0,k);
            pto   = ijktop(0,0,j);
            if ((phere==taskid) && (pto==proc_to))
               for (i=0; i<(nx/2+1); ++i)
               {
                  iq_to_i1[0][ijktoindex(i,j,k)] = index1;
                  ++index1;
               }

            /* unpacking scheme */
            phere = ijktop(0,0,j);
            pfrom = ijktop(0,0,k);
            if ((phere==taskid) && (pfrom==proc_from))
               for (i=0; i<(nx/2+1); ++i)
               {
                  iq_to_i2[0][ijktoindex(i,k,j)] = index2;
                  ++index2;
               }
         }
      }
      i1_start[0][np] = index1;
      i2_start[0][np] = index2;

   }
   else
   {
      iq_to_i1 = new int*[6];
      iq_to_i1[0] = new int[(nx/2+1)*nq1];
      iq_to_i1[1] = new int[ny*nq2];
      iq_to_i1[2] = new int[nz*nq3];
      iq_to_i1[3] = new int[ny*nq2];
      iq_to_i1[4] = new int[(nx/2+1)*nq1];
      iq_to_i1[5] = new int[nz*nq3];

      iq_to_i2 = new int*[6];
      iq_to_i2[0] = new int[ny*nq2];
      iq_to_i2[1] = new int[nz*nq3];
      iq_to_i2[2] = new int[ny*nq2];
      iq_to_i2[3] = new int[(nx/2+1)*nq1];
      iq_to_i2[4] = new int[nz*nq3];
      iq_to_i2[5] = new int[(nx/2+1)*nq1];

      i1_start = new int*[6];
      for (i=0; i<6; ++i) 
         i1_start[i] = new int[np+1];

      i2_start = new int*[6];
      for (i=0; i<6; ++i) 
         i2_start[i] = new int[np+1];

      /**************************************************/
      /* map1to2 mapping - done - tranpose operation #0 */
      /**************************************************/
      index1 = 0;
      index2 = 0;
      for (it=0; it<np; ++it)
      {
         proc_to   = (taskid+it)%np;
         proc_from = (taskid-it+np)%np;
         i1_start[0][it] = index1;
         i2_start[0][it] = index2;
         for (k=0; k<nz;       ++k)
         for (j=0; j<ny;       ++j)
         for (i=0; i<(nx/2+1); ++i)
         {
            phere = ijktop2(i,j,k);
            pto   = ijktop1(i,j,k);

            /* packing scheme */
            if ((phere==taskid) && (pto==proc_to))
            {
               iq_to_i1[0][ijktoindex2t(i,j,k)] = index1;
               ++index1;
            }
            /* unpacking scheme */
            if ((pto==taskid) && (phere==proc_from))
            {
               iq_to_i2[0][ijktoindex1(i,j,k)] = index2;
               ++index2;
            }
         }
      }
      i1_start[0][np] = index1;
      i2_start[0][np] = index2;

      /**************************************************/
      /* map2to3 mapping - done - tranpose operation #1 */
      /**************************************************/
      index1 = 0;
      index2 = 0;
      for (it=0; it<np; ++it)
      {
         proc_to   = (taskid+it)%np;
         proc_from = (taskid-it+np)%np;
         i1_start[1][it] = index1;
         i2_start[1][it] = index2;
         for (k=0; k<nz;       ++k)
         for (j=0; j<ny;       ++j)
         for (i=0; i<(nx/2+1); ++i)
         {
            phere = ijktop1(i,j,k);
            pto   = ijktop(i,j,k);

            /* packing scheme */
            if ((phere==taskid) && (pto==proc_to))
            {
               iq_to_i1[1][ijktoindex1(i,j,k)] = index1;
               ++index1;
            }
            /* unpacking scheme */
            if ((pto==taskid) && (phere==proc_from))
            {
               iq_to_i2[1][ijktoindex(i,j,k)] = index2;
               ++index2;
            }
         }
      }
      i1_start[1][np] = index1;
      i2_start[1][np] = index2;

      /**************************************************/
      /* map3to2 mapping - done - tranpose operation #2 */
      /**************************************************/
      index1 = 0;
      index2 = 0;
      for (it=0; it<np; ++it)
      {
         proc_to   = (taskid+it)%np;
         proc_from = (taskid-it+np)%np;
         i1_start[2][it] = index1;
         i2_start[2][it] = index2;
         for (k=0; k<nz;       ++k)
         for (j=0; j<ny;       ++j)
         for (i=0; i<(nx/2+1); ++i)
         {
            phere = ijktop(i,j,k);
            pto   = ijktop1(i,j,k);

            /* packing scheme */
            if ((phere==taskid) && (pto==proc_to))
            {
               iq_to_i1[2][ijktoindex(i,j,k)] = index1;
               ++index1;
            }
            /* unpacking scheme */
            if ((pto==taskid) && (phere==proc_from))
            {
               iq_to_i2[2][ijktoindex1(i,j,k)] = index2;
               ++index2;
            }
         }
      }
      i1_start[2][np] = index1;
      i2_start[2][np] = index2;

      /**************************************************/
      /* map2to1 mapping - done - tranpose operation #3 */
      /**************************************************/
      index1 = 0;
      index2 = 0;
      for (it=0; it<np; ++it)
      {
         proc_to   = (taskid+it)%np;
         proc_from = (taskid-it+np)%np;
         i1_start[3][it] = index1;
         i2_start[3][it] = index2;
         for (k=0; k<nz;       ++k)
         for (j=0; j<ny;       ++j)
         for (i=0; i<(nx/2+1); ++i)
         {
            phere = ijktop1(i,j,k);
            pto   = ijktop2(i,j,k);

            /* packing scheme */
            if ((phere==taskid) && (pto==proc_to))
            {
               iq_to_i1[3][ijktoindex1(i,j,k)] = index1;
               ++index1;
            }
            /* unpacking scheme */
            if ((pto==taskid) && (phere==proc_from))
            {
               iq_to_i2[3][ijktoindex2t(i,j,k)] = index2;
               ++index2;
            }
         }
      }
      i1_start[3][np] = index1;
      i2_start[3][np] = index2;

      /**************************************************/
      /* map1to3 mapping - done - tranpose operation #4 */
      /**************************************************/
      index1 = 0;
      index2 = 0;
      for (it=0; it<np; ++it)
      {
         proc_to   = (taskid+it)%np;
         proc_from = (taskid-it+np)%np;
         i1_start[4][it] = index1;
         i2_start[4][it] = index2;
         for (k=0; k<nz;       ++k)
         for (j=0; j<ny;       ++j)
         for (i=0; i<(nx/2+1); ++i)
         {
            phere = ijktop2(i,j,k);
            pto   = ijktop(i,j,k);

            /* packing scheme */
            if ((phere==taskid) && (pto==proc_to))
            {
               iq_to_i1[4][ijktoindex2t(i,j,k)] = index1;
               ++index1;
            }
            /* unpacking scheme */
            if ((pto==taskid) && (phere==proc_from))
            {
               iq_to_i2[4][ijktoindex(i,j,k)] = index2;
               ++index2;
            }
         }
      }
      i1_start[4][np] = index1;
      i2_start[4][np] = index2;

      /**************************************************/
      /* map3to1 mapping - done - tranpose operation #5 */
      /**************************************************/
      index1 = 0;
      index2 = 0;
      for (it=0; it<np; ++it)
      {
         proc_to   = (taskid+it)%np;
         proc_from = (taskid-it+np)%np;
         i1_start[5][it] = index1;
         i2_start[5][it] = index2;
         for (k=0; k<nz;       ++k)
         for (j=0; j<ny;       ++j)
         for (i=0; i<(nx/2+1); ++i)
         {
            phere = ijktop(i,j,k);
            pto   = ijktop2(i,j,k);

            /* packing scheme */
            if ((phere==taskid) && (pto==proc_to))
            {
               iq_to_i1[5][ijktoindex(i,j,k)] = index1;
               ++index1;
            }
            /* unpacking scheme */
            if ((pto==taskid) && (phere==proc_from))
            {
               iq_to_i2[5][ijktoindex2t(i,j,k)] = index2;
               ++index2;
            }
         }
      }
      i1_start[5][np] = index1;
      i2_start[5][np] = index2;

   }

   /* setup timereverse indexes */
   zplane_size = timereverse_size();
   t_iq_to_i1 = new int[zplane_size];
   t_iq_to_i2 = new int[zplane_size];
   t_i1_start = new int[np+1];
   t_i2_start = new int[np+1];
   nyh = ny/2;
   nzh = nz/2;
   index1 = 0;
   index2 = 0;
   for (it=0; it<np; ++it)
   {
      proc_to   = (taskid+it)%np;
      proc_from = (taskid-it+np)%np;
      t_i1_start[it] = index1;
      t_i2_start[it] = index2;

      /* k=(0,0,k3) */
      for (k=1; k<nzh; ++k)
      {
         k1 = k;
         k2 = -k;
         if (k1<0) k1 += nz;
         if (k2<0) k2 += nz;
         phere = ijktop(0,0,k1);
         pto   = ijktop(0,0,k2);

         /* packing scheme */
         if ((phere==taskid) && (pto==proc_to))
         {
            t_iq_to_i1[index1] = ijktoindex(0,0,k1);
            ++index1;
         }
         /* unpacking scheme */
         if ((pto==taskid) && (phere==proc_from))
         {
            t_iq_to_i2[index2] = ijktoindex(0,0,k2);
            ++index2;
         }
      }

       /* k=(0,k2,k3) */
      for (k=(-nzh+1); k<nzh; ++k)
      for (j=1; j<nyh; ++j)
      {
         j1 = j;  
         k1 = k;
         if (j1<0) j1 += ny;
         if (k1<0) k1 += nz;
         j2 = -j; 
         k2 = -k;
         if (j2<0) j2 += ny;
         if (k2<0) k2 += nz;
         phere = ijktop(0,j1,k1);
         pto   = ijktop(0,j2,k2);

          /* packing scheme */
         if ((phere==taskid) && (pto==proc_to))
         {
            t_iq_to_i1[index1] = ijktoindex(0,j1,k1);
            ++index1;
         }
         /* unpacking scheme */
         if ((pto==taskid) && (phere==proc_from))
         {
            t_iq_to_i2[index2] = ijktoindex(0,j2,k2);
            ++index2;
         }
      }
   }
   t_i1_start[np] = index1;
   t_i2_start[np] = index2;


   /* setup ffts */
   tmpx = new double[2*(2*nx+15)];
   tmpy = new double[2*(2*ny+15)];
   tmpz = new double[2*(2*nz+15)];
   drffti_(&nx,tmpx);
   dcffti_(&ny,tmpy);
   dcffti_(&nz,tmpz);

}

/********************************
 *                              *
 *          Destructors         *
 *                              *
 ********************************/
d3db::~d3db()
{
   int i;

   if (maptype==1)
   {
      delete [] iq_to_i1[0];
      delete [] iq_to_i1;

      delete [] iq_to_i2[0];
      delete [] iq_to_i2;

      delete [] i1_start[0];
      delete [] i1_start;
      delete [] i2_start[0];
      delete [] i2_start;
   }
   else
   {
      for (i=0; i<6; ++i)
      {
         delete [] iq_to_i1[i];
         delete [] iq_to_i2[i];
         delete [] i1_start[i];
         delete [] i2_start[i];
      }
      delete [] iq_to_i1;
      delete [] iq_to_i2;
      delete [] i1_start;
      delete [] i2_start;
   }

   delete [] t_iq_to_i1;
   delete [] t_iq_to_i2;
   delete [] t_i1_start;
   delete [] t_i2_start;

   delete [] tmpx;
   delete [] tmpy;
   delete [] tmpz;
}



/********************************
 *                              *
 *         d3db::r_alloc        *
 *                              *
 ********************************/
double * d3db::r_alloc()
{
   double *ptr = new double[n2ft3d];
   return ptr;
}

/********************************
 *                              *
 *         d3db::r_nalloc       *
 *                              *
 ********************************/
double * d3db::r_nalloc(const int nn)
{
   double *ptr = new double[n2ft3d*nn];
   return ptr;
}

/********************************
 *                              *
 *         d3db::r_dealloc      *
 *                              *
 ********************************/
void d3db::r_dealloc(double *ptr)
{
   delete [] ptr;
}

/********************************
 *                              *
 *         d3db::r_zero         *
 *                              *
 ********************************/
void d3db::r_zero(double *ptr)
{
   int i;
   int m = n2ft3d%7;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr[i] = 0.0;
   if (n2ft3d<7) 
      return;

   for (i=m; i<n2ft3d; i+=7)
   {
      ptr[i]   = 0.0;
      ptr[i+1] = 0.0;
      ptr[i+2] = 0.0;
      ptr[i+3] = 0.0;
      ptr[i+4] = 0.0;
      ptr[i+5] = 0.0;
      ptr[i+6] = 0.0;
   }
   return;
}

/********************************
 *                              *
 *        d3db::rr_copy         *
 *                              *
 ********************************/
void d3db::rr_copy(const double *ptr1, double *ptr2)
{
   int i;
   int m = n2ft3d%7;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr2[i] = ptr1[i];
   if (n2ft3d<7) 
      return;

   for (i=m; i<n2ft3d; i+=7)
   {
      ptr2[i]   = ptr1[i];
      ptr2[i+1] = ptr1[i+1];
      ptr2[i+2] = ptr1[i+2];
      ptr2[i+3] = ptr1[i+3];
      ptr2[i+4] = ptr1[i+4];
      ptr2[i+5] = ptr1[i+5];
      ptr2[i+6] = ptr1[i+6];
   }
   return;
}

/********************************
 *                              *
 *        d3db::tt_copy         *
 *                              *
 ********************************/
void d3db::tt_copy(const double *ptr1, double *ptr2)
{
   int i;
   int m = nfft3d%7;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr2[i] = ptr1[i];
   if (nfft3d<7)
      return;

   for (i=m; i<nfft3d; i+=7)
   {
      ptr2[i]   = ptr1[i];
      ptr2[i+1] = ptr1[i+1];
      ptr2[i+2] = ptr1[i+2];
      ptr2[i+3] = ptr1[i+3];
      ptr2[i+4] = ptr1[i+4];
      ptr2[i+5] = ptr1[i+5];
      ptr2[i+6] = ptr1[i+6];
   }
   return;
}



/********************************
 *                              *
 *         d3db::rr_SMul         *
 *                              *
 ********************************/
void d3db::rr_SMul(const double da, const double *ptr1, double *ptr2)
{
   int i;
   int m = n2ft3d%5;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr2[i] = da*ptr1[i];
   if (n2ft3d<5)
      return;
   for (i=m; i<n2ft3d; i+=5)
   {
      ptr2[i]   = da*ptr1[i];
      ptr2[i+1] = da*ptr1[i+1];
      ptr2[i+2] = da*ptr1[i+2];
      ptr2[i+3] = da*ptr1[i+3];
      ptr2[i+4] = da*ptr1[i+4];
   }
   return;
}

/********************************
 *                              *
 *        d3db::rrr_SMulAdd     *
 *                              *
 ********************************/
void d3db::rrr_SMulAdd(const double da, const double *ptr1, const double *ptr2, double *ptr3)
{
   int i;
   int m = n2ft3d%5;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr3[i] = da*ptr1[i] + ptr2[i];
   if (n2ft3d<5)
      return;
   for (i=m; i<n2ft3d; i+=5)
   {
      ptr3[i]   = da*ptr1[i]   + ptr2[i];
      ptr3[i+1] = da*ptr1[i+1] + ptr2[i];
      ptr3[i+2] = da*ptr1[i+2] + ptr2[i];
      ptr3[i+3] = da*ptr1[i+3] + ptr2[i];
      ptr3[i+4] = da*ptr1[i+4] + ptr2[i];
   }
   return;
}




/********************************
 *                              *
 *          d3db::r_SMul        *
 *                              *
 ********************************/
void d3db::r_SMul(const double da, double *ptr2)
{
   int i;
   int m = n2ft3d%5;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr2[i] *= da;
   if (n2ft3d<5)
      return;
   for (i=m; i<n2ft3d; i+=5)
   {
      ptr2[i]   *= da;
      ptr2[i+1] *= da;
      ptr2[i+2] *= da;
      ptr2[i+3] *= da;
      ptr2[i+4] *= da;
   }
   return;
}

/********************************
 *                              *
 *          d3db::r_abs         *
 *                              *
 ********************************/
void d3db::r_abs(double *ptr2)
{
   int i;
   int m = n2ft3d%5;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr2[i] = std::abs(ptr2[i]);
   if (n2ft3d<5)
      return;
   for (i=m; i<n2ft3d; i+=5)
   {
      ptr2[i]   = std::abs(ptr2[i]);
      ptr2[i+1] = std::abs(ptr2[i+1]);
      ptr2[i+2] = std::abs(ptr2[i+2]);
      ptr2[i+3] = std::abs(ptr2[i+3]);
      ptr2[i+4] = std::abs(ptr2[i+4]);
   }
   return;
}

/********************************
 *                              *
 *          d3db::r_sqr         *
 *                              *
 ********************************/
void d3db::r_sqr(double *ptr2)
{
   int i;
   int m = n2ft3d%5;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr2[i] *= ptr2[i];
   if (n2ft3d<5)
      return;
   for (i=m; i<n2ft3d; i+=5)
   {
      ptr2[i]   *= ptr2[i];
      ptr2[i+1] *= ptr2[i+1];
      ptr2[i+2] *= ptr2[i+2];
      ptr2[i+3] *= ptr2[i+3];
      ptr2[i+4] *= ptr2[i+4];
   }      
   return;
}


/********************************
 *                              *
 *         d3db::r_dsum        *
 *                              *
 ********************************/
double d3db::r_dsum(const double *ptr)
{
   int i;
   int m = n2ft3d%5;
   double sum=0.0;
   if (m>0)
      for (i=0; i<m; ++i)
         sum += ptr[i];
   if (n2ft3d<5)
      return sum;
   for (i=m; i<n2ft3d; i+=5)
   {
      sum += ptr[i]
          +  ptr[i+1]
          +  ptr[i+2]
          +  ptr[i+3]
          +  ptr[i+4];
   }
   return sum;
}

/********************************
 *                              *
 *         d3db::rrr_Sum        *
 *                              *
 ********************************/
void d3db::rrr_Sum(const double *ptr1, const double *ptr2, double *ptr3)
{
   int i;
   int m = n2ft3d%5;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr3[i] = ptr1[i] + ptr2[i];
   if (n2ft3d<5)
      return;
   for (i=m; i<n2ft3d; i+=5)
   {
      ptr3[i]   = ptr1[i]   + ptr2[i];
      ptr3[i+1] = ptr1[i+1] + ptr2[i+1];
      ptr3[i+2] = ptr1[i+2] + ptr2[i+2]; 
      ptr3[i+3] = ptr1[i+3] + ptr2[i+3];
      ptr3[i+4] = ptr1[i+4] + ptr2[i+4];
   }
   return;
}

/********************************
 *                              *
 *         d3db::rr_Sum         *
 *                              *
 ********************************/
void d3db::rr_Sum(const double *ptr2, double *ptr3)
{
   int i;
   int m = n2ft3d%5;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr3[i] += ptr2[i];
   if (n2ft3d<5)
      return;
   for (i=m; i<n2ft3d; i+=5)
   {
      ptr3[i]   += ptr2[i];
      ptr3[i+1] += ptr2[i+1];
      ptr3[i+2] += ptr2[i+2]; 
      ptr3[i+3] += ptr2[i+3];
      ptr3[i+4] += ptr2[i+4];
   }
   return;
}



/********************************
 *                              *
 *         d3db::rrr_Mul        *
 *                              *
 ********************************/
void d3db::rrr_Mul(const double *ptr1, const double *ptr2, double *ptr3)
{
   int i;
   int m = n2ft3d%5;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr3[i] = ptr1[i] * ptr2[i];
   if (n2ft3d<5)
      return;
   for (i=m; i<n2ft3d; i+=5)
   {
      ptr3[i]   = ptr1[i]   * ptr2[i];
      ptr3[i+1] = ptr1[i+1] * ptr2[i+1];
      ptr3[i+2] = ptr1[i+2] * ptr2[i+2];
      ptr3[i+3] = ptr1[i+3] * ptr2[i+3];
      ptr3[i+4] = ptr1[i+4] * ptr2[i+4];
   }
   return;
}



/********************************
 *                              *
 *         d3db::rr_dot        *
 *                              *
 ********************************/
double d3db::rr_dot(const double *ptr1, const double *ptr2) 
{
   int i;
   int m = n2ft3d%5;
   double sum=0.0;
   if (m>0)
      for (i=0; i<m; ++i)
         sum += ptr1[i]*ptr2[i];
   if (n2ft3d<5)
      return sum;
   for (i=m; i<n2ft3d; i+=5)
   {
      sum += ptr1[i]   * ptr2[i] 
          +  ptr1[i+1] * ptr2[i+1] 
          +  ptr1[i+2] * ptr2[i+2] 
          +  ptr1[i+3] * ptr2[i+3] 
          +  ptr1[i+4] * ptr2[i+4];
   }
   return parall->SumAll(1,sum);
}


/********************************
 *                              *
 *         d3db::t_alloc        *
 *                              *
 ********************************/
double * d3db::t_alloc()
{
   double *ptr = new double[nfft3d];
   return ptr;
}

/********************************
 *                              *
 *         d3db::t_dealloc      *
 *                              *
 ********************************/
void d3db::t_dealloc(double *ptr)
{
   delete [] ptr;
}


/********************************
 *                              *
 *         d3db::c_read         *
 *                              *
 ********************************/
void d3db::c_read(const int iunit, double *a, const int jcol)
{
   int jstart,jend,fillcolumn,index,ii,jj,p_to,p_here;
   int taskid   = parall->taskid();
   int taskid_j = parall->taskid_j();
   int np_j     = parall->np_j();
   
   if (jcol<0) 
   {
      jstart = 0;
      jend   = np_j-1;
      fillcolumn = 1;
   }
   else
   {
      jstart = jend = jcol;
      fillcolumn = (taskid_j==jcol);
   }

   /**********************
    **** slab mapping ****
    **********************/
   if (maptype==1)
   { 
      double *tmp = new double[(nx+2)*ny];
      int   bsize = (nx+2)*ny;
   
      /**** master node reads from file and distributes ****/
      if (taskid==MASTER)
         for (int k=0; k<nz; ++k)
         {
            dread(iunit,tmp,bsize);

            index = 2*ijktoindex(0,0,k);
            ii    = ijktop(0,0,k);
            for (jj=jstart; jj<=jend; ++jj)
            {
               p_to = parall->convert_taskid_ij(ii,jj);
               if (p_to==MASTER)
                  for (int k=0; k<bsize; ++k) a[index+k] = tmp[k];
               else
                  parall->dsend(0,9,p_to,bsize,tmp);
            }
         }
      
      /**** not master node ****/
      else if (fillcolumn) 
         for (int k=0; k<nz; ++k)
         {
            index = 2*ijktoindex(0,0,k);
            ii    = ijktop(0,0,k);
            p_here = parall->convert_taskid_ij(ii,taskid_j);
            if (p_here==taskid)
            { 
               parall->dreceive(0,9,MASTER,bsize,tmp);
               for (int k=0; k<bsize; ++k) a[index+k] = tmp[k];
            }
         }
	  
      delete [] tmp;
   }
   
   /*************************
    **** hilbert mapping ****
    *************************/
   else
   {
      double *tmp = new double[nx+2];
      int bsize = (nx+2);
         
      /**** master node reads from file and distributes ****/
      if (taskid==MASTER)
         for (int k=0; k<nz; ++k)
         for (int j=0; j<ny; ++j)
         {

            dread(iunit,tmp,bsize);

            index  = ijktoindex2(0,j,k);
            ii     = ijktop2(0,j,k);
            for (int jj=jstart; jj<=jend; ++jj)
            {
               p_to = parall->convert_taskid_ij(ii,jj);

               if (p_to==MASTER) 
                  for (int k=0; k<bsize; ++k) a[index+k] = tmp[k];
               else
                  parall->dsend(0,9,p_to,bsize,tmp);
            }
         }
      
      /**** not master node ****/
      else if (fillcolumn) 
         for (int k=0; k<nz; ++k)
         for (int j=0; j<ny; ++j)
         {
            index  = ijktoindex2(0,j,k);
            ii     = ijktop2(0,j,k);
            p_here = parall->convert_taskid_ij(ii,taskid_j);
            if (p_here==taskid) 
            {
               parall->dreceive(0,9,MASTER,bsize,tmp);
               for (int k=0; k<bsize; ++k) a[index+k] = tmp[k];
            }
         }

      delete [] tmp;

      double *tmp1 = new double[2*nfft3d];
      double *tmp2 = new double[2*nfft3d];
      c_transpose_ijk(4,a,tmp1,tmp2);
      delete [] tmp2;
      delete [] tmp1;
   }

}


/********************************
 *                              *
 *         d3db::c_write        *
 *                              *
 ********************************/
void d3db::c_write(const int iunit, double *a, const int jcol)
{
   int index,ii,jj,p_from,p_here;
   int taskid   = parall->taskid();
   int taskid_j = parall->taskid_j();
   int np_j     = parall->np_j();

   /**********************
    **** slab mapping ****
    **********************/
   if (maptype==1)
   {
      double *tmp = new double[(nx+2)*ny];
      int   bsize = (nx+2)*ny;

      /**** master node gathers and write to file ****/
      if (taskid==MASTER)
         for (int k=0; k<nz; ++k)
         {
            ii    = ijktop(0,0,k);
            p_from = parall->convert_taskid_ij(ii,jcol);
            if (p_from==MASTER) 
            {
               index = 2*ijktoindex(0,0,k);
               for (int k=0; k<bsize; ++k) a[index+k] = tmp[k];
               dwrite(iunit,tmp,bsize);
            }
            else
            {
                parall->dreceive(0,9,p_from,bsize,tmp);
            }
            dwrite(iunit,tmp,bsize);

         }

      /**** not master node ****/
      else
         for (int k=0; k<nz; ++k)
         {
            index = 2*ijktoindex(0,0,k);
            ii    = ijktop(0,0,k);
            p_here = parall->convert_taskid_ij(ii,taskid_j);
            if (p_here==taskid)
            {
               for (int k=0; k<bsize; ++k) tmp[k] = a[index+k];
               parall->dsend(0,9,MASTER,bsize,tmp);
            }
         }

      delete [] tmp;
   }

   /*************************
    **** hilbert mapping ****
    *************************/
   else
   {  
      double *tmp1 = new double[2*nfft3d];
      double *tmp2 = new double[2*nfft3d];
      c_transpose_ijk(5,a,tmp1,tmp2);
      delete [] tmp2;
      delete [] tmp1;

      double *tmp = new double[nx+2];
      int bsize = (nx+2);
         
      /**** master node write to file and fetches from other nodes ****/
      if (taskid==MASTER)
         for (int k=0; k<nz; ++k)
         for (int j=0; j<ny; ++j)
         {
            ii     = ijktop2(0,j,k);
            p_from = parall->convert_taskid_ij(ii,jcol);
            if (p_from==MASTER) 
            {
               index  = ijktoindex2(0,j,k);
               for (int k=0; k<bsize; ++k) tmp[k] = a[index+k];
            }
            else
            {
               parall->dreceive(0,9,p_from,bsize,tmp);
            }
            dwrite(iunit,tmp,bsize);
         }
      
      /**** not master node ****/
      else
         for (int k=0; k<nz; ++k)
         for (int j=0; j<ny; ++j)
         {  
            ii     = ijktop2(0,j,k);
            p_here = parall->convert_taskid_ij(ii,taskid_j);
            if (p_here==taskid)
            {  
               index  = ijktoindex2(0,j,k);
               for (int k=0; k<bsize; ++k) tmp[k] = a[index+k];
               parall->dsend(0,9,MASTER,bsize,tmp);
            }
         }
      
      delete [] tmp;
      
   }
}



static void cshift1_fftb(const int n1,const int n2, const int n3, const int n4, double *a)
{
   int i,j,indx;
   indx = 1;
   for (j=0; j<(n2*n3*n4); ++j)
   {
      for (i=2; i<=n1; ++i) 
      {
         a[indx+i-2] = a[indx+i-1];
      }
      indx += (n1+2);
   }
}
static void cshift_fftf(const int n1,const int n2, const int n3, const int n4, double *a)
{
   int i,j,indx;
   indx = 1;
   for (j=0; j<(n2*n3*n4); ++j)
   {
      for (i=n1; i>=2; --i) 
      {
         a[indx+i-1] = a[indx+i-2];
      }
      a[indx+1-1]    = 0.0;
      a[indx+n1+1-1] = 0.0;
      indx += (n1+2);
   }
}
static void cshift_fftf_ab(const int n1,const int n2, const int n3, const int n4, double *a, double *b)
{
   int i,j,indx;
   indx = 0;
   for (j=0; j<(n2*n3*n4); ++j)
   {
      for (i=n1; i>0; --i) 
      {
         b[indx+i] = a[indx+i-1];
      }
      b[indx+1]    = 0.0;
      b[indx+n1+1] = 0.0;
      indx += (n1+2);
   }
}
static void zeroend_fftb(const int n1,const int n2, const int n3, const int n4, double *a)
{
   int i,indx;
   indx = n1+1;
   for (i=0; i<(n2*n3*n4); ++i)
   {
      a[indx-1]   = 0.0;
      a[indx+1-1] = 0.0;
      indx += (n1+2);
   }
}


/********************************
 *                              *
 *         d3db::cr_fft3d       *
 *                              *
 ********************************/
void d3db::cr_fft3d(double *a)
{
   int i,j,k,jj,kk,q,indx,indx0,nxh,nxh2,nxhy,nxhy2,nxhz,nxhz2;
   double *tmp2,*tmp3;
   
   nxh  = nx/2+1;
   nxhy = nxh*ny;
   nxhz = nxh*nz;
   nxh2  = nx+2;
   nxhy2 = nxh2*ny;
   nxhz2 = nxh2*nz;

   tmp2 = new double[2*nfft3d];
   tmp3 = new double[2*nfft3d];

   /**********************
    **** slab mapping ****
    **********************/
   if (maptype==1)
   {
      /***************************************************
       ***     do fft along kz dimension               ***
       ***     A(kx,nz,ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
       ***************************************************/
       indx0=0;
       for (q=0; q<nq; ++q)
       {
          for (i=0; i<nxh; ++i)
          {
             kk   = 0; 
             indx = 2*i+indx0;
             for (k=0; k<nz; ++k)
             {
                tmp2[kk]   = a[indx];
                tmp2[kk+1] = a[indx+1];
                kk   += 2;
                indx += nxh2;
             }

             dcfftb_(&nz,tmp2,tmpz);

             kk   = 0; 
             indx = 2*i+indx0;
             for (k=0; k<nz; ++k)
             {
                a[indx]   = tmp2[kk];
                a[indx+1] = tmp2[kk+1];
                kk   += 2;
                indx += nxh2;
             }
          }
          indx0 += nxhz2;
       }

      /***********************************************
       ***         Do a transpose of A             ***
       ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
       ************************************************/
       c_transpose_jk(a,tmp2,tmp3);

      /*************************************************
       ***        do fft along ky dimension          ***
       ***    A(kx,ny,nz) <- fft1d^(-1)[A(kx,ky,nz)] ***
       *************************************************/
       indx0=0;
       for (q=0; q<nq; ++q)
       {
          for (i=0; i<nxh; ++i)
          {
             jj   = 0;
             indx = 2*i+indx0;
             for (j=0; j<ny; ++j)
             {
                tmp2[jj]   = a[indx];
                tmp2[jj+1] = a[indx+1];
                jj   += 2;
                indx += nxh2;
             }

             dcfftb_(&ny,tmp2,tmpy);

             jj   = 0; 
             indx = 2*i+indx0;
             for (j=0; j<ny; ++j)
             {
                a[indx]   = tmp2[jj];
                a[indx+1] = tmp2[jj+1];
                jj   += 2;
                indx += nxh2;
             }
          }
          indx0 += nxhy2;
       }


      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
       cshift1_fftb(nx,ny,nq,1,a);
       indx = 0;
       for (q=0; q<nq; ++q)
       for (j=0; j<ny; ++j)
       {
          drfftb_(&nx,&a[indx],tmpx);
          indx += nxh2;
       }
       zeroend_fftb(nx,ny,nq,1,a);


   }
   /*************************
    **** hilbert mapping ****
    *************************/
   else
   {

      /************************************************
       ***     do fft along kz dimension            ***
       ***   A(nz,kx,ky) <- fft1d^(-1)[A(kz,kx,ky)] ***
       ************************************************/
       indx = 0;
       for (q=0; q<nq3; ++q)
       {
          dcfftb_(&nz,&a[indx],tmpz);
          indx += 2*nz;
       }
       c_transpose_ijk(2,a,tmp2,tmp3);

      /************************************************
       ***     do fft along ky dimension            ***
       ***   A(ny,nz,kx) <- fft1d^(-1)[A(ky,nz,kx)] ***
       ************************************************/
       indx = 0;
       for (q=0; q<nq2; ++q)
       {
          dcfftb_(&ny,&a[indx],tmpy);
          indx += (2*ny);
       }
       c_transpose_ijk(3,a,tmp2,tmp3);

      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
       cshift1_fftb(nx,nq1,1,1,a);
       indx = 0;
       for (q=0; q<nq1; ++q)
       {
          drfftb_(&nx,&a[indx],tmpx);
          indx += nxh2;
       }
       zeroend_fftb(nx,nq1,1,1,a);

   }

   delete [] tmp3;
   delete [] tmp2;
}

/********************************
 *                              *
 *         d3db::rc_fft3d       *
 *                              *
 ********************************/
void d3db::rc_fft3d(double *a)
{
   int i,j,k,jj,kk,q,indx,indx0,nxh,nxh2,nxhy,nxhy2,nxhz,nxhz2;
   double *tmp2,*tmp3;
   
   nxh  = nx/2+1;
   nxhy = nxh*ny;
   nxhz = nxh*nz;
   nxh2  = nx+2;
   nxhy2 = nxh2*ny;
   nxhz2 = nxh2*nz;

   tmp2 = new double[2*nfft3d];
   tmp3 = new double[2*nfft3d];

   /**********************
    **** slab mapping ****
    **********************/
   if (maptype==1)
   {


      /********************************************
       ***     do fft along nx dimension        ***
       ***   A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]  ***
       ********************************************/
       indx = 0;
       for (q=0; q<nq; ++q)
       for (j=0; j<ny; ++j)
       {
          drfftf_(&nx,&a[indx],tmpx);
          indx += nxh2;
       }
       cshift_fftf(nx,ny,nq,1,a);



      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
       ********************************************/
       indx0=0;
       for (q=0; q<nq; ++q)
       {
          for (i=0; i<nxh; ++i)
          {
             jj   = 0;
             indx = 2*i+indx0;
             for (j=0; j<ny; ++j)
             {
                tmp2[jj]   = a[indx];
                tmp2[jj+1] = a[indx+1];
                jj   += 2;
                indx += nxh2;
             }

             dcfftf_(&ny,tmp2,tmpy);

             jj   = 0; 
             indx = 2*i+indx0;
             for (j=0; j<ny; ++j)
             {
                a[indx]   = tmp2[jj];
                a[indx+1] = tmp2[jj+1];
                jj   += 2;
                indx += nxh2;
             }
          }
          indx0 += nxhy2;
       }


      /********************************************
       ***         Do a transpose of A          ***
       ***      A(ky,nz,ky) <- A(kx,ky,nz)      ***
       ********************************************/
       c_transpose_jk(a,tmp2,tmp3);



      /********************************************
       ***     do fft along nz dimension        ***
       ***   A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]  ***
       ********************************************/
       indx0=0;
       for (q=0; q<nq; ++q)
       {
          for (i=0; i<nxh; ++i)
          {
             kk   = 0; 
             indx = 2*i+indx0;
             for (k=0; k<nz; ++k)
             {
                tmp2[kk]   = a[indx];
                tmp2[kk+1] = a[indx+1];
                kk   += 2;
                indx += nxh2;
             }

             dcfftf_(&nz,tmp2,tmpz);

             kk   = 0; 
             indx = 2*i+indx0;
             for (k=0; k<nz; ++k)
             {
                a[indx]   = tmp2[kk];
                a[indx+1] = tmp2[kk+1];
                kk   += 2;
                indx += nxh2;
             }
          }
          indx0 += nxhz2;
       }


   }
   /*************************
    **** hilbert mapping ****
    *************************/
   else
   {
      /********************************************
       ***     do fft along nx dimension        ***
       ***   A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]  ***
       ********************************************/
       indx = 0;
       for (q=0; q<nq1; ++q)
       {
          drfftf_(&nx,&a[indx],tmpx);
          indx += nxh2;
       }
       cshift_fftf(nx,nq1,1,1,a);
       c_transpose_ijk(0,a,tmp2,tmp3);

      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
       ********************************************/
       indx = 0;
       for (q=0; q<nq2; ++q)
       {
          dcfftf_(&ny,&a[indx],tmpy);
          indx += (2*ny);
       }
       c_transpose_ijk(1,a,tmp2,tmp3);

      /********************************************
       ***     do fft along nz dimension        ***
       ***   A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]  ***
       ********************************************/
       indx = 0;
       for (q=0; q<nq3; ++q)
       {
          dcfftf_(&nz,&a[indx],tmpz);
          indx += 2*nz;
       }
   }

   delete [] tmp3;
   delete [] tmp2;
}





/********************************
 *                              *
 *         d3db::t_read         *
 *                              *
 ********************************/
void d3db::t_read(const int iunit, double *a, const int jcol)
{
   int jstart,jend,fillcolumn,index,ii,jj,p_to,p_here;
   int taskid   = parall->taskid();
   int taskid_j = parall->taskid_j();
   int np_j     = parall->np_j();
   
   if (jcol<0) 
   {
      jstart = 0;
      jend   = np_j-1;
      fillcolumn = 1;
   }
   else
   {
      jstart = jend = jcol;
      fillcolumn = (taskid_j==jcol);
   }

   /**********************
    **** slab mapping ****
    **********************/
   if (maptype==1)
   { 
      double *tmp = new double[(nx/2+1)*ny];
      int   bsize = (nx/2+1)*ny;
   
      /**** master node reads from file and distributes ****/
      if (taskid==MASTER)
         for (int k=0; k<nz; ++k)
         {
            dread(iunit,tmp,bsize);

            index = ijktoindex(0,0,k);
            ii    = ijktop(0,0,k);
            for (jj=jstart; jj<=jend; ++jj)
            {
               p_to = parall->convert_taskid_ij(ii,jj);

               if (p_to==MASTER)
                  for (int k=0; k<bsize; ++k) a[index+k] = tmp[k];
               else
                  parall->dsend(0,9,p_to,bsize,tmp);
            }
         }
      
      /**** not master node ****/
      else if (fillcolumn) 
         for (int k=0; k<nz; ++k)
         {
            index = ijktoindex(0,0,k);
            ii    = ijktop(0,0,k);
            p_here = parall->convert_taskid_ij(ii,taskid_j);
            if (p_here==taskid)
            { 
               parall->dreceive(0,9,MASTER,bsize,tmp);
               for (int k=0; k<bsize; ++k) a[index+k] = tmp[k];
            }
         }
	  
      delete [] tmp;
   }
   
   /*************************
    **** hilbert mapping ****
    *************************/
   else
   {
      double *tmp = new double[nx/2+1];
      int bsize = (nx/2+1);
         
      /**** master node reads from file and distributes ****/
      if (taskid==MASTER)
         for (int k=0; k<nz; ++k)
         for (int j=0; j<ny; ++j)
         {

            dread(iunit,tmp,bsize);

            index  = ijktoindex2t(0,j,k);
            ii     = ijktop2(0,j,k);
            for (int jj=jstart; jj<=jend; ++jj)
            {
               p_to = parall->convert_taskid_ij(ii,jj);

               if (p_to==MASTER) 
                  for (int k=0; k<bsize; ++k) a[index+k] = tmp[k];
               else
                  parall->dsend(0,9,p_to,bsize,tmp);
            }
         }
      
      /**** not master node ****/
      else if (fillcolumn) 
         for (int k=0; k<nz; ++k)
         for (int j=0; j<ny; ++j)
         {
            index  = ijktoindex2t(0,j,k);
            ii     = ijktop2(0,j,k);
            p_here = parall->convert_taskid_ij(ii,taskid_j);
            if (p_here==taskid) 
            {
               parall->dreceive(0,9,MASTER,bsize,tmp);
               for (int k=0; k<bsize; ++k) a[index+k] = tmp[k];
            }
         }
      delete [] tmp;

      double *tmp1 = new double[2*nfft3d];
      double *tmp2 = new double[2*nfft3d];
      t_transpose_ijk(4,a,tmp1,tmp2);
      delete [] tmp2;
      delete [] tmp1;
   }
}


/********************************
 *                              *
 *         d3db::t_write        *
 *                              *
 ********************************/
void d3db::t_write(const int iunit, double *a, const int jcol)
{
   int index,ii,jj,p_from,p_here;
   int taskid   = parall->taskid();
   int taskid_j = parall->taskid_j();
   int np_j     = parall->np_j();

   /**********************
    **** slab mapping ****
    **********************/
   if (maptype==1)
   {
      double *tmp = new double[(nx+2)*ny];
      int   bsize = (nx/2+1)*ny;

      /**** master node gathers and write to file ****/
      if (taskid==MASTER)
         for (int k=0; k<nz; ++k)
         {
            ii    = ijktop(0,0,k);
            p_from = parall->convert_taskid_ij(ii,jcol);
            if (p_from==MASTER) 
            {
               index = ijktoindex(0,0,k);
               for (int k=0; k<bsize; ++k) a[index+k] = tmp[k];
               dwrite(iunit,tmp,bsize);
            }
            else
            {
                parall->dreceive(0,9,p_from,bsize,tmp);
            }
            dwrite(iunit,tmp,bsize);

         }

      /**** not master node ****/
      else
         for (int k=0; k<nz; ++k)
         {
            index = ijktoindex(0,0,k);
            ii    = ijktop(0,0,k);
            p_here = parall->convert_taskid_ij(ii,taskid_j);
            if (p_here==taskid)
            {
               for (int k=0; k<bsize; ++k) tmp[k] = a[index+k];
               parall->dsend(0,9,MASTER,bsize,tmp);
            }
         }

      delete [] tmp;
   }

   /*************************
    **** hilbert mapping ****
    *************************/
   else
   {  
      double *tmp1 = new double[2*nfft3d];
      double *tmp2 = new double[2*nfft3d];
      t_transpose_ijk(5,a,tmp1,tmp2);
      delete [] tmp2;
      delete [] tmp1;

      double *tmp = new double[nx/2+1];
      int bsize = (nx/2+1);
         
      /**** master node write to file and fetches from other nodes ****/
      if (taskid==MASTER)
         for (int k=0; k<nz; ++k)
         for (int j=0; j<ny; ++j)
         {
            ii     = ijktop2(0,j,k);
            p_from = parall->convert_taskid_ij(ii,jcol);
            if (p_from==MASTER) 
            {
               index  = ijktoindex2t(0,j,k);
               for (int k=0; k<bsize; ++k) tmp[k] = a[index+k];
            }
            else
            {
               parall->dreceive(0,9,p_from,bsize,tmp);
            }
            dwrite(iunit,tmp,bsize);
         }
      
      /**** not master node ****/
      else
         for (int k=0; k<nz; ++k)
         for (int j=0; j<ny; ++j)
         {  
            ii     = ijktop2(0,j,k);
            p_here = parall->convert_taskid_ij(ii,taskid_j);
            if (p_here==taskid)
            {  
               index  = ijktoindex2t(0,j,k);
               for (int k=0; k<bsize; ++k) tmp[k] = a[index+k];
               parall->dsend(0,9,MASTER,bsize,tmp);
            }
         }
      
      delete [] tmp;
      
   }
}




/********************************
 *                              *
 *    d3db::c_transpose_jk      *
 *                              *
 ********************************/
void d3db::c_transpose_jk(double *a, double *tmp1, double *tmp2)
{
   int it,proc_from,proc_to;
   int one=1;
   int msglen;

   parall->astart(1,np);

   c_bindexcopy(nfft3d,iq_to_i1[0],a,tmp1);

   /* it = 0, transpose data on same thread */
   msglen = 2*(i2_start[0][1] - i2_start[0][0]);
   dcopy_(&msglen,
          &(tmp1[2*i1_start[0][0]]),&one,
          &(tmp2[2*i2_start[0][0]]),&one);

   /* receive packed array data */
   for (it=1; it<np; ++it)
   {
      /* synchronous receive of tmp */
      proc_from = (taskid-it+np)%np;
      msglen = 2*(i2_start[0][it+1] - i2_start[0][it]);
      if (msglen>0)
         parall->adreceive(1,1,proc_from,msglen,&tmp2[2*i2_start[0][it]]);
   }
   for (it=1; it<np; ++it)
   {
      proc_to = (taskid+it)%np;
      msglen = 2*(i1_start[0][it+1] - i1_start[0][it]);
      if (msglen>0)
         parall->dsend(1,1,proc_to,msglen,&tmp1[2*i1_start[0][it]]);
   }
   parall->aend(1);
   c_aindexcopy(nfft3d,iq_to_i2[0],tmp2,a);
}

/********************************
 *                              *
 *    d3db::t_transpose_jk      *
 *                              *
 ********************************/
void d3db::t_transpose_jk(double *a, double *tmp1, double *tmp2)
{
   int it,proc_from,proc_to;
   int one=1;
   int msglen;

   parall->astart(1,np);

   t_bindexcopy(nfft3d,iq_to_i1[0],a,tmp1);

   /* it = 0, transpose data on same thread */
   msglen = (i2_start[0][1] - i2_start[0][0]);
   dcopy_(&msglen,
          &(tmp1[i1_start[0][0]]),&one,
          &(tmp2[i2_start[0][0]]),&one);

   /* receive packed array data */
   for (it=1; it<np; ++it)
   {
      /* synchronous receive of tmp */
      proc_from = (taskid-it+np)%np;
      msglen = (i2_start[0][it+1] - i2_start[0][it]);
      if (msglen>0)
         parall->adreceive(1,1,proc_from,msglen,&tmp2[i2_start[0][it]]);
   }
   for (it=1; it<np; ++it)
   {
      proc_to = (taskid+it)%np;
      msglen = (i1_start[0][it+1] - i1_start[0][it]);
      if (msglen>0)
         parall->dsend(1,1,proc_to,msglen,&tmp1[i1_start[0][it]]);
   }
   parall->aend(1);
   c_aindexcopy(nfft3d,iq_to_i2[0],tmp2,a);
}




/********************************
 *                              *
 *    d3db::c_transpose_ijk     *
 *                              *
 ********************************/
void d3db::c_transpose_ijk(const int op,double *a,double *tmp1,double *tmp2)
{
   int nnfft3d,it,proc_from,proc_to;
   int one=1;
   int msglen;

   parall->astart(1,np);

   /* pack a array */
   if ((op==0)||(op==4)) nnfft3d = (nx/2+1)*nq1;
   if ((op==1)||(op==3)) nnfft3d = (ny)    *nq2;
   if ((op==2)||(op==5)) nnfft3d = (nz)    *nq3;
   c_bindexcopy(nnfft3d,iq_to_i1[op],a,tmp1);

   /* it = 0, transpose data on same thread */
   msglen = 2*(i2_start[op][1] - i2_start[op][0]);
   dcopy_(&msglen,
          &(tmp1[2*i1_start[op][0]]),&one,
          &(tmp2[2*i2_start[op][0]]),&one);


   /* receive packed array data */
   for (it=1; it<np; ++it)
   {
      /* synchronous receive of tmp */
      proc_from = (taskid-it+np)%np;
      msglen = 2*(i2_start[op][it+1] - i2_start[op][it]);
      if (msglen>0)
         parall->adreceive(1,1,proc_from,msglen,&tmp2[2*i2_start[op][it]]);
   }
   for (it=1; it<np; ++it)
   {
      proc_to = (taskid+it)%np;
      msglen = 2*(i1_start[op][it+1] - i1_start[op][it]);
      if (msglen>0)
         parall->dsend(1,1,proc_to,msglen,&tmp1[2*i1_start[op][it]]);
   }

   /* wait for completion of mp_send, also do a sync */
   parall->aend(1);

   /* unpack a array */
   if ((op==3)||(op==5)) nnfft3d = (nx/2+1)*nq1;
   if ((op==0)||(op==2)) nnfft3d = (ny)    *nq2;
   if ((op==1)||(op==4)) nnfft3d = (nz)    *nq3;
   c_aindexcopy(nnfft3d,iq_to_i2[op],tmp2,a);

}


/********************************
 *                              *
 *    d3db::t_transpose_ijk     *
 *                              *
 ********************************/
void d3db::t_transpose_ijk(const int op,double *a,double *tmp1,double *tmp2)
{
   int nnfft3d,it,proc_from,proc_to;
   int one=1;
   int msglen;

   parall->astart(1,np);

   /* pack a array */
   if ((op==0)||(op==4)) nnfft3d = (nx/2+1)*nq1;
   if ((op==1)||(op==3)) nnfft3d = (ny)    *nq2;
   if ((op==2)||(op==5)) nnfft3d = (nz)    *nq3;
   t_bindexcopy(nnfft3d,iq_to_i1[op],a,tmp1);


   /* it = 0, transpose data on same thread */
   msglen = (i2_start[op][1] - i2_start[op][0]);
   dcopy_(&msglen,
          &(tmp1[i1_start[op][0]]),&one,
          &(tmp2[i2_start[op][0]]),&one);


   /* receive packed array data */
   for (it=1; it<np; ++it)
   {
      /* synchronous receive of tmp */
      proc_from = (taskid-it+np)%np;
      msglen = (i2_start[op][it+1] - i2_start[op][it]);
      if (msglen>0)
         parall->adreceive(1,1,proc_from,msglen,&tmp2[i2_start[op][it]]);
   }
   for (it=1; it<np; ++it)
   {
      proc_to = (taskid+it)%np;
      msglen = (i1_start[op][it+1] - i1_start[op][it]);
      if (msglen>0)
         parall->dsend(1,1,proc_to,msglen,&tmp1[i1_start[op][it]]);
   }

   /* wait for completion of mp_send, also do a sync */
   parall->aend(1);

   /* unpack a array */
   if ((op==3)||(op==5)) nnfft3d = (nx/2+1)*nq1;
   if ((op==0)||(op==2)) nnfft3d = (ny)    *nq2;
   if ((op==1)||(op==4)) nnfft3d = (nz)    *nq3;
   t_aindexcopy(nnfft3d,iq_to_i2[op],tmp2,a);
}


/********************************
 *                              *
 *    d3db::timereverse_size    *
 *                              *
 ********************************/
int d3db::timereverse_size()
{
   int it,indx1,indx2,i2,i3,j2,j3,k2,k3;
   int proc_to,proc_from,phere,pto,nyh,nzh;
   int sz;

   nzh = nz/2;
   nyh = ny/2;
   indx1 = 0;
   indx2 = 0;
   for (it=0; it<np; ++it)
   {
      proc_to   = (taskid+it)%np;
      proc_from = (taskid-it+np)%np;

      /**** K=(0,0,k3)  ****/
      for (k3=1; k3<nzh; ++k3)
      {
         i3 =  k3;
         j3 = -k3;
         if (i3<0) i3 = i3 + nz;
         if (j3<0) j3 = j3 + nz;

         phere = ijktop(0,0,i3);
         pto    =ijktop(0,0,j3);

         /* packing scheme */
         if ((phere==taskid)&&(pto==proc_to)) ++indx1;

         /* unpacking scheme */
         if ((pto==taskid)&&(phere==proc_from)) ++indx2;
      }
   
      /**** K=(0,k2,k3)  ****/
      for (k3=(-nzh+1); k3<nzh; ++k3)
      for (k2=1; k2<nyh; ++k2)
      {
         i2 =  k2;
         i3 =  k3;
         j2 = -k2;
         j3 = -k3;
         if (i2<0) i2 = i2 + ny;
         if (i3<0) i3 = i3 + nz;
         if (j2<0) j2 = j2 + ny;
         if (j3<0) j3 = j3 + nz;

         phere = ijktop(0,i2,i3);
         pto    =ijktop(0,j2,j3);

         /* packing scheme */
         if ((phere==taskid)&&(pto==proc_to)) ++indx1;

         /* unpacking scheme */
         if ((pto==taskid)&&(phere==proc_from)) ++indx2;
      }
   }

   sz = indx1;
   if (sz < indx2) sz = indx2;

   return sz;
}


/********************************
 *                              *
 *      d3db::c_timereverse     *
 *                              *
 ********************************/
void d3db::c_timereverse(double *a, double *tmp1, double *tmp2)
{
   int nnfft3d,indx,it,proc_from,proc_to;
   int one=1;
   int msglen;

   parall->astart(1,np);

   indx    = t_i1_start[0];
   nnfft3d = (t_i1_start[np] - t_i1_start[0] + 0);
   c_aindexcopy(nnfft3d,&t_iq_to_i1[indx],a,&tmp1[indx]);

   /* it = 0, transpose data on same thread */
   msglen = 2*(t_i2_start[1] - t_i2_start[0]);
   dcopy_(&msglen,
          &(tmp1[2*t_i1_start[0]]),&one,
          &(tmp2[2*t_i2_start[0]]),&one);

   /* receive packed array data */
   for (it=1; it<np; ++it)
   {
      /* synchronous receive of tmp */
      proc_from = (taskid-it+np)%np;
      msglen = 2*(t_i2_start[it+1] - t_i2_start[it]);
      if (msglen>0)
         parall->adreceive(1,1,proc_from,msglen,&tmp2[2*t_i2_start[it]]);
   }
   for (it=1; it<np; ++it)
   {
      proc_to = (taskid+it)%np;
      msglen = 2*(t_i1_start[it+1] - t_i1_start[it]);
      if (msglen>0)
         parall->dsend(1,1,proc_to,msglen,&tmp1[2*t_i1_start[it]]);
   }
   parall->aend(1);

   indx    = t_i2_start[0];
   nnfft3d = (t_i2_start[np] - t_i2_start[0] + 0);
   c_bindexcopy_conjg(nnfft3d,&t_iq_to_i2[indx],&tmp2[indx],a);
}


/********************************
 *                              *
 *         d3db::c_setpw        *
 *                              *
 ********************************/
void d3db::c_setpw(const int filling[], const double *cvalue, double *a)
{
   int i = filling[0];
   int j = filling[1];
   int k = filling[2];

   int indx = ijktoindex(i,j,k);
   int p    = ijktop(i,j,k);
   if (p==parall->taskid_i())
   {
      a[2*indx]   = cvalue[0];
      a[2*indx+1] = cvalue[1];
   }

   /* set the conjugate on the i==0 plane */
   if ((i==0) &&  (j!=0) && (k!=0))
   {
      int jc = j; if (jc>(ny/2)) jc -= ny;
      jc = -jc;   if (jc<0)   jc += ny;

      int kc = k; if (kc>(nz/2)) kc -= nz;
      kc = -kc;   if (kc<0)   kc += nz;
      indx = ijktoindex(i,jc,kc);
      p    = ijktop(i,jc,kc);
      if (p==parall->taskid_i())
      {
         a[2*indx]   = cvalue[0];
         a[2*indx+1] = -cvalue[1];
      }
   }
}

/********************************
 *                              *
 *        d3db::c_addrandom     *
 *                              *
 ********************************/

void d3db::c_addrandom(double *a)
{
   double fac = 1.0/sqrt(1.0*nfft3d);
   for (auto i=0; i<n2ft3d; ++i)
      a[i] += fac*(0.50-util_random(0));
}

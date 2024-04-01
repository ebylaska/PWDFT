/* c3db.cpp
   Author - Eric Bylaska

        this class is used for defining 3d parallel maps
*/

/**
 * @class c3db
 * @brief Container class for distributed 3D blocks operations.
 */

#include "compressed_io.hpp"
#include "fft.h"
#include "util.hpp"
#include <cmath>
#include <vector>

#include "blas.h"

//#include "gdevice.hpp"
#include "nwpw_timing.hpp"

#include "c3db.hpp"

#include "iofmt.hpp"
#include <cstring>
//#include <math.h>

#define mytaskid 1

namespace pwdft {

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/
/**
 * @brief Constructor for the `c3db` class.
 *
 * This constructor initializes an instance of the `c3db` class with the given parameters.
 *
 * @param inparall A pointer to a `Parallel` object.
 * @param inmaptype An integer specifying the mapping type.
 * @param nx The number of grid points in the x-direction.
 * @param ny The number of grid points in the y-direction.
 * @param nz The number of grid points in the z-direction.
 */
c3db::c3db(Parallel *inparall, const int inmaptype, const int nx, const int ny, const int nz, const int nbrillq)
     : Mapping3c(inmaptype, inparall->np_i(), inparall->taskid_i(), nx, ny, nz)
{
   int index1, index2, proc_to, proc_from;
   int nyh, nzh;
   int phere, pto, pfrom;
   p_nbrillq0 = nbrillq;

   parall = inparall;

   p_iq_to_i1 = new (std::nothrow) int **[nbrillq+1]; 
   p_iq_to_i2 = new (std::nothrow) int **[nbrillq+1]; 
   p_iz_to_i2 = new (std::nothrow) int **[nbrillq+1]; 
   p_iz_to_i2_count = new (std::nothrow) int *[nbrillq+1]; 
   p_i1_start = new (std::nothrow) int **[nbrillq+1]; 
   p_i2_start = new (std::nothrow) int **[nbrillq+1]; 

   p_jq_to_i1 = new (std::nothrow) int **[nbrillq+1]; 
   p_jq_to_i2 = new (std::nothrow) int **[nbrillq+1]; 
   p_jz_to_i2 = new (std::nothrow) int **[nbrillq+1]; 
   p_j1_start = new (std::nothrow) int **[nbrillq+1]; 
   p_j2_start = new (std::nothrow) int **[nbrillq+1]; 



   if (maptype==1)
   {
      iq_to_i1 = new (std::nothrow) int *[1];
      iq_to_i1[0] = new (std::nothrow) int[nx * ny * nq]();
      iq_to_i2 = new (std::nothrow) int *[1];
      iq_to_i2[0] = new (std::nothrow) int[nx * ny * nq]();
      i1_start = new (std::nothrow) int *[1];
      i1_start[0] = new (std::nothrow) int[nz+1]();
      i2_start = new (std::nothrow) int *[1];
      i2_start[0] = new (std::nothrow) int[nz+1]();
      index1 = 0;
      index2 = 0;
      for (auto it = 0; it < np; ++it) {
        proc_to = (taskid + it) % np;
        proc_from = (taskid - it + np) % np;
        i1_start[0][it] = index1;
        i2_start[0][it] = index2;
     
        for (auto k = 0; k < nz; ++k)
          for (auto j = 0; j < ny; ++j) 
          {
             /* packing scheme */
             phere = cijktop(0, 0, k);
             pto = cijktop(0, 0, j);
             if ((phere == taskid) && (pto == proc_to))
                for (auto i = 0; i < (nx); ++i) {
                   iq_to_i1[0][cijktoindex(i, j, k)] = index1;
                   ++index1;
                }
            
             /* unpacking scheme */
             phere = cijktop(0, 0, j);
             pfrom = cijktop(0, 0, k);
             if ((phere == taskid) && (pfrom == proc_from))
                for (auto i = 0; i < (nx); ++i) 
                {
                   iq_to_i2[0][cijktoindex(i, k, j)] = index2;
                   ++index2;
                }
          }
      }
      i1_start[0][np] = index1;
      i2_start[0][np] = index2;
     
      /* allocate ptranspose indexes */
      for (auto nb = 0; nb < 2; ++nb) 
      {
         p_iz_to_i2_count[nb] = new (std::nothrow) int[6];
         p_iq_to_i1[nb] = new (std::nothrow) int *[1]();
         p_iq_to_i1[nb][0] = new (std::nothrow) int[(nx) * ny * nq]();
         p_iq_to_i2[nb] = new (std::nothrow) int *[1]();
         p_iq_to_i2[nb][0] = new (std::nothrow) int[(nx) * ny * nq]();
         p_iz_to_i2[nb] = new (std::nothrow) int *[1]();
         p_iz_to_i2[nb][0] = new (std::nothrow) int[(nx) * ny * nq]();
         p_i1_start[nb] = new (std::nothrow) int *[1]();
         p_i1_start[nb][0] = new (std::nothrow) int[nz + 1]();
         p_i2_start[nb] = new (std::nothrow) int *[1]();
         p_i2_start[nb][0] = new (std::nothrow) int[nz + 1]();
        
         p_jq_to_i1[nb] = new (std::nothrow) int *[1]();
         p_jq_to_i1[nb][0] = new (std::nothrow) int[(nx) * ny * nq]();
         p_jq_to_i2[nb] = new (std::nothrow) int *[1]();
         p_jq_to_i2[nb][0] = new (std::nothrow) int[(nx) * ny * nq]();
         p_jz_to_i2[nb] = new (std::nothrow) int *[1]();
         p_jz_to_i2[nb][0] = new (std::nothrow) int[(nx) * ny * nq]();
         p_j1_start[nb] = new (std::nothrow) int *[1]();
         p_j1_start[nb][0] = new (std::nothrow) int[nz + 1]();
         p_j2_start[nb] = new (std::nothrow) int *[1]();
         p_j2_start[nb][0] = new (std::nothrow) int[nz + 1]();
      }

   } 
   else 
   {
      iq_to_i1 = new (std::nothrow) int *[6]();
      iq_to_i1[0] = new (std::nothrow) int[(nx) * nq1]();
      iq_to_i1[1] = new (std::nothrow) int[ny * nq2]();
      iq_to_i1[2] = new (std::nothrow) int[nz * nq3]();
     
      iq_to_i1[3] = new (std::nothrow) int[ny * nq2]();
      iq_to_i1[4] = new (std::nothrow) int[(nx) * nq1]();
      iq_to_i1[5] = new (std::nothrow) int[nz * nq3]();
     
      iq_to_i2 = new (std::nothrow) int *[6]();
      iq_to_i2[0] = new (std::nothrow) int[ny * nq2]();
      iq_to_i2[1] = new (std::nothrow) int[nz * nq3]();
      iq_to_i2[2] = new (std::nothrow) int[ny * nq2]();
     
      iq_to_i2[3] = new (std::nothrow) int[(nx) * nq1]();
      iq_to_i2[4] = new (std::nothrow) int[nz * nq3]();
      iq_to_i2[5] = new (std::nothrow) int[(nx) * nq1]();
     
      i1_start = new (std::nothrow) int *[6]();
      for (auto i = 0; i < 6; ++i)
         i1_start[i] = new (std::nothrow) int[np + 1]();
     
      i2_start = new (std::nothrow) int *[6]();
      for (auto i = 0; i < 6; ++i)
         i2_start[i] = new (std::nothrow) int[np + 1]();
     
      /***********************************************************************/
      /* map1to2 mapping - done - tranpose operation #0  (j,k,i) <-- (i,j,k) */
      /***********************************************************************/
      index1 = 0;
      index2 = 0;
      for (auto it = 0; it < np; ++it) 
      {
         proc_to = (taskid + it) % np;
         proc_from = (taskid - it + np) % np;
         i1_start[0][it] = index1;
         i2_start[0][it] = index2;
         for (auto k = 0; k < nz; ++k)
           for (auto j = 0; j < ny; ++j)
             for (auto i = 0; i < (nx); ++i) 
             {
                phere = cijktop2(i, j, k);
                pto = cijktop1(i, j, k);
               
                /* packing scheme */
                if ((phere == taskid) && (pto == proc_to)) 
                {
                   iq_to_i1[0][cijktoindex2t(i, j, k)] = index1;
                   ++index1;
                }
                /* unpacking scheme */
                if ((pto == taskid) && (phere == proc_from)) 
                {
                   iq_to_i2[0][cijktoindex1(i, j, k)] = index2;
                   ++index2;
                }
             }
      }
      i1_start[0][np] = index1;
      i2_start[0][np] = index2;
     
      /***********************************************************************/
      /* map2to3 mapping - done - tranpose operation #1 (k,i,j) <-- (j,k,i)  */
      /***********************************************************************/
      index1 = 0;
      index2 = 0;
      for (auto it = 0; it < np; ++it) 
      {
         proc_to = (taskid + it) % np;
         proc_from = (taskid - it + np) % np;
         i1_start[1][it] = index1;
         i2_start[1][it] = index2;
         for (auto k = 0; k < nz; ++k)
           for (auto j = 0; j < ny; ++j)
             for (auto i = 0; i < (nx); ++i) 
             {
                phere = cijktop1(i, j, k);
                pto = cijktop(i, j, k);
               
                /* packing scheme */
                if ((phere == taskid) && (pto == proc_to)) 
                {
                   iq_to_i1[1][cijktoindex1(i, j, k)] = index1;
                   ++index1;
                }
                /* unpacking scheme */
                if ((pto == taskid) && (phere == proc_from)) 
                {
                   iq_to_i2[1][cijktoindex(i, j, k)] = index2;
                   ++index2;
                }
             }
      }
      i1_start[1][np] = index1;
      i2_start[1][np] = index2;
     
      /***********************************************************************/
      /* map3to2 mapping - done - tranpose operation #2  (j,k,i) <-- (k,i,j) */
      /***********************************************************************/
      index1 = 0;
      index2 = 0;
      for (auto it = 0; it < np; ++it) 
      {
         proc_to = (taskid + it) % np;
         proc_from = (taskid - it + np) % np;
         i1_start[2][it] = index1;
         i2_start[2][it] = index2;
         for (auto k = 0; k < nz; ++k)
           for (auto j = 0; j < ny; ++j)
             for (auto i = 0; i < (nx); ++i) 
             {
                phere = cijktop(i, j, k);
                pto = cijktop1(i, j, k);
               
                /* packing scheme */
                if ((phere == taskid) && (pto == proc_to)) 
                {
                   iq_to_i1[2][cijktoindex(i, j, k)] = index1;
                   ++index1;
                }
                /* unpacking scheme */
                if ((pto == taskid) && (phere == proc_from)) 
                {
                   iq_to_i2[2][cijktoindex1(i, j, k)] = index2;
                   ++index2;
                }
             }
      }
      i1_start[2][np] = index1;
      i2_start[2][np] = index2;
     
      /***********************************************************************/
      /* map2to1 mapping - done - tranpose operation #3  (i,j,k) <-- (j,k,i) */
      /***********************************************************************/
      index1 = 0;
      index2 = 0;
      for (auto it = 0; it < np; ++it) 
      {
         proc_to = (taskid + it) % np;
         proc_from = (taskid - it + np) % np;
         i1_start[3][it] = index1;
         i2_start[3][it] = index2;
         for (auto k = 0; k < nz; ++k)
           for (auto j = 0; j < ny; ++j)
             for (auto i = 0; i < (nx); ++i) 
             {
                phere = cijktop1(i, j, k);
                pto = cijktop2(i, j, k);
               
                /* packing scheme */
                if ((phere == taskid) && (pto == proc_to)) 
                {
                   iq_to_i1[3][cijktoindex1(i, j, k)] = index1;
                   ++index1;
                }
                /* unpacking scheme */
                if ((pto == taskid) && (phere == proc_from)) 
                {
                   iq_to_i2[3][cijktoindex2t(i, j, k)] = index2;
                   ++index2;
                }
             }
      }
      i1_start[3][np] = index1;
      i2_start[3][np] = index2;
     
      /**********************************************************************/
      /* map1to3 mapping - done - tranpose operation #4 (k,i,j) <-- (i,j,k) */
      /**********************************************************************/
      index1 = 0;
      index2 = 0;
      for (auto it = 0; it < np; ++it) 
      {
         proc_to = (taskid + it) % np;
         proc_from = (taskid - it + np) % np;
         i1_start[4][it] = index1;
         i2_start[4][it] = index2;
         for (auto k = 0; k < nz; ++k)
           for (auto j = 0; j < ny; ++j)
             for (auto i = 0; i < (nx); ++i) 
             {
                phere = cijktop2(i, j, k);
                pto = cijktop(i, j, k);
               
                /* packing scheme */
                if ((phere == taskid) && (pto == proc_to)) 
                {
                   iq_to_i1[4][cijktoindex2t(i, j, k)] = index1;
                   ++index1;
                }
                /* unpacking scheme */
                if ((pto == taskid) && (phere == proc_from)) 
                {
                   iq_to_i2[4][cijktoindex(i, j, k)] = index2;
                   ++index2;
                }
             }
      }
      i1_start[4][np] = index1;
      i2_start[4][np] = index2;
     
      /**********************************************************************/
      /* map3to1 mapping - done - tranpose operation #5 (i,j,k) <-- (k,i,j) */
      /**********************************************************************/
      index1 = 0;
      index2 = 0;
      for (auto it = 0; it < np; ++it) 
      {
         proc_to = (taskid + it) % np;
         proc_from = (taskid - it + np) % np;
         i1_start[5][it] = index1;
         i2_start[5][it] = index2;
         for (auto k = 0; k < nz; ++k)
           for (auto j = 0; j < ny; ++j)
             for (auto i = 0; i < (nx); ++i) 
             {
                phere = cijktop(i, j, k);
                pto = cijktop2(i, j, k);
               
                /* packing scheme */
                if ((phere == taskid) && (pto == proc_to)) 
                {
                   iq_to_i1[5][cijktoindex(i, j, k)] = index1;
                   ++index1;
                }
                /* unpacking scheme */
                if ((pto == taskid) && (phere == proc_from)) 
                {
                   iq_to_i2[5][cijktoindex2t(i, j, k)] = index2;
                   ++index2;
                }
             }
      }
      i1_start[5][np] = index1;
      i2_start[5][np] = index2;
     
      /* allocate ptranspose indexes */
      for (auto nb=0; nb<=nbrillq; ++nb) 
      {
         p_iq_to_i1[nb] = new (std::nothrow) int *[6];
         p_iq_to_i1[nb][0] = new (std::nothrow) int[nx * nq1]();
         p_iq_to_i1[nb][1] = new (std::nothrow) int[ny * nq2]();
         p_iq_to_i1[nb][2] = new (std::nothrow) int[nz * nq3]();
         p_iq_to_i1[nb][3] = new (std::nothrow) int[ny * nq2]();
         p_iq_to_i1[nb][4] = new (std::nothrow) int[nx * nq1]();
         p_iq_to_i1[nb][5] = new (std::nothrow) int[nz * nq3]();
        
         p_iq_to_i2[nb] = new (std::nothrow) int *[6];
         p_iq_to_i2[nb][0] = new (std::nothrow) int[ny * nq2]();
         p_iq_to_i2[nb][1] = new (std::nothrow) int[nz * nq3]();
         p_iq_to_i2[nb][2] = new (std::nothrow) int[ny * nq2]();
         p_iq_to_i2[nb][3] = new (std::nothrow) int[nx * nq1]();
         p_iq_to_i2[nb][4] = new (std::nothrow) int[nz * nq3]();
         p_iq_to_i2[nb][5] = new (std::nothrow) int[nx * nq1]();
        
         p_iz_to_i2_count[nb] = new (std::nothrow) int[6];
         p_iz_to_i2[nb] = new (std::nothrow) int *[6];
         p_iz_to_i2[nb][0] = new (std::nothrow) int[ny * nq2]();
         p_iz_to_i2[nb][1] = new (std::nothrow) int[nz * nq3]();
         p_iz_to_i2[nb][2] = new (std::nothrow) int[ny * nq2]();
         p_iz_to_i2[nb][3] = new (std::nothrow) int[nx * nq1]();
         p_iz_to_i2[nb][4] = new (std::nothrow) int[nz * nq3]();
         p_iz_to_i2[nb][5] = new (std::nothrow) int[nx * nq1]();
        
         p_i1_start[nb] = new (std::nothrow) int *[6];
         for (auto i=0; i<6; ++i)
            p_i1_start[nb][i] = new (std::nothrow) int[np + 1]();
        
         p_i2_start[nb] = new (std::nothrow) int *[6];
         for (auto i = 0; i<6; ++i)
            p_i2_start[nb][i] = new (std::nothrow) int[np + 1]();
      }
   }
     

   /* setup ffts */
   c3db_tmp1 = new (std::nothrow) double[2*nfft3d]();
   c3db_tmp2 = new (std::nothrow) double[2*nfft3d]();


   /* setup ffts */
   tmpx = new (std::nothrow) double[2*(2*nx+15)]();
   tmpy = new (std::nothrow) double[2*(2*ny+15)]();
   tmpz = new (std::nothrow) double[2*(2*nz+15)]();
   dcffti_(&nx,tmpx);
   dcffti_(&ny,tmpy);
   dcffti_(&nz,tmpz);

/*if (defined NWPW_SYCL) || (defined NWPW_CUDA) || (defined NWPW_HIP) */
#if (defined NWPW_CUDA) || (defined NWPW_HIP)
   if (maptype==1) 
     fft_tag = mygdevice.batch_fft_init(nx,ny,nz,ny*nq,(nx)*nq,(nx)*nq);
   else
     fft_tag = mygdevice.batch_fft_init(nx,ny,nz,nq1,nq2,nq3);
#endif
}


/********************************
 *                              *
 *          Destructors         *
 *                              *
 ********************************/
c3db::~c3db() 
{
   int i, nb;

#if (defined NWPW_CUDA) || (defined NWPW_HIP)
   if (mygdevice.has_gpu())
      mygdevice.batch_fft_end(fft_tag);
#endif

   if (maptype == 1) {
      delete[] iq_to_i1[0];
      delete[] iq_to_i1;
     
      delete[] iq_to_i2[0];
      delete[] iq_to_i2;
     
      delete[] i1_start[0];
      delete[] i1_start;
      delete[] i2_start[0];
      delete[] i2_start;
      for (nb=0; nb<=p_nbrillq0; ++nb) 
      {
         delete[] p_iq_to_i1[nb][0];
         delete[] p_iq_to_i1[nb];
         delete[] p_iq_to_i2[nb][0];
         delete[] p_iq_to_i2[nb];
         delete[] p_iz_to_i2[nb][0];
         delete[] p_iz_to_i2[nb];
         delete[] p_i1_start[nb][0];
         delete[] p_i1_start[nb];
         delete[] p_i2_start[nb][0];
         delete[] p_i2_start[nb];
         delete[] p_jq_to_i1[nb][0];
         delete[] p_jq_to_i1[nb];
         delete[] p_jq_to_i2[nb][0];
         delete[] p_jq_to_i2[nb];
         delete[] p_jz_to_i2[nb][0];
         delete[] p_jz_to_i2[nb];
         delete[] p_j1_start[nb][0];
         delete[] p_j1_start[nb];
         delete[] p_j2_start[nb][0];
         delete[] p_j2_start[nb];
      }
   } 
   else 
   {
      this->r_transpose_ijk_end();
     
      for (i = 0; i < 6; ++i) {
         delete[] iq_to_i1[i];
         delete[] iq_to_i2[i];
         delete[] i1_start[i];
         delete[] i2_start[i];
      }
      delete[] iq_to_i1;
      delete[] iq_to_i2;
      delete[] i1_start;
      delete[] i2_start;
     
      for (nb=0; nb<=p_nbrillq0; ++nb) {
         for (i = 0; i < 6; ++i) {
            delete[] p_iq_to_i1[nb][i];
            delete[] p_iq_to_i2[nb][i];
            delete[] p_iz_to_i2[nb][i];
            delete[] p_i1_start[nb][i];
            delete[] p_i2_start[nb][i];
         }
         delete[] p_iq_to_i1[nb];
         delete[] p_iq_to_i2[nb];
         delete[] p_iz_to_i2[nb];
         delete[] p_i1_start[nb];
         delete[] p_i2_start[nb];
         delete[] p_iz_to_i2_count[nb];
      }
   }
   delete [] p_iq_to_i1;
   delete [] p_iq_to_i2;
   delete [] p_iz_to_i2; 
   delete [] p_iz_to_i2_count;
   delete [] p_i1_start; 
   delete [] p_i2_start;

   delete [] p_jq_to_i1;
   delete [] p_jq_to_i2;
   delete [] p_jz_to_i2;
   delete [] p_j1_start;
   delete [] p_j2_start;
 
   delete [] tmpx;
   delete [] tmpy;
   delete [] tmpz;

   delete [] c3db_tmp1;
   delete [] c3db_tmp2;
 
   //#endif
}

/******************************************
 *                                        *
 *      c3db::c_ptranspose_jk_init        *
 *                                        *
 ******************************************/
void c3db::c_ptranspose_jk_init(const int nb, bool *zero_arow3) 
{
   int index1, index2, index3;
   int jndex1, jndex2, jndex3;
   int proc_to, proc_from, phere, pto, pfrom;
   bool iszero_ii, iszero_jj;
 
   index1 = 0;
   index2 = 0;
   index3 = 0;
   jndex1 = 0;
   jndex2 = 0;
   jndex3 = 0;
   for (auto it = 0; it < np; ++it) {
     proc_to = (taskid + it) % np;
     proc_from = (taskid - it + np) % np;
     p_i1_start[nb][0][it] = index1;
     p_i2_start[nb][0][it] = index2;
 
     p_j1_start[nb][0][it] = jndex1;
     p_j2_start[nb][0][it] = jndex2;
 
     for (auto k = 0; k < nz; ++k)
       for (auto j = 0; j < ny; ++j) 
       {
          /* packing scheme */
          phere = cijktop(0, 0, k);
          pto = cijktop(0, 0, j);
          if ((phere == taskid) && (pto == proc_to))
            for (auto i = 0; i < nx; ++i) 
            {
               iszero_ii = (zero_arow3[i + (nx) * k]);
               iszero_jj = (zero_arow3[i + (nx) * j]);
               if (!iszero_ii) 
               {
                  p_iq_to_i1[nb][0][index1] = cijktoindex(i, j, k);
                  ++index1;
               }
               if (!iszero_jj) 
               {
                  p_jq_to_i1[nb][0][jndex1] = cijktoindex(i, j, k);
                  ++jndex1;
               }
            }
         
          /* unpacking scheme */
          phere = cijktop(0, 0, j);
          pfrom = cijktop(0, 0, k);
          if ((phere == taskid) && (pfrom == proc_from))
            for (auto i = 0; i < (nx); ++i) 
            {
               iszero_ii = (zero_arow3[i + (nx) * k]);
               iszero_jj = (zero_arow3[i + (nx) * j]);
               if (!iszero_ii) 
               {
                  p_iq_to_i2[nb][0][index2] = cijktoindex(i, k, j);
                  ++index2;
               } 
               else 
               {
                  p_iz_to_i2[nb][0][index3] = cijktoindex(i, k, j);
                  ++index3;
               }
               if (!iszero_jj) 
               {
                  p_jq_to_i2[nb][0][jndex2] = cijktoindex(i, k, j);
                  ++jndex2;
               } 
               else 
               {
                  p_jz_to_i2[nb][0][jndex3] = cijktoindex(i, k, j);
                  ++jndex3;
               }
            }
       }
   }
   p_i1_start[nb][0][np] = index1;
   p_i2_start[nb][0][np] = index2;
   p_j1_start[nb][0][np] = jndex1;
   p_j2_start[nb][0][np] = jndex2;
}


/******************************************
 *                                        *
 *      c3db::c_ptranspose_ijk_init       *
 *                                        *
 ******************************************/
void c3db::c_ptranspose_ijk_init(const int nb, bool *zero_arow2, bool *zero_arow3) 
{
   int index1, index2, index3, proc_to, proc_from, phere, pto;
   bool iszero;
 
   /**************************************************/
   /* map1to2 mapping - done - tranpose operation #0 */
   /*  (ny,nz,nx)  <-- (nx,ny,nz)                    */
   /*   use zero_arow2                               */
   /**************************************************/
   index1 = 0;
   index2 = 0;
   index3 = 0;
   for (auto it = 0; it < np; ++it) {
     proc_to = (taskid + it) % np;
     proc_from = (taskid - it + np) % np;
     p_i1_start[nb][0][it] = index1;
     p_i2_start[nb][0][it] = index2;
     for (auto k = 0; k < nz; ++k)
       for (auto j = 0; j < ny; ++j)
         for (auto i = 0; i < nx; ++i) 
         {
            iszero = (zero_arow2[i + k*nx]);
           
            // phere = int_mb(p_map1(1,id)+(j-1)+(k-1)*ny(id))
            // pto   = int_mb(p_map2(1,id)+(k-1)+(i-1)*nz(id))
           
            phere = cijktop2(i, j, k);
            pto = cijktop1(i, j, k);
           
            /* packing scheme */
            if ((phere == taskid) && (pto == proc_to)) 
            {
               if (!iszero) 
               {
                  p_iq_to_i1[nb][0][index1] = cijktoindex2t(i, j, k);
                  ++index1;
               }
            }
            /* unpacking scheme */
            if ((pto == taskid) && (phere == proc_from)) 
            {
               if (!iszero) 
               {
                  p_iq_to_i2[nb][0][index2] = cijktoindex1(i, j, k);
                  ++index2;
               } 
               else 
               {
                  p_iz_to_i2[nb][0][index3] = cijktoindex1(i, j, k);
                  ++index3;
               }
            }
         }
   }
   p_i1_start[nb][0][np] = index1;
   p_i2_start[nb][0][np] = index2;
   p_iz_to_i2_count[nb][0] = index3;
 
   /**************************************************/
   /* map2to3 mapping - done - tranpose operation #1 */
   /*     (nz,nx,ny)  <-- (ny,nz,nx)                 */
   /*     use zero_arow3                             */
   /**************************************************/
   index1 = 0;
   index2 = 0;
   index3 = 0;
   for (auto it = 0; it < np; ++it) 
   {
      proc_to = (taskid + it) % np;
      proc_from = (taskid - it + np) % np;
      p_i1_start[nb][1][it] = index1;
      p_i2_start[nb][1][it] = index2;
      for (auto k = 0; k < nz; ++k)
        for (auto j = 0; j < ny; ++j)
          for (auto i = 0; i < nx; ++i) 
          {
             iszero = (zero_arow3[i + j*nx]);
            
             // phere = int_mb(p_map2(1,id)+(k-1)+(i-1)*nz(id))
             // pto   = int_mb(p_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))
            
             phere = cijktop1(i, j, k);
             pto = cijktop(i, j, k);
            
             /* packing scheme */
             if ((phere == taskid) && (pto == proc_to)) {
               if (!iszero) {
                 p_iq_to_i1[nb][1][index1] = cijktoindex1(i, j, k);
                 ++index1;
               }
             }
             /* unpacking scheme */
             if ((pto == taskid) && (phere == proc_from)) {
               if (!iszero) {
                 p_iq_to_i2[nb][1][index2] = cijktoindex(i, j, k);
                 ++index2;
               } else {
                 p_iz_to_i2[nb][1][index3] = cijktoindex(i, j, k);
                 ++index3;
               }
             }
          }
   }
   p_i1_start[nb][1][np] = index1;
   p_i2_start[nb][1][np] = index2;
   p_iz_to_i2_count[nb][1] = index3;

 
   /**************************************************/
   /* map3to2 mapping - done - tranpose operation #2 */
   /*     (ny,nz,nx)  <-- (nz,nx,ny)                 */
   /*     use zero_arow3                             */
   /**************************************************/
   index1 = 0;
   index2 = 0;
   index3 = 0;
   for (auto it = 0; it < np; ++it) 
   {
      proc_to = (taskid + it) % np;
      proc_from = (taskid - it + np) % np;
      p_i1_start[nb][2][it] = index1;
      p_i2_start[nb][2][it] = index2;
      for (auto k = 0; k < nz; ++k)
        for (auto j = 0; j < ny; ++j)
          for (auto i = 0; i < nx; ++i) 
          {
             iszero = (zero_arow3[i + j*nx]);
            
             phere = cijktop(i, j, k);
             pto = cijktop1(i, j, k);
            
             /* packing scheme */
             if ((phere == taskid) && (pto == proc_to)) {
               if (!iszero) {
                 p_iq_to_i1[nb][2][index1] = cijktoindex(i, j, k);
                 ++index1;
               }
             }
             /* unpacking scheme */
             if ((pto == taskid) && (phere == proc_from)) 
             {
                if (!iszero) 
                {
                   p_iq_to_i2[nb][2][index2] = cijktoindex1(i, j, k);
                   ++index2;
                } 
                else 
                {
                   p_iz_to_i2[nb][2][index3] = cijktoindex1(i, j, k);
                   ++index3;
                }
             }
          }
   }
   p_i1_start[nb][2][np] = index1;
   p_i2_start[nb][2][np] = index2;
   p_iz_to_i2_count[nb][2] = index3;
 
   /**************************************************/
   /* map2to1 mapping - done - tranpose operation #3 */
   /*     (nx,ny,nz)  <-- (ny,nz,nx)                 */
   /*     use zero_arow2                             */
   /**************************************************/
   index1 = 0;
   index2 = 0;
   index3 = 0;
   for (auto it = 0; it < np; ++it) 
   {
      proc_to = (taskid + it) % np;
      proc_from = (taskid - it + np) % np;
      p_i1_start[nb][3][it] = index1;
      p_i2_start[nb][3][it] = index2;
      for (auto k = 0; k < nz; ++k)
        for (auto j = 0; j < ny; ++j)
          for (auto i = 0; i < nx; ++i) 
          {
             iszero = (zero_arow2[i + k*nx]);
            
            
             phere = cijktop1(i, j, k);
             pto = cijktop2(i, j, k);
            
             /* packing scheme */
             if ((phere == taskid) && (pto == proc_to)) 
             {
                if (!iszero) 
                {
                   p_iq_to_i1[nb][3][index1] = cijktoindex1(i, j, k);
                   ++index1;
                }
             }
             /* unpacking scheme */
             if ((pto == taskid) && (phere == proc_from)) 
             {
                if (!iszero) 
                {
                   p_iq_to_i2[nb][3][index2] = cijktoindex2t(i, j, k);
                   ++index2;
                } 
                else 
                {
                   p_iz_to_i2[nb][3][index3] = cijktoindex2t(i, j, k);
                   ++index3;
                }
             }
          }
   }
   p_i1_start[nb][3][np] = index1;
   p_i2_start[nb][3][np] = index2;
   p_iz_to_i2_count[nb][3] = index3;

 
   /**************************************************/
   /* map1to3 mapping - done - tranpose operation #4 */
   /**************************************************/
   index1 = 0;
   index2 = 0;
   index3 = 0;
   for (auto it = 0; it < np; ++it) 
   {
      proc_to = (taskid + it) % np;
      proc_from = (taskid - it + np) % np;
      p_i1_start[nb][4][it] = index1;
      p_i2_start[nb][4][it] = index2;
      for (auto k = 0; k < nz; ++k)
        for (auto j = 0; j < ny; ++j)
          for (auto i = 0; i < nx; ++i) 
          {
            
             phere = cijktop2(i, j, k);
             pto = cijktop(i, j, k);
            
             /* packing scheme */
             if ((phere == taskid) && (pto == proc_to)) 
             {
                p_iq_to_i1[nb][4][index1] = cijktoindex2t(i, j, k);
                ++index1;
             }
             /* unpacking scheme */
             if ((pto == taskid) && (phere == proc_from)) 
             {
                p_iq_to_i2[nb][4][index2] = cijktoindex(i, j, k);
                ++index2;
             }
          }
   }
   p_i1_start[nb][4][np] = index1;
   p_i2_start[nb][4][np] = index2;
   p_iz_to_i2_count[nb][4] = index3;

 
   /**************************************************/
   /* map3to1 mapping - done - tranpose operation #5 */
   /**************************************************/
   index1 = 0;
   index2 = 0;
   index3 = 0;
   for (auto it = 0; it < np; ++it) 
   {
      proc_to = (taskid + it) % np;
      proc_from = (taskid - it + np) % np;
      p_i1_start[nb][5][it] = index1;
      p_i2_start[nb][5][it] = index2;
      for (auto k = 0; k < nz; ++k)
        for (auto j = 0; j < ny; ++j)
          for (auto i = 0; i < nx; ++i) 
          {
             // phere = int_mb(p_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))
             // pto   = int_mb(p_map1(1,id)+(j-1)+(k-1)*ny(id))
            
             phere = cijktop(i, j, k);
             pto = cijktop2(i, j, k);
            
             /* packing scheme */
             if ((phere == taskid) && (pto == proc_to)) 
             {
                p_iq_to_i1[nb][5][index1] = cijktoindex(i, j, k);
                ++index1;
             }
             /* unpacking scheme */
             if ((pto == taskid) && (phere == proc_from)) 
             {
                p_iq_to_i2[nb][5][index2] = cijktoindex2t(i, j, k);
                ++index2;
             }
          }
   }
   p_i1_start[nb][5][np] = index1;
   p_i2_start[nb][5][np] = index2;
   p_iz_to_i2_count[nb][5] = index3;
}




/********************************
 *                              *
 *         c3db::c_alloc        *
 *                              *
 ********************************/
/**
 * @brief Allocate memory for a real-space array and initialize all elements to zero.
 *
 * This function allocates memory for a real-space array and initializes all its
 * elements to zero. It is used to create real-space arrays for various calculations.
 *
 * @return A pointer to the allocated real-space array.
 */
double *c3db::c_alloc()
{
   double *ptr = new (std::nothrow) double[2*nfft3d]();
   return ptr;
}


/********************************
 *                              *
 *         c3db::c_nalloc       *
 *                              *
 ********************************/
/**
 * @brief Allocate a 1D array with extended size for real-space data.
 *
 * This function allocates a 1D array for storing real-space data with extended size.
 *
 * @param nn The number of elements to allocate.
 * @return A pointer to the allocated array or nullptr if allocation fails.
 */
double *c3db::c_nalloc(const int nn)
{
   double *ptr = new (std::nothrow) double[2*nfft3d * nn]();
   return ptr;
}

/********************************
 *                              *
 *         c3db::c_dealloc      *
 *                              *
 ********************************/
void c3db::c_dealloc(double *ptr) { delete[] ptr; }
 



/********************************
 *                              *
 *         c3db::r_alloc        *
 *                              *
 ********************************/
/**
 * @brief Allocate memory for a real-space array and initialize all elements to zero.
 *
 * This function allocates memory for a real-space array and initializes all its
 * elements to zero. It is used to create real-space arrays for various calculations.
 *
 * @return A pointer to the allocated real-space array.
 */
double *c3db::r_alloc() 
{
   double *ptr = new (std::nothrow) double[nfft3d]();
   return ptr;
}


/********************************
 *                              *
 *         c3db::r_nalloc       *
 *                              *
 ********************************/
/**
 * @brief Allocate a 1D array with extended size for real-space data.
 *
 * This function allocates a 1D array for storing real-space data with extended size.
 *
 * @param nn The number of elements to allocate.
 * @return A pointer to the allocated array or nullptr if allocation fails.
 */
double *c3db::r_nalloc(const int nn) 
{
   double *ptr = new (std::nothrow) double[nfft3d * nn]();
   return ptr;
}


/********************************
 *                              *
 *         c3db::r_dealloc      *
 *                              *
 ********************************/
void c3db::r_dealloc(double *ptr) { delete[] ptr; }
 


/********************************
 *                              *
 *         c3db::r_zero         *
 *                              *
 ********************************/
/**
 * @brief Fill a double array with zeros.
 *
 * This function fills a given double array `ptr` with zeros. The size of the array to be filled is determined by the value of `n2ft3d`.
 *
 * @param ptr A pointer to the double array to be filled with zeros.
 *
 * @return None.
 *
 * @note The function efficiently fills the array with zeros using memory operations. The size of the operation is determined by the value of `n2ft3d`, which specifies the number of zeros to fill.
 */
void c3db::r_zero(double *ptr) 
{
   std::memset(ptr, 0, nfft3d * sizeof(double));
}

/********************************
 *                              *
 *         c3db::r_nzero        *
 *                              *
 ********************************/
/**
 * @brief Fill a double array with zeros.
 *
 * This function fills a given double array `ptr` with zeros. The number of zeros to be filled is determined by the value of `n` multiplied by `n2ft3d`.
 *
 * @param n The number of repetitions for filling zeros.
 * @param ptr A pointer to the double array to be filled with zeros.
 *
 * @return None.
 *
 * @note The function efficiently fills the array with zeros using memory operations. The size of the operation is determined by the product of `n` and `n2ft3d`, which specifies the number of zeros to fill.
 */
void c3db::r_nzero(int n, double *ptr) 
{
   std::memset(ptr, 0, n * nfft3d * sizeof(double));
}


/********************************
 *                              *
 *        c3db::rr_copy         *
 *                              *
 ********************************/
/**
 * @brief Copy the elements of one double array to another.
 *
 * This function efficiently copies the elements from the source double array `ptr1` to the destination double array `ptr2`.
 * The size of the arrays determines the number of elements to be copied, which is specified by the value of `n2ft3d`.
 *
 * @param ptr1 A pointer to the source double array containing the elements to be copied.
 * @param ptr2 A pointer to the destination double array where the copied elements will be stored.
 *
 * @return None.
 *
 * @note The function uses memory copying to efficiently transfer the elements from `ptr1` to `ptr2`. The size of the operation
 *       is determined by the value of `n2ft3d`, which specifies the number of elements to copy.
 */
void c3db::rr_copy(const double *ptr1, double *ptr2) 
{
   std::memcpy(ptr2, ptr1, nfft3d * sizeof(double));
}


/********************************
 *                              *
 *        c3db::rr_SMul         *
 *                              *
 ********************************/
 /**
 * @brief Multiply each element of a double array by a scalar value.
 *
 * This function multiplies each element of the source double array `ptr1` by the scalar value `da` and stores
 * the result in the destination double array `ptr2`. The size of the arrays and the scalar value `da` determine
 * the number of elements to process.
 *
 * @param da The scalar value to multiply with each element of `ptr1`.
 * @param ptr1 A pointer to the source double array containing the elements to be multiplied.
 * @param ptr2 A pointer to the destination double array where the multiplied elements will be stored.
 *
 * @return None.
 *
 * @note The function efficiently performs element-wise multiplication of `ptr1` with the scalar `da` and stores
 *       the result in `ptr2`. The size of the operation is determined by the value of `n2ft3d_map`, which specifies
 *       the number of elements to process.
 */
void c3db::rr_SMul(const double da, const double *ptr1, double *ptr2) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr2[i] = da * ptr1[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr2[i]   = da * ptr1[i];
      ptr2[i+1] = da * ptr1[i+1];
      ptr2[i+2] = da * ptr1[i+2];
      ptr2[i+3] = da * ptr1[i+3];
      ptr2[i+4] = da * ptr1[i+4];
   }
   return;
}

void c3db::cc_SMul(const double da, const double *ptr1, double *ptr2)
{
   int m = n2ft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr2[i] = da*ptr1[i];
   if (n2ft3d_map < 5)
      return;
   for (auto i=m; i<n2ft3d_map; i += 5)
   {
      ptr2[i]   = da * ptr1[i];
      ptr2[i+1] = da * ptr1[i+1];
      ptr2[i+2] = da * ptr1[i+2];
      ptr2[i+3] = da * ptr1[i+3];
      ptr2[i+4] = da * ptr1[i+4];
   }
   return;
}

void c3db::c_ZMul(const double da, const double db, double *ptr2)
{
   int m = nfft3d_map % 8;
   if (m > 0)
      for (auto i=0; i<m; ++i)
      {
         double c0 = ptr2[i]; double d0 = ptr2[i+1];
         ptr2[i]   = da*c0 - db*d0;
         ptr2[i+1] = da*d0 + db*c0;
      }
   if (nfft3d_map < 8)
      return;
   for (auto i=m; i<nfft3d_map; i += 8)
   {
      double c0 = ptr2[i];   double d0 = ptr2[i+1];
      double c1 = ptr2[i+2]; double d1 = ptr2[i+3];
      double c2 = ptr2[i+4]; double d2 = ptr2[i+5];
      double c3 = ptr2[i+6]; double d3 = ptr2[i+7];

      ptr2[i]   = da*c0 - db*d0;
      ptr2[i+1] = da*d0 + db*c0;

      ptr2[i+2] = da*c1 - db*d1;
      ptr2[i+3] = da*d1 + db*c1;

      ptr2[i+4] = da*c2 - db*d2;
      ptr2[i+5] = da*d2 + db*c2;

      ptr2[i+6] = da*c3 - db*d3;
      ptr2[i+7] = da*d3 + db*c3;
   }
   return;
}





/********************************
 *                              *
 *        c3db::rrr_SMulAdd     *
 *                              *
 ********************************/
 /**
 * @brief Perform element-wise scaling, addition, and addition of double arrays and store the result.
 *
 * This function performs an element-wise scaling of the input double array `ptr1` by the scalar `da`, then
 * adds the elements of the scaled array to the elements of the `ptr2` array. The final result is stored in
 * the `ptr3` array. It applies these operations to each element of the arrays individually. The function is
 * optimized for processing arrays with a size greater than or equal to 5 elements.
 *
 * @param da   A scalar value for scaling the elements of `ptr1`.
 * @param ptr1 A pointer to the input double array to be scaled.
 * @param ptr2 A pointer to the second input double array for addition.
 * @param ptr3 A pointer to the output double array where the result is stored.
 *
 * @return None.
 *
 * @note The function scales each element of `ptr1` by the scalar `da`, adds the elements of the scaled
 *       array to the elements of `ptr2`, and stores the result in `ptr3` for corresponding elements in
 *       the arrays. Loop unrolling is applied for efficiency when processing arrays with a size greater
 *       than or equal to 5 elements.
 */
void c3db::rrr_SMulAdd(const double da, const double *ptr1, const double *ptr2, double *ptr3) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr3[i] = da*ptr1[i] + ptr2[i];
   if (nfft3d_map < 5)
      return;
   for (auto i = m; i<nfft3d_map; i+=5) 
   {
      ptr3[i]   = da * ptr1[i]   + ptr2[i];
      ptr3[i+1] = da * ptr1[i+1] + ptr2[i+1];
      ptr3[i+2] = da * ptr1[i+2] + ptr2[i+2];
      ptr3[i+3] = da * ptr1[i+3] + ptr2[i+3];
      ptr3[i+4] = da * ptr1[i+4] + ptr2[i+4];
   }
   return;
}

/********************************
 *                              *
 *     c3db::rrrrr_SumMulAdd    *
 *                              *
 ********************************/
/*
   ptr5 = (ptr1+ptr2)*ptr3 + pt4
*/
/**
 * @brief Perform element-wise addition, multiplication, and addition of double arrays and store the result.
 *
 * This function performs an element-wise addition of the input double arrays `ptr1` and `ptr2`, then
 * multiplies the result by the elements of the `ptr3` array, and finally adds the elements of the `ptr4`
 * array. The result is stored in the `ptr5` array. It applies these operations to each element of the
 * arrays individually. The function is optimized for processing arrays with a size greater than or equal
 * to 5 elements.
 *
 * @param ptr1 A pointer to the first input double array for addition.
 * @param ptr2 A pointer to the second input double array for addition.
 * @param ptr3 A pointer to the double array for multiplication.
 * @param ptr4 A pointer to the double array for addition.
 * @param ptr5 A pointer to the output double array where the result is stored.
 *
 * @return None.
 *
 * @note The function computes the element-wise sum of `ptr1` and `ptr2`, multiplies the result by the
 *       elements of `ptr3`, and then adds the elements of `ptr4` to the final result stored in `ptr5`
 *       for corresponding elements in the arrays. Loop unrolling is applied for efficiency when
 *       processing arrays with a size greater than or equal to 5 elements.
 */
void c3db::rrrrr_SumMulAdd(const double *ptr1, const double *ptr2,
                           const double *ptr3, const double *ptr4,
                           double *ptr5)
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr5[i] = (ptr1[i]+ptr2[i])*ptr3[i] + ptr4[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr5[i] = (ptr1[i]+ptr2[i])*ptr3[i] + ptr4[i];
      ptr5[i+1] = (ptr1[i+1]+ptr2[i+1])*ptr3[i+1] + ptr4[i+1];
      ptr5[i+2] = (ptr1[i+2]+ptr2[i+2])*ptr3[i+2] + ptr4[i+2];
      ptr5[i+3] = (ptr1[i+3]+ptr2[i+3])*ptr3[i+3] + ptr4[i+3];
      ptr5[i+4] = (ptr1[i+4]+ptr2[i+4])*ptr3[i+4] + ptr4[i+4];
   }
   return;
}

/********************************
 *                              *
 *          c3db::r_SMul        *
 *                              *
 ********************************/
void c3db::r_SMul(const double da, double *ptr2)
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr2[i] *= da;
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i+=5)
   {
      ptr2[i]   *= da;
      ptr2[i+1] *= da;
      ptr2[i+2] *= da;
      ptr2[i+3] *= da;
      ptr2[i+4] *= da;
   }
}


/********************************
 *                              *
 *          c3db::c_SMul        *
 *                              *
 ********************************/
void c3db::c_SMul(const double da, double *ptr2)
{
   int m = n2ft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr2[i] *= da;
   if (n2ft3d_map < 5)
      return;
   for (auto i=m; i<n2ft3d_map; i+=5)
   {
      ptr2[i]   *= da;
      ptr2[i+1] *= da;
      ptr2[i+2] *= da;
      ptr2[i+3] *= da;
      ptr2[i+4] *= da;
   }
}




/********************************
 *                              *
 *          c3db::r_abs         *
 *                              *
 ********************************/
void c3db::r_abs(double *ptr2) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr2[i] = std::abs(ptr2[i]);
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr2[i] = std::abs(ptr2[i]);
      ptr2[i+1] = std::abs(ptr2[i+1]);
      ptr2[i+2] = std::abs(ptr2[i+2]);
      ptr2[i+3] = std::abs(ptr2[i+3]);
      ptr2[i+4] = std::abs(ptr2[i+4]);
   }
}

/********************************
 *                              *
 *          c3db::r_sqr         *
 *                              *
 ********************************/
void c3db::r_sqr(double *ptr2) 
{
   int i;
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr2[i] *= ptr2[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr2[i] *= ptr2[i];
      ptr2[i+1] *= ptr2[i+1];
      ptr2[i+2] *= ptr2[i+2];
      ptr2[i+3] *= ptr2[i+3];
      ptr2[i+4] *= ptr2[i+4];
   }
}

/********************************
 *                              *
 *          c3db::rr_sqr        *
 *                              *
 ********************************/
void c3db::rr_sqr(const double *ptr2, double *ptr3) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr3[i] = ptr2[i]*ptr2[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr3[i] = ptr2[i]*ptr2[i];
      ptr3[i+1] = ptr2[i+1]*ptr2[i+1];
      ptr3[i+2] = ptr2[i+2]*ptr2[i+2];
      ptr3[i+3] = ptr2[i+3]*ptr2[i+3];
      ptr3[i+4] = ptr2[i+4]*ptr2[i+4];
   }
}

/********************************
 *                              *
 *        c3db::rr_addsqr       *
 *                              *
 ********************************/
void c3db::rr_addsqr(const double *ptr2, double *ptr3) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr3[i] += ptr2[i]*ptr2[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr3[i] += ptr2[i]*ptr2[i];
      ptr3[i+1] += ptr2[i+1]*ptr2[i+1];
      ptr3[i+2] += ptr2[i+2]*ptr2[i+2];
      ptr3[i+3] += ptr2[i+3]*ptr2[i+3];
      ptr3[i+4] += ptr2[i+4]*ptr2[i+4];
   }
}

/********************************
 *                              *
 *          c3db::r_sqrt        *
 *                              *
 ********************************/
void c3db::r_sqrt(double *ptr2) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr2[i] = sqrt(ptr2[i]);
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr2[i] = sqrt(ptr2[i]);
      ptr2[i+1] = sqrt(ptr2[i+1]);
      ptr2[i+2] = sqrt(ptr2[i+2]);
      ptr2[i+3] = sqrt(ptr2[i+3]);
      ptr2[i+4] = sqrt(ptr2[i+4]);
   }
}

/********************************
 *                              *
 *         c3db::r_dsum         *
 *                              *
 ********************************/
double c3db::r_dsum(const double *ptr) 
{
   int m = nfft3d_map % 5;
   double sum = 0.0;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         sum += ptr[i];
   if (nfft3d_map < 5)
      return sum;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      sum += ptr[i] + ptr[i+1] + ptr[i+2] + ptr[i+3] + ptr[i+4];
   }
   return parall->SumAll(1, sum);
}

/********************************
 *                              *
 *         c3db::rrr_Sum2Add    *
 *                              *
 ********************************/
void c3db::rrr_Sum2Add(const double *ptr1, const double *ptr2, double *ptr3) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr3[i] += ptr1[i] + ptr2[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr3[i] += ptr1[i] + ptr2[i];
      ptr3[i+1] += ptr1[i+1] + ptr2[i+1];
      ptr3[i+2] += ptr1[i+2] + ptr2[i+2];
      ptr3[i+3] += ptr1[i+3] + ptr2[i+3];
      ptr3[i+4] += ptr1[i+4] + ptr2[i+4];
   }
}

/********************************
 *                              *
 *         c3db::rrrr_Sum       *
 *                              *
 ********************************/
void c3db::rrrr_Sum(const double *ptr1, const double *ptr2, const double *ptr3, double *ptr4)
{
   int m = nfft3d_map%5;
   if (m>0)
   {
      for (auto i=0; i<m; ++i)
         ptr4[i] = ptr1[i] + ptr2[i] + ptr3[i];
   }
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i+=5)
   {
      ptr4[i]   = ptr1[i]   + ptr2[i]   + ptr3[i];
      ptr4[i+1] = ptr1[i+1] + ptr2[i+1] + ptr3[i+1];
      ptr4[i+2] = ptr1[i+2] + ptr2[i+2] + ptr3[i+2];
      ptr4[i+3] = ptr1[i+3] + ptr2[i+3] + ptr3[i+3];
      ptr4[i+4] = ptr1[i+4] + ptr2[i+4] + ptr3[i+4];
   }
}

/********************************
 *                              *
 *         c3db::rrr_Sum        *
 *                              *
 ********************************/
void c3db::rrr_Sum(const double *ptr1, const double *ptr2, double *ptr3)
{
   int m = nfft3d_map%5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr3[i] = ptr1[i] + ptr2[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i+=5)
   {
      ptr3[i] = ptr1[i] + ptr2[i];
      ptr3[i+1] = ptr1[i+1] + ptr2[i+1];
      ptr3[i+2] = ptr1[i+2] + ptr2[i+2];
      ptr3[i+3] = ptr1[i+3] + ptr2[i+3];
      ptr3[i+4] = ptr1[i+4] + ptr2[i+4];
   }
}


/********************************
 *                              *
 *         c3db::rrc_Sum        *
 *                              *
 ********************************/
void c3db::rrc_Sum(const double *ptr1, const double *ptr2, double *ptr3)
{     
   int m = nfft3d_map%5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
      {
         ptr3[2*i]   = ptr1[i] + ptr2[i];    ptr3[2*i+1] = 0.0;
      }
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i+=5)
   {  
      ptr3[2*i]     = ptr1[i] + ptr2[i];     ptr3[2*i+1]     = 0.0;
      ptr3[2*(i+1)] = ptr1[i+1] + ptr2[i+1]; ptr3[2*(i+1)+1] = 0.0;
      ptr3[2*(i+2)] = ptr1[i+2] + ptr2[i+2]; ptr3[2*(i+2)+1] = 0.0;
      ptr3[2*(i+3)] = ptr1[i+3] + ptr2[i+3]; ptr3[2*(i+3)+1] = 0.0;
      ptr3[2*(i+4)] = ptr1[i+4] + ptr2[i+4]; ptr3[2*(i+4)+1] = 0.0;
   }  
}    



/********************************
 *                              *
 *         c3db::rr_Sum         *
 *                              *
 ********************************/
void c3db::rr_Sum(const double *ptr2, double *ptr3) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr3[i] += ptr2[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr3[i] += ptr2[i];
      ptr3[i+1] += ptr2[i+1];
      ptr3[i+2] += ptr2[i+2];
      ptr3[i+3] += ptr2[i+3];
      ptr3[i+4] += ptr2[i+4];
   }
}

/********************************
 *                              *
 *         c3db::rcc_Sum        *
 *                              *
 ********************************/
void c3db::rcc_Sum(const double *ptr1, const double *ptr2, double *ptr3)
{     
   int m = nfft3d_map%5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
      {
         ptr3[2*i]   = ptr1[i] + ptr2[2*i];  
         ptr3[2*i+1] = ptr2[2*i+1];  
      }
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i+=5)
   {  
      ptr3[2*i]       = ptr1[i]   + ptr2[2*i];   
      ptr3[2*i+1]     =             ptr2[2*i+1];   
      ptr3[2*(i+1)]   = ptr1[i+1] + ptr2[2*(i+1)]; 
      ptr3[2*(i+1)+1] =             ptr2[2*(i+1)+1]; 
      ptr3[2*(i+2)]   = ptr1[i+2] + ptr2[2*(i+2)]; 
      ptr3[2*(i+2)+1] =             ptr2[2*(i+2)+1]; 
      ptr3[2*(i+3)]   = ptr1[i+3] + ptr2[2*(i+3)]; 
      ptr3[2*(i+3)+1] =             ptr2[2*(i+3)+1]; 
      ptr3[2*(i+4)]   = ptr1[i+4] + ptr2[2*(i+4)]; 
      ptr3[2*(i+4)+1] =             ptr2[2*(i+4)+1]; 
   }
}     


/********************************
 *                              *
 *         c3db::rrr_Minus      *
 *                              *
 ********************************/
void c3db::rrr_Minus(const double *ptr1, const double *ptr2, double *ptr3) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
        ptr3[i] = ptr1[i] - ptr2[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr3[i] = ptr1[i] - ptr2[i];
      ptr3[i+1] = ptr1[i+1] - ptr2[i+1];
      ptr3[i+2] = ptr1[i+2] - ptr2[i+2];
      ptr3[i+3] = ptr1[i+3] - ptr2[i+3];
      ptr3[i+4] = ptr1[i+4] - ptr2[i+4];
   }
}

/********************************
 *                              *
 *         c3db::arrr_Minus     *
 *                              *
 ********************************/
void c3db::arrr_Minus(const double a, const double *ptr1, const double *ptr2, double *ptr3)
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr3[i] = a*(ptr1[i] - ptr2[i]);
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr3[i] = a*(ptr1[i] - ptr2[i]);
      ptr3[i+1] = a*(ptr1[i+1] - ptr2[i+1]);
      ptr3[i+2] = a*(ptr1[i+2] - ptr2[i+2]);
      ptr3[i+3] = a*(ptr1[i+3] - ptr2[i+3]);
      ptr3[i+4] = a*(ptr1[i+4] - ptr2[i+4]);
   }
}

/********************************
 *                              *
 *         c3db::rr_Minus       *
 *                              *
 ********************************/
 /**
 * @brief Perform element-wise subtraction of two double arrays and store the result.
 *
 * This function performs an element-wise subtraction between the input double array
 * `ptr2` and the output double array `ptr3`, and stores the result in `ptr3`. It subtracts
 * each element `ptr2[i]` from the corresponding element `ptr3[i]` and assigns the result
 * to `ptr3[i]`. The function is optimized for processing arrays with a size greater than
 * or equal to 5 elements.
 *
 * @param ptr2 A pointer to the input double array for subtraction.
 * @param ptr3 A pointer to the output double array where the result is stored.
 *
 * @return None.
 *
 * @note The function computes the element-wise subtraction of `ptr2` from `ptr3` for
 *       corresponding elements in the arrays. The result is stored in the `ptr3` array.
 *       Loop unrolling is applied for efficiency when processing arrays with a size
 *       greater than or equal to 5 elements.
 */
void c3db::rr_Minus(const double *ptr2, double *ptr3) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr3[i] -= ptr2[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr3[i] -= ptr2[i];
      ptr3[i+1] -= ptr2[i+1];
      ptr3[i+2] -= ptr2[i+2];
      ptr3[i+3] -= ptr2[i+3];
      ptr3[i+4] -= ptr2[i+4];
   }
}

/********************************
 *                              *
 *         c3db::rrr_Mul        *
 *                              *
 ********************************/
 /**
 * @brief Perform element-wise multiplication between two double arrays and store the result.
 *
 * This function performs an element-wise multiplication between the input double arrays
 * `ptr1` and `ptr2`, and stores the result in the output double array `ptr3`. It multiplies
 * each element `ptr1[i]` with the corresponding element `ptr2[i]` and assigns the result
 * to `ptr3[i]`. The function is optimized for processing arrays with a size greater than
 * or equal to 5 elements.
 *
 * @param ptr1 A pointer to the first input double array for multiplication.
 * @param ptr2 A pointer to the second input double array for multiplication.
 * @param ptr3 A pointer to the output double array where the result is stored.
 *
 * @return None.
 *
 * @note The function computes the element-wise multiplication of `ptr1` and `ptr2` for
 *       corresponding elements in the arrays. The result is stored in the `ptr3` array.
 *       Loop unrolling is applied for efficiency when processing arrays with a size
 *       greater than or equal to 5 elements.
 */
void c3db::rrr_Mul(const double *ptr1, const double *ptr2, double *ptr3) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr3[i] = ptr1[i]*ptr2[i];
   if (nfft3d_map<5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr3[i] = ptr1[i]*ptr2[i];
      ptr3[i+1] = ptr1[i+1]*ptr2[i+1];
      ptr3[i+2] = ptr1[i+2]*ptr2[i+2];
      ptr3[i+3] = ptr1[i+3]*ptr2[i+3];
      ptr3[i+4] = ptr1[i+4]*ptr2[i+4];
   }
   return;
}

/********************************
 *                              *
 *         c3db::rr_Mul         *
 *                              *
 ********************************/
 /**
 * @brief Perform an element-wise multiplication with another double array.
 *
 * This function performs an element-wise multiplication on the input double arrays
 * `ptr1` and `ptr3`. It multiplies each element `ptr3[i]` by the corresponding
 * element `ptr1[i]` and stores the result in `ptr3`. The function is optimized for
 * processing arrays with a size greater than or equal to 5 elements.
 *
 * @param ptr1 A pointer to the input double array used for multiplication.
 * @param ptr3 A pointer to the output double array where the result is stored.
 *
 * @return None.
 *
 * @note The function computes the element-wise multiplication of `ptr3` by `ptr1` for
 *       corresponding elements in the arrays. The result is stored in the `ptr3`
 *       array. Loop unrolling is applied for efficiency when processing arrays with
 *       a size greater than or equal to 5 elements.
 */
void c3db::rr_Mul(const double *ptr1, double *ptr3) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr3[i] *= ptr1[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr3[i] *= ptr1[i];
      ptr3[i+1] *= ptr1[i+1];
      ptr3[i+2] *= ptr1[i+2];
      ptr3[i+3] *= ptr1[i+3];
      ptr3[i+4] *= ptr1[i+4];
   }
   return;
}


/********************************
 *                              *
 *         c3db::cc_Mul         *
 *                              *
 ********************************/
void c3db::cc_Mul(const double *ptr1, double *ptr3)
{      
   int m = nfft3d_map % 5; 
   if (m > 0)
      for (auto i=0; i<m; ++i) 
      {
         int i0r = 2*i;  int i0i = 2*i+1;
         double x0 = ptr3[i0r]*ptr1[i0r] - ptr3[i0i]*ptr1[i0i];
         double y0 = ptr3[i0r]*ptr1[i0i] + ptr3[i0i]*ptr1[i0r];
         ptr3[i0r] = x0;
         ptr3[i0i] = y0;
      }
   if (nfft3d_map < 5) 
      return;
   for (auto i=m; i<nfft3d_map; i += 5)
   {
      int i1=i+1; int i2=i+2; int i3=i+3; int i4=i+4;
      int i0r = 2*i;  int i0i = 2*i+1;
      int i1r = 2*i1; int i1i = 2*i1+1;
      int i2r = 2*i2; int i2i = 2*i2+1;
      int i3r = 2*i3; int i3i = 2*i3+1;
      int i4r = 2*i4; int i4i = 2*i4+1;

      double x0 = ptr3[i0r]*ptr1[i0r] - ptr3[i0i]*ptr1[i0i];
      double y0 = ptr3[i0r]*ptr1[i0i] + ptr3[i0i]*ptr1[i0r];

      double x1 = ptr3[i1r]*ptr1[i1r] - ptr3[i1i]*ptr1[i1i];
      double y1 = ptr3[i1r]*ptr1[i1i] + ptr3[i1i]*ptr1[i1r];

      double x2 = ptr3[i2r]*ptr1[i2r] - ptr3[i2i]*ptr1[i2i];
      double y2 = ptr3[i2r]*ptr1[i2i] + ptr3[i2i]*ptr1[i2r];

      double x3 = ptr3[i3r]*ptr1[i3r] - ptr3[i3i]*ptr1[i3i];
      double y3 = ptr3[i3r]*ptr1[i3i] + ptr3[i3i]*ptr1[i3r];

      double x4 = ptr3[i4r]*ptr1[i4r] - ptr3[i4i]*ptr1[i4i];
      double y4 = ptr3[i4r]*ptr1[i4i] + ptr3[i4i]*ptr1[i4r];

      ptr3[i0r] = x0;
      ptr3[i0i] = y0;

      ptr3[i1r] = x1;
      ptr3[i1i] = y1;

      ptr3[i2r] = x2;
      ptr3[i2i] = y2;

      ptr3[i3r] = x3;
      ptr3[i3i] = y3;

      ptr3[i4r] = x4;
      ptr3[i4i] = y4;
   }
   return;
}  

/********************************
 *                              *
 *         c3db::bb_Mul         *
 *                              *
 ********************************/
//ptr3 = conjg(ptr3)*ptr1
void c3db::bb_Mul(const double *ptr1, double *ptr3)
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
      {
         int i0r = 2*i;  int i0i = 2*i+1;
         double x0  = ptr3[i0r]*ptr1[i0r] + ptr3[i0i]*ptr1[i0i];
         double y0  = ptr3[i0r]*ptr1[i0i] - ptr3[i0i]*ptr1[i0r];
         ptr3[i0r] = x0;
         ptr3[i0i] = y0;
      }
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5)
   {
      int i1=i+1; int i2=i+2; int i3=i+3; int i4=i+4;
      int i0r = 2*i;  int i0i = 2*i+1;
      int i1r = 2*i1; int i1i = 2*i1+1;
      int i2r = 2*i2; int i2i = 2*i2+1;
      int i3r = 2*i3; int i3i = 2*i3+1;
      int i4r = 2*i4; int i4i = 2*i4+1;

      double x0 = ptr3[i0r]*ptr1[i0r] + ptr3[i0i]*ptr1[i0i];
      double y0 = ptr3[i0r]*ptr1[i0i] - ptr3[i0i]*ptr1[i0r];

      double x1 = ptr3[i1r]*ptr1[i1r] + ptr3[i1i]*ptr1[i1i];
      double y1 = ptr3[i1r]*ptr1[i1i] - ptr3[i1i]*ptr1[i1r];

      double x2 = ptr3[i2r]*ptr1[i2r] + ptr3[i2i]*ptr1[i2i];
      double y2 = ptr3[i2r]*ptr1[i2i] - ptr3[i2i]*ptr1[i2r];

      double x3 = ptr3[i3r]*ptr1[i3r] + ptr3[i3i]*ptr1[i3i];
      double y3 = ptr3[i3r]*ptr1[i3i] - ptr3[i3i]*ptr1[i3r];

      double x4 = ptr3[i4r]*ptr1[i4r] + ptr3[i4i]*ptr1[i4i];
      double y4 = ptr3[i4r]*ptr1[i4i] - ptr3[i4i]*ptr1[i4r];

      ptr3[i0r] = x0;
      ptr3[i0i] = y0;

      ptr3[i1r] = x1;
      ptr3[i1i] = y1;

      ptr3[i2r] = x2;
      ptr3[i2i] = y2;

      ptr3[i3r] = x3;
      ptr3[i3i] = y3;

      ptr3[i4r] = x4;
      ptr3[i4i] = y4;
   }
   return;
} 




/********************************
 *                              *
 *     c3db::rrr_SqrMulAdd      *
 *                              *
 ********************************/
 /**
 * @brief Compute an element-wise operation involving squares and multiplication.
 *
 * This function performs an element-wise operation on the input arrays `ptr1`, `ptr2`,
 * and `ptr3`. It computes `(ptr1[i]^2) * ptr2[i] + ptr3[i]` for each element `i` in
 * the arrays and stores the result in the `ptr3` array. The function is optimized for
 * processing arrays with a size greater than or equal to 5 elements.
 *
 * @param ptr1 A pointer to the first input double array.
 * @param ptr2 A pointer to the second input double array.
 * @param ptr3 A pointer to the output double array where the result is stored.
 *
 * @return None.
 *
 * @note The function computes the specified element-wise operation for corresponding
 *       elements in the input arrays `ptr1`, `ptr2`, and `ptr3`. The result is stored
 *       in the `ptr3` array. Loop unrolling is applied for efficiency when processing
 *       arrays with a size greater than or equal to 5 elements.
 */
void c3db::rrr_SqrMulAdd(const double *ptr1, const double *ptr2, double *ptr3) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr3[i] += ptr1[i]*ptr1[i]*ptr2[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr3[i]   += (ptr1[i]*ptr1[i])*ptr2[i];
      ptr3[i+1] += (ptr1[i+1]*ptr1[i+1])*ptr2[i+1];
      ptr3[i+2] += (ptr1[i+2]*ptr1[i+2])*ptr2[i+2];
      ptr3[i+3] += (ptr1[i+3]*ptr1[i+3])*ptr2[i+3];
      ptr3[i+4] += (ptr1[i+4]*ptr1[i+4])*ptr2[i+4];
   }
   return;
}

/******************************************
 *                                        *
 *     c3db::rrrrrrr_Sqr3MulPlusMul2      *
 *                                        *
 ******************************************/
 /**
 * @brief Compute a complex expression involving multiple array operations.
 *
 * This function computes a complex expression involving multiple array operations:
 * `(ptr1[i]^2 + ptr2[i]^2 + ptr3[i]^2) * ptr4[i] + ptr5[i] * ptr6[i]`, where each
 * term in the expression corresponds to elements in different arrays. The result is
 * stored in the `ptr7` array. The function is optimized for processing arrays with
 * a size greater than or equal to 5 elements.
 *
 * @param ptr1 A pointer to the first input double array.
 * @param ptr2 A pointer to the second input double array.
 * @param ptr3 A pointer to the third input double array.
 * @param ptr4 A pointer to the fourth input double array.
 * @param ptr5 A pointer to the fifth input double array.
 * @param ptr6 A pointer to the sixth input double array.
 * @param ptr7 A pointer to the output double array where the result is stored.
 *
 * @return None.
 *
 * @note The function computes the specified complex expression element-wise for
 *       corresponding elements in the input arrays `ptr1`, `ptr2`, `ptr3`, `ptr4`,
 *       `ptr5`, and `ptr6`. The result is stored in the `ptr7` array. Loop unrolling
 *       is applied for efficiency when processing arrays with a size greater than or
 *       equal to 5 elements.
 */
void c3db::rrrrrrr_Sqr3MulPlusMul2(const double *ptr1, const double *ptr2,
                                   const double *ptr3, const double *ptr4,
                                   const double *ptr5, const double *ptr6,
                                   double *ptr7) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr7[i] = (ptr1[i]*ptr1[i] + ptr2[i]*ptr2[i] + ptr3[i]*ptr3[i])*ptr4[i] + ptr5[i]*ptr6[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr7[i] = (ptr1[i]*ptr1[i] + ptr2[i]*ptr2[i] + ptr3[i]*ptr3[i])* ptr4[i] + ptr5[i]*ptr6[i];
      ptr7[i+1] = (ptr1[i+1]*ptr1[i+1] + ptr2[i+1]*ptr2[i+1] + ptr3[i+1]*ptr3[i+1])*ptr4[i+1] + ptr5[i+1]*ptr6[i+1];
      ptr7[i+2] = (ptr1[i+1]*ptr1[i+2] + ptr2[i+2]*ptr2[i+2] + ptr3[i+2]*ptr3[i+2])*ptr4[i+2] + ptr5[i+2]*ptr6[i+2];
      ptr7[i+3] = (ptr1[i+1]*ptr1[i+3] + ptr2[i+3]*ptr2[i+3] + ptr3[i+3]*ptr3[i+3])*ptr4[i+3] + ptr5[i+3]*ptr6[i+3];
      ptr7[i+4] = (ptr1[i+1]*ptr1[i+4] + ptr2[i+4]*ptr2[i+4] + ptr3[i+4]*ptr3[i+4])*ptr4[i+4] + ptr5[i+4]*ptr6[i+4];
   }
   return;
}

/********************************
 *                              *
 *         c3db::rc_Mul         *
 *                              *
 ********************************/
 /**
 * @brief Element-wise complex multiplication of two complex double arrays.
 *
 * This function performs element-wise complex multiplication of two complex double
 * arrays, `ptr1` and `ptr3`. It multiplies the real and imaginary parts of each
 * complex number in `ptr3` by the corresponding real value in `ptr1`. The result is
 * updated in place in the `ptr3` array.
 *
 * @param ptr1 A pointer to the input double array containing the real values.
 * @param ptr3 A pointer to the complex double array where the result is updated.
 *
 * @return None.
 *
 * @note The function computes element-wise complex multiplication by multiplying both
 *       the real and imaginary parts of each complex number in `ptr3` by the
 *       corresponding real value in `ptr1`. Loop unrolling is applied for efficiency
 *       when processing arrays with a size greater than or equal to 5 complex numbers.
 */
void c3db::rc_Mul(const double *ptr1, double *ptr3) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i) 
      {
         ptr3[2*i]   *= ptr1[i];
         ptr3[2*i+1] *= ptr1[i];
      }
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr3[2*(i)]   *= ptr1[i];
      ptr3[2*(i)+1] *= ptr1[i];
     
      ptr3[2*(i+1)]   *= ptr1[i+1];
      ptr3[2*(i+1)+1] *= ptr1[i+1];
     
      ptr3[2*(i+2)]   *= ptr1[i+2];
      ptr3[2*(i+2)+1] *= ptr1[i+2];
     
      ptr3[2*(i+3)]   *= ptr1[i+3];
      ptr3[2*(i+3)+1] *= ptr1[i+3];
     
      ptr3[2*(i+4)]   *= ptr1[i+4];
      ptr3[2*(i+4)+1] *= ptr1[i+4];
   }
   return;
}

/********************************
 *                              *
 *         c3db::rrr_Mul2Add    *
 *                              *
 ********************************/
 /**
 * @brief Element-wise multiplication and addition of two double arrays.
 *
 * This function performs element-wise multiplication of two double arrays, `ptr1` and
 * `ptr2`, and adds the result to the corresponding elements of `ptr3`. It is designed
 * to optimize the operation with loop unrolling for efficiency.
 *
 * @param ptr1 A pointer to the first input double array for multiplication.
 * @param ptr2 A pointer to the second input double array for multiplication.
 * @param ptr3 A pointer to the output double array where the result is added.
 *
 * @return None.
 *
 * @note The function computes the element-wise multiplication of `ptr1` and `ptr2` and
 *       adds the result to the corresponding elements of `ptr3`. Loop unrolling is
 *       applied for efficiency in processing arrays with a size greater than or equal
 *       to 5 elements.
 */
void c3db::rrr_Mul2Add(const double *ptr1, const double *ptr2, double *ptr3) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr3[i] += ptr1[i] * ptr2[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i < nfft3d_map; i += 5) 
   {
      ptr3[i] += ptr1[i]*ptr2[i];
      ptr3[i+1] += ptr1[i+1]*ptr2[i+1];
      ptr3[i+2] += ptr1[i+2]*ptr2[i+2];
      ptr3[i+3] += ptr1[i+3]*ptr2[i+3];
      ptr3[i+4] += ptr1[i+4]*ptr2[i+4];
   }
   return;
}

/********************************
 *                              *
 *         c3db::rrr_Divide     *
 *                              *
 ********************************/
 /**
 * @brief Element-wise division of two double arrays with safety checks.
 *
 * This function performs element-wise division of two double arrays, `ptr1` and `ptr2`,
 * with appropriate safety checks to avoid division by zero. The result is stored in
 * `ptr3`. If the absolute value of an element in `ptr2` is greater than the defined
 * constant `ETA_DIV` (1.0e-9), the division is performed; otherwise, the corresponding
 * element in `ptr3` is set to 0.0.
 *
 * @param ptr1 A pointer to the first input double array for division.
 * @param ptr2 A pointer to the second input double array for division.
 * @param ptr3 A pointer to the output double array to store the result.
 *
 * @return None.
 *
 * @note The function computes the division element-wise and checks if the absolute
 *       value of each element in `ptr2` is greater than `ETA_DIV` before performing
 *       division. If the condition is met, it computes `ptr1[i] / ptr2[i]` for each
 *       element `i`, otherwise, it sets `ptr3[i]` to 0.0.
 */
#define ETA_DIV 1.0e-9
void c3db::rrr_Divide(const double *ptr1, const double *ptr2, double *ptr3)
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr3[i] = (std::abs(ptr2[i])>ETA_DIV) ? (ptr1[i]/ptr2[i]) : (0.0);
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5)
   {
      ptr3[i] = (std::abs(ptr2[i]) > ETA_DIV) ? (ptr1[i] / ptr2[i]) : (0.0);
      ptr3[i+1] = (std::abs(ptr2[i+1]) > ETA_DIV) ? (ptr1[i+1]/ptr2[i+1]) : (0.0);
      ptr3[i+2] = (std::abs(ptr2[i+2]) > ETA_DIV) ? (ptr1[i+2]/ptr2[i+2]) : (0.0);
      ptr3[i+3] = (std::abs(ptr2[i+3]) > ETA_DIV) ? (ptr1[i+3]/ptr2[i+3]) : (0.0);
      ptr3[i+4] = (std::abs(ptr2[i+4]) > ETA_DIV) ? (ptr1[i+4]/ptr2[i+4]) : (0.0);
   }
   return;
}

/********************************
 *                              *
 *         c3db::rr_Divide      *
 *                              *
 ********************************/
 /**
 * @brief Element-wise division of two double arrays with safety checks.
 *
 * This function performs element-wise division of two double arrays, `ptr3` and `ptr2`,
 * with appropriate safety checks to avoid division by zero. The result is stored in
 * `ptr3`. If the absolute value of an element in `ptr2` is greater than `ETA_DIV`, the
 * division is performed; otherwise, the corresponding element in `ptr3` is set to 0.0.
 *
 * @param ptr2 A pointer to the input double array for division.
 * @param ptr3 A pointer to the output double array to store the result.
 *
 * @return None.
 *
 * @note The function computes the division element-wise and checks if the absolute
 *       value of each element in `ptr2` is greater than `ETA_DIV` before performing
 *       division. If the condition is met, it computes `ptr3[i] / ptr2[i]` for each
 *       element `i`, otherwise, it sets `ptr3[i]` to 0.0.
 */
void c3db::rr_Divide(const double *ptr2, double *ptr3)
{
   int m = nfft3d_map%5;
   if (m > 0)
      for (auto i = 0; i < m; ++i)
         ptr3[i] = (std::abs(ptr2[i]) > ETA_DIV) ? (ptr3[i]/ptr2[i]) : (0.0);
   if (nfft3d_map<5)
      return;
   for (auto i = m; i < nfft3d_map; i += 5)
   {
      ptr3[i] = (std::abs(ptr2[i]) > ETA_DIV) ? (ptr3[i] / ptr2[i]) : (0.0);
      ptr3[i+1] = (std::abs(ptr2[i+1]) > ETA_DIV) ? (ptr3[i+1]/ptr2[i+1]) : (0.0);
      ptr3[i+2] = (std::abs(ptr2[i+2]) > ETA_DIV) ? (ptr3[i+2]/ptr2[i+2]) : (0.0);
      ptr3[i+3] = (std::abs(ptr2[i+3]) > ETA_DIV) ? (ptr3[i+3]/ptr2[i+3]) : (0.0);
      ptr3[i+4] = (std::abs(ptr2[i+4]) > ETA_DIV) ? (ptr3[i+4]/ptr2[i+4]) : (0.0);
   }
   return;
}

/********************************
 *                              *
 *         c3db::rr_screen0     *
 *                              *
 ********************************/
/**
 * @brief Calculate the screening factor screen0 = (1.0/epsilon - 1.0).
 *
 * This function computes the screening factor `screen0` as (1.0/epsilon - 1.0), where
 * `epsilon` is a double array provided as `ptr2`. The result is stored in `ptr3`. The
 * calculation is performed element-wise with appropriate handling for values close to
 * zero to avoid division by zero.
 *
 * @param ptr2 A pointer to the input double array `epsilon`.
 * @param ptr3 A pointer to the output double array `screen0`.
 *
 * @return None.
 *
 * @note The function computes `screen0` element-wise, and for each element of `epsilon`,
 *       it checks if the absolute value is greater than `ETA_DIV` to avoid division by
 *       zero. If the absolute value is greater, it calculates (1.0/epsilon - 1.0),
 *       otherwise, it sets the corresponding element of `screen0` to 0.0.
 */
void c3db::rr_screen0(const double *ptr2, double *ptr3)
{
   int m = nfft3d_map%5;
   if (m>0)
      for (auto i=0; i<m; ++i)
         ptr3[i] = (std::abs(ptr2[i]) > ETA_DIV) ? (1.0/ptr2[i]-1.0) : (0.0);
   if (nfft3d_map<5)
      return;
   for (auto i=m; i<nfft3d_map; i+=5)
   {
      ptr3[i]   = (std::abs(ptr2[i])   > ETA_DIV) ? (1.0/ptr2[i]-1.0) : (0.0);
      ptr3[i+1] = (std::abs(ptr2[i+1]) > ETA_DIV) ? (1.0/ptr2[i+1]-1.0) : (0.0);
      ptr3[i+2] = (std::abs(ptr2[i+2]) > ETA_DIV) ? (1.0/ptr2[i+2]-1.0) : (0.0);
      ptr3[i+3] = (std::abs(ptr2[i+3]) > ETA_DIV) ? (1.0/ptr2[i+3]-1.0) : (0.0);
      ptr3[i+4] = (std::abs(ptr2[i+4]) > ETA_DIV) ? (1.0/ptr2[i+4]-1.0) : (0.0);
   }
   return;
}



/********************************
 *                              *
 *         c3db::rr_daxpy       *
 *                              *
 ********************************/
 /**
 * @brief Perform a scaled vector addition of two double arrays.
 *
 * This function performs a scaled vector addition of two double arrays, `ptr1` and
 * `ptr2`, with a scaling factor `alpha`. The operation updates the values of `ptr2`
 * in place. It uses loop unrolling for efficiency when applicable.
 *
 * @param alpha The scaling factor applied to `ptr1` before addition.
 * @param ptr1 A pointer to the first double array.
 * @param ptr2 A pointer to the second double array, which is updated with the result.
 *
 * @return None.
 */
void c3db::rr_daxpy(const double alpha, const double *ptr1, double *ptr2) 
{
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         ptr2[i] += alpha * ptr1[i];
   if (nfft3d_map < 5)
      return;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      ptr2[i]   += alpha*ptr1[i];
      ptr2[i+1] += alpha*ptr1[i+1];
      ptr2[i+2] += alpha*ptr1[i+2];
      ptr2[i+3] += alpha*ptr1[i+3];
      ptr2[i+4] += alpha*ptr1[i+4];
   }
}

/********************************
 *                              *
 *         c3db::rr_dot         *
 *                              *
 ********************************/
 /**
 * @brief Compute the dot product of two double arrays.
 *
 * This function calculates the dot product of two double arrays, `ptr1` and `ptr2`,
 * with appropriate parallelization. It uses a loop unrolling technique to optimize
 * the calculation for efficiency. The result is a double value representing the dot
 * product of the arrays.
 *
 * @param ptr1 A pointer to the first double array.
 * @param ptr2 A pointer to the second double array.
 *
 * @return The dot product of the two input arrays as a double value.
 */
double c3db::rr_dot(const double *ptr1, const double *ptr2) 
{
   double sum = 0.0;
 
   int m = nfft3d_map % 5;
   if (m > 0)
      for (auto i=0; i<m; ++i)
         sum += ptr1[i]*ptr2[i];
   if (nfft3d_map < 5)
      return sum;
   for (auto i=m; i<nfft3d_map; i += 5) 
   {
      sum += ptr1[i]*ptr2[i] 
           + ptr1[i+1]*ptr2[i+1] 
           + ptr1[i+2]*ptr2[i+2] 
           + ptr1[i+3]*ptr2[i+3] 
           + ptr1[i+4]*ptr2[i+4];
   }
 
   return parall->SumAll(1,sum);
}

/********************************
 *                              *
 *         c3db::nrr_vdot       *
 *                              *
 ********************************/
void c3db::nrr_vdot(const int n, const double *ptr1, const double *ptr2, double *v) 
{
   int m = nfft3d_map % 5;
   for (auto k = 0; k < n; ++k)
      v[k] = 0.0;
   if (m > 0)
      for (auto i = 0; i < m; ++i)
         for (auto k = 0; k < n; ++k)
            v[k] += ptr1[n*i + k] * ptr2[i];
   if (nfft3d_map >= 5) {
      for (auto i=m; i < nfft3d_map; i += 5)
         for (auto k=0; k<n; ++k) 
         {
            v[k] += ptr1[n*i + k] * ptr2[i] +
                    ptr1[n*(i+1) + k]*ptr2[i+1] +
                    ptr1[n*(i+2) + k]*ptr2[i+2] +
                    ptr1[n*(i+3) + k]*ptr2[i+3] +
                    ptr1[n*(i+4) + k]*ptr2[i+4];
         }
     parall->Vector_SumAll(1, n, v);
   }
}

/********************************
 *                              *
 *         c3db::cc_dot         *
 *                              *
 ********************************/
double c3db::cc_dot(const double *ptr1, const double *ptr2) 
{
   double sum = 0.0;

  /*************************
   ****   slab mapping  ****
   *************************/
   if (maptype==1)
   {
      //***** kx!=0 plane, so double count *****
      for (auto q=0; q<nq; ++q)
      for (auto j=0; j<ny; ++j)
      for (auto i=1; i<(nx/2+1); ++i)
      {
        auto index = q*(nx/2+1)*ny + j*(nx/2+1) + i;
        sum += (ptr1[2*index]  *ptr2[2*index] + ptr1[2*index+1]*ptr2[2*index+1]);
      }
      sum *= 2.0;

      //***** kx==0 plane, so single count *****
      for (auto q=0; q<nq; ++q)
      for (auto j=0; j<ny; ++j)
      {
         auto index = q*(nx/2+1)*ny + j*(nx/2+1);
         sum += (ptr1[2*index]  *ptr2[2*index] + ptr1[2*index+1]*ptr2[2*index+1]);
      }
   }

  /*************************
   **** hilbert mapping ****
   *************************/
   else
   {
      int taskid_i = parall->taskid_i();

      //***** kx!=0 plane, so double count *****
      for (auto index=0; index<nfft3d; ++index)
      {
         sum += (ptr1[2*index]*ptr2[2*index] + ptr1[2*index+1]*ptr2[2*index+1]);
      }
      sum *= 2.0;

      //***** kx==0 plane, so single count *****
      for (auto k=0; k<nz; ++k)
      for (auto j=0; j<ny; ++j)
      {
         auto index = cijktoindex(0, j, k);
         auto p     = cijktop(0, j, k);
         if (p==taskid_i)
         {
            sum -= (ptr1[2*index]*ptr2[2*index] + ptr1[2*index+1]*ptr2[2*index+1]);
         }
      }

   }

   return parall->SumAll(1, sum);
}




/********************************
 *                              *
 *         c3db::c_read         *
 *                              *
 ********************************/
 /**
 * @brief Read data from an input unit with support for parallel computing.
 *
 * This function is responsible for reading data from a specified input unit, taking
 * into account parallel computing. It supports two mapping types: "slab mapping"
 * and "hilbert mapping." The function uses various parameters and parallelization
 * techniques to handle data retrieval and distribution efficiently.
 *
 * @param iunit An integer specifying the input unit to read data from.
 * @param a A pointer to a double array where the read data will be stored.
 * @param jcol An integer specifying the column index for parallelization.
 *
 * @return None.
 *
 * @note The behavior of this function depends on the mapping type set in the 'maptype'
 *       variable, the parallelization parameters, and the distribution of data from
 *       the input unit to the specified array.
 */
void c3db::c_read(const int iunit, double *a, const int jcol, const int kcol) 
{
   bool fillcolumn,fillzone;
   int jstart, jend,  kstart, kend;
   int taskid = parall->taskid();
   int taskid_j = parall->taskid_j();
   int taskid_k = parall->taskid_k();
   int np_j = parall->np_j();
   int np_k = parall->np_k();

   if (jcol < 0) 
   {
      jstart = 0;
      jend = np_j - 1;
      fillcolumn = true;
   } 
   else 
   {
      jstart = jend = jcol;
      fillcolumn = (taskid_j == jcol);
   }
   if (kcol < 0) 
   {
      kstart = 0;
      kend = np_k - 1;
      fillzone = true;
   } 
   else 
   {
      kstart = kend = kcol;
      fillzone = (taskid_k == kcol);
   }

 
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype==1) 
   {
      int bsize = 2*(nx)*ny;
      double *tmp = new (std::nothrow) double[bsize]();
 
      /**** master node reads from file and distributes ****/
      if (taskid == MASTER)
      {
         for (auto k=0; k<nz; ++k) 
         {
            dread(iunit, tmp, bsize);
        
            int index = 2*cijktoindex(0, 0, k);
            int ii    = cijktop(0, 0, k);
            for (auto kk=kstart; kk<=kend; ++kk) 
            for (auto jj=jstart; jj<=jend; ++jj) 
            {
               int p_to = parall->convert_taskid_ijk(ii,jj,kk);
               if (p_to==MASTER)
                  std::memcpy(a+index,tmp,bsize*sizeof(double));
               else
                  parall->dsend(0, 9,p_to,bsize,tmp);
            }
         }
      }
 
      /**** not master node ****/
      else if (fillcolumn && fillzone)
      {
         for (auto k=0; k<nz; ++k) 
         {
            int index = 2*cijktoindex(0, 0, k);
            int ii = cijktop(0, 0, k);
            int p_here = parall->convert_taskid_ijk(ii,taskid_j,taskid_k);
            if (p_here == taskid) 
            {
               parall->dreceive(0, 9, MASTER, bsize, tmp);
               std::memcpy(a+index,tmp,bsize*sizeof(double));
            }
         }
      }
      delete[] tmp;
   }
 
   /************************************
    **** hilbert or hcurve  mapping ****
    ************************************/
   else 
   {
      int bsize = 2*(nx);
      double tmp[bsize];
     
      /**** master node reads from file and distributes ****/
      if (taskid==MASTER)
      {
         for (auto k=0; k<nz; ++k)
         for (auto j=0; j<ny; ++j) 
         {
            dread(iunit, tmp, bsize);
        
            int index = 2*cijktoindex2(0, j, k);
            int ii    = cijktop2(0, j, k);
            for (int kk = kstart; kk <= kend; ++kk) 
            for (int jj = jstart; jj <= jend; ++jj) 
            {
               int p_to = parall->convert_taskid_ijk(ii, jj, kk);
               if (p_to == MASTER)
                  std::memcpy(a+index,tmp,bsize*sizeof(double));
               else
                  parall->dsend(0, 9, p_to, bsize, tmp);
            }
         }
      }
      /**** not master node ****/
      else if (fillcolumn && fillzone)
      {
         for (auto k=0; k<nz; ++k)
         for (auto j=0; j<ny; ++j) 
         {
            int index = 2*cijktoindex2(0, j, k);
            int ii = cijktop2(0, j, k);
            int p_here = parall->convert_taskid_ijk(ii, taskid_j, taskid_k);
            if (p_here == taskid) 
            {
               parall->dreceive(0, 9, MASTER, bsize, tmp);
               std::memcpy(a+index,tmp,bsize*sizeof(double));
            }
         }
      }
     
      if (fillcolumn && fillzone)
      {
         double *tmp1 = c3db::c3db_tmp1;
         double *tmp2 = c3db::c3db_tmp2;
         c_transpose_ijk(4, a, tmp1, tmp2);
      }
   }
}



/********************************
 *                              *
 *         c3db::r_read         *
 *                              *
 ********************************/
 /**
 * @brief Read data from an input unit with support for parallel computing.
 *
 * This function is responsible for reading data from a specified input unit, taking
 * into account parallel computing. It supports two mapping types: "slab mapping"
 * and "hilbert mapping." The function uses various parameters and parallelization
 * techniques to handle data retrieval and distribution efficiently.
 *
 * @param iunit An integer specifying the input unit to read data from.
 * @param a A pointer to a double array where the read data will be stored.
 * @param jcol An integer specifying the column index for parallelization.
 *
 * @return None.
 *
 * @note The behavior of this function depends on the mapping type set in the 'maptype'
 *       variable, the parallelization parameters, and the distribution of data from
 *       the input unit to the specified array.
 */
void c3db::r_read(const int iunit, double *a, const int jcol, const int kcol, const bool dotrans)
{
   bool fillcolumn,fillzone;
   int jstart, jend,  kstart, kend;
   int taskid = parall->taskid();
   int taskid_j = parall->taskid_j();
   int taskid_k = parall->taskid_k();
   int np_j = parall->np_j();
   int np_k = parall->np_k();

   if (jcol < 0)
   {
      jstart = 0;
      jend = np_j - 1;
      fillcolumn = true;
   }
   else
   {
      jstart = jend = jcol;
      fillcolumn = (taskid_j == jcol);
   }
   if (kcol < 0)
   {
      kstart = 0;
      kend = np_k - 1;
      fillzone = true;
   }
   else
   {
      kstart = kend = kcol;
      fillzone = (taskid_k == kcol);
   }

   /**********************
    **** slab mapping ****
    **********************/
   if (maptype==1)
   {
      int bsize = (nx)*ny;
      double *tmp = new (std::nothrow) double[bsize]();

      /**** master node reads from file and distributes ****/
      if (taskid == MASTER)
      {
         for (auto k=0; k<nz; ++k)
         {
            dread(iunit, tmp, bsize);

            int index = cijktoindex(0, 0, k);
            int ii    = cijktop(0, 0, k);
            for (auto kk=kstart; kk<=kend; ++kk)
            for (auto jj=jstart; jj<=jend; ++jj)
            {
               int p_to = parall->convert_taskid_ijk(ii,jj,kk);
               if (p_to==MASTER)
                  std::memcpy(a+index,tmp,bsize*sizeof(double));
               else
                  parall->dsend(0, 9,p_to,bsize,tmp);
            }
         }
      }

      /**** not master node ****/
      else if (fillcolumn && fillzone)
      {
         for (auto k=0; k<nz; ++k)
         {
            int index = cijktoindex(0, 0, k);
            int ii = cijktop(0, 0, k);
            int p_here = parall->convert_taskid_ijk(ii,taskid_j,taskid_k);
            if (p_here == taskid)
            {
               parall->dreceive(0, 9, MASTER, bsize, tmp);
               std::memcpy(a+index,tmp,bsize*sizeof(double));
            }
         }
      }
      delete[] tmp;
   }


   /*************************
    **** hilbert mapping ****
    *************************/
   else
   {
      int bsize = (nx);
      double tmp[nx];

      /**** master node reads from file and distributes ****/
      if (taskid==MASTER)
      {
         for (auto k=0; k<nz; ++k)
         for (auto j=0; j<ny; ++j)
         {
            dread(iunit, tmp, bsize);

            int index = cijktoindex2(0, j, k);
            int ii    = cijktop2(0, j, k);
            for (int kk = kstart; kk <= kend; ++kk)
            for (int jj = jstart; jj <= jend; ++jj)
            {
               int p_to = parall->convert_taskid_ijk(ii, jj, kk);
               if (p_to==MASTER)
                  std::memcpy(a+index,tmp,bsize*sizeof(double));
               else
                  parall->dsend(0, 9, p_to, bsize, tmp);
            }
         }
      }
      /**** not master node ****/
      else if (fillcolumn && fillzone)
      {
         for (auto k=0; k<nz; ++k)
         for (auto j=0; j<ny; ++j)
         {
            int index = cijktoindex2(0, j, k);
            int ii = cijktop2(0, j, k);
            int p_here = parall->convert_taskid_ijk(ii, taskid_j, taskid_k);
            if (p_here==taskid)
            {
               parall->dreceive(0, 9, MASTER, bsize, tmp);
               std::memcpy(a+index,tmp,bsize*sizeof(double));
            }
         }
      }

      if (fillcolumn && fillzone && dotrans)
      {
         double *tmp1 = c3db::c3db_tmp1;
         double *tmp2 = c3db::c3db_tmp2;
         r_transpose_ijk(4, a, tmp1, tmp2);
      }
   }
}




/********************************
 *                              *
 *         c3db::c_write        *
 *                              *
 ********************************/
 /**
 * @brief Write data to an output unit with support for parallel computing.
 *
 * This function is responsible for writing data to a specified output unit, taking
 * into account parallel computing. It supports two mapping types: "slab mapping"
 * and "hilbert mapping." The function uses various parameters and parallelization
 * techniques to handle data transfer and writing efficiently.
 *
 * @param iunit An integer specifying the output unit to write data to.
 * @param a A pointer to a double array containing the data to be written.
 * @param jcol An integer specifying the column index for parallelization.
 *
 * @return None.
 *
 * @note The behavior of this function depends on the mapping type set in the 'maptype'
 *       variable and the parallelization parameters. It supports both "slab mapping"
 *       and "hilbert mapping" for data distribution.
 */
void c3db::c_write(const int iunit, double *a, const int jcol, const int kcol) 
{
   int taskid = parall->taskid();
   int taskid_j = parall->taskid_j();
   int taskid_k = parall->taskid_k();
   int np_j = parall->np_j();
   int np_k = parall->np_k();
   int idum[1] = {1};


   /**********************
    **** slab mapping ****
    **********************/
   if (maptype==1) 
   {
      int bsize = 2*nx*ny;
      double *tmp = new (std::nothrow) double[bsize]();
      
      /**** master node gathers and write to file ****/
      if (taskid==MASTER)
      {
         for (auto k=0; k<nz; ++k) 
         {
            int ii = cijktop(0, 0, k);
            int p_from = parall->convert_taskid_ijk(ii,jcol,kcol);
            if (p_from == MASTER) 
            {
               int index = 2*cijktoindex(0, 0, k);
               std::memcpy(tmp,a+index,bsize*sizeof(double));
            } 
            else 
            {
               parall->isend(0, 7, p_from, 1, idum);
               parall->dreceive(0, 9, p_from, bsize, tmp);
            }
            dwrite(iunit, tmp, bsize);
         }
      }
      
      /**** not master node ****/
      else
      {
         for (auto k=0; k<nz; ++k) 
         {
            int index = 2 * cijktoindex(0, 0, k);
            int ii = cijktop(0, 0, k);
            int p_here = parall->convert_taskid_ijk(ii, jcol, kcol);
            if (p_here==taskid) 
            {
               std::memcpy(tmp,a+index,bsize*sizeof(double));
               parall->ireceive(0, 7, MASTER, 1, idum);
               parall->dsend(0, 9, MASTER, bsize, tmp);
            }
         }
      }
      delete[] tmp;
   }

  /*************************
   **** hilbert mapping ****
   *************************/
   else 
   {


      if ((taskid_j==jcol) && (taskid_k==kcol))
      {
         // double *tmp1 = new (std::nothrow) double[2*nfft3d];
         // double *tmp2 = new (std::nothrow) double[2*nfft3d];
         //double *tmp1 = new (std::nothrow) double[2 * nfft3d]();
         //double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
         double *tmp1 = c3db::c3db_tmp1;
         double *tmp2 = c3db::c3db_tmp2;
         c_transpose_ijk(5, a, tmp1, tmp2);
     
         // delete [] tmp2;
         // delete [] tmp1;
      }
     
      int bsize = 2*nx;
      double tmp[bsize];
     
      /**** master node write to file and fetches from other nodes ****/
      if (taskid == MASTER)
      {
         for (auto k=0; k<nz; ++k)
         for (auto j=0; j<ny; ++j) 
         {

            int ii = cijktop2(0, j, k);
            int p_from = parall->convert_taskid_ijk(ii, jcol, kcol);
     
            if (p_from==MASTER) 
            {
               int index = 2*cijktoindex2(0, j, k);
               std::memcpy(tmp,a+index,bsize*sizeof(double));
            } 
            else 
            {
               parall->isend(0, 7, p_from, 1, idum);
               parall->dreceive(0, 9, p_from, bsize, tmp);
            }

            dwrite(iunit, tmp, bsize);
         }
      } 
      /**** not master node ****/
      else
      {
         for (auto k=0; k<nz; ++k)
         for (auto j=0; j<ny; ++j) 
         {
            int ii = cijktop2(0, j, k);
            int p_here = parall->convert_taskid_ijk(ii, jcol, kcol);
            if (p_here==taskid) 
            {
               int index = 2*cijktoindex2(0, j, k);
               std::memcpy(tmp,a+index,bsize*sizeof(double));
               parall->ireceive(0, 7, MASTER, 1, idum);
               parall->dsend(0, 9, MASTER, bsize, tmp);
            }
         }
      }
      // delete [] tmp;
   }
}

/********************************
 *                              *
 *         c3db::r_write        *
 *                              *
 ********************************/
 /**
 * @brief Write data to an output unit with support for parallel computing.
 *
 * This function is responsible for writing data to a specified output unit, taking
 * into account parallel computing. It supports two mapping types: "slab mapping"
 * and "hilbert mapping." The function uses various parameters and parallelization
 * techniques to handle data transfer and writing efficiently.
 *
 * @param iunit An integer specifying the output unit to write data to.
 * @param a A pointer to a double array containing the data to be written.
 * @param jcol An integer specifying the column index for parallelization.
 *    
 * @return None.
 * 
 * @note The behavior of this function depends on the mapping type set in the 'maptype'
 *       variable and the parallelization parameters. It supports both "slab mapping"
 *       and "hilbert mapping" for data distribution.
 */
void c3db::r_write(const int iunit, double *a, const int jcol, const int kcol, const bool dotrans)
{  
   int taskid = parall->taskid();
   int taskid_j = parall->taskid_j();
   int taskid_k = parall->taskid_k();
   int np_j = parall->np_j();
   int np_k = parall->np_k();
   int idum[1] = {1};

   /**********************
    **** slab mapping ****
    **********************/
   if (maptype==1) 
   {     
      int bsize = nx*ny;
      double *tmp = new (std::nothrow) double[bsize]();
      
      /**** master node gathers and write to file ****/
      if (taskid==MASTER)
      {
         for (auto k=0; k<nz; ++k)
         {
            int ii = cijktop(0, 0, k);
            int p_from = parall->convert_taskid_ijk(ii,jcol,kcol);
            if (p_from == MASTER)
            {
               int index = cijktoindex(0, 0, k);
               std::memcpy(tmp,a+index,bsize*sizeof(double));
            }
            else 
            {  
               parall->isend(0, 7, p_from, 1, idum);
               parall->dreceive(0, 9, p_from, bsize, tmp);
            }
            dwrite(iunit, tmp, bsize);
         }     
      }     
            
      /**** not master node ****/
      else
      {
         for (auto k=0; k<nz; ++k)
         {
            int index = cijktoindex(0, 0, k);
            int ii = cijktop(0, 0, k);
            int p_here = parall->convert_taskid_ijk(ii, jcol, kcol);
            if (p_here==taskid)
            {
               std::memcpy(tmp,a+index,bsize*sizeof(double));
               parall->ireceive(0, 7, MASTER, 1, idum);
               parall->dsend(0, 9, MASTER, bsize, tmp);
            }
         }
      }
      delete[] tmp;
   }

  /*************************
   **** hilbert mapping ****
   *************************/
   else
   {
      if ((taskid_j==jcol) && (taskid_k==kcol) && dotrans)
      {
         // double *tmp1 = new (std::nothrow) double[2*nfft3d];
         // double *tmp2 = new (std::nothrow) double[2*nfft3d];
         //double *tmp1 = new (std::nothrow) double[2 * nfft3d]();
         //double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
         double *tmp1 = c3db::c3db_tmp1;
         double *tmp2 = c3db::c3db_tmp2;
         r_transpose_ijk(5, a, tmp1, tmp2);

         // delete [] tmp2;
         // delete [] tmp1;
      }

      int bsize = nx;
      double tmp[bsize];

      /**** master node write to file and fetches from other nodes ****/
      if (taskid == MASTER)
      {
         for (auto k=0; k<nz; ++k)
         for (auto j=0; j<ny; ++j)
         {
            int ii = cijktop2(0, j, k);
            int p_from = parall->convert_taskid_ijk(ii, jcol, kcol);
            if (p_from==MASTER)
            {
               int index = cijktoindex2(0, j, k);
               std::memcpy(tmp,a+index,bsize*sizeof(double));
            }
            else
            {
               parall->isend(0, 7, p_from, 1, idum);
               parall->dreceive(0, 9, p_from, bsize, tmp);
            }
            dwrite(iunit, tmp, bsize);
         }
      }
      /**** not master node ****/
      else
      {
         for (auto k=0; k<nz; ++k)
         for (auto j=0; j<ny; ++j)
         {
            int ii = cijktop2(0, j, k);
            int p_here = parall->convert_taskid_ijk(ii, jcol, kcol);
            if (p_here==taskid)
            {
               int index = cijktoindex2(0, j, k);
               std::memcpy(tmp,a+index,bsize*sizeof(double));
               parall->ireceive(0, 7, MASTER, 1, idum);
               parall->dsend(0, 9, MASTER, bsize, tmp);
            }
         }
      }
   }
}




/********************************
 *                              *
 *     c3db::c_write_buffer     *
 *                              *
 ********************************/
void c3db::c_write_buffer(const int iunit, double *a, const int jcol, const int kcol)
{
   int taskid = parall->taskid();
   int taskid_j = parall->taskid_j();
   int taskid_k = parall->taskid_k();

   int buff_count = 2*(nx)*ny*nz;
   double *buffer = new (std::nothrow) double[buff_count]();
   std::memset(buffer,0,buff_count*sizeof(double));
 
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1)
   {
      int bsize = 2*nx*ny;
      for (auto k=0; k<nz; ++k)
      {
         int ii = cijktop(0, 0, k);
         int p_here = parall->convert_taskid_ijk(ii, jcol, kcol);
         if (p_here==taskid)
         {
            int index = 2*cijktoindex(0, 0, k);
            std::memcpy(buffer,a+index,bsize*sizeof(double));
         }
      }
   }

   /*************************
    **** hilbert mapping ****
    *************************/
   else
   {
      int bsize = 2*nx;
      if ((taskid_j==jcol) && (taskid_k==kcol))
      {
         double *tmp1 = c3db::c3db_tmp1;
         double *tmp2 = c3db::c3db_tmp2;
         c_transpose_ijk(5, a, tmp1, tmp2);
      }
      for (int k = 0; k < nz; ++k)
      for (int j = 0; j < ny; ++j)
      {
         int ii = cijktop2(0, j, k);
         int p_here = parall->convert_taskid_ijk(ii, jcol, kcol);
         if (p_here==taskid)
         {
            int index = 2*cijktoindex2(0, j, k);
            std::memcpy(buffer,a+index,bsize*sizeof(double));
         }
      }
   }

   double *buffer2 = new (std::nothrow) double[buff_count]();
   parall->Reduce_Values(1,MASTER,buff_count,buffer,buffer2);
   if (parall->is_master())
      dwrite(iunit,buffer2,buff_count);
   delete [] buffer2;
   delete [] buffer;
}





/********************************
 *                              *
 *         c3db::t_read         *
 *                              *
 ********************************/
 /**
 * @brief Read data from an input unit with support for parallel computing.
 *
 * This function is responsible for reading data from a specified input unit, taking
 * into account parallel computing. It supports two mapping types: "slab mapping"
 * and "hilbert mapping." The function uses various parameters and parallelization
 * techniques to handle data retrieval and distribution efficiently.
 *
 * @param iunit An integer specifying the input unit to read data from.
 * @param a A pointer to a double array where the read data will be stored.
 * @param jcol An integer specifying the column index for parallelization.
 *
 * @return None.
 *
 * @note The behavior of this function depends on the mapping type set in the 'maptype'
 *       variable, the parallelization parameters, and the distribution of data from
 *       the input unit to the specified array.
 */
void c3db::t_read(const int iunit, double *a, const int jcol, const int kcol) 
{
   bool fillcolumn,fillzone;
   int jstart,jend,kstart,kend;
   int taskid = parall->taskid();
   int taskid_j = parall->taskid_j();
   int taskid_k = parall->taskid_k();
   int np_j = parall->np_j();
   int np_k = parall->np_k();
 
   if (jcol<0) 
   {
      jstart = 0;
      jend = np_j - 1;
      fillcolumn = true;
   } 
   else 
   {
      jstart = jend = jcol;
      fillcolumn = (taskid_j == jcol);
   }
   if (kcol<0) 
   {
      kstart = 0;
      kend = np_k - 1;
      fillzone = true;
   } 
   else 
   {
      kstart = kend = kcol;
      fillzone = (taskid_k == kcol);
   }
 
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype==1) 
   {
      int bsize = (nx)*ny;
      double *tmp = new (std::nothrow) double[bsize]();
 
      /**** master node reads from file and distributes ****/
      if (taskid == MASTER)
      {
         for (int k=0; k<nz; ++k) 
         {
            dread(iunit, tmp, bsize);
        
            int index = cijktoindex(0, 0, k);
            int ii = cijktop(0, 0, k);
            for (auto kk=kstart; kk<=kend; ++kk) 
            for (auto jj=jstart; jj<=jend; ++jj) 
            {
               int p_to = parall->convert_taskid_ijk(ii,jj,kk);
               if (p_to == MASTER)
                  std::memcpy(a+index,tmp,bsize*sizeof(double));
               else
                  parall->dsend(0, 9, p_to, bsize, tmp);
            }
         }
      }
 
      /**** not master node ****/
      else if (fillcolumn && fillzone)
      {
         for (int k = 0; k < nz; ++k) 
         {
            int index = cijktoindex(0, 0, k);
            int ii = cijktop(0, 0, k);
            int p_here = parall->convert_taskid_ijk(ii, taskid_j, taskid_k);
            if (p_here == taskid) 
            {
               parall->dreceive(0, 9, MASTER, bsize, tmp);
               std::memcpy(a+index,tmp,bsize*sizeof(double));
            }
         }
      }
      delete[] tmp;
   }
 
   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      int bsize = (nx);
      double tmp[bsize];
     
      /**** master node reads from file and distributes ****/
      if (taskid==MASTER)
      {
         for (auto k=0; k<nz; ++k)
         for (auto j=0; j<ny; ++j) 
         {
            dread(iunit, tmp, bsize);
        
            int index = cijktoindex2(0, j, k);
            int ii = cijktop2(0, j, k);
            for (auto kk=kstart; kk<=kend; ++kk) 
            for (auto jj=jstart; jj<=jend; ++jj) 
            {
               int p_to = parall->convert_taskid_ijk(ii, jj, kk);
               if (p_to == MASTER)
                  std::memcpy(a+index,tmp,bsize*sizeof(double));
               else
                  parall->dsend(0, 9, p_to, bsize, tmp);
            }
         }
      }
     
      /**** not master node ****/
      else if (fillcolumn && fillzone)
      {
         for (auto k=0; k<nz; ++k)
         for (auto j=0; j<ny; ++j) 
         {
            int index = cijktoindex2(0, j, k);
            int ii = cijktop2(0, j, k);
            int p_here = parall->convert_taskid_ijk(ii, taskid_j, taskid_k);
            if (p_here==taskid) 
            {
               parall->dreceive(0, 9, MASTER, bsize, tmp);
               std::memcpy(a+index,tmp,bsize*sizeof(double));
            }
         }
      }
     
      double *tmp1 = c3db::c3db_tmp1;
      double *tmp2 = c3db::c3db_tmp2;
      r_transpose_ijk(4, a, tmp1, tmp2);
   }
}


/********************************
 *                              *
 *         c3db::t_write        *
 *                              *
 ********************************/
 /**
 * @brief Write data to an output unit with support for parallel computing.
 *
 * This function is responsible for writing data to a specified output unit, taking
 * into account parallel computing. It supports two mapping types: "slab mapping"
 * and "hilbert mapping." The function uses various parameters and parallelization
 * techniques to handle data transfer and writing efficiently.
 *
 * @param iunit An integer specifying the output unit to write data to.
 * @param a A pointer to a double array containing the data to be written.
 * @param jcol An integer specifying the column index for parallelization.
 *
 * @return None.
 *
 * @note The behavior of this function depends on the mapping type set in the 'maptype'
 *       variable and the parallelization parameters. It supports both "slab mapping"
 *       and "hilbert mapping" for data distribution.
 */
void c3db::t_write(const int iunit, double *a, const int jcol, const int kcol) 
{
   int taskid = parall->taskid();
   int taskid_j = parall->taskid_j();
   int taskid_k = parall->taskid_k();
   int idum[1] = {1};

   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      int bsize = nx*ny;
      double *tmp = new (std::nothrow) double[bsize]();
      
      /**** master node gathers and write to file ****/
      if (taskid==MASTER)
      {
         for (auto k=0; k<nz; ++k) 
         {
            int ii = cijktop(0, 0, k);
            int p_from = parall->convert_taskid_ijk(ii, jcol, kcol);
            if (p_from==MASTER) 
            {
               int index = cijktoindex(0, 0, k);
               std::memcpy(tmp,a+index,bsize*sizeof(double));
            } 
            else 
            {
               parall->isend(0, 7, p_from, 1, idum);
               parall->dreceive(0, 9, p_from, bsize, tmp);
            }
            dwrite(iunit, tmp, bsize);
         }
      }
      
      /**** not master node ****/
      else
      {
         for (auto k=0; k<nz; ++k) 
         {
            int index = cijktoindex(0, 0, k);
            int ii = cijktop(0, 0, k);
            int p_here = parall->convert_taskid_ijk(ii, jcol, kcol);
            if (p_here == taskid) 
            {
               std::memcpy(tmp,a+index,bsize*sizeof(double));
               parall->ireceive(0, 7, MASTER, 1, idum);
               parall->dsend(0, 9, MASTER, bsize, tmp);
            }
         }
      }
      
      delete[] tmp;
   }

  /*************************
   **** hilbert mapping ****
   *************************/
   else 
   {
      if ((taskid_j==jcol) && (taskid_k==kcol))
      {
         double *tmp1 = c3db::c3db_tmp1;
         double *tmp2 = c3db::c3db_tmp2;
         r_transpose_ijk(5, a, tmp1, tmp2);
      }
     
      int bsize = nx;
      double tmp[bsize];
     
      /**** master node write to file and fetches from other nodes ****/
      if (taskid==MASTER)
      {
         for (auto k=0; k<nz; ++k)
         for (auto j=0; j<ny; ++j) 
         {
            int ii = cijktop2(0, j, k);
            int p_from = parall->convert_taskid_ijk(ii, jcol, kcol);
            if (p_from==MASTER) 
            {
               int index = cijktoindex2(0, j, k);
               std::memcpy(tmp,a+index,bsize*sizeof(double));
            } 
            else 
            {
               parall->isend(0, 7, p_from, 1, idum);
               parall->dreceive(0, 9, p_from, bsize, tmp);
            }
            dwrite(iunit, tmp, bsize);
         }
      } 
      /**** not master node ****/
      else
      {
         for (auto k=0; k<nz; ++k)
         for (auto j=0; j<ny; ++j) 
         {
            int ii = cijktop2(0, j, k);
            int p_here = parall->convert_taskid_ijk(ii, jcol, kcol);
            if (p_here==taskid) 
            {
               int index = cijktoindex2(0, j, k);
               std::memcpy(tmp,a+index,bsize*sizeof(double));
               parall->ireceive(0, 7, MASTER, 1, idum);
               parall->dsend(0, 9, MASTER, bsize, tmp);
            }
         }
      }
   }
}


/********************************
 *                              *
 *     c3db::t_write_buffer     *
 *                              *
 ********************************/
void c3db::t_write_buffer(const int iunit, double *a, const int jcol, const int kcol)
{
   int taskid = parall->taskid();
   int taskid_j = parall->taskid_j();
   int taskid_k = parall->taskid_k();

   int buff_count = (nx)*ny*nz;
   double *buffer = new (std::nothrow) double[buff_count]();
   std::memset(buffer,0,buff_count*sizeof(double));

   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1)
   {
      int bsize = nx*ny;
      for (auto k=0; k<nz; ++k)
      {
         int ii = cijktop(0, 0, k);
         int p_here = parall->convert_taskid_ijk(ii, jcol, kcol);
         if (p_here==taskid)
         {
            int index = cijktoindex(0, 0, k);
            std::memcpy(buffer+k*nx*ny,a+index,nx*ny*sizeof(double));
         }
      }
   }

   /*************************
    **** hilbert mapping ****
    *************************/
   else
   {
      if ((taskid_j==jcol) && (taskid_k==kcol))
      {
         double *tmp1 = c3db::c3db_tmp1;
         double *tmp2 = c3db::c3db_tmp2;
         r_transpose_ijk(5, a, tmp1, tmp2);
      }
      for (auto k=0; k<nz; ++k)
      for (auto j=0; j<ny; ++j)
      {
         int ii = cijktop2(0, j, k);
         int p_here = parall->convert_taskid_ijk(ii, jcol, kcol);
         if (p_here==taskid)
         {
            int index = cijktoindex2(0, j, k);
            std::memcpy(buffer + j*nx + k*nx*ny, a+index, nx*sizeof(double));
         }
      }
   }

   double *buffer2 = new (std::nothrow) double[buff_count]();
   parall->Reduce_Values(1,MASTER,buff_count,buffer,buffer2);
   if (parall->is_master())
      dwrite(iunit,buffer2,buff_count);
   delete [] buffer2;
   delete [] buffer;
}



/**********************************
 *                                *
 * c3db::c_write_buffer_max_final *
 *                                *
 **********************************/
void c3db::c_write_buffer_max_final(const int iunit, int &buff_count, double *buffer)
{
   if (taskid == MASTER)
   {
      if (buff_count > 0) 
         dwrite(iunit,buffer,buff_count);
      buff_count = 0;
   }
}

/********************************
 *                              *
 *   c3db::c_write_buffer_max   *
 *                              *
 ********************************/
void c3db::c_write_buffer_max(const int iunit, double *a, const int jcol, const int kcol,
                              const int buff_max, int &buff_count, double *buffer)
{
   int index, ii, jj, p_from, p_here;
   int taskid = parall->taskid();
   int taskid_j = parall->taskid_j();
   int taskid_k = parall->taskid_j();

   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1)
   {
      int bsize = 2*nx*ny;
      double *tmp = new (std::nothrow) double[bsize]();

      /**** master node gathers and write to file ****/
      if (taskid == MASTER)
      {
         if ((buff_max - buff_count) < (bsize*nz))
         {
            dwrite(iunit,buffer,buff_count);
            buff_count = 0;
         }

         for (auto k=0; k<nz; ++k)
         {
            int ii = cijktop(0, 0, k);
            int p_from = parall->convert_taskid_ijk(ii, jcol, kcol);
            if (p_from==MASTER)
            {
               int index = 2 * cijktoindex(0, 0, k);
               std::memcpy(tmp,a+index,bsize*sizeof(double));
            }
            else
            {
               parall->dreceive(0, 9, p_from, bsize, tmp);
            }

            std::memcpy(buffer+buff_count,tmp,bsize*sizeof(double));
            buff_count += bsize;
         }
      }
      /**** not master node ****/
      else
      {
         for (auto k=0; k<nz; ++k)
         {
            int index = 2 * cijktoindex(0, 0, k);
            int ii = cijktop(0, 0, k);
            int p_here = parall->convert_taskid_ijk(ii, jcol, kcol);
            if (p_here==taskid)
            {
               std::memcpy(tmp,a+index,bsize*sizeof(double));
               parall->dsend(0, 9, MASTER, bsize, tmp);
            }
         }
      }
      delete[] tmp;
   }

  /*************************
   **** hilbert mapping ****
   *************************/
   else
   {
      if ((taskid_j==jcol) && (taskid_k==kcol))
      {
         double *tmp1 = c3db::c3db_tmp1;
         double *tmp2 = c3db::c3db_tmp2;
         c_transpose_ijk(5, a, tmp1, tmp2);
      }

      int bsize = 2*nx;
      double tmp[nx];

      /**** master node write to file and fetches from other nodes ****/
      if (taskid==MASTER)
      {
         if ((buff_max - buff_count) < (bsize*ny*nz))
         {
            dwrite(iunit, buffer, buff_count);
            buff_count = 0;
         }

         for (auto k=0; k<nz; ++k)
         for (auto j=0; j<ny; ++j)
         {
            ii = cijktop2(0, j, k);
            p_from = parall->convert_taskid_ijk(ii, jcol, kcol);
            if (p_from == MASTER)
            {
               int index = 2*cijktoindex2(0, j, k);
               std::memcpy(tmp,a+index,bsize*sizeof(double));
            }
            else
            {
               parall->dreceive(0, 9, p_from, bsize, tmp);
            }
            std::memcpy(buffer+buff_count,tmp,bsize*sizeof(double));
            buff_count += bsize;
         }
      }

      /**** not master node ****/
      else
      {
         for (auto k=0; k<nz; ++k)
         for (auto j=0; j<ny; ++j)
         {
            int ii = cijktop2(0, j, k);
            int p_here = parall->convert_taskid_ijk(ii, jcol, kcol);
            if (p_here==taskid)
            {
               int index = 2*cijktoindex2(0, j, k);
               std::memcpy(tmp,a+index,bsize*sizeof(double));
               parall->dsend(0, 9, MASTER, bsize, tmp);
            }
         }
      }
   }
}


/***********************************************
 *                                            *
 *            c3db::r_formatwrite             *
 *                                            *
 **********************************************/
/**
 * @brief Format and write a double array to a string stream.
 *
 * This function formats and writes a double array `a` to a string stream. The formatted data is organized as a table, where each row represents a 1D slice of the array. The data is formatted with specified width and precision and then written to the string stream.
 *
 * @param a A pointer to the double array to be formatted and written.
 *
 * @return A string containing the formatted data.
 *
 * @note The function performs data formatting and writing based on task-specific conditions, such as the task's ID and parallelization. It efficiently manages the data exchange between the master node and other nodes in a parallel environment.
 */
std::string c3db::r_formatwrite(double *a) 
{
   std::stringstream stream;
 
   int taskid = parall->taskid();
 
   double tmp[nx];
 
   /************************************
    **** slab and hilbert  mappings ****
    ************************************/
 
   /**** master node gathers and write to file ****/
   if (taskid == MASTER) {
     for (auto k = 0; k < nz; ++k)
       for (auto j = 0; j < ny; ++j) {
         for (auto i = 0; i < nx; ++i) {
           int index = cijktoindex2(i, j, k);
           int p_from = cijktop2(i, j, k);
           if (p_from == MASTER)
             tmp[i] = a[index];
           else
             parall->dreceive(0, 189, p_from, 1, tmp + index);
         }
         for (auto i = 0; i < nx; i += 6) {
           for (auto i1 = i; i1 < std::min(i + 6, nx); ++i1)
             stream << Efmt(13, 5) << tmp[i1];
           stream << std::endl;
         }
         // stream << std::endl;
       }
   }
   /**** not master node ****/
   else {
     for (auto k = 0; k < nz; ++k)
       for (auto j = 0; j < ny; ++j) {
         for (auto i = 0; i < nx; ++i) {
           int index = cijktoindex2(i, j, k);
           int p_here = cijktop2(i, j, k);
           if (p_here == taskid)
             parall->dsend(0, 189, MASTER, 1, a + index);
         }
       }
   }
 
   return stream.str();
}

/***********************************************
 *                                            *
 *         c3db::r_formatwrite_reverse        *
 *                                            *
 **********************************************/
/**
 * @brief Format and write a double array to a string stream in reverse order.
 *
 * This function is similar to `c3db::r_formatwrite`, but it formats and writes the given double array `a` in reverse order (along the k dimension) to a string stream. The formatted data is written in rows and columns, with each row containing up to 6 elements from the array. The output is intended to be used for debugging or saving data to a file.
 *
 * @param a A pointer to the double array to be formatted and written in reverse order.
 *
 * @return A string containing the formatted data in reverse order along the k dimension.
 *
 * @note The function performs a data gathering operation among parallel tasks to construct the formatted string. The size and structure of the output string are determined by the dimensions of the array and the number of tasks involved in the parallel computation.
 */
std::string c3db::r_formatwrite_reverse(double *a) 
{
   std::stringstream stream;
 
   int taskid = parall->taskid();
   double tmp[nz];
 
   /************************************
    **** slab and hilbert  mappings ****
    ************************************/
 
   /**** master node gathers and write to file ****/
   if (taskid == MASTER) {
     for (auto i = 0; i < nx; ++i)
       for (auto j = 0; j < ny; ++j) {
         for (auto k = 0; k < nz; ++k) {
           int index = cijktoindex2(i, j, k);
           int p_from = cijktop2(i, j, k);
           if (p_from == MASTER)
             tmp[k] = a[index];
           else
             parall->dreceive(0, 189, p_from, 1, tmp + k);
         }
         for (auto k = 0; k < nz; k += 6) {
           for (auto k1 = k; k1 < std::min(k + 6, nz); ++k1)
             stream << Efmt(13, 5) << tmp[k1];
           stream << std::endl;
         }
         // stream << std::endl;
       }
   }
   /**** not master node ****/
   else {
     for (auto i = 0; i < nx; ++i)
       for (auto j = 0; j < ny; ++j) {
         for (auto k = 0; k < nz; ++k) {
           int index = cijktoindex2(i, j, k);
           int p_here = cijktop2(i, j, k);
           if (p_here == taskid)
             parall->dsend(0, 189, MASTER, 1, a + index);
         }
       }
   }
 
   return stream.str();
}

void c3db::cshift1_fftb(const int n1, const int n2, const int n3, const int n4,
                        double *a) {
  int i, j, indx;
  indx = 1;
  for (j = 0; j < (n2 * n3 * n4); ++j) {
    for (i = 2; i <= n1; ++i) {
      a[indx + i - 2] = a[indx + i - 1];
    }
    indx += (n1 + 2);
  }
}
void c3db::cshift_fftf(const int n1, const int n2, const int n3, const int n4,
                       double *a) {
  int i, j, indx;
  indx = 1;
  for (j = 0; j < (n2 * n3 * n4); ++j) {
    for (i = n1; i >= 2; --i) {
      a[indx + i - 1] = a[indx + i - 2];
    }
    a[indx + 1 - 1] = 0.0;
    a[indx + n1 + 1 - 1] = 0.0;
    indx += (n1 + 2);
  }
}

static void cshift_fftf_ab(const int n1, const int n2, const int n3, const int n4, double *a, double *b) 
{
   int i, j, indx;
   indx = 0;
   for (j = 0; j < (n2 * n3 * n4); ++j) 
   {
      for (i = n1; i > 0; --i) 
      {
         b[indx + i] = a[indx + i - 1];
      }
      b[indx + 1] = 0.0;
      b[indx + n1 + 1] = 0.0;
      indx += (n1 + 2);
   }
}

void c3db::zeroend_fftb(const int n1, const int n2, const int n3, const int n4, double *a) 
{
   int i, indx;
   indx = n1 + 1;
   for (i = 0; i < (n2 * n3 * n4); ++i) 
   {
      a[indx - 1] = 0.0;
      a[indx + 1 - 1] = 0.0;
      indx += (n1 + 2);
   }
}

/********************************
 *                              *
 *         c3db::cr_fft3d       *
 *                              *
 ********************************/
void c3db::cr_fft3d(double *a) 
{

   nwpw_timing_function ftime(1);
   int i, j, k, jj, kk, q, indx, indx0, nxy, nxz, nn, shift;
   double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
   double *tmp3 = new (std::nothrow) double[2 * nfft3d]();
 
   nxy = nx * ny;
   nxz = nx * nz;
 
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      std::memset(tmp2,0,2*nfft3d*sizeof(double));

      /***************************************************
       ***     do fft along kz dimension               ***
       ***     A(kx,nz,ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
       ***************************************************/
      indx0 = 0;
      nn = 0;
      for (q = 0; q < nq; ++q) 
      {
         for (i = 0; i < nx; ++i) 
         {
            kk = 0;
            indx = 2 * i + indx0;
            shift = 2 * nz * nn;
            for (k = 0; k < nz; ++k) 
            {
               tmp2[kk   + shift] = a[indx];
               tmp2[kk+1 + shift] = a[indx + 1];
               kk += 2;
               indx += nx;
            }
            ++nn;
         }
         indx0 += nxz;
      }
     
      mygdevice.batch_cfftz_tmpz(fft_tag,false, nz, nn, 2*nfft3d, tmp2, tmpz);
     
      indx0 = 0;
      nn = 0;
      for (q = 0; q < nq; ++q) 
      {
         for (i = 0; i < nx; ++i) 
         {
            kk = 0;
            indx = 2 * i + indx0;
            shift = 2 * nz * nn;
            for (k = 0; k < nz; ++k) 
            {
               a[indx]   = tmp2[kk   + shift];
               a[indx+1] = tmp2[kk+1 + shift];
               kk += 2;
               indx += nx;
            }
            ++nn;
         }
         indx0 += nxz;
      }
     
      /***********************************************
       ***         Do a transpose of A             ***
       ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
       ************************************************/
      c_transpose_jk(a, tmp2, tmp3);

      /*************************************************
       ***        do fft along ky dimension          ***
       ***    A(kx,ny,nz) <- fft1d^(-1)[A(kx,ky,nz)] ***
       *************************************************/
      indx0 = 0;
      nn = 0;
      for (q = 0; q < nq; ++q) 
      {
         for (i = 0; i < nx; ++i) 
         {
            jj = 0;
            indx = 2 * i + indx0;
            shift = 2 * ny * nn;
            for (j = 0; j < ny; ++j) 
            {
               tmp2[jj + shift] = a[indx];
               tmp2[jj + 1 + shift] = a[indx + 1];
               jj += 2;
               indx += nx;
            }
            ++nn;
         }
         indx0 += nxy;
      }
     
      mygdevice.batch_cffty_tmpy(fft_tag,false, ny, nn, 2*nfft3d, tmp2, tmpy);
     
      indx0 = 0;
      nn = 0;
      for (q = 0; q < nq; ++q) 
      {
         for (i = 0; i < nx; ++i) 
         {
            jj = 0;
            indx = 2 * i + indx0;
            shift = 2 * ny * nn;
            for (j = 0; j < ny; ++j) 
            {
               a[indx]   = tmp2[jj   + shift];
               a[indx+1] = tmp2[jj+1 + shift];
               jj += 2;
               indx += nx;
            }
            ++nn;
         }
         indx0 += nxy;
      }

      mygdevice.batch_cfftx_tmpx(fft_tag,false, nx, ny * nq, 2*nfft3d, a, tmpx);
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
      mygdevice.batch_cfftz_tmpz(fft_tag,false, nz, nq3, 2*nfft3d, a, tmpz);
      
      c_transpose_ijk(2, a, tmp2, tmp3);
     
      /************************************************
       ***     do fft along ky dimension            ***
       ***   A(ny,nz,kx) <- fft1d^(-1)[A(ky,nz,kx)] ***
       ************************************************/
      mygdevice.batch_cffty_tmpy(fft_tag,false, ny, nq2, 2*nfft3d, a, tmpy);
     
      c_transpose_ijk(3, a, tmp2, tmp3);
     
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      mygdevice.batch_cfftx_tmpx(fft_tag,false, nx, nq1, 2*nfft3d, a, tmpx);
     
   }
 
   delete[] tmp3;
   delete[] tmp2;
}

/********************************
 *                              *
 *         c3db::rc_fft3d       *
 *                              *
 ********************************/
void c3db::rc_fft3d(double *a) 
{

   nwpw_timing_function ftime(1);
   int i, j, k, jj, kk, q, indx, indx0, nxy, nxz, nn, shift;
   double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
   double *tmp3 = new (std::nothrow) double[2 * nfft3d]();
 
   nxy = nx * ny;
   nxz = nx * nz;
 
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      /********************************************
       ***     do fft along nx dimension        ***
       ***   A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]  ***
       ********************************************/
      mygdevice.batch_cfftx_tmpx(fft_tag,true, nx, ny*nq, 2*nfft3d, a, tmpx);

      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
       ********************************************/
      indx0 = 0;
      nn = 0;
      for (q = 0; q < nq; ++q) 
      {
         for (i = 0; i < nx; ++i) 
         {
            jj = 0;
            indx = 2 * i + indx0;
            shift = 2 * ny * nn;
            for (j = 0; j < ny; ++j) 
            {
               tmp2[jj   + shift] = a[indx];
               tmp2[jj+1 + shift] = a[indx + 1];
               jj += 2;
               indx += nx;
            }
            ++nn;
         }
         indx0 += nxy;
      }
     
      mygdevice.batch_cffty_tmpy(fft_tag,true, ny, nn, 2*nfft3d, tmp2, tmpy);
     
      indx0 = 0;
      nn = 0;
      for (q = 0; q < nq; ++q) 
      {
         for (i = 0; i < nx; ++i) 
         {
            jj = 0;
            indx = 2 * i + indx0;
            shift = 2 * ny * nn;
            for (j = 0; j < ny; ++j) 
            {
               a[indx]     = tmp2[jj   + shift];
               a[indx + 1] = tmp2[jj+1 + shift];
               jj += 2;
               indx += nx;
            }
            ++nn;
         }
         indx0 += nxy;
      }
     
     
      /********************************************
       ***         Do a transpose of A          ***
       ***      A(ky,nz,ky) <- A(kx,ky,nz)      ***
       ********************************************/
      c_transpose_jk(a, tmp2, tmp3);
     
      /********************************************
       ***     do fft along nz dimension        ***
       ***   A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]  ***
       ********************************************/
      indx0 = 0;
      nn = 0;
      for (q = 0; q < nq; ++q) 
      {
         for (i = 0; i < nx; ++i) 
         {
            kk = 0;
            indx = 2 * i + indx0;
            shift = 2 * nz * nn;
            for (k = 0; k < nz; ++k) 
            {
               tmp2[kk   + shift] = a[indx];
               tmp2[kk+1 + shift] = a[indx + 1];
               kk += 2;
               indx += nx;
            }
            ++nn;
         }
         indx0 += nxz;
      }
     
      mygdevice.batch_cfftz_tmpz(fft_tag,true, nz, nn, 2*n2ft3d, tmp2, tmpz);
     
      indx0 = 0;
      nn = 0;
      for (q = 0; q < nq; ++q) 
      {
         for (i = 0; i < nx; ++i) 
         {
            kk = 0;
            indx = 2 * i + indx0;
            shift = 2 * nz * nn;
            for (k = 0; k < nz; ++k) 
            {
               a[indx]   = tmp2[kk   + shift];
               a[indx+1] = tmp2[kk+1 + shift];
               kk += 2;
               indx += nx;
            }
            ++nn;
         }
         indx0 += nxz;
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
      mygdevice.batch_cfftx_tmpx(fft_tag,true, nx, nq1, 2*nfft3d, a, tmpx);
      
      c_transpose_ijk(0, a, tmp2, tmp3);
     
      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
       ********************************************/
      mygdevice.batch_cffty_tmpy(fft_tag,true, ny, nq2, 2*nfft3d, a, tmpy);
     
      c_transpose_ijk(1, a, tmp2, tmp3);
     
      /********************************************
       ***     do fft along nz dimension        ***
       ***   A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]  ***
       ********************************************/
      mygdevice.batch_cfftz_tmpz(fft_tag,true, nz, nq3, 2*n2ft3d, a, tmpz);
   }
 
   delete[] tmp3;
   delete[] tmp2;
}


/**************************************
 *                                    *
 *    c3db::c_ptranspose1_jk_start    *
 *                                    *
 **************************************/
void c3db::c_ptranspose1_jk_start(const int nb, double *a, double *tmp1,
                                  double *tmp2, const int request_indx,
                                  const int msgtype) 
{
   int msglen;
   int n1 = p_i1_start[nb][0][np];
 
   c_aindexcopy(n1, p_iq_to_i1[nb][0], a, tmp1);
 
   /* it = 0, transpose data on same thread */
   msglen = 2*(p_i2_start[nb][0][1] - p_i2_start[nb][0][0]);
   std::memcpy(tmp2 + 2*p_i2_start[nb][0][0], tmp1 + 2*p_i1_start[nb][0][0], msglen*sizeof(double));
 
   /* receive packed array data */
   for (auto it=1; it<np; ++it) 
   {
      /* synchronous receive of tmp */
      int proc_from = (taskid - it + np) % np;
      msglen = 2*(p_i2_start[nb][0][it+1] - p_i2_start[nb][0][it]);
      if (msglen > 0)
         parall->adreceive(request_indx, msgtype, proc_from, msglen, tmp2 + 2*p_i2_start[nb][0][it]);
     // parall->adreceive(request_indx,msgtype,proc_from,msglen,&tmp2[2*p_i2_start[nb][0][it]]);
   }
   for (auto it=1; it<np; ++it) 
   {
      int proc_to = (taskid + it) % np;
      msglen = 2*(p_i1_start[nb][0][it+1] - p_i1_start[nb][0][it]);
      if (msglen > 0)
         parall->adsend(request_indx, msgtype, proc_to, msglen, tmp1 + 2*p_i1_start[nb][0][it]);
     // parall->adsend(request_indx,msgtype,proc_to,msglen,&tmp1[2*p_i1_start[nb][0][it]]);
   }
}

/**************************************
 *                                    *
 *    c3db::c_ptranspose1_jk_end      *
 *                                    *
 **************************************/
void c3db::c_ptranspose1_jk_end(const int nb, double *a, double *tmp2, const int request_indx) 
{
   parall->awaitall(request_indx);

   int n2 = p_i2_start[nb][0][np];
   c_bindexcopy(n2, p_iq_to_i2[nb][0], tmp2, a);
   c_bindexzero(nfft3d - n2, p_iz_to_i2[nb][0], a);
}

/**************************************
 *                                    *
 *    c3db::c_ptranspose2_jk_start    *
 *                                    *
 **************************************/
void c3db::c_ptranspose2_jk_start(const int nb, double *a, double *tmp1,
                                  double *tmp2, const int request_indx,
                                  const int msgtype) 
{
   int msglen;
   int n1 = p_j1_start[nb][0][np];
 
   c_aindexcopy(n1, p_jq_to_i1[nb][0], a, tmp1);
 
   /* it = 0, transpose data on same thread */
   msglen = 2*(p_j2_start[nb][0][1] - p_j2_start[nb][0][0]);
   std::memcpy(tmp2 + 2 * p_j2_start[nb][0][0], tmp1 + 2 * p_j1_start[nb][0][0],
               msglen * sizeof(double));
 
   /* receive packed array data */
   for (auto it=1; it<np; ++it) 
   {
      /* synchronous receive of tmp */
      int proc_from = (taskid - it + np) % np;
      msglen = 2 * (p_j2_start[nb][0][it+1] - p_j2_start[nb][0][it]);
      if (msglen > 0)
         parall->adreceive(request_indx, msgtype, proc_from, msglen, tmp2 + 2*p_j2_start[nb][0][it]);
   }
   for (auto it = 1; it < np; ++it) 
   {
      int proc_to = (taskid + it) % np;
      msglen = 2 * (p_j1_start[nb][0][it + 1] - p_j1_start[nb][0][it]);
      if (msglen > 0)
         parall->adsend(request_indx, msgtype, proc_to, msglen, tmp1 + 2*p_j1_start[nb][0][it]);
   }
}

/**************************************
 *                                    *
 *    c3db::c_ptranspose2_jk_end      *
 *                                    *
 **************************************/
void c3db::c_ptranspose2_jk_end(const int nb, double *a, double *tmp2,
                                const int request_indx) 
{
   parall->awaitall(request_indx);

   int n2 = p_j2_start[nb][0][np];
   c_bindexcopy(n2, p_jq_to_i2[nb][0], tmp2, a);
   c_bindexzero(nfft3d - n2, p_jz_to_i2[nb][0], a);
}

/**************************************
 *                                    *
 *    c3db::c_ptranspose_ijk_start    *
 *                                    *
 **************************************/
void c3db::c_ptranspose_ijk_start(const int nb, const int op, double *a,
                                  double *tmp1, double *tmp2,
                                  const int request_indx, const int msgtype) 
{
   int nnfft3d, msglen;
 
   int n1 = p_i1_start[nb][op][np];
 
   /* pack a array - tmp1->tmp2 */
   c_aindexcopy(n1, p_iq_to_i1[nb][op], a, tmp1);
 
   /* it = 0, transpose data on same thread  - tmp2->tmp1*/
   msglen = 2*(p_i2_start[nb][op][1] - p_i2_start[nb][op][0]);
   std::memcpy(tmp2 + 2 * p_i2_start[nb][op][0],
               tmp1 + 2 * p_i1_start[nb][op][0], msglen * sizeof(double));
 
   /* receive packed array data */
   for (auto it=1; it<np; ++it) 
   {
      /* synchronous receive of tmp */
      int proc_from = (taskid - it + np) % np;
      msglen = 2*(p_i2_start[nb][op][it+1] - p_i2_start[nb][op][it]);
      if (msglen > 0)
        parall->adreceive(request_indx, msgtype, proc_from, msglen,
                          tmp2 + 2 * p_i2_start[nb][op][it]);
      // parall->adreceive(request_indx,msgtype,proc_from,msglen,&tmp2[2*p_i2_start[nb][op][it]]);
   }
   for (auto it=1; it<np; ++it) 
   {
      int proc_to = (taskid + it) % np;
      msglen = 2*(p_i1_start[nb][op][it+1] - p_i1_start[nb][op][it]);
      if (msglen > 0)
        parall->adsend(request_indx, msgtype, proc_to, msglen,
                       tmp1 + 2 * p_i1_start[nb][op][it]);
      // parall->adsend(request_indx,msgtype,proc_to,msglen,&tmp1[2*p_i1_start[nb][op][it]]);
   }
}


/**************************************
 *                                    *
 *    c3db::c_ptranspose_ijk_end      *
 *                                    *
 **************************************/
void c3db::c_ptranspose_ijk_end(const int nb, const int op, double *a,
                                double *tmp2, const int request_indx) 
{
   int n2 = p_i2_start[nb][op][np];
   int n3 = p_iz_to_i2_count[nb][op];

   /* wait for completion of mp_send, also do a sync */
   parall->awaitall(request_indx);

   /* unpack a array */
   c_bindexcopy(n2, p_iq_to_i2[nb][op], tmp2, a);
   c_bindexzero(n3, p_iz_to_i2[nb][op], a);
}

/********************************
 *                              *
 *    c3db::c_ptranspose1_jk    *
 *                              *
 ********************************/
void c3db::c_ptranspose1_jk(const int nb, double *a, double *tmp1, double *tmp2) 
{
   int msglen;
 
   int n1 = p_i1_start[nb][0][np];
   int n2 = p_i2_start[nb][0][np];
 
   parall->astart(1, np);
 
   c_aindexcopy(n1, p_iq_to_i1[nb][0], a, tmp1);
 
   /* it = 0, transpose data on same thread */
   msglen = 2*(p_i2_start[nb][0][1] - p_i2_start[nb][0][0]);
   // int one=1;
   // DCOPY_PWDFT(msglen,&(tmp1[2*p_i1_start[nb][0][0]]),one,&(tmp2[2*p_i2_start[nb][0][0]]),one);
   std::memcpy(tmp2 + 2*p_i2_start[nb][0][0],tmp1+2*p_i1_start[nb][0][0],msglen*sizeof(double));
 
   /* receive packed array data */
   for (auto it=1; it<np; ++it) 
   {
      /* synchronous receive of tmp */
      int proc_from = (taskid - it + np) % np;
      msglen = 2*(p_i2_start[nb][0][it+1] - p_i2_start[nb][0][it]);
      if (msglen > 0)
         parall->adreceive(1, 1, proc_from, msglen, &tmp2[2 * p_i2_start[nb][0][it]]);
   }

   for (auto it=1; it<np; ++it) 
   {
      int proc_to = (taskid + it) % np;
      msglen = 2 * (p_i1_start[nb][0][it + 1] - p_i1_start[nb][0][it]);
      if (msglen > 0)
         parall->dsend(1, 1, proc_to, msglen, &tmp1[2 * p_i1_start[nb][0][it]]);
   }
   parall->aend(1);
 
   c_bindexcopy(n2, p_iq_to_i2[nb][0], tmp2, a);
   c_bindexzero(nfft3d - n2, p_iz_to_i2[nb][0], a);
}

/********************************
 *                              *
 *    c3db::c_ptranspose2_jk    *
 *                              *
 ********************************/
void c3db::c_ptranspose2_jk(const int nb, double *a, double *tmp1, double *tmp2) 
{
   int msglen;
 
   int n1 = p_j1_start[nb][0][np];
   int n2 = p_j2_start[nb][0][np];
 
   parall->astart(1,np);
 
   c_aindexcopy(n1,p_jq_to_i1[nb][0],a,tmp1);
 
   /* it = 0, transpose data on same thread */
   msglen = 2*(p_j2_start[nb][0][1]-p_j2_start[nb][0][0]);
   std::memcpy(tmp2+2*p_j2_start[nb][0][0], tmp1+2*p_j1_start[nb][0][0], msglen*sizeof(double));
 
   /* receive packed array data */
   for (auto it=1; it<np; ++it) 
   {
      /* synchronous receive of tmp */
      int proc_from = (taskid - it + np) % np;
      msglen = 2*(p_j2_start[nb][0][it + 1] - p_j2_start[nb][0][it]);
      if (msglen > 0)
         parall->adreceive(1,1,proc_from,msglen,&tmp2[2*p_j2_start[nb][0][it]]);
   }
   for (auto it=1; it<np; ++it) 
   {
      int proc_to = (taskid + it) % np;
      msglen = 2*(p_j1_start[nb][0][it+1] - p_j1_start[nb][0][it]);
      if (msglen > 0)
         parall->dsend(1,1,proc_to,msglen,&tmp1[2*p_j1_start[nb][0][it]]);
   }

   parall->aend(1);
 
   c_bindexcopy(n2, p_jq_to_i2[nb][0],tmp2,a);
   c_bindexzero(nfft3d-n2,p_jz_to_i2[nb][0],a);
}


/********************************
 *                              *
 *    c3db::c_transpose_jk      *
 *                              *
 ********************************/
void c3db::c_transpose_jk(double *a, double *tmp1, double *tmp2)
{
   int msglen;

   parall->astart(1, np);

   c_bindexcopy(nfft3d,iq_to_i1[0],a,tmp1);

   /* it = 0, transpose data on same thread */
   msglen = 2*(i2_start[0][1]-i2_start[0][0]);
   // int one=1;
   // DCOPY_PWDFT(msglen,&(tmp1[2*i1_start[0][0]]),one,&(tmp2[2*i2_start[0][0]]),one);
   std::memcpy(tmp2+2*i2_start[0][0],tmp1+2*i1_start[0][0],msglen*sizeof(double));

   /* receive packed array data */
   for (auto it=1; it<np; ++it)
   {
      /* synchronous receive of tmp */
      auto proc_from = (taskid - it + np) % np;
      msglen = 2*(i2_start[0][it+1] - i2_start[0][it]);
      if (msglen > 0)
         parall->adreceive(1,1,proc_from,msglen,&tmp2[2*i2_start[0][it]]);
   }
   for (auto it=1; it<np; ++it)
   {
      auto proc_to = (taskid + it) % np;
      msglen = 2*(i1_start[0][it+1] - i1_start[0][it]);
      if (msglen > 0)
         parall->dsend(1,1,proc_to,msglen,&tmp1[2*i1_start[0][it]]);
   }
   parall->aend(1);
   c_aindexcopy(nfft3d,iq_to_i2[0],tmp2,a);
}


/********************************
 *                              *
 *    c3db::c_ptranspose_ijk    *
 *                              *
 ********************************/
void c3db::c_ptranspose_ijk(const int nb, const int op, double *a, double *tmp1, double *tmp2) 
{
   int nnfft3d, msglen;
 
   int n1 = p_i1_start[nb][op][np];
   int n2 = p_i2_start[nb][op][np];
   int n3 = p_iz_to_i2_count[nb][op];
 
   parall->astart(1, np);
 
   /* pack a array */
   c_aindexcopy(n1, p_iq_to_i1[nb][op], a, tmp1);
 
   /* it = 0, transpose data on same thread */
   msglen = 2*(p_i2_start[nb][op][1] - p_i2_start[nb][op][0]);
   std::memcpy(tmp2 + 2*p_i2_start[nb][op][0], tmp1 + 2*p_i1_start[nb][op][0], msglen*sizeof(double));
 
   /* receive packed array data */
   for (auto it=1; it<np; ++it) 
   {
      /* synchronous receive of tmp */
      int proc_from = (taskid - it + np) % np;
      msglen = 2*(p_i2_start[nb][op][it + 1] - p_i2_start[nb][op][it]);
      if (msglen > 0)
         parall->adreceive(1,1, proc_from, msglen, tmp2+2*p_i2_start[nb][op][it]);
   }
   for (auto it=1; it<np; ++it) 
   {
      int proc_to = (taskid + it) % np;
      msglen = 2*(p_i1_start[nb][op][it + 1] - p_i1_start[nb][op][it]);
      if (msglen > 0)
         parall->dsend(1,1,proc_to, msglen, tmp1+2*p_i1_start[nb][op][it]);
   }
 
   /* wait for completion of mp_send, also do a sync */
   parall->aend(1);
 
   /* unpack a array */
   c_bindexcopy(n2,p_iq_to_i2[nb][op],tmp2,a);
   c_bindexzero(n3,p_iz_to_i2[nb][op],a);
}

/********************************
 *                              *
 *    c3db::c_transpose_ijk     *
 *                              *
 ********************************/
void c3db::c_transpose_ijk(const int op, double *a, double *tmp1, double *tmp2) 
{
   int nnfft3d,msglen;
 
   parall->astart(1,np);
 
   /* pack a array */
   if ((op == 0) || (op == 4)) nnfft3d = (nx)*nq1;
   if ((op == 1) || (op == 3)) nnfft3d = (ny)*nq2;
   if ((op == 2) || (op == 5)) nnfft3d = (nz)*nq3;
   c_bindexcopy(nnfft3d,iq_to_i1[op],a,tmp1);
 
   /* it = 0, transpose data on same thread */
   msglen = 2*(i2_start[op][1] - i2_start[op][0]);
   std::memcpy(tmp2 + 2*i2_start[op][0], tmp1 + 2*i1_start[op][0], msglen*sizeof(double));
 
   /* receive packed array data */
   for (auto it=1; it<np; ++it) 
   {
      /* synchronous receive of tmp */
      int proc_from = (taskid - it + np) % np;
      msglen = 2*(i2_start[op][it+1]-i2_start[op][it]);
      if (msglen > 0)
         parall->adreceive(1, 1, proc_from, msglen, tmp2 + 2*i2_start[op][it]);
   }
   for (auto it=1; it<np; ++it) 
   {
      int proc_to = (taskid + it) % np;
      msglen = 2*(i1_start[op][it+1]-i1_start[op][it]);
      if (msglen > 0)
         parall->dsend(1, 1, proc_to, msglen, tmp1 + 2*i1_start[op][it]);
   }
 
   /* wait for completion of mp_send, also do a sync */
   parall->aend(1);
 
   /* unpack a array */
   if ((op == 3) || (op == 5)) nnfft3d = (nx)*nq1;
   if ((op == 0) || (op == 2)) nnfft3d = (ny)*nq2;
   if ((op == 1) || (op == 4)) nnfft3d = (nz)*nq3;
   c_aindexcopy(nnfft3d,iq_to_i2[op],tmp2,a);
}


/********************************
 *                              *
 *         c3db::c_setpw        *
 *                              *
 ********************************/
/**
 * @brief Set complex values in an array based on a specific pattern.
 *
 * This function sets complex values in the array 'a' based on the specified 'filling' indices and 'cvalue'.
 *
 * @param filling   An integer array representing the indices (i, j, k) where the complex value should be set.
 * @param cvalue    A complex value (double[2]) to be assigned to the array 'a' at the specified indices.
 * @param a         The complex array where values will be set.
 *
 * @note The function assumes that memory for 'a' has been allocated appropriately by the caller.
 *
 * @details
 * - The function sets complex values in the array 'a' based on the indices provided in the 'filling' array and the 'cvalue'.
 * - The 'filling' array should contain the indices (i, j, k) specifying where the complex value should be assigned in 'a'.
 * - The complex value 'cvalue' is assigned to the array 'a' at the specified indices.
 * - The function also sets the conjugate value in a specific case where i == 0, j != 0, and k != 0.
 *
 * @param i         The 'i' index indicating the position where the complex value should be set.
 * @param j         The 'j' index indicating the position where the complex value should be set.
 * @param k         The 'k' index indicating the position where the complex value should be set.
 * @param indx      The calculated index in the array 'a' based on the provided (i, j, k) indices.
 * @param p         The task ID associated with the calculated (i, j, k) indices.
 *
 */
void c3db::c_setpw(const int filling[], const double *cvalue, double *a) 
{
   int i = filling[0];
   int j = filling[1];
   int k = filling[2];
 
   int indx = cijktoindex(i, j, k);
   int p = cijktop(i, j, k);
   if (p == parall->taskid_i()) 
   {
      a[2 * indx] = cvalue[0];
      a[2 * indx + 1] = cvalue[1];
   }
 
   /* set the conjugate on the i==0 plane */
   if ((i == 0) && (j != 0) && (k != 0)) 
   {
      int jc = j;
      if (jc > (ny / 2))
        jc -= ny;
      jc = -jc;
      if (jc < 0)
        jc += ny;
     
      int kc = k;
      if (kc > (nz / 2))
        kc -= nz;
      kc = -kc;
      if (kc < 0)
        kc += nz;
      indx = cijktoindex(i, jc, kc);
      p = cijktop(i, jc, kc);
      if (p == parall->taskid_i()) 
      {
         a[2 * indx] = cvalue[0];
         a[2 * indx + 1] = -cvalue[1];
      }
   }
}

/********************************
 *                              *
 *        c3db::c_addrandom     *
 *                              *
 ********************************/
/**
 * @brief Add random values to the elements of an array.
 *
 * This function adds random values to each element of the input array 'a'. The random values are generated between -0.25
 * and 0.25, scaled by a factor based on the size of the grid (nfft3d), and then added to the corresponding elements of 'a'.
 *
 * @param a  The array of elements to which random values will be added.
 *
 * @note The function assumes that memory for 'a' has been allocated appropriately by the caller.
 *
 * @details
 * - The function uses a linear congruential pseudorandom number generator to generate random values.
 * - The scaling factor 'fac' is calculated as 1.0 divided by the square root of the total number of elements in the grid ('nfft3d').
 * - For each element of the array 'a', a random value between -0.25 and 0.25, scaled by 'fac', is generated and added.
 *
 * @par Pseudorandom Number Generator:
 * - The pseudorandom number generator used by this function is based on the 'util_random' function.
 * - It generates pseudorandom numbers in the range [0, 1).
 * - The generated random values are then transformed to the range [-0.25, 0.25].
 *
 */
void c3db::c_addrandom(double *a) 
{
   double fac = 1.0 / sqrt(1.0 * nfft3d);
   for (auto i = 0; i < n2ft3d; ++i)
      a[i] += fac * (0.50 - util_random(0));
}

/********************************
 *                              *
 *        c3db::r_setrandom     *
 *                              *
 ********************************/
/**
 * @brief Set the elements of an array to random values.
 *
 * This function sets the elements of the input array 'a' to random values between -0.25 and 0.25, scaled by a factor
 * based on the size of the grid (nfft3d).
 *
 * @param a  The array of elements to be filled with random values.
 *
 * @note The function assumes that memory for 'a' has been allocated appropriately by the caller.
 *
 * @details
 * - The function uses a linear congruential pseudorandom number generator to generate random values.
 * - The scaling factor 'fac' is calculated as 1.0 divided by the total number of elements in the grid ('nfft3d').
 * - Each element of the array 'a' is set to a random value between -0.25 and 0.25, scaled by 'fac'.
 *
 * @par Pseudorandom Number Generator:
 * - The pseudorandom number generator used by this function is based on the 'util_random' function.
 * - It generates pseudorandom numbers in the range [0, 1).
 * - The generated random values are then transformed to the range [-0.25, 0.25].
 *
 */
void c3db::r_setrandom(double *a) 
{
   double fac = 1.0 / (1.0 * nfft3d);
   for (auto i = 0; i < n2ft3d; ++i)
      a[i] = fac*(0.50 - util_random(0));
}

/********************************
 *                              *
 *        c3db::hr2r_expand     *
 *                              *
 ********************************/
// expands a grid that is ((nx/2+2),ny/2,nz/2) to (nx,ny,nz)
/**
 * @brief Expand a grid from ((nx/2+2), ny/2, nz/2) to (nx, ny, nz).
 *
 * This function performs a grid expansion operation from the input array 'ah' to the output array 'a'. It is designed
 * for a specific mapping type (maptype) and grid sizes.
 *
 * @param ah  The input contracted grid array of size ((nx/2), ny/2, nz/2).
 * @param a   The output expanded grid array of size (nx, ny, nz).
 *
 * @note The function assumes that memory for 'ah' and 'a' has been allocated appropriately by the caller.
 * @note The function uses 'maptype' to determine the appropriate expansion method.
 *
 * @param maptype The mapping type used for expansion:
 * - If 'maptype' is 1, a 3D expansion is performed.
 * - If 'maptype' is not 1, a 2D expansion is performed.
 *
 * @details
 * - For maptype = 1 (3D expansion):
 *   - The input contracted grid 'ah' with dimensions ((nx/2+2), ny/2, nz/2) is expanded to 'a' with dimensions (nx, ny, nz).
 *
 * - For maptype != 1 (2D expansion):
 *   - The input contracted grid 'ah' with dimensions ((nx/2+2), ny/2) is expanded to 'a' with dimensions (nx, ny).
 *
 * @par Implementation Details:
 * - The function uses nested loops to iterate over the input grid and perform the expansion according to the specified maptype.
 * - The resulting expanded grid 'a' is filled with the expanded values.
 *
 */
void c3db::hr2r_expand(const double *ah, double *a) 
{
   std::memset(a, 0, n2ft3d * sizeof(double));
   if (maptype == 1) 
   {
      int nxh = nx/2;
      int nqh = nq/2;
      int nyh = ny/2;
      for (auto j=0; j<nyh; ++j)
      for (auto q=0; q<nqh; ++q)
      for (auto i=0; i<nxh; ++i)
         a[i+j*(nx) + q*(nx)*ny] = ah[i + j*(nxh) + q*(nxh)*nyh];
   } 
   else 
   {
      int nxh = nx/2;
      int nq1h = nq1/4;
      for (auto q=0; q<nq1h; ++q)
      for (auto i=0; i<nxh; ++i)
         a[i+q*(nx)] = ah[i + q*(nxh)];
   }
}

/********************************
 *                              *
 *      c3db::r2hr_contract     *
 *                              *
 ********************************/
// contracts to a grid that is (nx,ny,nz) --> ((nx)/2,ny/2,nz/2)
/**
 * @brief Contract a grid from (nx, ny, nz) to ((nx) / 2, ny / 2, nz / 2).
 *
 * This function performs a grid contraction operation from the input array 'a' to the output array 'ah'. It is designed
 * for a specific mapping type (maptype) and grid sizes.
 *
 * @param a   The input grid array of size (nx, ny, nz).
 * @param ah  The output contracted grid array of size ((nx) / 2, ny / 2, nz / 2).
 *
 * @note The function assumes that memory for 'a' and 'ah' has been allocated appropriately by the caller.
 * @note The function uses 'maptype' to determine the appropriate contraction method.
 *
 * @param maptype The mapping type used for contraction:
 * - If 'maptype' is 1, a 3D contraction is performed.
 * - If 'maptype' is not 1, a 2D contraction is performed.
 *
 * @details
 * - For maptype = 1 (3D contraction):
 *   - The input grid 'a' with dimensions (nx, ny, nz) is contracted to 'ah' with dimensions ((nx ) / 2, ny / 2, nz / 2).
 *
 * - For maptype != 1 (2D contraction):
 *   - The input grid 'a' with dimensions (nx, ny) is contracted to 'ah' with dimensions ((nx ) / 2, ny / 4).
 *
 * @par Implementation Details:
 * - The function uses nested loops to iterate over the input grid and perform the contraction according to the specified maptype.
 * - The resulting contracted grid 'ah' is filled with the contracted values.
 *
 */
void c3db::r2hr_contract(const double *a, double *ah) 
{
   std::memset(ah, 0, n2ft3d / 8 * sizeof(double));
   if (maptype == 1) 
   {
      int nxh = nx/2;
      int nyh = ny/2;
      int nqh = nq/2;
      for (auto q=0; q<nqh; ++q)
      for (auto j=0; j<nyh; ++j)
      for (auto i=0; i<nxh; ++i)
         ah[i + j*(nxh) + q*(nxh)*nyh] = a[i + j*(nx) + q*(nx)*ny];
   } 
   else 
   {
      int nxh = nx/2;
      int nq1h = nq1/4;
      for (auto q=0; q<nq1h; ++q)
      for (auto i=0; i<nxh; ++i)
         ah[i + q*(nxh)] = a[i + q*(nx)];
   }
}

/************************************************
 *           r_tranpose routines block          *
 ************************************************/
/******************************************
 *                                        *
 *      c3db::r_transpose_ijk_init        *
 *                                        *
 ******************************************/
void c3db::r_transpose_ijk_init() {
  if ((maptype != 1) && (!initialized_r_transpose)) {
    initialized_r_transpose = true;

    iq_to_ir1 = new (std::nothrow) int *[6]();
    iq_to_ir1[0] = new (std::nothrow) int[nx * nqr1]();
    iq_to_ir1[1] = new (std::nothrow) int[ny * nqr2]();
    iq_to_ir1[2] = new (std::nothrow) int[nz * nqr3]();

    iq_to_ir1[3] = new (std::nothrow) int[ny * nqr2]();
    iq_to_ir1[4] = new (std::nothrow) int[nx* nqr1]();
    iq_to_ir1[5] = new (std::nothrow) int[nz * nqr3]();

    iq_to_ir2 = new (std::nothrow) int *[6]();
    iq_to_ir2[0] = new (std::nothrow) int[ny * nqr2]();
    iq_to_ir2[1] = new (std::nothrow) int[nz * nqr3]();
    iq_to_ir2[2] = new (std::nothrow) int[ny * nqr2]();

    iq_to_ir2[3] = new (std::nothrow) int[nx * nq1]();
    iq_to_ir2[4] = new (std::nothrow) int[nz * nqr3]();
    iq_to_ir2[5] = new (std::nothrow) int[nx * nqr1]();

    ir1_start = new (std::nothrow) int *[6]();
    for (auto i = 0; i < 6; ++i)
      ir1_start[i] = new (std::nothrow) int[np + 1]();

    ir2_start = new (std::nothrow) int *[6]();
    for (auto i = 0; i < 6; ++i)
      ir2_start[i] = new (std::nothrow) int[np + 1]();

    /***********************************************************************/
    /* map1to2 mapping - done - tranpose operation #0  (j,k,i) <-- (i,j,k) */
    /***********************************************************************/
    int index1 = 0;
    int index2 = 0;
    for (auto it = 0; it < np; ++it) 
    {
       auto proc_to = (taskid + it) % np;
       auto proc_from = (taskid - it + np) % np;
       ir1_start[0][it] = index1;
       ir2_start[0][it] = index2;
       for (auto k = 0; k < nz; ++k)
         for (auto j = 0; j < ny; ++j)
           for (auto i = 0; i < nx; ++i) 
           {
              auto phere = cijkrtop2(i, j, k);
              auto pto = cijkrtop1(i, j, k);
             
              /* packing scheme */
              if ((phere == taskid) && (pto == proc_to)) 
              {
                 // int indx = ijktoindex2(i,j,k);
                 iq_to_ir1[0][cijkrtoindex2t(i, j, k)] = index1;
                 ++index1;
              }
              /* unpacking scheme */
              if ((pto == taskid) && (phere == proc_from)) 
              {
                 iq_to_ir2[0][cijkrtoindex1(i, j, k)] = index2;
                 ++index2;
              }
           }
    }
    ir1_start[0][np] = index1;
    ir2_start[0][np] = index2;

    /**********************************************************************/
    /* map2to3 mapping - done - tranpose operation #1 (k,i,j) <-- (j,k,i) */
    /**********************************************************************/
    index1 = 0;
    index2 = 0;
    for (auto it = 0; it < np; ++it) 
    {
       auto proc_to = (taskid + it) % np;
       auto proc_from = (taskid - it + np) % np;
       ir1_start[1][it] = index1;
       ir2_start[1][it] = index2;
       for (auto k = 0; k < nz; ++k)
         for (auto j = 0; j < ny; ++j)
           for (auto i = 0; i < nx; ++i) 
           {
              auto phere = cijkrtop1(i, j, k);
              auto pto = cijkrtop(i, j, k);
             
              /* packing scheme */
              if ((phere == taskid) && (pto == proc_to)) 
              {
                 iq_to_ir1[1][cijkrtoindex1(i, j, k)] = index1;
                 ++index1;
              }
              /* unpacking scheme */
              if ((pto == taskid) && (phere == proc_from)) 
              {
                 iq_to_ir2[1][cijkrtoindex(i, j, k)] = index2;
                 ++index2;
              }
           }
    }
    ir1_start[1][np] = index1;
    ir2_start[1][np] = index2;

    /**********************************************************************/
    /* map3to2 mapping - done - tranpose operation #2 (j,k,i) <-- (k,i,j) */
    /**********************************************************************/
    index1 = 0;
    index2 = 0;
    for (auto it = 0; it < np; ++it) 
    {
       auto proc_to = (taskid + it) % np;
       auto proc_from = (taskid - it + np) % np;
       ir1_start[2][it] = index1;
       ir2_start[2][it] = index2;
       for (auto k = 0; k < nz; ++k)
         for (auto j = 0; j < ny; ++j)
           for (auto i = 0; i < nx; ++i) 
           {
              auto phere = cijkrtop(i, j, k);
              auto pto = cijkrtop1(i, j, k);
             
              /* packing scheme */
              if ((phere == taskid) && (pto == proc_to)) 
              {
                 iq_to_ir1[2][cijkrtoindex(i, j, k)] = index1;
                 ++index1;
              }
              /* unpacking scheme */
              if ((pto == taskid) && (phere == proc_from)) 
              {
                 iq_to_ir2[2][cijkrtoindex1(i, j, k)] = index2;
                 ++index2;
              }
           }
    }
    ir1_start[2][np] = index1;
    ir2_start[2][np] = index2;


    /***********************************************************************/
    /* map2to1 mapping - done - tranpose operation #3  (i,j,k) <-- (j,k,i) */
    /***********************************************************************/
    index1 = 0;
    index2 = 0;
    for (auto it = 0; it < np; ++it) 
    {
       auto proc_to = (taskid + it) % np;
       auto proc_from = (taskid - it + np) % np;
       ir1_start[3][it] = index1;
       ir2_start[3][it] = index2;
       for (auto k = 0; k < nz; ++k)
         for (auto j = 0; j < ny; ++j)
           for (auto i = 0; i < nx; ++i) 
           {
              auto phere = cijkrtop1(i, j, k);
              auto pto = cijkrtop2(i, j, k);
             
              /* packing scheme */
              if ((phere == taskid) && (pto == proc_to)) 
              {
                 iq_to_ir1[3][cijkrtoindex1(i, j, k)] = index1;
                 ++index1;
              }
              /* unpacking scheme */
              if ((pto == taskid) && (phere == proc_from)) 
              {
                 iq_to_ir2[3][cijkrtoindex2t(i, j, k)] = index2;
                 ++index2;
              }
           }
    }
    ir1_start[3][np] = index1;
    ir2_start[3][np] = index2;

    /***********************************************************************/
    /* map1to3 mapping - done - tranpose operation #4  (k,i,j) <-- (i,j,k) */
    /***********************************************************************/
    index1 = 0;
    index2 = 0;
    for (auto it = 0; it < np; ++it) 
    {
       auto proc_to = (taskid + it) % np;
       auto proc_from = (taskid - it + np) % np;
       ir1_start[4][it] = index1;
       ir2_start[4][it] = index2;
       for (auto k = 0; k < nz; ++k)
         for (auto j = 0; j < ny; ++j)
           for (auto i = 0; i < nx; ++i) 
           {
              auto phere = cijkrtop2(i, j, k);
              auto pto = cijkrtop(i, j, k);
             
              /* packing scheme */
              if ((phere == taskid) && (pto == proc_to)) 
              {
                 iq_to_ir1[4][cijkrtoindex2t(i, j, k)] = index1;
                 ++index1;
              }
              /* unpacking scheme */
              if ((pto == taskid) && (phere == proc_from)) 
              {
                 iq_to_ir2[4][cijkrtoindex(i, j, k)] = index2;
                 ++index2;
              }
           }
    }
    ir1_start[4][np] = index1;
    ir2_start[4][np] = index2;

    /***********************************************************************/
    /* map3to1 mapping - done - tranpose operation #5  (i,j,k) <-- (k,i,j) */
    /***********************************************************************/
    index1 = 0;
    index2 = 0;
    for (auto it = 0; it < np; ++it) 
    {
       auto proc_to = (taskid + it) % np;
       auto proc_from = (taskid - it + np) % np;
       ir1_start[5][it] = index1;
       ir2_start[5][it] = index2;
       for (auto k = 0; k < nz; ++k)
         for (auto j = 0; j < ny; ++j)
           for (auto i = 0; i < nx; ++i) 
           {
              auto phere = cijkrtop(i, j, k);
              auto pto = cijkrtop2(i, j, k);
             
              /* packing scheme */
              if ((phere == taskid) && (pto == proc_to)) 
              {
                 iq_to_ir1[5][cijkrtoindex(i, j, k)] = index1;
                 ++index1;
              }
              /* unpacking scheme */
              if ((pto == taskid) && (phere == proc_from)) 
              {
                 iq_to_ir2[5][cijkrtoindex2t(i, j, k)] = index2;
                 ++index2;
              }
           }
    }
    ir1_start[5][np] = index1;
    ir2_start[5][np] = index2;
  }
}

/******************************************
 *                                        *
 *      c3db::r_transpose_ijk_end         *
 *                                        *
 ******************************************/
void c3db::r_transpose_ijk_end() {
  if (initialized_r_transpose) {
    initialized_r_transpose = false;
    {
      for (auto i = 0; i < 6; ++i) {
        delete[] iq_to_ir1[i];
        delete[] iq_to_ir2[i];
        delete[] ir1_start[i];
        delete[] ir2_start[i];
      }
      delete[] iq_to_ir1;
      delete[] iq_to_ir2;
      delete[] ir1_start;
      delete[] ir2_start;
    }
  }
}

/********************************
 *                              *
 *    c3db::r_transpose_jk      *
 *                              *
 ********************************/
/**
 * @brief Perform a transpose operation along the jk-plane on arrays 'a', 'tmp1', and 'tmp2'.
 *
 * This function is designed for parallel data communication and transposition in a scientific computing context.
 *
 * @param a The input array to be transposed.
 * @param tmp1 Temporary storage array for intermediate data.
 * @param tmp2 Temporary storage array for intermediate data.
 *
 * @note The function assumes that it has been initialized and is executed within a parallel computing environment.
 * @note The function performs parallel data communication and transposition along the jk-plane.
 *
 * @par Parallel Communication:
 * - The function is designed for parallel processing with 'np' processes. It uses 'parall' for communication and synchronization.
 * - Parallel communication operations include data packing, in-thread transposition, data receiving, and data sending.
 *
 * @par Memory Management:
 * - The function assumes that memory for 'a', 'tmp1', and 'tmp2' has been allocated appropriately by the caller.
 *
 * @author [Your Name]
 * @date [Date]
 */
void c3db::r_transpose_jk(double *a, double *tmp1, double *tmp2)
{
   int msglen;

   parall->astart(1,np);

   c_bindexcopy(nfft3d, iq_to_i1[0], a, tmp1);

   /* it = 0, transpose data on same thread */
   msglen = 2*(i2_start[0][1] - i2_start[0][0]);
   std::memcpy(tmp2+2*i2_start[0][0],tmp1+2*i1_start[0][0],msglen*sizeof(double));

   /* receive packed array data */
   for (auto it=1; it<np; ++it)
   {
      /* synchronous receive of tmp */
      auto proc_from = (taskid - it + np) % np;
      msglen = 2 * (i2_start[0][it + 1] - i2_start[0][it]);
      if (msglen > 0)
         parall->adreceive(1,1,proc_from, msglen, &tmp2[i2_start[0][it]]);
   }
   for (auto it=1; it<np; ++it)
   {
      auto proc_to = (taskid+it)%np;
      msglen = 2*(i1_start[0][it+1] - i1_start[0][it]);
      if (msglen > 0)
         parall->dsend(1,1,proc_to,msglen,&tmp1[i1_start[0][it]]);
   }
   parall->aend(1);
   c_aindexcopy(nfft3d,iq_to_i2[0],tmp2,a);
}

/********************************
 *                              *
 *    c3db::r_transpose_ijk     *
 *                              *
 ********************************/
/**
 * @brief Perform a transpose operation on arrays 'a', 'tmp1', and 'tmp2' based on the specified operation 'op'.
 *
 * This function is designed for parallel data communication and transposition in a scientific computing context.
 *
 * @param op The operation code that specifies the type of transpose operation to perform.
 *           Valid values are:
 *           - 0: Transpose along the x-axis
 *           - 1: Transpose along the y-axis
 *           - 2: Transpose along the z-axis
 *           - 3: Transpose along the x-axis (with additional operations)
 *           - 4: Transpose along the y-axis (with additional operations)
 *           - 5: Transpose along the z-axis (with additional operations)
 * @param a The input array to be transposed.
 * @param tmp1 Temporary storage array for intermediate data.
 * @param tmp2 Temporary storage array for intermediate data.
 *
 * @note The function assumes that it has been initialized using 'r_transpose_ijk_init' prior to invocation.
 * @note The function performs parallel data communication and transposition based on the specified 'op'.
 *
 * @par Parallel Communication:
 * - The function is designed for parallel processing with 'np' processes. It uses 'parall' for communication and synchronization.
 * - Parallel communication operations include data packing, in-thread transposition, data receiving, and data sending.
 *
 * @par Memory Management:
 * - The function assumes that memory for 'a', 'tmp1', and 'tmp2' has been allocated appropriately by the caller.
 *
 */
void c3db::r_transpose_ijk(const int op, double *a, double *tmp1, double *tmp2) 
{
   if (!initialized_r_transpose)
      this->r_transpose_ijk_init();
 
   int nnfft3d, it, proc_from, proc_to;
   int msglen;
 
   parall->astart(1, np);
 
   /* pack a array */
   if ((op == 0) || (op == 4)) nnfft3d = (nx)*nqr1;
   if ((op == 1) || (op == 3)) nnfft3d = (ny)*nqr2;
   if ((op == 2) || (op == 5)) nnfft3d = (nz)*nqr3;
   t_bindexcopy(nnfft3d, iq_to_ir1[op], a, tmp1);
 
   /* it = 0, transpose data on same thread */
   msglen = (ir2_start[op][1] - ir2_start[op][0]);
   // int one=1;
   // DCOPY_PWDFT(msglen,&(tmp1[i1_start[op][0]]),one,&(tmp2[i2_start[op][0]]),one);
   std::memcpy(tmp2 + ir2_start[op][0], tmp1 + ir1_start[op][0], msglen*sizeof(double));
 
   /* receive packed array data */
   for (it = 1; it < np; ++it) 
   {
      /* synchronous receive of tmp */
      proc_from = (taskid - it + np) % np;
      msglen = (ir2_start[op][it + 1] - ir2_start[op][it]);
      if (msglen > 0)
         parall->adreceive(1, 1, proc_from, msglen, &tmp2[ir2_start[op][it]]);
   }
   for (it = 1; it < np; ++it) {
     proc_to = (taskid + it) % np;
     msglen = (ir1_start[op][it + 1] - ir1_start[op][it]);
     if (msglen > 0)
       parall->dsend(1, 1, proc_to, msglen, &tmp1[ir1_start[op][it]]);
   }
 
   /* wait for completion of mp_send, also do a sync */
   parall->aend(1);
 
   /* unpack a array */
   if ((op == 3) || (op == 5)) nnfft3d = (nx) * nqr1;
   if ((op == 0) || (op == 2)) nnfft3d = (ny)*nqr2;
   if ((op == 1) || (op == 4)) nnfft3d = (nz)*nqr3;
   t_aindexcopy(nnfft3d, iq_to_ir2[op], tmp2, a);
}

/**************************************
 *                                    *
 *     c3db::rrrr_periodic_gradient   *
 *                                    *
 **************************************/
// computes the periodic gradient on a (n1,n2,n3) grid.
#define one_over_60 1.66666666666667e-2
void c3db::rrrr_periodic_gradient(const double *rho, double *drho1, double *drho2, double *drho3)
{

   if (maptype == 1)
   {
      double *a = this->r_alloc();
      double *tmp2 = c3db::c3db_tmp1;
      double *tmp3 = c3db::c3db_tmp2;
      std::memcpy(a,rho,n2ft3d*sizeof(double));
      this->r_transpose_jk(a,tmp2,tmp3);

      // drho1 gradient
      for (auto q=0; q<nq; ++q)
      {
         int im3,im2,im1,ip1,ip2,ip3;
         const double *f = rho + q * (nx)*nz;
         double *df1 = drho1 + q*(nx)*nz;
         for (auto i=0; i<nx; ++i)
         {
            im3 = i-3; if (im3<0) im3 += nx;
            im2 = i-2; if (im2<0) im2 += nx;
            im1 = i-1; if (im1<0) im1 += nx;
            ip1 = i+1; if (ip1>=nx) ip1 -= nx;
            ip2 = i+2; if (ip2>=nx) ip2 -= nx;
            ip3 = i+3; if (ip3>=nx) ip3 -= nx;

            for (auto j = 0; j < ny; ++j)
               df1[i+ j*nx] = one_over_60*(  -1.0*f[im3+j*nx] +  9.0*f[im2+j*nx]
                                           - 45.0*f[im1+j*nx] + 45.0*f[ip1+j*nx]
                                           -  9.0*f[ip2+j*nx] +  1.0*f[ip3+j*nx]);
         }
      }

      // drho2 gradient
      for (auto q=0; q<nq; ++q)
      {
         int jm3,jm2,jm1,jp1,jp2,jp3;
         const double *f = a + q*(nx)*ny;
         double *df2 = drho2 + q*(nx)*ny;
         for (auto j=0; j<ny; ++j)
         {
            jm3 = j-3; if (jm3<0) jm3 += ny;
            jm2 = j-2; if (jm2<0) jm2 += ny;
            jm1 = j-1; if (jm1<0) jm1 += ny;
            jp1 = j+1; if (jp1>=ny) jp1 -= ny;
            jp2 = j+2; if (jp2>=ny) jp2 -= ny;
            jp3 = j+3; if (jp3>=ny) jp3 -= ny;
            for (auto i = 0; i<nx; ++i)
               df2[i+j*nx] = one_over_60*( -1.0*f[i+jm3*nx] +  9.0*f[i+jm2*nx]
                                         - 45.0*f[i+jm1*nx] + 45.0*f[i+jp1*nx]
                                         -  9.0*f[i+jp2*nx] +  1.0*f[i+jp3*nx]);
         }
      }
      r_transpose_jk(drho2,tmp2,tmp3);

      // drho3 gradient
      for (auto q = 0; q < nq; ++q)
      {
         int km3,km2,km1,kp1,kp2,kp3;
         const double *f = rho + q*(nx)*nz;
         double *df3 = drho3 + q*(nx)*nz;
         for (auto k=0; k<nz; ++k)
         {
            km3 = k-3; if (km3<0) km3 += nz;
            km2 = k-2; if (km2<0) km2 += nz;
            km1 = k-1; if (km1<0) km1 += nz;
            kp1 = k+1; if (kp1>=nz) kp1 -= nz;
            kp2 = k+2; if (kp2>=nz) kp2 -= nz;
            kp3 = k+3; if (kp3>=nz) kp3 -= nz;
            for (auto i=0; i<nx; ++i)
               df3[i+k*nx] = one_over_60 * ( -1.0*f[i+km3*nx] +  9.0*f[i+km2*nx]
                                           - 45.0*f[i+km1*nx] + 45.0*f[i+kp1*nx]
                                           -  9.0*f[i+kp2*nx] +  1.0*f[i+kp3*nx]);
         }
      }
      this->r_dealloc(a);
   }
   else
   {
      double *a = this->r_alloc();
      double *tmp2 = new (std::nothrow) double[nrft3d]();
      double *tmp3 = new (std::nothrow) double[nrft3d]();

      // drho1 gradient
      for (auto q=0; q<nqr1; ++q)
      {
         int im3,im2,im1,ip1,ip2,ip3;
         const double *f = rho + q*(nx);
         double *df1 = drho1 + q*(nx);
         for (auto i=0; i<nx; ++i)
         {
            im3 = i-3; if (im3 < 0) im3 += nx;
            im2 = i-2; if (im2 < 0) im2 += nx;
            im1 = i-1; if (im1 < 0) im1 += nx;
            ip1 = i+1; if (ip1 >= nx) ip1 -= nx;
            ip2 = i+2; if (ip2 >= nx) ip2 -= nx;
            ip3 = i+3; if (ip3 >= nx) ip3 -= nx;
            df1[i] = one_over_60*( -1.0*f[im3] + 9.0*f[im2] - 45.0*f[im1]
                                 + 45.0*f[ip1] - 9.0*f[ip2] +  1.0*f[ip3]);
         }
      }

      // drho2 gradient
      std::memcpy(a,rho,nrft3d*sizeof(double));
      this->r_transpose_ijk(0,a,tmp2,tmp3); //#(i,j,k) --> (j,k,i)

      for (auto q=0; q<nqr2; ++q)
      {
         int jm3,jm2,jm1,jp1,jp2,jp3;
         const double *f = a + q*(ny);
         double *df2 = drho2 + q*(ny);
         for (auto j=0; j<ny; ++j)
         {
            jm3 = j-3; if (jm3<0) jm3 += ny;
            jm2 = j-2; if (jm2<0) jm2 += ny;
            jm1 = j-1; if (jm1<0) jm1 += ny;
            jp1 = j+1; if (jp1>=ny) jp1 -= ny;
            jp2 = j+2; if (jp2>=ny) jp2 -= ny;
            jp3 = j+3; if (jp3>=ny) jp3 -= ny;
            df2[j] = one_over_60*( -1.0*f[jm3] + 9.0*f[jm2] - 45.0*f[jm1]
                                 + 45.0*f[jp1] - 9.0*f[jp2] +  1.0*f[jp3]);
         }
      }
      this->r_transpose_ijk(3,drho2,tmp2,tmp3); //#(j,k,i) --> (i,j,k)

      // drho3 gradient
      std::memcpy(a,rho,nrft3d*sizeof(double));
      this->r_transpose_ijk(4,a,tmp2, tmp3);
      for (auto q=0; q<nqr3; ++q)
      {
         int km3,km2,km1,kp1,kp2,kp3;
         const double *f = a + q*(nz);
         double *df3 = drho3 + q*(nz);
         for (auto k=0; k<nz; ++k)
         {
            km3 = k-3; if (km3<0) km3 += nz;
            km2 = k-2; if (km2<0) km2 += nz;
            km1 = k-1; if (km1<0) km1 += nz;
            kp1 = k+1; if (kp1>=nz) kp1 -= nz;
            kp2 = k+2; if (kp2>=nz) kp2 -= nz;
            kp3 = k+3; if (kp3>=nz) kp3 -= nz;
            df3[k] = one_over_60*( -1.0*f[km3] + 9.0*f[km2] - 45.0*f[km1]
                                 + 45.0*f[kp1] - 9.0*f[kp2] +  1.0*f[kp3]);
         }
      }
      this->r_transpose_ijk(5,drho3,tmp2,tmp3);

      delete[] tmp3;
      delete[] tmp2;
      this->r_dealloc(a);
   }
}


/**************************************
 *                                    *
 *    c3db::rrrr_periodic_laplacian   *
 *                                    *
 **************************************/
// computes the periodic gradient on a (nx,ny,nz) grid.
#define one_over_180 5.55555555555556e-3
void c3db::rrrr_periodic_laplacian(const double *rho, double *grxx,
                                   double *gryy, double *grzz) {

  if (maptype == 1) 
  {
     double *a = this->r_alloc();
     double *tmp2 = c3db::c3db_tmp1;
     double *tmp3 = c3db::c3db_tmp2;
     std::memcpy(a, rho, n2ft3d * sizeof(double));
     this->r_transpose_jk(a, tmp2, tmp3);
    
     // xx gradient
     for (auto q = 0; q < nq; ++q) {
       int im3, im2, im1, ip1, ip2, ip3;
       const double *f = rho + q*nx*nz;
       double *ddf = grxx + q*nx*nz;
       for (auto i = 0; i < nx; ++i) {
         im3 = i - 3; if (im3 < 0) im3 += nx;
         im2 = i - 2; if (im2 < 0) im2 += nx;
         im1 = i - 1; if (im1 < 0) im1 += nx;
         ip1 = i + 1; if (ip1 >= nx) ip1 -= nx;
         ip2 = i + 2; if (ip2 >= nx) ip2 -= nx;
         ip3 = i + 3; if (ip3 >= nx) ip3 -= nx;
    
         for (auto j = 0; j < ny; ++j)
           ddf[i + j * nx] = one_over_180 * (2.0 * f[im3 + j * nx] - 27.0 * f[im2 + j * nx] +
                                           270.0 * f[im1 + j * nx] - 490.0 * f[i + j * nx] +
                                           270.0 * f[ip1 + j * nx] - 27.0 * f[ip2 + j * nx] + 2.0 * f[ip3 + j * nx]);
       }
     }
    
     // yy gradient
     for (auto q = 0; q < nq; ++q) {
       int jm3, jm2, jm1, jp1, jp2, jp3;
       const double *f = a + q*nx*ny;
       double *ddf = gryy + q*nx*ny;
       for (auto j = 0; j < ny; ++j) {
         jm3 = j - 3; if (jm3 < 0) jm3 += ny;
         jm2 = j - 2; if (jm2 < 0) jm2 += ny;
         jm1 = j - 1; if (jm1 < 0) jm1 += ny;
         jp1 = j + 1; if (jp1 >= ny) jp1 -= ny;
         jp2 = j + 2; if (jp2 >= ny) jp2 -= ny;
         jp3 = j + 3; if (jp3 >= ny) jp3 -= ny;
         for (auto i = 0; i < nx; ++i)
           ddf[i + j * nx] = one_over_180 * (2.0 * f[i + jm3 * nx] - 27.0 * f[i + jm2 * nx] +
                                           270.0 * f[i + jm1 * nx] - 490.0 * f[i + j * nx] +
                                           270.0 * f[i + jp1 * nx] - 27.0 * f[i + jp2 * nx] + 2.0 * f[i + jp3 * nx]);
       }
     }
     this->r_transpose_jk(gryy, tmp2, tmp3);
    
     // zz gradient
     for (auto q = 0; q < nq; ++q) {
       int km3, km2, km1, kp1, kp2, kp3;
       const double *f = rho + q*nx*nz;
       double *ddf = grzz + q*nx*nz;
       for (auto k = 0; k < nz; ++k) {
         km3 = k - 3; if (km3 < 0) km3 += nz;
         km2 = k - 2; if (km2 < 0) km2 += nz;
         km1 = k - 1; if (km1 < 0) km1 += nz;
         kp1 = k + 1; if (kp1 >= nz) kp1 -= nz;
         kp2 = k + 2; if (kp2 >= nz) kp2 -= nz;
         kp3 = k + 3; if (kp3 >= nz) kp3 -= nz;
         for (auto i = 0; i < nx; ++i)
           ddf[i + k * nx] = one_over_180 * (2.0 * f[i + km3 * nx] - 27.0 * f[i + km2 * nx] +
                                           270.0 * f[i + km1 * nx] - 490.0 * f[i + k * nx] +
                                           270.0 * f[i + kp1 * nx] - 27.0 * f[i + kp2 * nx] + 2.0 * f[i + kp3 * nx]);
       }
     }
    
     this->r_dealloc(a);
  } 
  else 
  {
    double *a = this->r_alloc();
    double *tmp2 = new (std::nothrow) double[nrft3d]();
    double *tmp3 = new (std::nothrow) double[nrft3d]();

    // xx gradient
    for (auto q = 0; q < nqr1; ++q) 
    {
       int im3, im2, im1, ip1, ip2, ip3;
       const double *f = rho + q*nx;
       double *ddf = grxx + q*nx;
       for (auto i = 0; i < nx; ++i) 
       {
          im3 = i - 3; if (im3 < 0) im3 += nx;
          im2 = i - 2; if (im2 < 0) im2 += nx;
          im1 = i - 1; if (im1 < 0) im1 += nx;
          ip1 = i + 1; if (ip1 >= nx) ip1 -= nx;
          ip2 = i + 2; if (ip2 >= nx) ip2 -= nx;
          ip3 = i + 3; if (ip3 >= nx) ip3 -= nx;
          ddf[i] = one_over_180 * (2.0 * f[im3] - 27.0 * f[im2] + 270.0 * f[im1] - 490.0 * f[i] +
                                 270.0 * f[ip1] - 27.0 * f[ip2] + 2.0 * f[ip3]);
       }
    }

    // yy gradient
    std::memcpy(a, rho, nrft3d * sizeof(double));
    this->r_transpose_ijk(0, a, tmp2, tmp3); //#(i,j,k) --> (j,k,i)
    for (auto q = 0; q < nqr2; ++q) {
      int im3, im2, im1, ip1, ip2, ip3;
      const double *f = a + q * (ny);
      double *ddf = gryy + q * (ny);
      for (auto i = 0; i < ny; ++i) {
        im3 = i-3; if (im3 < 0) im3 += ny;
        im2 = i-2; if (im2 < 0) im2 += ny;
        im1 = i-1; if (im1 < 0) im1 += ny;
        ip1 = i+1; if (ip1 >= ny) ip1 -= ny;
        ip2 = i+2; if (ip2 >= ny) ip2 -= ny;
        ip3 = i+3; if (ip3 >= ny) ip3 -= ny;
        ddf[i] = one_over_180 *
                 (2.0 * f[im3] - 27.0 * f[im2] + 270.0 * f[im1] - 490.0 * f[i] +
                  270.0 * f[ip1] - 27.0 * f[ip2] + 2.0 * f[ip3]);
      }
    }
    this->r_transpose_ijk(3, gryy, tmp2, tmp3); //#(j,k,i) --> (i,j,k)

    // zz gradient
    std::memcpy(a, rho, nrft3d * sizeof(double));
    this->r_transpose_ijk(4, a, tmp2, tmp3);
    for (auto q = 0; q < nqr3; ++q) {
      int im3, im2, im1, ip1, ip2, ip3;
      const double *f = a + q * (nz);
      double *ddf = grzz + q * (nz);
      for (auto i = 0; i < nz; ++i) {
        im3 = i-3; if (im3 < 0) im3 += nz;
        im2 = i-2; if (im2 < 0) im2 += nz;
        im1 = i-1; if (im1 < 0) im1 += nz;
        ip1 = i+1; if (ip1 >= nz) ip1 -= nz;
        ip2 = i+2; if (ip2 >= nz) ip2 -= nz;
        ip3 = i+3; if (ip3 >= nz) ip3 -= nz;
        ddf[i] = one_over_180 *
                 (2.0 * f[im3] - 27.0 * f[im2] + 270.0 * f[im1] - 490.0 * f[i] +
                  270.0 * f[ip1] - 27.0 * f[ip2] + 2.0 * f[ip3]);
      }
    }
    this->r_transpose_ijk(5, grzz, tmp2, tmp3);

    delete[] tmp3;
    delete[] tmp2;
    this->r_dealloc(a);
  }
}

/************************************************
 *           r_tranpose routines block          *
 ************************************************/


/************************************************
 *       Gaussian filters routines block        *
 ************************************************/

 #define C_PERIODIC_GAUSSIAN_FILTER(nx, nfilter, coeff, a, b) \
    do { \
        for (auto i=0; i<(nx); ++i) { \
            double sum = 0.0; \
            for (auto r=0; r<(nfilter); ++r) { \
                auto ii = i+r-((nfilter)-1)/2; \
                if (ii < 0)   ii += (nx); \
                if (ii >= (nx)) ii -= (nx); \
                sum += (coeff)[r]*(a)[ii]; \
            } \
            (b)[i] = sum; \
        } \
    } while (0)

#define C_GAUSSIAN_FILTER(nx, nfilter, coeff, a, b) \
    do { \
        for (auto i=0; i<(nx); ++i) { \
            double sum = 0.0; \
            for (auto r=0; r<(nfilter); ++r) { \
                auto ii = i+r-((nfilter)-1)/2; \
                if ((ii >= 0) && (ii < (nx))) \
                   sum += (coeff)[r]*(a)[ii]; \
            } \
            (b)[i] = sum; \
        } \
    } while (0)


#define C_SMOOTH_PERIODIC_2D(nx, ny, nfilter, coeff, a, b) \
    do { \
        for (auto j=0; j<(ny); ++j) \
        for (auto i=0; i<(nx); ++i) \
        { \
            double sum = 0.0; \
            for (auto r=0; r<(nfilter); ++r) \
            { \
                auto ii = i+r-((nfilter)-1)/2; \
                if (ii < 0)     ii += (nx); \
                if (ii >= (nx)) ii -= (nx); \
                sum += (coeff)[r]*(a)[ii+j*(nx)]; \
            } \
            (b)[i+j*(nx)] = sum; \
        } \
    } while(0)

#define C_SMOOTH_PERIODIC_2D_TRANS(nx, ny, nfilter, coeff, a, b) \
    do { \
        for (auto i=0; i<(nx); ++i) \
        for (auto j=0; j<(ny); ++j) \
        { \
            double sum = 0.0; \
            for (auto r=0; r<(nfilter); ++r) \
            { \
                auto jj = j+r-((nfilter)-1)/2; \
                if (jj < 0)     jj += (ny); \
                if (jj >= (ny)) jj -= (ny); \
                sum += (coeff)[r]*(a)[i+jj*(nx)]; \
            } \
            (b)[i+j*(nx)] = sum; \
        } \
    } while(0)


/**********************************************
 *                                            *
 *     c3db::rr_periodic_gaussian_filter      *
 *                                            *
 **********************************************/
// computes the 3D Gaussian filter function
/**
 * @brief Compute a 3D Gaussian filter for periodic data.
 *
 * This function computes a 3D Gaussian filter for periodic data with the specified sigma value.
 *
 * @param sigma The standard deviation of the Gaussian filter.
 * @param a The input real-space data.
 * @param b The output filtered data.
 */
void c3db::rr_periodic_gaussian_filter(const double sigma, const double *a, double *b)
{

   // Compute the filter size based on the sigma value
   int nfilter = 2*(int)(2.5*sigma) + 1;

   // Compute the filter coefficients
   std::vector<double> coeff(nfilter);
   double sum = 0.0;
   for (auto i=0; i<nfilter; ++i)
   {
      coeff[i] = std::exp(-0.5*std::pow((i-(nfilter-1)/2)/sigma,2));
      sum += coeff[i];
   }
   for (auto i=0; i<nfilter; ++i)
      coeff[i] /= sum;

   if (maptype==1)
   {
      double *c    = this->r_alloc();
      double *tmp2 = this->r_alloc();
      double *tmp3 = this->r_alloc();

      // Apply the filter along the x axis
      for (auto q=0; q<nq; ++q)
      {
         const double *f  = a + q*(nx)*nz;
         double       *f1 = b + q*(nx)*nz;
         C_SMOOTH_PERIODIC_2D(nx,ny,nfilter,coeff,f,f1);
      }

      // Apply the filter along the y axis
      for (auto q=0; q<nq; ++q)
      {
         const double *f  = b + q*(nx)*nz;
         double       *f1 = c + q*(nx)*nz;
         C_SMOOTH_PERIODIC_2D_TRANS(nx,ny,nfilter,coeff,f,f1);
      }

      // transpose y and z
      this->c_transpose_jk(c,tmp2,tmp3);

      // Apply the filter along the z axis
      for (auto q=0; q<nq; ++q)
      {
         const double *f  = c + q*(nx)*ny;
         double       *f1 = b + q*(nx)*ny;
         C_SMOOTH_PERIODIC_2D_TRANS(nx,nz,nfilter,coeff,f,f1);
      }

      // transpose z and y
      this->c_transpose_jk(b,tmp2,tmp3);

      this->r_dealloc(tmp3);
      this->r_dealloc(tmp2);
      this->r_dealloc(c);

   }
   else
   {
      double *c = this->r_alloc();
      double *tmp2 = this->r_alloc();
      double *tmp3 = this->r_alloc();
      //double *tmp2 = new (std::nothrow) double[nrft3d]();
      ////double *tmp3 = new (std::nothrow) double[nrft3d]();

      // Apply the filter along the x axis
      for (auto q=0; q<nqr1; ++q)
      {
         const double *f = a + q*(nx);
         double      *f1 = b + q*(nx);
         C_PERIODIC_GAUSSIAN_FILTER(nx,nfilter,coeff,f,f1);
      }

      // tranpose operation #0  (j,k,i) <-- (i,j,k)
      // rotate b(ny,nz,nx) <-- b(nx,ny,nz)
      this->r_transpose_ijk(0,b,tmp2,tmp3);

      // Apply the filter along the y axis
      for (auto q=0; q<nqr2; ++q)
      {
         const double *f = b + q*(ny);
         double      *f1 = c + q*(ny);
         C_PERIODIC_GAUSSIAN_FILTER(ny,nfilter,coeff,f,f1);
      }

      // tranpose operation #1 (k,i,j) <-- (j,k,i)
      // rotate c(nz,nx,ny) <-- c(ny,nz,nx)
      this->r_transpose_ijk(1,c,tmp2,tmp3);

      // Apply the filter along the z axis
      for (auto q=0; q<nqr3; ++q)
      {
         const double *f = c + q*(nz);
         double      *f1 = b + q*(nz);
         C_PERIODIC_GAUSSIAN_FILTER(nz,nfilter,coeff,f,f1);
      }

      // tranpose operation #5 (i,j,k) <-- (k,i,j)
      // rotate b(nx,ny,nz) <-- b(nz,nx,ny)
      this->r_transpose_ijk(5,b,tmp2,tmp3);

      //delete [] tmp3;
      //delete [] tmp2;
      this->r_dealloc(tmp3);
      this->r_dealloc(tmp2);
      this->r_dealloc(c);
   }
}

/************************************************
 *       End Gaussian filters routines block    *
 ************************************************/



} // namespace pwdft

/* d3db.C
   Author - Eric Bylaska

        this class is used for defining 3d parallel maps
*/

#include "compressed_io.hpp"
#include "fft.h"
#include "util.hpp"
#include <cmath>
#include <vector>

#include "blas.h"

//#include "gdevice.hpp"
#include "nwpw_timing.hpp"

#include "d3db.hpp"

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

d3db::d3db(Parallel *inparall, const int inmaptype, const int nx, const int ny, const int nz)
    : Mapping3(inmaptype, inparall->np_i(), inparall->taskid_i(), nx, ny, nz) 
{
   int index1, index2, proc_to, proc_from;
   int nyh, nzh;
   int phere, pto, pfrom;
 
   parall = inparall;

 
   if (maptype==1) 
   {
     iq_to_i1 = new (std::nothrow) int *[1];
     iq_to_i1[0] = new (std::nothrow) int[(nx/2+1) * ny * nq]();
     iq_to_i2 = new (std::nothrow) int *[1];
     iq_to_i2[0] = new (std::nothrow) int[(nx/2+1) * ny * nq]();
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
         for (auto j = 0; j < ny; ++j) {
           /* packing scheme */
           phere = ijktop(0, 0, k);
           pto = ijktop(0, 0, j);
           if ((phere == taskid) && (pto == proc_to))
             for (auto i = 0; i < (nx/2+1); ++i) {
               iq_to_i1[0][ijktoindex(i, j, k)] = index1;
               ++index1;
             }
 
           /* unpacking scheme */
           phere = ijktop(0, 0, j);
           pfrom = ijktop(0, 0, k);
           if ((phere == taskid) && (pfrom == proc_from))
             for (auto i = 0; i < (nx/2+1); ++i) {
               iq_to_i2[0][ijktoindex(i, k, j)] = index2;
               ++index2;
             }
         }
     }
     i1_start[0][np] = index1;
     i2_start[0][np] = index2;
 
     /* allocate ptranspose indexes */
     for (auto nb = 0; nb < 2; ++nb) {
       p_iq_to_i1[nb] = new (std::nothrow) int *[1]();
       p_iq_to_i1[nb][0] = new (std::nothrow) int[(nx / 2 + 1) * ny * nq]();
       p_iq_to_i2[nb] = new (std::nothrow) int *[1]();
       p_iq_to_i2[nb][0] = new (std::nothrow) int[(nx / 2 + 1) * ny * nq]();
       p_iz_to_i2[nb] = new (std::nothrow) int *[1]();
       p_iz_to_i2[nb][0] = new (std::nothrow) int[(nx / 2 + 1) * ny * nq]();
       p_i1_start[nb] = new (std::nothrow) int *[1]();
       p_i1_start[nb][0] = new (std::nothrow) int[nz + 1]();
       p_i2_start[nb] = new (std::nothrow) int *[1]();
       p_i2_start[nb][0] = new (std::nothrow) int[nz + 1]();
 
       p_jq_to_i1[nb] = new (std::nothrow) int *[1]();
       p_jq_to_i1[nb][0] = new (std::nothrow) int[(nx / 2 + 1) * ny * nq]();
       p_jq_to_i2[nb] = new (std::nothrow) int *[1]();
       p_jq_to_i2[nb][0] = new (std::nothrow) int[(nx / 2 + 1) * ny * nq]();
       p_jz_to_i2[nb] = new (std::nothrow) int *[1]();
       p_jz_to_i2[nb][0] = new (std::nothrow) int[(nx / 2 + 1) * ny * nq]();
       p_j1_start[nb] = new (std::nothrow) int *[1]();
       p_j1_start[nb][0] = new (std::nothrow) int[nz + 1]();
       p_j2_start[nb] = new (std::nothrow) int *[1]();
       p_j2_start[nb][0] = new (std::nothrow) int[nz + 1]();
     }
 
   } else {
     iq_to_i1 = new (std::nothrow) int *[6]();
     iq_to_i1[0] = new (std::nothrow) int[(nx / 2 + 1) * nq1]();
     iq_to_i1[1] = new (std::nothrow) int[ny * nq2]();
     iq_to_i1[2] = new (std::nothrow) int[nz * nq3]();
 
     iq_to_i1[3] = new (std::nothrow) int[ny * nq2]();
     iq_to_i1[4] = new (std::nothrow) int[(nx / 2 + 1) * nq1]();
     iq_to_i1[5] = new (std::nothrow) int[nz * nq3]();
 
     iq_to_i2 = new (std::nothrow) int *[6]();
     iq_to_i2[0] = new (std::nothrow) int[ny * nq2]();
     iq_to_i2[1] = new (std::nothrow) int[nz * nq3]();
     iq_to_i2[2] = new (std::nothrow) int[ny * nq2]();
 
     iq_to_i2[3] = new (std::nothrow) int[(nx / 2 + 1) * nq1]();
     iq_to_i2[4] = new (std::nothrow) int[nz * nq3]();
     iq_to_i2[5] = new (std::nothrow) int[(nx / 2 + 1) * nq1]();
 
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
     for (auto it = 0; it < np; ++it) {
       proc_to = (taskid + it) % np;
       proc_from = (taskid - it + np) % np;
       i1_start[0][it] = index1;
       i2_start[0][it] = index2;
       for (auto k = 0; k < nz; ++k)
         for (auto j = 0; j < ny; ++j)
           for (auto i = 0; i < (nx / 2 + 1); ++i) {
             phere = ijktop2(i, j, k);
             pto = ijktop1(i, j, k);
 
             /* packing scheme */
             if ((phere == taskid) && (pto == proc_to)) {
               iq_to_i1[0][ijktoindex2t(i, j, k)] = index1;
               ++index1;
             }
             /* unpacking scheme */
             if ((pto == taskid) && (phere == proc_from)) {
               iq_to_i2[0][ijktoindex1(i, j, k)] = index2;
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
     for (auto it = 0; it < np; ++it) {
       proc_to = (taskid + it) % np;
       proc_from = (taskid - it + np) % np;
       i1_start[1][it] = index1;
       i2_start[1][it] = index2;
       for (auto k = 0; k < nz; ++k)
         for (auto j = 0; j < ny; ++j)
           for (auto i = 0; i < (nx / 2 + 1); ++i) {
             phere = ijktop1(i, j, k);
             pto = ijktop(i, j, k);
 
             /* packing scheme */
             if ((phere == taskid) && (pto == proc_to)) {
               iq_to_i1[1][ijktoindex1(i, j, k)] = index1;
               ++index1;
             }
             /* unpacking scheme */
             if ((pto == taskid) && (phere == proc_from)) {
               iq_to_i2[1][ijktoindex(i, j, k)] = index2;
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
     for (auto it = 0; it < np; ++it) {
       proc_to = (taskid + it) % np;
       proc_from = (taskid - it + np) % np;
       i1_start[2][it] = index1;
       i2_start[2][it] = index2;
       for (auto k = 0; k < nz; ++k)
         for (auto j = 0; j < ny; ++j)
           for (auto i = 0; i < (nx / 2 + 1); ++i) {
             phere = ijktop(i, j, k);
             pto = ijktop1(i, j, k);
 
             /* packing scheme */
             if ((phere == taskid) && (pto == proc_to)) {
               iq_to_i1[2][ijktoindex(i, j, k)] = index1;
               ++index1;
             }
             /* unpacking scheme */
             if ((pto == taskid) && (phere == proc_from)) {
               iq_to_i2[2][ijktoindex1(i, j, k)] = index2;
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
     for (auto it = 0; it < np; ++it) {
       proc_to = (taskid + it) % np;
       proc_from = (taskid - it + np) % np;
       i1_start[3][it] = index1;
       i2_start[3][it] = index2;
       for (auto k = 0; k < nz; ++k)
         for (auto j = 0; j < ny; ++j)
           for (auto i = 0; i < (nx / 2 + 1); ++i) {
             phere = ijktop1(i, j, k);
             pto = ijktop2(i, j, k);
 
             /* packing scheme */
             if ((phere == taskid) && (pto == proc_to)) {
               iq_to_i1[3][ijktoindex1(i, j, k)] = index1;
               ++index1;
             }
             /* unpacking scheme */
             if ((pto == taskid) && (phere == proc_from)) {
               iq_to_i2[3][ijktoindex2t(i, j, k)] = index2;
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
     for (auto it = 0; it < np; ++it) {
       proc_to = (taskid + it) % np;
       proc_from = (taskid - it + np) % np;
       i1_start[4][it] = index1;
       i2_start[4][it] = index2;
       for (auto k = 0; k < nz; ++k)
         for (auto j = 0; j < ny; ++j)
           for (auto i = 0; i < (nx/2 + 1); ++i) {
             phere = ijktop2(i, j, k);
             pto = ijktop(i, j, k);
 
             /* packing scheme */
             if ((phere == taskid) && (pto == proc_to)) {
               iq_to_i1[4][ijktoindex2t(i, j, k)] = index1;
               ++index1;
             }
             /* unpacking scheme */
             if ((pto == taskid) && (phere == proc_from)) {
               iq_to_i2[4][ijktoindex(i, j, k)] = index2;
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
     for (auto it = 0; it < np; ++it) {
       proc_to = (taskid + it) % np;
       proc_from = (taskid - it + np) % np;
       i1_start[5][it] = index1;
       i2_start[5][it] = index2;
       for (auto k = 0; k < nz; ++k)
         for (auto j = 0; j < ny; ++j)
           for (auto i = 0; i < (nx / 2 + 1); ++i) {
             phere = ijktop(i, j, k);
             pto = ijktop2(i, j, k);
 
             /* packing scheme */
             if ((phere == taskid) && (pto == proc_to)) {
               iq_to_i1[5][ijktoindex(i, j, k)] = index1;
               ++index1;
             }
             /* unpacking scheme */
             if ((pto == taskid) && (phere == proc_from)) {
               iq_to_i2[5][ijktoindex2t(i, j, k)] = index2;
               ++index2;
             }
           }
     }
     i1_start[5][np] = index1;
     i2_start[5][np] = index2;
 
     /* allocate ptranspose indexes */
     for (auto nb = 0; nb < 2; ++nb) {
       p_iq_to_i1[nb] = new (std::nothrow) int *[6];
       p_iq_to_i1[nb][0] = new (std::nothrow) int[(nx / 2 + 1) * nq1]();
       p_iq_to_i1[nb][1] = new (std::nothrow) int[ny * nq2]();
       p_iq_to_i1[nb][2] = new (std::nothrow) int[nz * nq3]();
       p_iq_to_i1[nb][3] = new (std::nothrow) int[ny * nq2]();
       p_iq_to_i1[nb][4] = new (std::nothrow) int[(nx / 2 + 1) * nq1]();
       p_iq_to_i1[nb][5] = new (std::nothrow) int[nz * nq3]();
 
       p_iq_to_i2[nb] = new (std::nothrow) int *[6];
       p_iq_to_i2[nb][0] = new (std::nothrow) int[ny * nq2]();
       p_iq_to_i2[nb][1] = new (std::nothrow) int[nz * nq3]();
       p_iq_to_i2[nb][2] = new (std::nothrow) int[ny * nq2]();
       p_iq_to_i2[nb][3] = new (std::nothrow) int[(nx / 2 + 1) * nq1]();
       p_iq_to_i2[nb][4] = new (std::nothrow) int[nz * nq3]();
       p_iq_to_i2[nb][5] = new (std::nothrow) int[(nx / 2 + 1) * nq1]();
 
       p_iz_to_i2[nb] = new (std::nothrow) int *[6];
       p_iz_to_i2[nb][0] = new (std::nothrow) int[ny * nq2]();
       p_iz_to_i2[nb][1] = new (std::nothrow) int[nz * nq3]();
       p_iz_to_i2[nb][2] = new (std::nothrow) int[ny * nq2]();
       p_iz_to_i2[nb][3] = new (std::nothrow) int[(nx / 2 + 1) * nq1]();
       p_iz_to_i2[nb][4] = new (std::nothrow) int[nz * nq3]();
       p_iz_to_i2[nb][5] = new (std::nothrow) int[(nx / 2 + 1) * nq1]();
 
       p_i1_start[nb] = new (std::nothrow) int *[6];
       for (auto i = 0; i < 6; ++i)
         p_i1_start[nb][i] = new (std::nothrow) int[np + 1]();
 
       p_i2_start[nb] = new (std::nothrow) int *[6];
       for (auto i = 0; i < 6; ++i)
         p_i2_start[nb][i] = new (std::nothrow) int[np + 1]();
     }
   }
 
   /* setup timereverse indexes */
   zplane_size = timereverse_size() + 1;
   t_iq_to_i1 = new (std::nothrow) int[zplane_size]();
   t_iq_to_i2 = new (std::nothrow) int[zplane_size]();
   t_i1_start = new (std::nothrow) int[np + 1]();
   t_i2_start = new (std::nothrow) int[np + 1]();
   nyh = ny / 2;
   nzh = nz / 2;
   index1 = 0;
   index2 = 0;
   for (auto it = 0; it < np; ++it) {
     proc_to = (taskid + it) % np;
     proc_from = (taskid - it + np) % np;
     t_i1_start[it] = index1;
     t_i2_start[it] = index2;
 
     /* k=(0,0,k3) */
     for (auto k = 1; k < nzh; ++k) {
       auto k1 = k;
       auto k2 = -k;
       if (k1 < 0)
         k1 += nz;
       if (k2 < 0)
         k2 += nz;
       phere = ijktop(0, 0, k1);
       pto = ijktop(0, 0, k2);
 
       /* packing scheme */
       if ((phere == taskid) && (pto == proc_to)) {
         t_iq_to_i1[index1] = ijktoindex(0, 0, k1);
         ++index1;
       }
       /* unpacking scheme */
       if ((pto == taskid) && (phere == proc_from)) {
         t_iq_to_i2[index2] = ijktoindex(0, 0, k2);
         ++index2;
       }
     }
 
     /* k=(0,k2,k3) */
     for (auto k = (-nzh + 1); k < nzh; ++k)
       for (auto j = 1; j < nyh; ++j) {
         auto j1 = j;
         auto k1 = k;
         if (j1 < 0)
           j1 += ny;
         if (k1 < 0)
           k1 += nz;
         auto j2 = -j;
         auto k2 = -k;
         if (j2 < 0)
           j2 += ny;
         if (k2 < 0)
           k2 += nz;
         phere = ijktop(0, j1, k1);
         pto = ijktop(0, j2, k2);
 
         /* packing scheme */
         if ((phere == taskid) && (pto == proc_to)) {
           t_iq_to_i1[index1] = ijktoindex(0, j1, k1);
           ++index1;
         }
         /* unpacking scheme */
         if ((pto == taskid) && (phere == proc_from)) {
           t_iq_to_i2[index2] = ijktoindex(0, j2, k2);
           ++index2;
         }
       }
   }
   t_i1_start[np] = index1;
   t_i2_start[np] = index2;
 
   /* setup ffts */
   tmpx = new (std::nothrow) double[2*(2*nx+15)]();
   tmpy = new (std::nothrow) double[2*(2*ny+15)]();
   tmpz = new (std::nothrow) double[2*(2*nz+15)]();
   drffti_(&nx,tmpx);
   dcffti_(&ny,tmpy);
   dcffti_(&nz,tmpz);

#if (defined NWPW_SYCL) || (defined NWPW_CUDA) || (defined NWPW_HIP)
   if (maptype==1) 
      fft_tag = mygdevice.batch_fft_init(nx,ny,nz,ny*nq,(nx/2+1)*nq,(nx/2+1)*nq);
   else
     fft_tag = mygdevice.batch_fft_init(nx,ny,nz,nq1,nq2,nq3);
#endif
}

/********************************
 *                              *
 *          Destructors         *
 *                              *
 ********************************/
d3db::~d3db() {
   int i, nb;

#if (defined NWPW_SYCL) || (defined NWPW_CUDA) || (defined NWPW_HIP)
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
      for (nb = 0; nb < 2; ++nb) {
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
   } else {
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
     
      for (nb = 0; nb < 2; ++nb) {
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
      }
   }
 
   delete[] t_iq_to_i1;
   delete[] t_iq_to_i2;
   delete[] t_i1_start;
   delete[] t_i2_start;
 
   delete[] tmpx;
   delete[] tmpy;
   delete[] tmpz;
 
   //#endif
}

/******************************************
 *                                        *
 *      d3db::c_ptranspose_jk_init        *
 *                                        *
 ******************************************/
void d3db::c_ptranspose_jk_init(const int nb, bool *zero_arow3) {
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
      for (auto j = 0; j < ny; ++j) {
        /* packing scheme */
        phere = ijktop(0, 0, k);
        pto = ijktop(0, 0, j);
        if ((phere == taskid) && (pto == proc_to))
          for (auto i = 0; i < (nx / 2 + 1); ++i) {
            // iszero_ii = (zero_arow3[i + (nx/2+1)*(k-1)]==1);
            // iszero_jj = (zero_arow3[i + (nx/2+1)*(j-1)]==1);
            iszero_ii = (zero_arow3[i + (nx / 2 + 1) * k]);
            iszero_jj = (zero_arow3[i + (nx / 2 + 1) * j]);
            if (!iszero_ii) {
              p_iq_to_i1[nb][0][index1] = ijktoindex(i, j, k);
              ++index1;
            }
            if (!iszero_jj) {
              p_jq_to_i1[nb][0][jndex1] = ijktoindex(i, j, k);
              ++jndex1;
            }
          }

        /* unpacking scheme */
        phere = ijktop(0, 0, j);
        pfrom = ijktop(0, 0, k);
        if ((phere == taskid) && (pfrom == proc_from))
          for (auto i = 0; i < (nx / 2 + 1); ++i) {
            // iszero_ii = (zero_arow3[i + (nx/2+1)*(k-1)]==1);
            // iszero_jj = (zero_arow3[i + (nx/2+1)*(j-1)]==1);
            iszero_ii = (zero_arow3[i + (nx / 2 + 1) * k]);
            iszero_jj = (zero_arow3[i + (nx / 2 + 1) * j]);
            if (!iszero_ii) {
              p_iq_to_i2[nb][0][index2] = ijktoindex(i, k, j);
              ++index2;
            } else {
              p_iz_to_i2[nb][0][index3] = ijktoindex(i, k, j);
              ++index3;
            }
            if (!iszero_jj) {
              p_jq_to_i2[nb][0][jndex2] = ijktoindex(i, k, j);
              ++jndex2;
            } else {
              p_jz_to_i2[nb][0][jndex3] = ijktoindex(i, k, j);
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
 *      d3db::c_ptranspose_ijk_init       *
 *                                        *
 ******************************************/
void d3db::c_ptranspose_ijk_init(const int nb, bool *zero_arow2,
                                 bool *zero_arow3) {
  int index1, index2, index3, proc_to, proc_from, phere, pto;
  bool iszero;

  /**************************************************/
  /* map1to2 mapping - done - tranpose operation #0 */
  /*  (ny,nz,nx/2+1)  <-- (nx/2+1,ny,nz)            */
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
        for (auto i = 0; i < (nx / 2 + 1); ++i) {
          iszero = (zero_arow2[i + k * (nx / 2 + 1)]);

          // phere = int_mb(p_map1(1,id)+(j-1)+(k-1)*ny(id))
          // pto   = int_mb(p_map2(1,id)+(k-1)+(i-1)*nz(id))

          phere = ijktop2(i, j, k);
          pto = ijktop1(i, j, k);

          /* packing scheme */
          if ((phere == taskid) && (pto == proc_to)) {
            if (!iszero) {
              p_iq_to_i1[nb][0][index1] = ijktoindex2t(i, j, k);
              ++index1;
            }
          }
          /* unpacking scheme */
          if ((pto == taskid) && (phere == proc_from)) {
            if (!iszero) {
              p_iq_to_i2[nb][0][index2] = ijktoindex1(i, j, k);
              ++index2;
            } else {
              p_iz_to_i2[nb][0][index3] = ijktoindex1(i, j, k);
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
  /*     (nz,nx/2+1,ny)  <-- (ny,nz,nx/2+1)         */
  /*     use zero_arow3                             */
  /**************************************************/
  index1 = 0;
  index2 = 0;
  index3 = 0;
  for (auto it = 0; it < np; ++it) {
    proc_to = (taskid + it) % np;
    proc_from = (taskid - it + np) % np;
    p_i1_start[nb][1][it] = index1;
    p_i2_start[nb][1][it] = index2;
    for (auto k = 0; k < nz; ++k)
      for (auto j = 0; j < ny; ++j)
        for (auto i = 0; i < (nx / 2 + 1); ++i) {
          iszero = (zero_arow3[i + j * (nx / 2 + 1)]);

          // phere = int_mb(p_map2(1,id)+(k-1)+(i-1)*nz(id))
          // pto   = int_mb(p_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))

          phere = ijktop1(i, j, k);
          pto = ijktop(i, j, k);

          /* packing scheme */
          if ((phere == taskid) && (pto == proc_to)) {
            if (!iszero) {
              p_iq_to_i1[nb][1][index1] = ijktoindex1(i, j, k);
              ++index1;
            }
          }
          /* unpacking scheme */
          if ((pto == taskid) && (phere == proc_from)) {
            if (!iszero) {
              p_iq_to_i2[nb][1][index2] = ijktoindex(i, j, k);
              ++index2;
            } else {
              p_iz_to_i2[nb][1][index3] = ijktoindex(i, j, k);
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
  /*     (ny,nz,nx/2+1)  <-- (nz,nx/2+1,ny)         */
  /*     use zero_arow3                             */
  /**************************************************/
  index1 = 0;
  index2 = 0;
  index3 = 0;
  for (auto it = 0; it < np; ++it) {
    proc_to = (taskid + it) % np;
    proc_from = (taskid - it + np) % np;
    p_i1_start[nb][2][it] = index1;
    p_i2_start[nb][2][it] = index2;
    for (auto k = 0; k < nz; ++k)
      for (auto j = 0; j < ny; ++j)
        for (auto i = 0; i < (nx / 2 + 1); ++i) {
          iszero = (zero_arow3[i + j * (nx / 2 + 1)]);

          // phere = int_mb(p_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))
          // pto   = int_mb(p_map2(1,id)+(k-1)+(i-1)*nz(id))

          phere = ijktop(i, j, k);
          pto = ijktop1(i, j, k);

          /* packing scheme */
          if ((phere == taskid) && (pto == proc_to)) {
            if (!iszero) {
              p_iq_to_i1[nb][2][index1] = ijktoindex(i, j, k);
              ++index1;
            }
          }
          /* unpacking scheme */
          if ((pto == taskid) && (phere == proc_from)) {
            if (!iszero) {
              p_iq_to_i2[nb][2][index2] = ijktoindex1(i, j, k);
              ++index2;
            } else {
              p_iz_to_i2[nb][2][index3] = ijktoindex1(i, j, k);
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
  /*     (nx/2+1,ny,nz)  <-- (ny,nz,nx/2+1)         */
  /*     use zero_arow2                             */
  /**************************************************/
  index1 = 0;
  index2 = 0;
  index3 = 0;
  for (auto it = 0; it < np; ++it) {
    proc_to = (taskid + it) % np;
    proc_from = (taskid - it + np) % np;
    p_i1_start[nb][3][it] = index1;
    p_i2_start[nb][3][it] = index2;
    for (auto k = 0; k < nz; ++k)
      for (auto j = 0; j < ny; ++j)
        for (auto i = 0; i < (nx / 2 + 1); ++i) {
          iszero = (zero_arow2[i + k * (nx / 2 + 1)]);

          // phere = int_mb(p_map2(1,id)+(k-1)+(i-1)*nz(id))
          // pto   = int_mb(p_map1(1,id)+(j-1)+(k-1)*ny(id))

          phere = ijktop1(i, j, k);
          pto = ijktop2(i, j, k);

          /* packing scheme */
          if ((phere == taskid) && (pto == proc_to)) {
            if (!iszero) {
              p_iq_to_i1[nb][3][index1] = ijktoindex1(i, j, k);
              ++index1;
            }
          }
          /* unpacking scheme */
          if ((pto == taskid) && (phere == proc_from)) {
            if (!iszero) {
              p_iq_to_i2[nb][3][index2] = ijktoindex2t(i, j, k);
              ++index2;
            } else {
              p_iz_to_i2[nb][3][index3] = ijktoindex2t(i, j, k);
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
  for (auto it = 0; it < np; ++it) {
    proc_to = (taskid + it) % np;
    proc_from = (taskid - it + np) % np;
    p_i1_start[nb][4][it] = index1;
    p_i2_start[nb][4][it] = index2;
    for (auto k = 0; k < nz; ++k)
      for (auto j = 0; j < ny; ++j)
        for (auto i = 0; i < (nx / 2 + 1); ++i) {
          // phere = int_mb(p_map1(1,id)+(j-1)+(k-1)*ny(id))
          // pto   = int_mb(p_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))

          phere = ijktop2(i, j, k);
          pto = ijktop(i, j, k);

          /* packing scheme */
          if ((phere == taskid) && (pto == proc_to)) {
            p_iq_to_i1[nb][4][index1] = ijktoindex2t(i, j, k);
            ++index1;
          }
          /* unpacking scheme */
          if ((pto == taskid) && (phere == proc_from)) {
            p_iq_to_i2[nb][4][index2] = ijktoindex(i, j, k);
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
  for (auto it = 0; it < np; ++it) {
    proc_to = (taskid + it) % np;
    proc_from = (taskid - it + np) % np;
    p_i1_start[nb][5][it] = index1;
    p_i2_start[nb][5][it] = index2;
    for (auto k = 0; k < nz; ++k)
      for (auto j = 0; j < ny; ++j)
        for (auto i = 0; i < (nx / 2 + 1); ++i) {
          // phere = int_mb(p_map3(1,id)+(i-1)+(j-1)*(nx(id)/2+1))
          // pto   = int_mb(p_map1(1,id)+(j-1)+(k-1)*ny(id))

          phere = ijktop(i, j, k);
          pto = ijktop2(i, j, k);

          /* packing scheme */
          if ((phere == taskid) && (pto == proc_to)) {
            p_iq_to_i1[nb][5][index1] = ijktoindex(i, j, k);
            ++index1;
          }
          /* unpacking scheme */
          if ((pto == taskid) && (phere == proc_from)) {
            p_iq_to_i2[nb][5][index2] = ijktoindex2t(i, j, k);
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
 *         d3db::r_alloc        *
 *                              *
 ********************************/
double *d3db::r_alloc() {
  double *ptr = new (std::nothrow) double[n2ft3d]();
  return ptr;
}

/********************************
 *                              *
 *         d3db::r_nalloc       *
 *                              *
 ********************************/
double *d3db::r_nalloc(const int nn) {
  double *ptr = new (std::nothrow) double[n2ft3d * nn]();
  return ptr;
}

/********************************
 *                              *
 *         d3db::r_dealloc      *
 *                              *
 ********************************/
void d3db::r_dealloc(double *ptr) { delete[] ptr; }

/********************************
 *                              *
 *         d3db::r_zero         *
 *                              *
 ********************************/
void d3db::r_zero(double *ptr) {
  std::memset(ptr, 0, n2ft3d * sizeof(double));
  /*int i;
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
  */
}

/********************************
 *                              *
 *         d3db::r_nzero        *
 *                              *
 ********************************/
void d3db::r_nzero(int n, double *ptr) {
  std::memset(ptr, 0, n * n2ft3d * sizeof(double));
  /*int i;
  int m = (n*n2ft3d)%7;
  if (m>0)
     for (i=0; i<m; ++i)
        ptr[i] = 0.0;
  if ((n*n2ft3d)<7)
     return;

  for (i=m; i<(n*n2ft3d); i+=7)
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
  */
}

/********************************
 *                              *
 *         d3db::t_nzero        *
 *                              *
 ********************************/
void d3db::t_nzero(int n, double *ptr) {
  std::memset(ptr, 0, n * nfft3d * sizeof(double));
  /*
  int i;
  int m = (n*nfft3d)%7;
  if (m>0)
     for (i=0; i<m; ++i)
        ptr[i] = 0.0;
  if ((n*nfft3d)<7)
     return;

  for (i=m; i<(n*nfft3d); i+=7)
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
  */
}

/********************************
 *                              *
 *        d3db::rr_copy         *
 *                              *
 ********************************/
void d3db::rr_copy(const double *ptr1, double *ptr2) {
  std::memcpy(ptr2, ptr1, n2ft3d * sizeof(double));
  /*
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
  */
}

/********************************
 *                              *
 *        d3db::tt_copy         *
 *                              *
 ********************************/
void d3db::tt_copy(const double *ptr1, double *ptr2) {
  std::memcpy(ptr2, ptr1, nfft3d_map * sizeof(double));
  /*
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
  */
}

/********************************
 *                              *
 *        d3db::rr_SMul         *
 *                              *
 ********************************/
void d3db::rr_SMul(const double da, const double *ptr1, double *ptr2) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr2[i] = da * ptr1[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr2[i] = da * ptr1[i];
    ptr2[i + 1] = da * ptr1[i + 1];
    ptr2[i + 2] = da * ptr1[i + 2];
    ptr2[i + 3] = da * ptr1[i + 3];
    ptr2[i + 4] = da * ptr1[i + 4];
  }
  return;
}

/********************************
 *                              *
 *        d3db::rrr_SMulAdd     *
 *                              *
 ********************************/
void d3db::rrr_SMulAdd(const double da, const double *ptr1, const double *ptr2,
                       double *ptr3) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr3[i] = da * ptr1[i] + ptr2[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr3[i] = da * ptr1[i] + ptr2[i];
    ptr3[i + 1] = da * ptr1[i + 1] + ptr2[i + 1];
    ptr3[i + 2] = da * ptr1[i + 2] + ptr2[i + 2];
    ptr3[i + 3] = da * ptr1[i + 3] + ptr2[i + 3];
    ptr3[i + 4] = da * ptr1[i + 4] + ptr2[i + 4];
  }
  return;
}

/********************************
 *                              *
 *     d3db::rrrrr_SumMulAdd    *
 *                              *
 ********************************/
/*
   ptr5 = (ptr1+ptr2)*ptr3 + pt4
*/
void d3db::rrrrr_SumMulAdd(const double *ptr1, const double *ptr2,
                           const double *ptr3, const double *ptr4,
                           double *ptr5) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr5[i] = (ptr1[i] + ptr2[i]) * ptr3[i] + ptr4[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr5[i] = (ptr1[i] + ptr2[i]) * ptr3[i] + ptr4[i];
    ptr5[i + 1] = (ptr1[i + 1] + ptr2[i + 1]) * ptr3[i + 1] + ptr4[i + 1];
    ptr5[i + 2] = (ptr1[i + 2] + ptr2[i + 2]) * ptr3[i + 2] + ptr4[i + 2];
    ptr5[i + 3] = (ptr1[i + 3] + ptr2[i + 3]) * ptr3[i + 3] + ptr4[i + 3];
    ptr5[i + 4] = (ptr1[i + 4] + ptr2[i + 4]) * ptr3[i + 4] + ptr4[i + 4];
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
   int m = n2ft3d_map % 5;
   if (m > 0)
      for (i = 0; i < m; ++i)
         ptr2[i] *= da;
   if (n2ft3d_map < 5)
      return;
   for (i=m; i<n2ft3d_map; i+=5)
   {
      ptr2[i] *= da;
      ptr2[i + 1] *= da;
      ptr2[i + 2] *= da;
      ptr2[i + 3] *= da;
      ptr2[i + 4] *= da;
   }
}

/********************************
 *                              *
 *       d3db::r_zero_ends      *
 *                              *
 ********************************/
void d3db::r_zero_ends(double *A)
{
   int index;
 
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      for (auto q=0; q<nq; ++q)
         for (auto k=0; k<nz; ++k)
         {
            // index = nx + (k-1)*(nx+2) + (q-1)*(nx+2)*ny;
            index = nx + k * (nx+2) + q*(nx+2)*ny;
            A[index] = 0.0;
            A[index+1] = 0.0;
         }
   }
   /*************************
    **** hilbert mapping ****
    *************************/
   else
   {
      for (auto q=0; q<nq1; ++q) 
      {
         index = nx + q*(nx+2);
         A[index]   = 0.0;
         A[index+1] = 0.0;
      }
      if (n2ft3d_map<n2ft3d)
         std::memset(A+n2ft3d_map,0,(n2ft3d-n2ft3d_map)*sizeof(double));
   }
}

/********************************
 *                              *
 *       d3db::r_zero_mends     *
 *                              *
 ********************************/
void d3db::r_zero_mends(double *A)
{
   int index;
 
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) {
      for (auto q = 0; q < nq; ++q)
         for (auto k = 0; k < nz; ++k)
         {
            // index = nx + (k-1)*(nx+2) + (q-1)*(nx+2)*ny;
            index = nx + k * (nx + 2) + q * (nx + 2) * ny;
            A[index] = 0.0;
            A[index + 1] = 0.0;
         }
   }
   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      for (auto q=0; q<nq1; ++q)
      {
         index = nx + q*(nx+2);
         A[index]   = 0.0;
         A[index+1] = 0.0;
      }
      if (n2ft3d_map<n2ft3d)
         std::memset(A+n2ft3d_map,0,(n2ft3d-n2ft3d_map)*sizeof(double));
   }
 
   /* grid points in coordination space */
   int nzh = nz / 2;
   int nyh = ny / 2;
   int nxh = nx / 2;
   for (auto k3 = (-nzh); k3 < nzh; ++k3)
      for (auto k2 = (-nyh); k2 < nyh; ++k2)
         for (auto k1 = (-nxh); k1 < nxh; ++k1) 
         {
            int i = k1 + nxh;
            int j = k2 + nyh;
            int k = k3 + nzh;
            int indx = ijktoindex2(i, j, k);
            int p = ijktop2(i, j, k);
            
            if (p == parall->taskid_i())
               if ((k1 == (-nxh))  || (k1 == (nxh-1)) || (k2 == (-nyh)) ||
                   (k2 == (nyh-1)) || (k3 == (-nzh))  || (k3 == (nzh-1)))
                  A[indx] = 0.0;
         }
}

/********************************
 *                              *
 *          d3db::r_abs         *
 *                              *
 ********************************/
void d3db::r_abs(double *ptr2) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr2[i] = std::abs(ptr2[i]);
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr2[i] = std::abs(ptr2[i]);
    ptr2[i + 1] = std::abs(ptr2[i + 1]);
    ptr2[i + 2] = std::abs(ptr2[i + 2]);
    ptr2[i + 3] = std::abs(ptr2[i + 3]);
    ptr2[i + 4] = std::abs(ptr2[i + 4]);
  }
}

/********************************
 *                              *
 *          d3db::r_sqr         *
 *                              *
 ********************************/
void d3db::r_sqr(double *ptr2) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr2[i] *= ptr2[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr2[i] *= ptr2[i];
    ptr2[i + 1] *= ptr2[i + 1];
    ptr2[i + 2] *= ptr2[i + 2];
    ptr2[i + 3] *= ptr2[i + 3];
    ptr2[i + 4] *= ptr2[i + 4];
  }
}

/********************************
 *                              *
 *          d3db::rr_sqr        *
 *                              *
 ********************************/
void d3db::rr_sqr(const double *ptr2, double *ptr3) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr3[i] = ptr2[i] * ptr2[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr3[i] = ptr2[i] * ptr2[i];
    ptr3[i + 1] = ptr2[i + 1] * ptr2[i + 1];
    ptr3[i + 2] = ptr2[i + 2] * ptr2[i + 2];
    ptr3[i + 3] = ptr2[i + 3] * ptr2[i + 3];
    ptr3[i + 4] = ptr2[i + 4] * ptr2[i + 4];
  }
}

/********************************
 *                              *
 *        d3db::rr_addsqr       *
 *                              *
 ********************************/
void d3db::rr_addsqr(const double *ptr2, double *ptr3) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr3[i] += ptr2[i] * ptr2[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr3[i] += ptr2[i] * ptr2[i];
    ptr3[i + 1] += ptr2[i + 1] * ptr2[i + 1];
    ptr3[i + 2] += ptr2[i + 2] * ptr2[i + 2];
    ptr3[i + 3] += ptr2[i + 3] * ptr2[i + 3];
    ptr3[i + 4] += ptr2[i + 4] * ptr2[i + 4];
  }
}

/********************************
 *                              *
 *          d3db::r_sqrt        *
 *                              *
 ********************************/
void d3db::r_sqrt(double *ptr2) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr2[i] = sqrt(ptr2[i]);
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr2[i] = sqrt(ptr2[i]);
    ptr2[i + 1] = sqrt(ptr2[i + 1]);
    ptr2[i + 2] = sqrt(ptr2[i + 2]);
    ptr2[i + 3] = sqrt(ptr2[i + 3]);
    ptr2[i + 4] = sqrt(ptr2[i + 4]);
  }
}

/********************************
 *                              *
 *         d3db::r_dsum         *
 *                              *
 ********************************/
double d3db::r_dsum(const double *ptr) {
  int i;
  int m = n2ft3d_map % 5;
  double sum = 0.0;
  if (m > 0)
    for (i = 0; i < m; ++i)
      sum += ptr[i];
  if (n2ft3d_map < 5)
    return sum;
  for (i = m; i < n2ft3d_map; i += 5) {
    sum += ptr[i] + ptr[i + 1] + ptr[i + 2] + ptr[i + 3] + ptr[i + 4];
  }
  return parall->SumAll(1, sum);
}

/********************************
 *                              *
 *         d3db::rrr_Sum2Add    *
 *                              *
 ********************************/
void d3db::rrr_Sum2Add(const double *ptr1, const double *ptr2, double *ptr3) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr3[i] += ptr1[i] + ptr2[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr3[i] += ptr1[i] + ptr2[i];
    ptr3[i + 1] += ptr1[i + 1] + ptr2[i + 1];
    ptr3[i + 2] += ptr1[i + 2] + ptr2[i + 2];
    ptr3[i + 3] += ptr1[i + 3] + ptr2[i + 3];
    ptr3[i + 4] += ptr1[i + 4] + ptr2[i + 4];
  }
}

/********************************
 *                              *
 *         d3db::rrrr_Sum       *
 *                              *
 ********************************/
void d3db::rrrr_Sum(const double *ptr1, const double *ptr2, const double *ptr3, double *ptr4) 
{
   int i;
   int m = n2ft3d_map%5;
   if (m>0)
   {
      for (i=0; i<m; ++i)
         ptr4[i] = ptr1[i] + ptr2[i] + ptr3[i];
   }
   if (n2ft3d_map < 5)
      return;
   for (i=m; i<n2ft3d_map; i+=5)
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
 *         d3db::rrr_Sum        *
 *                              *
 ********************************/
void d3db::rrr_Sum(const double *ptr1, const double *ptr2, double *ptr3) 
{
   int i;
   int m = n2ft3d_map%5;
   if (m > 0)
      for (i = 0; i<m; ++i)
         ptr3[i] = ptr1[i] + ptr2[i];
   if (n2ft3d_map < 5)
      return;
   for (i=m; i<n2ft3d_map; i+=5)
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
 *         d3db::rr_Sum         *
 *                              *
 ********************************/
void d3db::rr_Sum(const double *ptr2, double *ptr3) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr3[i] += ptr2[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr3[i] += ptr2[i];
    ptr3[i + 1] += ptr2[i + 1];
    ptr3[i + 2] += ptr2[i + 2];
    ptr3[i + 3] += ptr2[i + 3];
    ptr3[i + 4] += ptr2[i + 4];
  }
}

/********************************
 *                              *
 *         d3db::rrr_Minus      *
 *                              *
 ********************************/
void d3db::rrr_Minus(const double *ptr1, const double *ptr2, double *ptr3) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr3[i] = ptr1[i] - ptr2[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr3[i] = ptr1[i] - ptr2[i];
    ptr3[i + 1] = ptr1[i + 1] - ptr2[i + 1];
    ptr3[i + 2] = ptr1[i + 2] - ptr2[i + 2];
    ptr3[i + 3] = ptr1[i + 3] - ptr2[i + 3];
    ptr3[i + 4] = ptr1[i + 4] - ptr2[i + 4];
  }
}

/********************************
 *                              *
 *         d3db::arrr_Minus     *
 *                              *
 ********************************/
void d3db::arrr_Minus(const double a, const double *ptr1, const double *ptr2,
                      double *ptr3) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr3[i] = a * (ptr1[i] - ptr2[i]);
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr3[i] = a * (ptr1[i] - ptr2[i]);
    ptr3[i + 1] = a * (ptr1[i + 1] - ptr2[i + 1]);
    ptr3[i + 2] = a * (ptr1[i + 2] - ptr2[i + 2]);
    ptr3[i + 3] = a * (ptr1[i + 3] - ptr2[i + 3]);
    ptr3[i + 4] = a * (ptr1[i + 4] - ptr2[i + 4]);
  }
}

/********************************
 *                              *
 *         d3db::rr_Minus       *
 *                              *
 ********************************/
void d3db::rr_Minus(const double *ptr2, double *ptr3) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr3[i] -= ptr2[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr3[i] -= ptr2[i];
    ptr3[i + 1] -= ptr2[i + 1];
    ptr3[i + 2] -= ptr2[i + 2];
    ptr3[i + 3] -= ptr2[i + 3];
    ptr3[i + 4] -= ptr2[i + 4];
  }
}

/********************************
 *                              *
 *         d3db::rrr_Mul        *
 *                              *
 ********************************/
void d3db::rrr_Mul(const double *ptr1, const double *ptr2, double *ptr3) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr3[i] = ptr1[i] * ptr2[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr3[i] = ptr1[i] * ptr2[i];
    ptr3[i + 1] = ptr1[i + 1] * ptr2[i + 1];
    ptr3[i + 2] = ptr1[i + 2] * ptr2[i + 2];
    ptr3[i + 3] = ptr1[i + 3] * ptr2[i + 3];
    ptr3[i + 4] = ptr1[i + 4] * ptr2[i + 4];
  }
  return;
}

/********************************
 *                              *
 *         d3db::rr_Mul         *
 *                              *
 ********************************/
void d3db::rr_Mul(const double *ptr1, double *ptr3) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr3[i] *= ptr1[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr3[i] *= ptr1[i];
    ptr3[i + 1] *= ptr1[i + 1];
    ptr3[i + 2] *= ptr1[i + 2];
    ptr3[i + 3] *= ptr1[i + 3];
    ptr3[i + 4] *= ptr1[i + 4];
  }
  return;
}

/********************************
 *                              *
 *     d3db::rrr_SqrMulAdd      *
 *                              *
 ********************************/
void d3db::rrr_SqrMulAdd(const double *ptr1, const double *ptr2, double *ptr3) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr3[i] += ptr1[i] * ptr1[i] * ptr2[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr3[i] += (ptr1[i] * ptr1[i]) * ptr2[i];
    ptr3[i + 1] += (ptr1[i + 1] * ptr1[i + 1]) * ptr2[i + 1];
    ptr3[i + 2] += (ptr1[i + 2] * ptr1[i + 2]) * ptr2[i + 2];
    ptr3[i + 3] += (ptr1[i + 3] * ptr1[i + 3]) * ptr2[i + 3];
    ptr3[i + 4] += (ptr1[i + 4] * ptr1[i + 4]) * ptr2[i + 4];
  }
  return;
}

/******************************************
 *                                        *
 *     d3db::rrrrrrr_Sqr3MulPlusMul2      *
 *                                        *
 ******************************************/
void d3db::rrrrrrr_Sqr3MulPlusMul2(const double *ptr1, const double *ptr2,
                                   const double *ptr3, const double *ptr4,
                                   const double *ptr5, const double *ptr6,
                                   double *ptr7) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr7[i] = (ptr1[i] * ptr1[i] + ptr2[i] * ptr2[i] + ptr3[i] * ptr3[i]) *
                    ptr4[i] +
                ptr5[i] * ptr6[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr7[i] =
        (ptr1[i] * ptr1[i] + ptr2[i] * ptr2[i] + ptr3[i] * ptr3[i]) * ptr4[i] +
        ptr5[i] * ptr6[i];
    ptr7[i + 1] = (ptr1[i + 1] * ptr1[i + 1] + ptr2[i + 1] * ptr2[i + 1] +
                   ptr3[i + 1] * ptr3[i + 1]) *
                      ptr4[i + 1] +
                  ptr5[i + 1] * ptr6[i + 1];
    ptr7[i + 2] = (ptr1[i + 1] * ptr1[i + 2] + ptr2[i + 2] * ptr2[i + 2] +
                   ptr3[i + 2] * ptr3[i + 2]) *
                      ptr4[i + 2] +
                  ptr5[i + 2] * ptr6[i + 2];
    ptr7[i + 3] = (ptr1[i + 1] * ptr1[i + 3] + ptr2[i + 3] * ptr2[i + 3] +
                   ptr3[i + 3] * ptr3[i + 3]) *
                      ptr4[i + 3] +
                  ptr5[i + 3] * ptr6[i + 3];
    ptr7[i + 4] = (ptr1[i + 1] * ptr1[i + 4] + ptr2[i + 4] * ptr2[i + 4] +
                   ptr3[i + 4] * ptr3[i + 4]) *
                      ptr4[i + 4] +
                  ptr5[i + 4] * ptr6[i + 4];
  }
  return;
}

/********************************
 *                              *
 *         d3db::tc_Mul         *
 *                              *
 ********************************/
void d3db::tc_Mul(const double *ptr1, double *ptr3) {
  int i;
  int m = nfft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i) {
      ptr3[2 * i] *= ptr1[i];
      ptr3[2 * i + 1] *= ptr1[i];
    }
  if (nfft3d_map < 5)
    return;
  for (i = m; i < nfft3d_map; i += 5) {
    ptr3[2 * (i)] *= ptr1[i];
    ptr3[2 * (i) + 1] *= ptr1[i];

    ptr3[2 * (i + 1)] *= ptr1[i + 1];
    ptr3[2 * (i + 1) + 1] *= ptr1[i + 1];

    ptr3[2 * (i + 2)] *= ptr1[i + 2];
    ptr3[2 * (i + 2) + 1] *= ptr1[i + 2];

    ptr3[2 * (i + 3)] *= ptr1[i + 3];
    ptr3[2 * (i + 3) + 1] *= ptr1[i + 3];

    ptr3[2 * (i + 4)] *= ptr1[i + 4];
    ptr3[2 * (i + 4) + 1] *= ptr1[i + 4];
  }
  return;
}

/********************************
 *                              *
 *         d3db::rrr_Mul2Add    *
 *                              *
 ********************************/
void d3db::rrr_Mul2Add(const double *ptr1, const double *ptr2, double *ptr3) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr3[i] += ptr1[i] * ptr2[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr3[i] += ptr1[i] * ptr2[i];
    ptr3[i + 1] += ptr1[i + 1] * ptr2[i + 1];
    ptr3[i + 2] += ptr1[i + 2] * ptr2[i + 2];
    ptr3[i + 3] += ptr1[i + 3] * ptr2[i + 3];
    ptr3[i + 4] += ptr1[i + 4] * ptr2[i + 4];
  }
  return;
}

/********************************
 *                              *
 *         d3db::rrr_Divide     *
 *                              *
 ********************************/
#define ETA_DIV 1.0e-9
void d3db::rrr_Divide(const double *ptr1, const double *ptr2, double *ptr3) 
{
   int i;
   int m = n2ft3d_map % 5;
   if (m > 0)
     for (i = 0; i < m; ++i)
       ptr3[i] = (std::abs(ptr2[i])>ETA_DIV) ? (ptr1[i]/ptr2[i]) : (0.0);
   if (n2ft3d_map < 5)
     return;
   for (i = m; i < n2ft3d_map; i += 5) 
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
 *         d3db::rr_Divide      *
 *                              *
 ********************************/
void d3db::rr_Divide(const double *ptr2, double *ptr3) 
{
   int i;
   int m = n2ft3d_map%5;
   if (m > 0)
      for (i = 0; i < m; ++i)
         ptr3[i] = (std::abs(ptr2[i]) > ETA_DIV) ? (ptr3[i]/ptr2[i]) : (0.0);
   if (n2ft3d_map<5)
      return;
   for (i = m; i < n2ft3d_map; i += 5) 
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
 *         d3db::rr_screen0     *
 *                              *
 ********************************/
// This subroutine is used to calculate screen0 = (1.0/epsilon-1.0)
void d3db::rr_screen0(const double *ptr2, double *ptr3)
{
   int i;
   int m = n2ft3d_map%5;
   if (m>0)
      for (i=0; i<m; ++i)
         ptr3[i] = (std::abs(ptr2[i]) > ETA_DIV) ? (1.0/ptr2[i]-1.0) : (0.0);
   if (n2ft3d_map<5)
      return;
   for (i=m; i<n2ft3d_map; i+=5)
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
 *         d3db::rr_daxpy       *
 *                              *
 ********************************/
void d3db::rr_daxpy(const double alpha, const double *ptr1, double *ptr2) {
  int i;
  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      ptr2[i] += alpha * ptr1[i];
  if (n2ft3d_map < 5)
    return;
  for (i = m; i < n2ft3d_map; i += 5) {
    ptr2[i] += alpha * ptr1[i];
    ptr2[i + 1] += alpha * ptr1[i + 1];
    ptr2[i + 2] += alpha * ptr1[i + 2];
    ptr2[i + 3] += alpha * ptr1[i + 3];
    ptr2[i + 4] += alpha * ptr1[i + 4];
  }
}

/********************************
 *                              *
 *         d3db::rr_dot         *
 *                              *
 ********************************/
double d3db::rr_dot(const double *ptr1, const double *ptr2) {
  int i;
  double sum = 0.0;

  int m = n2ft3d_map % 5;
  if (m > 0)
    for (i = 0; i < m; ++i)
      sum += ptr1[i] * ptr2[i];
  if (n2ft3d_map < 5)
    return sum;
  for (i = m; i < n2ft3d_map; i += 5) {
    sum += ptr1[i] * ptr2[i] + ptr1[i + 1] * ptr2[i + 1] +
           ptr1[i + 2] * ptr2[i + 2] + ptr1[i + 3] * ptr2[i + 3] +
           ptr1[i + 4] * ptr2[i + 4];
  }

  return parall->SumAll(1, sum);
}

/********************************
 *                              *
 *         d3db::nrr_vdot       *
 *                              *
 ********************************/
void d3db::nrr_vdot(const int n, const double *ptr1, const double *ptr2,
                    double *v) {
  int m = n2ft3d_map % 5;
  for (auto k = 0; k < n; ++k)
    v[k] = 0.0;
  if (m > 0)
    for (auto i = 0; i < m; ++i)
      for (auto k = 0; k < n; ++k)
        v[k] += ptr1[n * i + k] * ptr2[i];
  if (n2ft3d_map >= 5) {
    for (auto i = m; i < n2ft3d_map; i += 5)
      for (auto k = 0; k < n; ++k) {
        v[k] += ptr1[n * i + k] * ptr2[i] +
                ptr1[n * (i + 1) + k] * ptr2[i + 1] +
                ptr1[n * (i + 2) + k] * ptr2[i + 2] +
                ptr1[n * (i + 3) + k] * ptr2[i + 3] +
                ptr1[n * (i + 4) + k] * ptr2[i + 4];
      }
    parall->Vector_SumAll(1, n, v);
  }
}

/********************************
 *                              *
 *         d3db::t_alloc        *
 *                              *
 ********************************/
double *d3db::t_alloc() {
  double *ptr = new (std::nothrow) double[nfft3d]();
  return ptr;
}

/********************************
 *                              *
 *         d3db::t_dealloc      *
 *                              *
 ********************************/
void d3db::t_dealloc(double *ptr) { delete[] ptr; }

/********************************
 *                              *
 *         d3db::c_read         *
 *                              *
 ********************************/
void d3db::c_read(const int iunit, double *a, const int jcol) {
  int jstart, jend, fillcolumn, index, ii, jj, p_to, p_here;
  int taskid = parall->taskid();
  int taskid_j = parall->taskid_j();
  int np_j = parall->np_j();

  if (jcol < 0) {
    jstart = 0;
    jend = np_j - 1;
    fillcolumn = 1;
  } else {
    jstart = jend = jcol;
    fillcolumn = (taskid_j == jcol);
  }

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) {
    double *tmp = new (std::nothrow) double[(nx + 2) * ny]();
    // double tmp[(nx+2)*ny];
    int bsize = (nx + 2) * ny;

    /**** master node reads from file and distributes ****/
    if (taskid == MASTER)
      for (int k = 0; k < nz; ++k) {
        dread(iunit, tmp, bsize);

        index = 2 * ijktoindex(0, 0, k);
        ii = ijktop(0, 0, k);
        for (jj = jstart; jj <= jend; ++jj) {
          p_to = parall->convert_taskid_ij(ii, jj);
          if (p_to == MASTER)
            for (int k = 0; k < bsize; ++k)
              a[index + k] = tmp[k];
          else
            parall->dsend(0, 9, p_to, bsize, tmp);
        }
      }

    /**** not master node ****/
    else if (fillcolumn)
      for (int k = 0; k < nz; ++k) {
        index = 2 * ijktoindex(0, 0, k);
        ii = ijktop(0, 0, k);
        p_here = parall->convert_taskid_ij(ii, taskid_j);
        if (p_here == taskid) {
          parall->dreceive(0, 9, MASTER, bsize, tmp);
          for (int k = 0; k < bsize; ++k)
            a[index + k] = tmp[k];
        }
      }

    delete[] tmp;
  }

  /*************************
   **** hilbert mapping ****
   *************************/
  else {
    // double *tmp = new (std::nothrow) double[nx+2]();
    double tmp[nx + 2];
    int bsize = (nx + 2);

    /**** master node reads from file and distributes ****/
    if (taskid == MASTER)
      for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j) {

          dread(iunit, tmp, bsize);

          index = ijktoindex2(0, j, k);
          ii = ijktop2(0, j, k);
          for (int jj = jstart; jj <= jend; ++jj) {
            p_to = parall->convert_taskid_ij(ii, jj);

            if (p_to == MASTER)
              for (int k = 0; k < bsize; ++k)
                a[index + k] = tmp[k];
            else
              parall->dsend(0, 9, p_to, bsize, tmp);
          }
        }

    /**** not master node ****/
    else if (fillcolumn)
      for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j) {
          index = ijktoindex2(0, j, k);
          ii = ijktop2(0, j, k);
          p_here = parall->convert_taskid_ij(ii, taskid_j);
          if (p_here == taskid) {
            parall->dreceive(0, 9, MASTER, bsize, tmp);
            for (int k = 0; k < bsize; ++k)
              a[index + k] = tmp[k];
          }
        }

    // double tmp1[2*nfft3d];
    // double tmp2[2*nfft3d];
    double *tmp1 = new (std::nothrow) double[2 * nfft3d]();
    double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
    c_transpose_ijk(4, a, tmp1, tmp2);
    delete[] tmp2;
    delete[] tmp1;
  }
}

/********************************
 *                              *
 *         d3db::c_write        *
 *                              *
 ********************************/
void d3db::c_write(const int iunit, double *a, const int jcol) {
  int index, ii, jj, p_from, p_here;
  int taskid = parall->taskid();
  int taskid_j = parall->taskid_j();
  int np_j = parall->np_j();

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) {
    double *tmp = new (std::nothrow) double[(nx + 2) * ny]();
    // double tmp[(nx+2)*ny];
    int bsize = (nx + 2) * ny;

    /**** master node gathers and write to file ****/
    if (taskid == MASTER)
      for (int k = 0; k < nz; ++k) {
        ii = ijktop(0, 0, k);
        p_from = parall->convert_taskid_ij(ii, jcol);
        if (p_from == MASTER) {
          index = 2 * ijktoindex(0, 0, k);
          for (int k = 0; k < bsize; ++k)
            tmp[k] = a[index + k];
        } else {
          parall->dreceive(0, 9, p_from, bsize, tmp);
        }
        dwrite(iunit, tmp, bsize);
      }

    /**** not master node ****/
    else
      for (int k = 0; k < nz; ++k) {
        index = 2 * ijktoindex(0, 0, k);
        ii = ijktop(0, 0, k);
        p_here = parall->convert_taskid_ij(ii, taskid_j);
        if (p_here == taskid) {
          for (int k = 0; k < bsize; ++k)
            tmp[k] = a[index + k];
          parall->dsend(0, 9, MASTER, bsize, tmp);
        }
      }

    delete[] tmp;
  }

  /*************************
   **** hilbert mapping ****
   *************************/
  else {
    // double *tmp1 = new (std::nothrow) double[2*nfft3d];
    // double *tmp2 = new (std::nothrow) double[2*nfft3d];
    {
      double *tmp1 = new (std::nothrow) double[2 * nfft3d]();
      double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
      c_transpose_ijk(5, a, tmp1, tmp2);
      delete[] tmp2;
      delete[] tmp1;
    }
    // delete [] tmp2;
    // delete [] tmp1;

    // double *tmp = new (std::nothrow) double[nx+2];
    double tmp[nx + 2];
    int bsize = (nx + 2);

    /**** master node write to file and fetches from other nodes ****/
    if (taskid == MASTER)
      for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j) {
          ii = ijktop2(0, j, k);
          p_from = parall->convert_taskid_ij(ii, jcol);
          if (p_from == MASTER) {
            index = ijktoindex2(0, j, k);
            for (int k = 0; k < bsize; ++k)
              tmp[k] = a[index + k];
          } else {
            parall->dreceive(0, 9, p_from, bsize, tmp);
          }
          dwrite(iunit, tmp, bsize);
        }

    /**** not master node ****/
    else
      for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j) {
          ii = ijktop2(0, j, k);
          p_here = parall->convert_taskid_ij(ii, taskid_j);
          if (p_here == taskid) {
            index = ijktoindex2(0, j, k);
            for (int k = 0; k < bsize; ++k)
              tmp[k] = a[index + k];
            parall->dsend(0, 9, MASTER, bsize, tmp);
          }
        }

    // delete [] tmp;
  }
}

/***********************************************
 *                                            *
 *            d3db::r_formatwrite             *
 *                                            *
 **********************************************/
std::string d3db::r_formatwrite(double *a) {
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
          int index = ijktoindex2(i, j, k);
          int p_from = ijktop2(i, j, k);
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
          int index = ijktoindex2(i, j, k);
          int p_here = ijktop2(i, j, k);
          if (p_here == taskid)
            parall->dsend(0, 189, MASTER, 1, a + index);
        }
      }
  }

  return stream.str();
}

/***********************************************
 *                                            *
 *         d3db::r_formatwrite_reverse        *
 *                                            *
 **********************************************/
std::string d3db::r_formatwrite_reverse(double *a) {
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
          int index = ijktoindex2(i, j, k);
          int p_from = ijktop2(i, j, k);
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
          int index = ijktoindex2(i, j, k);
          int p_here = ijktop2(i, j, k);
          if (p_here == taskid)
            parall->dsend(0, 189, MASTER, 1, a + index);
        }
      }
  }

  return stream.str();
}

void d3db::cshift1_fftb(const int n1, const int n2, const int n3, const int n4,
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
void d3db::cshift_fftf(const int n1, const int n2, const int n3, const int n4,
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
static void cshift_fftf_ab(const int n1, const int n2, const int n3,
                           const int n4, double *a, double *b) {
  int i, j, indx;
  indx = 0;
  for (j = 0; j < (n2 * n3 * n4); ++j) {
    for (i = n1; i > 0; --i) {
      b[indx + i] = a[indx + i - 1];
    }
    b[indx + 1] = 0.0;
    b[indx + n1 + 1] = 0.0;
    indx += (n1 + 2);
  }
}
void d3db::zeroend_fftb(const int n1, const int n2, const int n3, const int n4,
                        double *a) {
  int i, indx;
  indx = n1 + 1;
  for (i = 0; i < (n2 * n3 * n4); ++i) {
    a[indx - 1] = 0.0;
    a[indx + 1 - 1] = 0.0;
    indx += (n1 + 2);
  }
}

/********************************
 *                              *
 *         d3db::cr_fft3d       *
 *                              *
 ********************************/
void d3db::cr_fft3d(double *a) {

  nwpw_timing_function ftime(1);
  int i, j, k, jj, kk, q, indx, indx0, nxh, nxh2, nxhy, nxhy2, nxhz, nxhz2, nn,
      shift;
  double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
  double *tmp3 = new (std::nothrow) double[2 * nfft3d]();

  nxh = nx / 2 + 1;
  nxhy = nxh * ny;
  nxhz = nxh * nz;
  nxh2 = nx + 2;
  nxhy2 = nxh2 * ny;
  nxhz2 = nxh2 * nz;

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) {
    /***************************************************
     ***     do fft along kz dimension               ***
     ***     A(kx,nz,ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
     ***************************************************/
    indx0 = 0;
    nn = 0;
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        kk = 0;
        indx = 2 * i + indx0;
        shift = 2 * nz * nn;
        for (k = 0; k < nz; ++k) {
          tmp2[kk + shift] = a[indx];
          tmp2[kk + 1 + shift] = a[indx + 1];
          kk += 2;
          indx += nxh2;
        }
        ++nn;
      }
      indx0 += nxhz2;
    }

    mygdevice.batch_cfftz_tmpz(fft_tag,false, nz, nn, n2ft3d, tmp2, tmpz);

    indx0 = 0;
    nn = 0;
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        kk = 0;
        indx = 2 * i + indx0;
        shift = 2 * nz * nn;
        for (k = 0; k < nz; ++k) {
          a[indx] = tmp2[kk + shift];
          a[indx + 1] = tmp2[kk + 1 + shift];
          kk += 2;
          indx += nxh2;
        }
        ++nn;
      }
      indx0 += nxhz2;
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
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        jj = 0;
        indx = 2 * i + indx0;
        shift = 2 * ny * nn;
        for (j = 0; j < ny; ++j) {
          tmp2[jj + shift] = a[indx];
          tmp2[jj + 1 + shift] = a[indx + 1];
          jj += 2;
          indx += nxh2;
        }
        ++nn;
      }
      indx0 += nxhy2;
    }

    mygdevice.batch_cffty_tmpy(fft_tag,false, ny, nn, n2ft3d, tmp2, tmpy);

    indx0 = 0;
    nn = 0;
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        jj = 0;
        indx = 2 * i + indx0;
        shift = 2 * ny * nn;
        for (j = 0; j < ny; ++j) {
          a[indx] = tmp2[jj + shift];
          a[indx + 1] = tmp2[jj + 1 + shift];
          jj += 2;
          indx += nxh2;
        }
        ++nn;
      }
      indx0 += nxhy2;
    }
    
    mygdevice.batch_cfftx_tmpx(fft_tag,false, nx, ny * nq, n2ft3d, a, tmpx);

  }
  /*************************
   **** hilbert mapping ****
   *************************/
  else {

    /************************************************
     ***     do fft along kz dimension            ***
     ***   A(nz,kx,ky) <- fft1d^(-1)[A(kz,kx,ky)] ***
     ************************************************/
    mygdevice.batch_cfftz_tmpz(fft_tag,false, nz, nq3, n2ft3d, a, tmpz);
    
    c_transpose_ijk(2, a, tmp2, tmp3);

    /************************************************
     ***     do fft along ky dimension            ***
     ***   A(ny,nz,kx) <- fft1d^(-1)[A(ky,nz,kx)] ***
     ************************************************/
    mygdevice.batch_cffty_tmpy(fft_tag,false, ny, nq2, n2ft3d, a, tmpy);
   
    c_transpose_ijk(3, a, tmp2, tmp3);

    /************************************************
     ***     do fft along kx dimension            ***
     ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
     ************************************************/
    mygdevice.batch_cfftx_tmpx(fft_tag,false, nx, nq1, n2ft3d, a, tmpx);
   
    zeroend_fftb(nx, nq1, 1, 1, a);
    if (n2ft3d_map < n2ft3d)
      std::memset(a + n2ft3d_map, 0, (n2ft3d - n2ft3d_map) * sizeof(double));
  }

  delete[] tmp3;
  delete[] tmp2;
}

/********************************
 *                              *
 *         d3db::rc_fft3d       *
 *                              *
 ********************************/
void d3db::rc_fft3d(double *a) {

  nwpw_timing_function ftime(1);
  int i, j, k, jj, kk, q, indx, indx0, nxh, nxh2, nxhy, nxhy2, nxhz, nxhz2, nn,
      shift;
  double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
  double *tmp3 = new (std::nothrow) double[2 * nfft3d]();

  nxh = nx / 2 + 1;
  nxhy = nxh * ny;
  nxhz = nxh * nz;
  nxh2 = nx + 2;
  nxhy2 = nxh2 * ny;
  nxhz2 = nxh2 * nz;

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) {

    /********************************************
     ***     do fft along nx dimension        ***
     ***   A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]  ***
     ********************************************/
    mygdevice.batch_cfftx_tmpx(fft_tag,true, nx, ny*nq, n2ft3d, a, tmpx);

    /********************************************
     ***     do fft along ny dimension        ***
     ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
     ********************************************/
    indx0 = 0;
    nn = 0;
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        jj = 0;
        indx = 2 * i + indx0;
        shift = 2 * ny * nn;
        for (j = 0; j < ny; ++j) {
          tmp2[jj + shift] = a[indx];
          tmp2[jj + 1 + shift] = a[indx + 1];
          jj += 2;
          indx += nxh2;
        }
        ++nn;
      }
      indx0 += nxhy2;
    }

    mygdevice.batch_cffty_tmpy(fft_tag,true, ny, nn, n2ft3d, tmp2, tmpy);

    indx0 = 0;
    nn = 0;
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        jj = 0;
        indx = 2 * i + indx0;
        shift = 2 * ny * nn;
        for (j = 0; j < ny; ++j) {
          a[indx] = tmp2[jj + shift];
          a[indx + 1] = tmp2[jj + 1 + shift];
          jj += 2;
          indx += nxh2;
        }
        ++nn;
      }
      indx0 += nxhy2;
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
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        kk = 0;
        indx = 2 * i + indx0;
        shift = 2 * nz * nn;
        for (k = 0; k < nz; ++k) {
          tmp2[kk + shift] = a[indx];
          tmp2[kk + 1 + shift] = a[indx + 1];
          kk += 2;
          indx += nxh2;
        }
        ++nn;
      }
      indx0 += nxhz2;
    }

    mygdevice.batch_cfftz_tmpz(fft_tag,true, nz, nn, n2ft3d, tmp2, tmpz);

    indx0 = 0;
    nn = 0;
    for (q = 0; q < nq; ++q) {
      for (i = 0; i < nxh; ++i) {
        kk = 0;
        indx = 2 * i + indx0;
        shift = 2 * nz * nn;
        for (k = 0; k < nz; ++k) {
          a[indx] = tmp2[kk + shift];
          a[indx + 1] = tmp2[kk + 1 + shift];
          kk += 2;
          indx += nxh2;
        }
        ++nn;
      }
      indx0 += nxhz2;
    }
    
  }
  /*************************
   **** hilbert mapping ****
   *************************/
  else {
    /********************************************
     ***     do fft along nx dimension        ***
     ***   A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]  ***
     ********************************************/
    mygdevice.batch_cfftx_tmpx(fft_tag,true, nx, nq1, n2ft3d, a, tmpx);
    
    c_transpose_ijk(0, a, tmp2, tmp3);

    /********************************************
     ***     do fft along ny dimension        ***
     ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
     ********************************************/
    mygdevice.batch_cffty_tmpy(fft_tag,true, ny, nq2, n2ft3d, a, tmpy);
   
    c_transpose_ijk(1, a, tmp2, tmp3);

    /********************************************
     ***     do fft along nz dimension        ***
     ***   A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]  ***
     ********************************************/
    mygdevice.batch_cfftz_tmpz(fft_tag,true, nz, nq3, n2ft3d, a, tmpz);
  }

  delete[] tmp3;
  delete[] tmp2;
}

/********************************
 *                              *
 *         d3db::t_read         *
 *                              *
 ********************************/
void d3db::t_read(const int iunit, double *a, const int jcol) {
  int jstart, jend, fillcolumn, index, ii, jj, p_to, p_here;
  int taskid = parall->taskid();
  int taskid_j = parall->taskid_j();
  int np_j = parall->np_j();

  if (jcol < 0) {
    jstart = 0;
    jend = np_j - 1;
    fillcolumn = 1;
  } else {
    jstart = jend = jcol;
    fillcolumn = (taskid_j == jcol);
  }

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) {
    double *tmp = new (std::nothrow) double[(nx / 2 + 1) * ny]();
    // double tmp[(nx/2+1)*ny];
    int bsize = (nx / 2 + 1) * ny;

    /**** master node reads from file and distributes ****/
    if (taskid == MASTER)
      for (int k = 0; k < nz; ++k) {
        dread(iunit, tmp, bsize);

        index = ijktoindex(0, 0, k);
        ii = ijktop(0, 0, k);
        for (jj = jstart; jj <= jend; ++jj) {
          p_to = parall->convert_taskid_ij(ii, jj);

          if (p_to == MASTER)
            for (int k = 0; k < bsize; ++k)
              a[index + k] = tmp[k];
          else
            parall->dsend(0, 9, p_to, bsize, tmp);
        }
      }

    /**** not master node ****/
    else if (fillcolumn)
      for (int k = 0; k < nz; ++k) {
        index = ijktoindex(0, 0, k);
        ii = ijktop(0, 0, k);
        p_here = parall->convert_taskid_ij(ii, taskid_j);
        if (p_here == taskid) {
          parall->dreceive(0, 9, MASTER, bsize, tmp);
          for (int k = 0; k < bsize; ++k)
            a[index + k] = tmp[k];
        }
      }

    delete[] tmp;
  }

  /*************************
   **** hilbert mapping ****
   *************************/
  else {
    // double *tmp = new (std::nothrow) double[nx/2+1]();
    double tmp[nx / 2 + 1];
    int bsize = (nx / 2 + 1);

    /**** master node reads from file and distributes ****/
    if (taskid == MASTER)
      for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j) {

          dread(iunit, tmp, bsize);

          index = ijktoindex2t(0, j, k);
          ii = ijktop2(0, j, k);
          for (int jj = jstart; jj <= jend; ++jj) {
            p_to = parall->convert_taskid_ij(ii, jj);

            if (p_to == MASTER)
              for (int k = 0; k < bsize; ++k)
                a[index + k] = tmp[k];
            else
              parall->dsend(0, 9, p_to, bsize, tmp);
          }
        }

    /**** not master node ****/
    else if (fillcolumn)
      for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j) {
          index = ijktoindex2t(0, j, k);
          ii = ijktop2(0, j, k);
          p_here = parall->convert_taskid_ij(ii, taskid_j);
          if (p_here == taskid) {
            parall->dreceive(0, 9, MASTER, bsize, tmp);
            for (int k = 0; k < bsize; ++k)
              a[index + k] = tmp[k];
          }
        }
    // delete [] tmp;

    // double tmp1[2*nfft3d];
    // double tmp2[2*nfft3d];
    double *tmp1 = new (std::nothrow) double[2 * nfft3d]();
    double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
    t_transpose_ijk(4, a, tmp1, tmp2);
    delete[] tmp2;
    delete[] tmp1;
  }
}

/********************************
 *                              *
 *         d3db::t_write        *
 *                              *
 ********************************/
void d3db::t_write(const int iunit, double *a, const int jcol) {
  int index, ii, jj, p_from, p_here;
  int taskid = parall->taskid();
  int taskid_j = parall->taskid_j();
  int np_j = parall->np_j();

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) {
    double *tmp = new (std::nothrow) double[(nx + 2) * ny]();
    // double tmp[(nx+2)*ny];
    int bsize = (nx / 2 + 1) * ny;

    /**** master node gathers and write to file ****/
    if (taskid == MASTER)
      for (int k = 0; k < nz; ++k) {
        ii = ijktop(0, 0, k);
        p_from = parall->convert_taskid_ij(ii, jcol);
        if (p_from == MASTER) {
          index = ijktoindex(0, 0, k);
          for (int k = 0; k < bsize; ++k)
            tmp[k] = a[index + k];
        } else {
          parall->dreceive(0, 9, p_from, bsize, tmp);
        }
        dwrite(iunit, tmp, bsize);
      }

    /**** not master node ****/
    else
      for (int k = 0; k < nz; ++k) {
        index = ijktoindex(0, 0, k);
        ii = ijktop(0, 0, k);
        p_here = parall->convert_taskid_ij(ii, taskid_j);
        if (p_here == taskid) {
          for (int k = 0; k < bsize; ++k)
            tmp[k] = a[index + k];
          parall->dsend(0, 9, MASTER, bsize, tmp);
        }
      }

    delete[] tmp;
  }

  /*************************
   **** hilbert mapping ****
   *************************/
  else {
    { // double tmp1[2*nfft3d];
      // double tmp2[2*nfft3d];
      double *tmp1 = new (std::nothrow) double[2 * nfft3d]();
      double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
      t_transpose_ijk(5, a, tmp1, tmp2);
      delete[] tmp2;
      delete[] tmp1;
    }

    // double *tmp = new (std::nothrow) double[nx/2+1];
    double tmp[nx / 2 + 1];
    int bsize = (nx / 2 + 1);

    /**** master node write to file and fetches from other nodes ****/
    if (taskid == MASTER)
      for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j) {
          ii = ijktop2(0, j, k);
          p_from = parall->convert_taskid_ij(ii, jcol);
          if (p_from == MASTER) {
            index = ijktoindex2t(0, j, k);
            for (int k = 0; k < bsize; ++k)
              tmp[k] = a[index + k];
          } else {
            parall->dreceive(0, 9, p_from, bsize, tmp);
          }
          dwrite(iunit, tmp, bsize);
        }

    /**** not master node ****/
    else
      for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j) {
          ii = ijktop2(0, j, k);
          p_here = parall->convert_taskid_ij(ii, taskid_j);
          if (p_here == taskid) {
            index = ijktoindex2t(0, j, k);
            for (int k = 0; k < bsize; ++k)
              tmp[k] = a[index + k];
            parall->dsend(0, 9, MASTER, bsize, tmp);
          }
        }

    // delete [] tmp;
  }
}

/**************************************
 *                                    *
 *    d3db::c_ptranspose1_jk_start    *
 *                                    *
 **************************************/
void d3db::c_ptranspose1_jk_start(const int nb, double *a, double *tmp1,
                                  double *tmp2, const int request_indx,
                                  const int msgtype) {
  int it, proc_from, proc_to;
  int msglen;

  int n1 = p_i1_start[nb][0][np];

  c_aindexcopy(n1, p_iq_to_i1[nb][0], a, tmp1);

  /* it = 0, transpose data on same thread */
  msglen = 2 * (p_i2_start[nb][0][1] - p_i2_start[nb][0][0]);
  // int one=1;
  // DCOPY_PWDFT(msglen,&(tmp1[2*p_i1_start[nb][0][0]]),one,&(tmp2[2*p_i2_start[nb][0][0]]),one);
  std::memcpy(tmp2 + 2 * p_i2_start[nb][0][0], tmp1 + 2 * p_i1_start[nb][0][0],
              msglen * sizeof(double));

  /* receive packed array data */
  for (it = 1; it < np; ++it) {
    /* synchronous receive of tmp */
    proc_from = (taskid - it + np) % np;
    msglen = 2 * (p_i2_start[nb][0][it + 1] - p_i2_start[nb][0][it]);
    if (msglen > 0)
      parall->adreceive(request_indx, msgtype, proc_from, msglen,
                        tmp2 + 2 * p_i2_start[nb][0][it]);
    // parall->adreceive(request_indx,msgtype,proc_from,msglen,&tmp2[2*p_i2_start[nb][0][it]]);
  }
  for (it = 1; it < np; ++it) {
    proc_to = (taskid + it) % np;
    msglen = 2 * (p_i1_start[nb][0][it + 1] - p_i1_start[nb][0][it]);
    if (msglen > 0)
      parall->adsend(request_indx, msgtype, proc_to, msglen,
                     tmp1 + 2 * p_i1_start[nb][0][it]);
    // parall->adsend(request_indx,msgtype,proc_to,msglen,&tmp1[2*p_i1_start[nb][0][it]]);
  }
}

/**************************************
 *                                    *
 *    d3db::c_ptranspose1_jk_end      *
 *                                    *
 **************************************/
void d3db::c_ptranspose1_jk_end(const int nb, double *a, double *tmp2,
                                const int request_indx) {
  parall->awaitall(request_indx);

  int n2 = p_i2_start[nb][0][np];
  c_bindexcopy(n2, p_iq_to_i2[nb][0], tmp2, a);
  c_bindexzero(nfft3d - n2, p_iz_to_i2[nb][0], a);
}

/**************************************
 *                                    *
 *    d3db::c_ptranspose2_jk_start    *
 *                                    *
 **************************************/
void d3db::c_ptranspose2_jk_start(const int nb, double *a, double *tmp1,
                                  double *tmp2, const int request_indx,
                                  const int msgtype) {
  int it, proc_from, proc_to;
  int msglen;

  int n1 = p_j1_start[nb][0][np];

  c_aindexcopy(n1, p_jq_to_i1[nb][0], a, tmp1);

  /* it = 0, transpose data on same thread */
  msglen = 2 * (p_j2_start[nb][0][1] - p_j2_start[nb][0][0]);
  // int one=1;
  // DCOPY_PWDFT(msglen,&(tmp1[2*p_j1_start[nb][0][0]]),one,&(tmp2[2*p_j2_start[nb][0][0]]),one);
  std::memcpy(tmp2 + 2 * p_j2_start[nb][0][0], tmp1 + 2 * p_j1_start[nb][0][0],
              msglen * sizeof(double));

  /* receive packed array data */
  for (it = 1; it < np; ++it) {
    /* synchronous receive of tmp */
    proc_from = (taskid - it + np) % np;
    msglen = 2 * (p_j2_start[nb][0][it + 1] - p_j2_start[nb][0][it]);
    if (msglen > 0)
      parall->adreceive(request_indx, msgtype, proc_from, msglen,
                        &tmp2[2 * p_j2_start[nb][0][it]]);
  }
  for (it = 1; it < np; ++it) {
    proc_to = (taskid + it) % np;
    msglen = 2 * (p_j1_start[nb][0][it + 1] - p_j1_start[nb][0][it]);
    if (msglen > 0)
      parall->adsend(request_indx, msgtype, proc_to, msglen,
                     &tmp1[2 * p_j1_start[nb][0][it]]);
  }
}

/**************************************
 *                                    *
 *    d3db::c_ptranspose2_jk_end      *
 *                                    *
 **************************************/
void d3db::c_ptranspose2_jk_end(const int nb, double *a, double *tmp2,
                                const int request_indx) {
  parall->awaitall(request_indx);

  int n2 = p_j2_start[nb][0][np];
  c_bindexcopy(n2, p_jq_to_i2[nb][0], tmp2, a);
  c_bindexzero(nfft3d - n2, p_jz_to_i2[nb][0], a);
}

/**************************************
 *                                    *
 *    d3db::c_ptranspose_ijk_start    *
 *                                    *
 **************************************/
void d3db::c_ptranspose_ijk_start(const int nb, const int op, double *a,
                                  double *tmp1, double *tmp2,
                                  const int request_indx, const int msgtype) {
  int nnfft3d, it, proc_from, proc_to;
  int msglen;

  int n1 = p_i1_start[nb][op][np];

  /* pack a array - tmp1->tmp2 */
  c_aindexcopy(n1, p_iq_to_i1[nb][op], a, tmp1);

  /* it = 0, transpose data on same thread  - tmp2->tmp1*/
  msglen = 2 * (p_i2_start[nb][op][1] - p_i2_start[nb][op][0]);
  // int one=1;
  // DCOPY_PWDFT(msglen,&(tmp1[2*p_i1_start[nb][op][0]]),one,&(tmp2[2*p_i2_start[nb][op][0]]),one);
  // std::memcpy(&(tmp2[2*p_i2_start[nb][op][0]]),&(tmp1[2*p_i1_start[nb][op][0]]),msglen*sizeof(double));
  std::memcpy(tmp2 + 2 * p_i2_start[nb][op][0],
              tmp1 + 2 * p_i1_start[nb][op][0], msglen * sizeof(double));

  /* receive packed array data */
  for (it = 1; it < np; ++it) {
    /* synchronous receive of tmp */
    proc_from = (taskid - it + np) % np;
    msglen = 2 * (p_i2_start[nb][op][it + 1] - p_i2_start[nb][op][it]);
    if (msglen > 0)
      parall->adreceive(request_indx, msgtype, proc_from, msglen,
                        tmp2 + 2 * p_i2_start[nb][op][it]);
    // parall->adreceive(request_indx,msgtype,proc_from,msglen,&tmp2[2*p_i2_start[nb][op][it]]);
  }
  for (it = 1; it < np; ++it) {
    proc_to = (taskid + it) % np;
    msglen = 2 * (p_i1_start[nb][op][it + 1] - p_i1_start[nb][op][it]);
    if (msglen > 0)
      parall->adsend(request_indx, msgtype, proc_to, msglen,
                     tmp1 + 2 * p_i1_start[nb][op][it]);
    // parall->adsend(request_indx,msgtype,proc_to,msglen,&tmp1[2*p_i1_start[nb][op][it]]);
  }
}

/**************************************
 *                                    *
 *    d3db::c_ptranspose_ijk_end      *
 *                                    *
 **************************************/
void d3db::c_ptranspose_ijk_end(const int nb, const int op, double *a,
                                double *tmp2, const int request_indx) {
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
 *    d3db::c_ptranspose1_jk    *
 *                              *
 ********************************/
void d3db::c_ptranspose1_jk(const int nb, double *a, double *tmp1,
                            double *tmp2) {
  int it, proc_from, proc_to;
  int msglen;

  int n1 = p_i1_start[nb][0][np];
  int n2 = p_i2_start[nb][0][np];

  parall->astart(1, np);

  c_aindexcopy(n1, p_iq_to_i1[nb][0], a, tmp1);

  /* it = 0, transpose data on same thread */
  msglen = 2 * (p_i2_start[nb][0][1] - p_i2_start[nb][0][0]);
  // int one=1;
  // DCOPY_PWDFT(msglen,&(tmp1[2*p_i1_start[nb][0][0]]),one,&(tmp2[2*p_i2_start[nb][0][0]]),one);
  std::memcpy(tmp2 + 2 * p_i2_start[nb][0][0], tmp1 + 2 * p_i1_start[nb][0][0],
              msglen * sizeof(double));

  /* receive packed array data */
  for (it = 1; it < np; ++it) {
    /* synchronous receive of tmp */
    proc_from = (taskid - it + np) % np;
    msglen = 2 * (p_i2_start[nb][0][it + 1] - p_i2_start[nb][0][it]);
    if (msglen > 0)
      parall->adreceive(1, 1, proc_from, msglen,
                        &tmp2[2 * p_i2_start[nb][0][it]]);
  }
  for (it = 1; it < np; ++it) {
    proc_to = (taskid + it) % np;
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
 *    d3db::c_ptranspose2_jk    *
 *                              *
 ********************************/
void d3db::c_ptranspose2_jk(const int nb, double *a, double *tmp1,
                            double *tmp2) {
  int it, proc_from, proc_to;
  int msglen;

  int n1 = p_j1_start[nb][0][np];
  int n2 = p_j2_start[nb][0][np];

  parall->astart(1, np);

  c_aindexcopy(n1, p_jq_to_i1[nb][0], a, tmp1);

  /* it = 0, transpose data on same thread */
  msglen = 2 * (p_j2_start[nb][0][1] - p_j2_start[nb][0][0]);
  // int one=1;
  // DCOPY_PWDFT(msglen,&(tmp1[2*p_j1_start[nb][0][0]]),one,&(tmp2[2*p_j2_start[nb][0][0]]),one);
  std::memcpy(tmp2 + 2 * p_j2_start[nb][0][0], tmp1 + 2 * p_j1_start[nb][0][0],
              msglen * sizeof(double));

  /* receive packed array data */
  for (it = 1; it < np; ++it) {
    /* synchronous receive of tmp */
    proc_from = (taskid - it + np) % np;
    msglen = 2 * (p_j2_start[nb][0][it + 1] - p_j2_start[nb][0][it]);
    if (msglen > 0)
      parall->adreceive(1, 1, proc_from, msglen,
                        &tmp2[2 * p_j2_start[nb][0][it]]);
  }
  for (it = 1; it < np; ++it) {
    proc_to = (taskid + it) % np;
    msglen = 2 * (p_j1_start[nb][0][it + 1] - p_j1_start[nb][0][it]);
    if (msglen > 0)
      parall->dsend(1, 1, proc_to, msglen, &tmp1[2 * p_j1_start[nb][0][it]]);
  }
  parall->aend(1);

  c_bindexcopy(n2, p_jq_to_i2[nb][0], tmp2, a);
  c_bindexzero(nfft3d - n2, p_jz_to_i2[nb][0], a);
}

/********************************
 *                              *
 *    d3db::c_transpose_jk      *
 *                              *
 ********************************/
void d3db::c_transpose_jk(double *a, double *tmp1, double *tmp2) 
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
 *    d3db::t_transpose_jk      *
 *                              *
 ********************************/
void d3db::t_transpose_jk(double *a, double *tmp1, double *tmp2) {
  int it, proc_from, proc_to;
  int msglen;

  parall->astart(1, np);

  t_bindexcopy(nfft3d,iq_to_i1[0],a,tmp1);

  /* it = 0, transpose data on same thread */
  msglen = (i2_start[0][1] - i2_start[0][0]);
  // int one=1;
  // DCOPY_PWDFT(msglen,&(tmp1[i1_start[0][0]]),one,&(tmp2[i2_start[0][0]]),one);
  std::memcpy(tmp2 + i2_start[0][0], tmp1 + i1_start[0][0],
              msglen * sizeof(double));

  /* receive packed array data */
  for (it = 1; it < np; ++it) {
    /* synchronous receive of tmp */
    proc_from = (taskid - it + np) % np;
    msglen = (i2_start[0][it + 1] - i2_start[0][it]);
    if (msglen > 0)
      parall->adreceive(1, 1, proc_from, msglen, &tmp2[i2_start[0][it]]);
  }
  for (it = 1; it < np; ++it) {
    proc_to = (taskid + it) % np;
    msglen = (i1_start[0][it + 1] - i1_start[0][it]);
    if (msglen > 0)
      parall->dsend(1, 1, proc_to, msglen, &tmp1[i1_start[0][it]]);
  }
  parall->aend(1);
  t_aindexcopy(nfft3d, iq_to_i2[0], tmp2, a);
}

/********************************
 *                              *
 *    d3db::c_ptranspose_ijk    *
 *                              *
 ********************************/
void d3db::c_ptranspose_ijk(const int nb, const int op, double *a, double *tmp1,
                            double *tmp2) {
  int nnfft3d, it, proc_from, proc_to;
  int msglen;

  int n1 = p_i1_start[nb][op][np];
  int n2 = p_i2_start[nb][op][np];
  int n3 = p_iz_to_i2_count[nb][op];

  parall->astart(1, np);

  /* pack a array */
  c_aindexcopy(n1, p_iq_to_i1[nb][op], a, tmp1);

  /* it = 0, transpose data on same thread */
  msglen = 2 * (p_i2_start[nb][op][1] - p_i2_start[nb][op][0]);
  // int one=1;
  // DCOPY_PWDFT(msglen,&(tmp1[2*p_i1_start[nb][op][0]]),one,&(tmp2[2*p_i2_start[nb][op][0]]),one);
  // std::memcpy(&(tmp2[2*p_i2_start[nb][op][0]]),&(tmp1[2*p_i1_start[nb][op][0]]),msglen*sizeof(double));
  std::memcpy(tmp2 + 2 * p_i2_start[nb][op][0],
              tmp1 + 2 * p_i1_start[nb][op][0], msglen * sizeof(double));

  /* receive packed array data */
  for (it = 1; it < np; ++it) {
    /* synchronous receive of tmp */
    proc_from = (taskid - it + np) % np;
    msglen = 2 * (p_i2_start[nb][op][it + 1] - p_i2_start[nb][op][it]);
    if (msglen > 0)
      parall->adreceive(1, 1, proc_from, msglen,
                        &tmp2[2 * p_i2_start[nb][op][it]]);
  }
  for (it = 1; it < np; ++it) {
    proc_to = (taskid + it) % np;
    msglen = 2 * (p_i1_start[nb][op][it + 1] - p_i1_start[nb][op][it]);
    if (msglen > 0)
      parall->dsend(1, 1, proc_to, msglen, &tmp1[2 * p_i1_start[nb][op][it]]);
  }

  /* wait for completion of mp_send, also do a sync */
  parall->aend(1);

  /* unpack a array */
  c_bindexcopy(n2, p_iq_to_i2[nb][op], tmp2, a);
  c_bindexzero(n3, p_iz_to_i2[nb][op], a);
}

/********************************
 *                              *
 *    d3db::c_transpose_ijk     *
 *                              *
 ********************************/
void d3db::c_transpose_ijk(const int op, double *a, double *tmp1,
                           double *tmp2) {
  int nnfft3d, it, proc_from, proc_to;
  int msglen;

  parall->astart(1, np);

  /* pack a array */
  if ((op == 0) || (op == 4))
    nnfft3d = (nx / 2 + 1) * nq1;
  if ((op == 1) || (op == 3))
    nnfft3d = (ny)*nq2;
  if ((op == 2) || (op == 5))
    nnfft3d = (nz)*nq3;
  c_bindexcopy(nnfft3d, iq_to_i1[op], a, tmp1);

  /* it = 0, transpose data on same thread */
  msglen = 2 * (i2_start[op][1] - i2_start[op][0]);
  // int one=1;
  // DCOPY_PWDFT(msglen,&(tmp1[2*i1_start[op][0]]),one,&(tmp2[2*i2_start[op][0]]),one);
  std::memcpy(tmp2 + 2 * i2_start[op][0], tmp1 + 2 * i1_start[op][0],
              msglen * sizeof(double));

  /* receive packed array data */
  for (it = 1; it < np; ++it) {
    /* synchronous receive of tmp */
    proc_from = (taskid - it + np) % np;
    msglen = 2 * (i2_start[op][it + 1] - i2_start[op][it]);
    if (msglen > 0)
      parall->adreceive(1, 1, proc_from, msglen, &tmp2[2 * i2_start[op][it]]);
  }
  for (it = 1; it < np; ++it) {
    proc_to = (taskid + it) % np;
    msglen = 2 * (i1_start[op][it + 1] - i1_start[op][it]);
    if (msglen > 0)
      parall->dsend(1, 1, proc_to, msglen, &tmp1[2 * i1_start[op][it]]);
  }

  /* wait for completion of mp_send, also do a sync */
  parall->aend(1);

  /* unpack a array */
  if ((op == 3) || (op == 5))
    nnfft3d = (nx / 2 + 1) * nq1;
  if ((op == 0) || (op == 2))
    nnfft3d = (ny)*nq2;
  if ((op == 1) || (op == 4))
    nnfft3d = (nz)*nq3;
  c_aindexcopy(nnfft3d, iq_to_i2[op], tmp2, a);
}

/********************************
 *                              *
 *    d3db::t_transpose_ijk     *
 *                              *
 ********************************/
void d3db::t_transpose_ijk(const int op, double *a, double *tmp1,
                           double *tmp2) {
  int nnfft3d, it, proc_from, proc_to;
  int msglen;

  parall->astart(1, np);

  /* pack a array */
  if ((op == 0) || (op == 4))
    nnfft3d = (nx / 2 + 1) * nq1;
  if ((op == 1) || (op == 3))
    nnfft3d = (ny)*nq2;
  if ((op == 2) || (op == 5))
    nnfft3d = (nz)*nq3;
  t_bindexcopy(nnfft3d, iq_to_i1[op], a, tmp1);

  /* it = 0, transpose data on same thread */
  msglen = (i2_start[op][1] - i2_start[op][0]);
  // int one=1;
  // DCOPY_PWDFT(msglen,&(tmp1[i1_start[op][0]]),one,&(tmp2[i2_start[op][0]]),one);
  std::memcpy(tmp2 + i2_start[op][0], tmp1 + i1_start[op][0],
              msglen * sizeof(double));

  /* receive packed array data */
  for (it = 1; it < np; ++it) {
    /* synchronous receive of tmp */
    proc_from = (taskid - it + np) % np;
    msglen = (i2_start[op][it + 1] - i2_start[op][it]);
    if (msglen > 0)
      parall->adreceive(1, 1, proc_from, msglen, &tmp2[i2_start[op][it]]);
  }
  for (it = 1; it < np; ++it) {
    proc_to = (taskid + it) % np;
    msglen = (i1_start[op][it + 1] - i1_start[op][it]);
    if (msglen > 0)
      parall->dsend(1, 1, proc_to, msglen, &tmp1[i1_start[op][it]]);
  }

  /* wait for completion of mp_send, also do a sync */
  parall->aend(1);

  /* unpack a array */
  if ((op == 3) || (op == 5))
    nnfft3d = (nx / 2 + 1) * nq1;
  if ((op == 0) || (op == 2))
    nnfft3d = (ny)*nq2;
  if ((op == 1) || (op == 4))
    nnfft3d = (nz)*nq3;
  t_aindexcopy(nnfft3d, iq_to_i2[op], tmp2, a);
}

/********************************
 *                              *
 *    d3db::timereverse_size    *
 *                              *
 ********************************/
int d3db::timereverse_size() {
  int it, indx1, indx2, i2, i3, j2, j3, k2, k3;
  int proc_to, proc_from, phere, pto, nyh, nzh;
  int sz;

  nzh = nz / 2;
  nyh = ny / 2;
  indx1 = 0;
  indx2 = 0;
  for (it = 0; it < np; ++it) {
    proc_to = (taskid + it) % np;
    proc_from = (taskid - it + np) % np;

    /**** K=(0,0,k3)  ****/
    for (k3 = 1; k3 < nzh; ++k3) {
      i3 = k3;
      j3 = -k3;
      if (i3 < 0)
        i3 = i3 + nz;
      if (j3 < 0)
        j3 = j3 + nz;

      phere = ijktop(0, 0, i3);
      pto = ijktop(0, 0, j3);

      /* packing scheme */
      if ((phere == taskid) && (pto == proc_to))
        ++indx1;

      /* unpacking scheme */
      if ((pto == taskid) && (phere == proc_from))
        ++indx2;
    }

    /**** K=(0,k2,k3)  ****/
    for (k3 = (-nzh + 1); k3 < nzh; ++k3)
      for (k2 = 1; k2 < nyh; ++k2) {
        i2 = k2;
        i3 = k3;
        j2 = -k2;
        j3 = -k3;
        if (i2 < 0)
          i2 = i2 + ny;
        if (i3 < 0)
          i3 = i3 + nz;
        if (j2 < 0)
          j2 = j2 + ny;
        if (j3 < 0)
          j3 = j3 + nz;

        phere = ijktop(0, i2, i3);
        pto = ijktop(0, j2, j3);

        /* packing scheme */
        if ((phere == taskid) && (pto == proc_to))
          ++indx1;

        /* unpacking scheme */
        if ((pto == taskid) && (phere == proc_from))
          ++indx2;
      }
  }

  sz = indx1;
  if (sz < indx2)
    sz = indx2;

  return sz;
}

/********************************
 *                              *
 *      d3db::t_timereverse     *
 *                              *
 ********************************/
void d3db::t_timereverse(double *a, double *tmp1, double *tmp2) {
  int nnfft3d, indx, it, proc_from, proc_to;
  int msglen;

  parall->astart(1, np);

  indx = t_i1_start[0];
  nnfft3d = (t_i1_start[np] - t_i1_start[0] + 0);
  t_aindexcopy(nnfft3d, t_iq_to_i1 + indx, a, tmp1 + indx);

  // !!! DEBUG WARNING possible issue with memcpy
  /* it = 0, transpose data on same thread */
  msglen = (t_i2_start[1] - t_i2_start[0]);
  // int one=1;
  // DCOPY_PWDFT(msglen,&(tmp1[2*t_i1_start[0]]),one,&(tmp2[2*t_i2_start[0]]),one);
  std::memcpy(tmp2 + t_i2_start[0], tmp1 + t_i1_start[0],
              msglen * sizeof(double));

  /* receive packed array data */
  for (it = 1; it < np; ++it) {
    /* asynchronous receive of tmp */
    proc_from = (taskid - it + np) % np;
    msglen = (t_i2_start[it + 1] - t_i2_start[it]);
    if (msglen > 0)
      parall->adreceive(1, 1, proc_from, msglen, &tmp2[2 * t_i2_start[it]]);
  }
  for (it = 1; it < np; ++it) {
    proc_to = (taskid + it) % np;
    msglen = (t_i1_start[it + 1] - t_i1_start[it]);
    if (msglen > 0)
      parall->dsend(1, 1, proc_to, msglen, &tmp1[2 * t_i1_start[it]]);
  }
  parall->aend(1);

  indx = t_i2_start[0];
  nnfft3d = (t_i2_start[np] - t_i2_start[0] + 0);
  t_bindexcopy(nnfft3d, &t_iq_to_i2[indx], &tmp2[indx], a);
}

/********************************
 *                              *
 *      d3db::c_timereverse     *
 *                              *
 ********************************/
void d3db::c_timereverse(double *a, double *tmp1, double *tmp2) {
  int nnfft3d, indx, it, proc_from, proc_to;
  int msglen;

  parall->astart(1, np);

  indx = t_i1_start[0];
  nnfft3d = (t_i1_start[np] - t_i1_start[0] + 0);
  c_aindexcopy(nnfft3d, t_iq_to_i1 + indx, a, tmp1 + indx);

  /* it = 0, transpose data on same thread */
  msglen = 2 * (t_i2_start[1] - t_i2_start[0]);
  // int one=1;
  // DCOPY_PWDFT(msglen,&(tmp1[2*t_i1_start[0]]),one,&(tmp2[2*t_i2_start[0]]),one);
  std::memcpy(tmp2 + 2 * t_i2_start[0], tmp1 + 2 * t_i1_start[0],
              msglen * sizeof(double));

  /* receive packed array data */
  for (it = 1; it < np; ++it) {
    /* synchronous receive of tmp */
    proc_from = (taskid - it + np) % np;
    msglen = 2 * (t_i2_start[it + 1] - t_i2_start[it]);
    if (msglen > 0)
      parall->adreceive(1, 1, proc_from, msglen, &tmp2[2 * t_i2_start[it]]);
  }
  for (it = 1; it < np; ++it) {
    proc_to = (taskid + it) % np;
    msglen = 2 * (t_i1_start[it + 1] - t_i1_start[it]);
    if (msglen > 0)
      parall->dsend(1, 1, proc_to, msglen, &tmp1[2 * t_i1_start[it]]);
  }
  parall->aend(1);

  indx = t_i2_start[0];
  nnfft3d = (t_i2_start[np] - t_i2_start[0] + 0);
  c_bindexcopy_conjg(nnfft3d, t_iq_to_i2 + indx, tmp2 + indx, a);
}

/**************************************
 *                                    *
 *      d3db::c_timereverse_start     *
 *                                    *
 **************************************/
void d3db::c_timereverse_start(double *a, double *tmp1, double *tmp2,
                               const int request_indx, const int msgtype) {
  int nnfft3d, indx, it, proc_from, proc_to;
  int msglen;

  indx = t_i1_start[0];
  nnfft3d = (t_i1_start[np] - t_i1_start[0] + 0);
  c_aindexcopy(nnfft3d, t_iq_to_i1 + indx, a, tmp1 + indx);

  /* it = 0, transpose data on same thread */
  msglen = 2 * (t_i2_start[1] - t_i2_start[0]);
  // int one=1;
  // DCOPY_PWDFT(msglen,&(tmp1[2*t_i1_start[0]]),one,&(tmp2[2*t_i2_start[0]]),one);
  std::memcpy(tmp2 + 2 * t_i2_start[0], tmp1 + 2 * t_i1_start[0],
              msglen * sizeof(double));

  /* receive packed array data */
  for (it = 1; it < np; ++it) {
    /* synchronous receive of tmp */
    proc_from = (taskid - it + np) % np;
    msglen = 2 * (t_i2_start[it + 1] - t_i2_start[it]);
    if (msglen > 0)
      parall->adreceive(request_indx, msgtype, proc_from, msglen,
                        &tmp2[2 * t_i2_start[it]]);
  }
  for (it = 1; it < np; ++it) {
    proc_to = (taskid + it) % np;
    msglen = 2 * (t_i1_start[it + 1] - t_i1_start[it]);
    if (msglen > 0)
      parall->adsend(request_indx, msgtype, proc_to, msglen,
                     &tmp1[2 * t_i1_start[it]]);
  }
}

/**************************************
 *                                    *
 *      d3db::c_timereverse_end       *
 *                                    *
 **************************************/
void d3db::c_timereverse_end(double *a, double *tmp1, double *tmp2,
                             const int request_indx) {
  int indx = t_i2_start[0];
  int nnfft3d = (t_i2_start[np] - t_i2_start[0] + 0);

  parall->awaitall(request_indx);
  c_bindexcopy_conjg(nnfft3d, t_iq_to_i2 + indx, tmp2 + indx, a);
}

/********************************
 *                              *
 *         d3db::c_setpw        *
 *                              *
 ********************************/
void d3db::c_setpw(const int filling[], const double *cvalue, double *a) {
  int i = filling[0];
  int j = filling[1];
  int k = filling[2];

  int indx = ijktoindex(i, j, k);
  int p = ijktop(i, j, k);
  if (p == parall->taskid_i()) {
    a[2 * indx] = cvalue[0];
    a[2 * indx + 1] = cvalue[1];
  }

  /* set the conjugate on the i==0 plane */
  if ((i == 0) && (j != 0) && (k != 0)) {
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
    indx = ijktoindex(i, jc, kc);
    p = ijktop(i, jc, kc);
    if (p == parall->taskid_i()) {
      a[2 * indx] = cvalue[0];
      a[2 * indx + 1] = -cvalue[1];
    }
  }
}

/********************************
 *                              *
 *        d3db::c_addrandom     *
 *                              *
 ********************************/

void d3db::c_addrandom(double *a) {
  double fac = 1.0 / sqrt(1.0 * nfft3d);
  for (auto i = 0; i < n2ft3d; ++i)
    a[i] += fac * (0.50 - util_random(0));
}

/********************************
 *                              *
 *        d3db::r_setrandom     *
 *                              *
 ********************************/
void d3db::r_setrandom(double *a) {
  double fac = 1.0 / (1.0 * nfft3d);
  for (auto i = 0; i < n2ft3d; ++i)
     a[i] = fac*(0.50 - util_random(0));
}

/********************************
 *                              *
 *        d3db::hr2r_expand     *
 *                              *
 ********************************/
// expands a grid that is ((nx/2+2),ny/2,nz/2) to (nx+2,ny,nz)
void d3db::hr2r_expand(const double *ah, double *a) {
  std::memset(a, 0, n2ft3d * sizeof(double));
  if (maptype == 1) 
  {
     int nxh = nx/2;
     int nqh = nq/2;
     int nyh = ny/2;
     for (auto j=0; j<nyh; ++j)
        for (auto q=0; q<nqh; ++q)
           for (auto i=0; i<nxh; ++i)
              a[i + j*(nx+2) + q*(nx+2)*ny] = ah[i + j*(nxh+2) + q*(nxh+2)*nyh];
  } 
  else 
  {
     int nxh = nx/2;
     int nq1h = nq1/4;
     for (auto q=0; q<nq1h; ++q)
        for (auto i=0; i<nxh; ++i)
           a[i + q*(nx+2)] = ah[i + q*(nxh+2)];
  }
}

/********************************
 *                              *
 *      d3db::r2hr_contract     *
 *                              *
 ********************************/
// contracts to a grid that is (nx,ny,nz) --> ((nx+2)/2,ny/2,nz/2)
void d3db::r2hr_contract(const double *a, double *ah) {
  std::memset(ah, 0, n2ft3d / 8 * sizeof(double));
  if (maptype == 1) 
  {
     int nxh = nx/2;
     int nyh = ny/2;
     int nqh = nq/2;
     for (auto q=0; q<nqh; ++q)
        for (auto j=0; j<nyh; ++j)
           for (auto i=0; i<nxh; ++i)
              ah[i + j*(nxh+2) + q*(nxh+2)*nyh] = a[i + j*(nx+2) + q*(nx+2)*ny];
  } 
  else 
  {
     int nxh = nx/2;
     int nq1h = nq1/4;
     for (auto q=0; q<nq1h; ++q)
        for (auto i=0; i<nxh; ++i)
           ah[i + q*(nxh+2)] = a[i + q*(nx+2)];
  }
}

/************************************************
 *           r_tranpose routines block          *
 ************************************************/
/******************************************
 *                                        *
 *      d3db::r_transpose_ijk_init        *
 *                                        *
 ******************************************/
void d3db::r_transpose_ijk_init() {
  if ((maptype != 1) && (!initialized_r_transpose)) {
    initialized_r_transpose = true;

    iq_to_ir1 = new (std::nothrow) int *[6]();
    iq_to_ir1[0] = new (std::nothrow) int[(nx + 2) * nqr1]();
    iq_to_ir1[1] = new (std::nothrow) int[ny * nqr2]();
    iq_to_ir1[2] = new (std::nothrow) int[nz * nqr3]();

    iq_to_ir1[3] = new (std::nothrow) int[ny * nqr2]();
    iq_to_ir1[4] = new (std::nothrow) int[(nx + 2) * nqr1]();
    iq_to_ir1[5] = new (std::nothrow) int[nz * nqr3]();

    iq_to_ir2 = new (std::nothrow) int *[6]();
    iq_to_ir2[0] = new (std::nothrow) int[ny * nqr2]();
    iq_to_ir2[1] = new (std::nothrow) int[nz * nqr3]();
    iq_to_ir2[2] = new (std::nothrow) int[ny * nqr2]();

    iq_to_ir2[3] = new (std::nothrow) int[(nx + 2) * nq1]();
    iq_to_ir2[4] = new (std::nothrow) int[nz * nqr3]();
    iq_to_ir2[5] = new (std::nothrow) int[(nx + 2) * nqr1]();

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
    for (auto it = 0; it < np; ++it) {
      auto proc_to = (taskid + it) % np;
      auto proc_from = (taskid - it + np) % np;
      ir1_start[0][it] = index1;
      ir2_start[0][it] = index2;
      for (auto k = 0; k < nz; ++k)
        for (auto j = 0; j < ny; ++j)
          for (auto i = 0; i < (nx + 2); ++i) {
            // int p    = ijktop2(i,j,k);
            auto phere = ijkrtop2(i, j, k);
            auto pto = ijkrtop1(i, j, k);

            /* packing scheme */
            if ((phere == taskid) && (pto == proc_to)) {
              // int indx = ijktoindex2(i,j,k);
              iq_to_ir1[0][ijkrtoindex2t(i, j, k)] = index1;
              ++index1;
            }
            /* unpacking scheme */
            if ((pto == taskid) && (phere == proc_from)) {
              iq_to_ir2[0][ijkrtoindex1(i, j, k)] = index2;
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
    for (auto it = 0; it < np; ++it) {
      auto proc_to = (taskid + it) % np;
      auto proc_from = (taskid - it + np) % np;
      ir1_start[1][it] = index1;
      ir2_start[1][it] = index2;
      for (auto k = 0; k < nz; ++k)
        for (auto j = 0; j < ny; ++j)
          for (auto i = 0; i < (nx + 2); ++i) {
            auto phere = ijkrtop1(i, j, k);
            auto pto = ijkrtop(i, j, k);

            /* packing scheme */
            if ((phere == taskid) && (pto == proc_to)) {
              iq_to_ir1[1][ijkrtoindex1(i, j, k)] = index1;
              ++index1;
            }
            /* unpacking scheme */
            if ((pto == taskid) && (phere == proc_from)) {
              iq_to_ir2[1][ijkrtoindex(i, j, k)] = index2;
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
    for (auto it = 0; it < np; ++it) {
      auto proc_to = (taskid + it) % np;
      auto proc_from = (taskid - it + np) % np;
      ir1_start[2][it] = index1;
      ir2_start[2][it] = index2;
      for (auto k = 0; k < nz; ++k)
        for (auto j = 0; j < ny; ++j)
          for (auto i = 0; i < (nx + 2); ++i) {
            auto phere = ijkrtop(i, j, k);
            auto pto = ijkrtop1(i, j, k);

            /* packing scheme */
            if ((phere == taskid) && (pto == proc_to)) {
              iq_to_ir1[2][ijkrtoindex(i, j, k)] = index1;
              ++index1;
            }
            /* unpacking scheme */
            if ((pto == taskid) && (phere == proc_from)) {
              iq_to_ir2[2][ijkrtoindex1(i, j, k)] = index2;
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
    for (auto it = 0; it < np; ++it) {
      auto proc_to = (taskid + it) % np;
      auto proc_from = (taskid - it + np) % np;
      ir1_start[3][it] = index1;
      ir2_start[3][it] = index2;
      for (auto k = 0; k < nz; ++k)
        for (auto j = 0; j < ny; ++j)
          for (auto i = 0; i < (nx + 2); ++i) {
            auto phere = ijkrtop1(i, j, k);
            auto pto = ijkrtop2(i, j, k);

            /* packing scheme */
            if ((phere == taskid) && (pto == proc_to)) {
              iq_to_ir1[3][ijkrtoindex1(i, j, k)] = index1;
              ++index1;
            }
            /* unpacking scheme */
            if ((pto == taskid) && (phere == proc_from)) {
              iq_to_ir2[3][ijkrtoindex2t(i, j, k)] = index2;
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
    for (auto it = 0; it < np; ++it) {
      auto proc_to = (taskid + it) % np;
      auto proc_from = (taskid - it + np) % np;
      ir1_start[4][it] = index1;
      ir2_start[4][it] = index2;
      for (auto k = 0; k < nz; ++k)
        for (auto j = 0; j < ny; ++j)
          for (auto i = 0; i < (nx + 2); ++i) {
            auto phere = ijkrtop2(i, j, k);
            auto pto = ijkrtop(i, j, k);

            /* packing scheme */
            if ((phere == taskid) && (pto == proc_to)) {
              iq_to_ir1[4][ijkrtoindex2t(i, j, k)] = index1;
              ++index1;
            }
            /* unpacking scheme */
            if ((pto == taskid) && (phere == proc_from)) {
              iq_to_ir2[4][ijkrtoindex(i, j, k)] = index2;
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
    for (auto it = 0; it < np; ++it) {
      auto proc_to = (taskid + it) % np;
      auto proc_from = (taskid - it + np) % np;
      ir1_start[5][it] = index1;
      ir2_start[5][it] = index2;
      for (auto k = 0; k < nz; ++k)
        for (auto j = 0; j < ny; ++j)
          for (auto i = 0; i < (nx + 2); ++i) {
            auto phere = ijkrtop(i, j, k);
            auto pto = ijkrtop2(i, j, k);

            /* packing scheme */
            if ((phere == taskid) && (pto == proc_to)) {
              iq_to_ir1[5][ijkrtoindex(i, j, k)] = index1;
              ++index1;
            }
            /* unpacking scheme */
            if ((pto == taskid) && (phere == proc_from)) {
              iq_to_ir2[5][ijkrtoindex2t(i, j, k)] = index2;
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
 *      d3db::r_transpose_ijk_end         *
 *                                        *
 ******************************************/
void d3db::r_transpose_ijk_end() {
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
 *    d3db::r_transpose_jk      *
 *                              *
 ********************************/
void d3db::r_transpose_jk(double *a, double *tmp1, double *tmp2) 
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
 *    d3db::r_transpose_ijk     *
 *                              *
 ********************************/
void d3db::r_transpose_ijk(const int op, double *a, double *tmp1,
                           double *tmp2) {
  if (!initialized_r_transpose)
    this->r_transpose_ijk_init();

  int nnfft3d, it, proc_from, proc_to;
  int msglen;

  parall->astart(1, np);

  /* pack a array */
  if ((op == 0) || (op == 4))
    nnfft3d = (nx + 2) * nqr1;
  if ((op == 1) || (op == 3))
    nnfft3d = (ny)*nqr2;
  if ((op == 2) || (op == 5))
    nnfft3d = (nz)*nqr3;
  t_bindexcopy(nnfft3d, iq_to_ir1[op], a, tmp1);

  /* it = 0, transpose data on same thread */
  msglen = (ir2_start[op][1] - ir2_start[op][0]);
  // int one=1;
  // DCOPY_PWDFT(msglen,&(tmp1[i1_start[op][0]]),one,&(tmp2[i2_start[op][0]]),one);
  std::memcpy(tmp2 + ir2_start[op][0], tmp1 + ir1_start[op][0],
              msglen * sizeof(double));

  /* receive packed array data */
  for (it = 1; it < np; ++it) {
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
  if ((op == 3) || (op == 5))
    nnfft3d = (nx + 2) * nqr1;
  if ((op == 0) || (op == 2))
    nnfft3d = (ny)*nqr2;
  if ((op == 1) || (op == 4))
    nnfft3d = (nz)*nqr3;
  t_aindexcopy(nnfft3d, iq_to_ir2[op], tmp2, a);
}

/**************************************
 *                                    *
 *     d3db::rrrr_periodic_gradient   *
 *                                    *
 **************************************/
// computes the periodic gradient on a (n1,n2,n3) grid.
#define one_over_60 1.66666666666667e-2
void d3db::rrrr_periodic_gradient(const double *rho, double *drho1, double *drho2, double *drho3) 
{
   int nx2 = nx + 2;
 
   if (maptype == 1) 
   {
      double *a = this->r_alloc();
      double *tmp2 = new (std::nothrow) double[2*nfft3d]();
      double *tmp3 = new (std::nothrow) double[2*nfft3d]();
      std::memcpy(a,rho,n2ft3d*sizeof(double));
      this->r_transpose_jk(a,tmp2,tmp3);
     
      // drho1 gradient
      for (auto q=0; q<nq; ++q) 
      {
         int im3,im2,im1,ip1,ip2,ip3;
         const double *f = rho + q * (nx2)*nz;
         double *df1 = drho1 + q*(nx2)*nz;
         for (auto i=0; i<nx; ++i) 
         {
            im3 = i-3; if (im3<0) im3 += nx;
            im2 = i-2; if (im2<0) im2 += nx;
            im1 = i-1; if (im1<0) im1 += nx;
            ip1 = i+1; if (ip1>=nx) ip1 -= nx;
            ip2 = i+2; if (ip2>=nx) ip2 -= nx;
            ip3 = i+3; if (ip3>=nx) ip3 -= nx;
           
            for (auto j = 0; j < ny; ++j)
               df1[i+ j*nx2] = one_over_60*(  -1.0*f[im3+j*nx2] +  9.0*f[im2+j*nx2] 
                                            - 45.0*f[im1+j*nx2] + 45.0*f[ip1+j*nx2] 
                                            -  9.0*f[ip2+j*nx2] +  1.0*f[ip3+j*nx2]);
         }
      }

      // drho2 gradient
      for (auto q=0; q<nq; ++q) 
      {
         int jm3,jm2,jm1,jp1,jp2,jp3;
         const double *f = a + q*(nx2)*ny;
         double *df2 = drho2 + q*(nx2)*ny;
         for (auto j=0; j<ny; ++j) 
         {
            jm3 = j-3; if (jm3<0) jm3 += ny;
            jm2 = j-2; if (jm2<0) jm2 += ny;
            jm1 = j-1; if (jm1<0) jm1 += ny;
            jp1 = j+1; if (jp1>=ny) jp1 -= ny;
            jp2 = j+2; if (jp2>=ny) jp2 -= ny;
            jp3 = j+3; if (jp3>=ny) jp3 -= ny;
            for (auto i = 0; i<nx; ++i)
               df2[i+j*nx2] = one_over_60*( -1.0*f[i+jm3*nx2] +  9.0*f[i+jm2*nx2] 
                                          - 45.0*f[i+jm1*nx2] + 45.0*f[i+jp1*nx2] 
                                          -  9.0*f[i+jp2*nx2] +  1.0*f[i+jp3*nx2]);
         }
      }
      r_transpose_jk(drho2,tmp2,tmp3);
     
      // drho3 gradient
      for (auto q = 0; q < nq; ++q) 
      {
         int km3,km2,km1,kp1,kp2,kp3;
         const double *f = rho + q*(nx2)*nz;
         double *df3 = drho3 + q*(nx2)*nz;
         for (auto k=0; k<nz; ++k) 
         {
            km3 = k-3; if (km3<0) km3 += nz;
            km2 = k-2; if (km2<0) km2 += nz;
            km1 = k-1; if (km1<0) km1 += nz;
            kp1 = k+1; if (kp1>=nz) kp1 -= nz;
            kp2 = k+2; if (kp2>=nz) kp2 -= nz;
            kp3 = k+3; if (kp3>=nz) kp3 -= nz;
            for (auto i=0; i<nx; ++i)
               df3[i+k*nx2] = one_over_60 * ( -1.0*f[i+km3*nx2] +  9.0*f[i+km2*nx2]
                                            - 45.0*f[i+km1*nx2] + 45.0*f[i+kp1*nx2]
                                            -  9.0*f[i+kp2*nx2] +  1.0*f[i+kp3*nx2]);
         }
      }
      this->r_dealloc(a);
      delete[] tmp3;
      delete[] tmp2;
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
         const double *f = rho + q*(nx+2);
         double *df1 = drho1 + q*(nx2);
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
   this->r_zero_ends(drho1);
   this->r_zero_ends(drho2);
   this->r_zero_ends(drho3);
}


/**************************************
 *                                    *
 *    d3db::rrrr_periodic_laplacian   *
 *                                    *
 **************************************/
// computes the periodic gradient on a (nx,ny,nz) grid.
#define one_over_180 5.55555555555556e-3
void d3db::rrrr_periodic_laplacian(const double *rho, double *grxx,
                                   double *gryy, double *grzz) {
  int nx2 = nx + 2;

  if (maptype == 1) {
    double *a = this->r_alloc();
    double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
    double *tmp3 = new (std::nothrow) double[2 * nfft3d]();
    std::memcpy(a, rho, n2ft3d * sizeof(double));
    this->r_transpose_jk(a, tmp2, tmp3);

    // xx gradient
    for (auto q = 0; q < nq; ++q) {
      int im3, im2, im1, ip1, ip2, ip3;
      const double *f = rho + q * (nx + 2) * nz;
      double *ddf = grxx + q * (nx + 2) * nz;
      for (auto i = 0; i < nx; ++i) {
        im3 = i - 3; if (im3 < 0) im3 += nx;
        im2 = i - 2; if (im2 < 0) im2 += nx;
        im1 = i - 1; if (im1 < 0) im1 += nx;
        ip1 = i + 1; if (ip1 >= nx) ip1 -= nx;
        ip2 = i + 2; if (ip2 >= nx) ip2 -= nx;
        ip3 = i + 3; if (ip3 >= nx) ip3 -= nx;

        for (auto j = 0; j < ny; ++j)
          ddf[i + j * nx2] =
              one_over_180 *
              (2.0 * f[im3 + j * nx2] - 27.0 * f[im2 + j * nx2] +
               270.0 * f[im1 + j * nx2] - 490.0 * f[i + j * nx2] +
               270.0 * f[ip1 + j * nx2] - 27.0 * f[ip2 + j * nx2] +
               2.0 * f[ip3 + j * nx2]);
      }
    }

    // yy gradient
    for (auto q = 0; q < nq; ++q) {
      int jm3, jm2, jm1, jp1, jp2, jp3;
      const double *f = a + q * (nx + 2) * ny;
      double *ddf = gryy + q * (nx + 2) * ny;
      for (auto j = 0; j < ny; ++j) {
        jm3 = j - 3; if (jm3 < 0) jm3 += ny;
        jm2 = j - 2; if (jm2 < 0) jm2 += ny;
        jm1 = j - 1; if (jm1 < 0) jm1 += ny;
        jp1 = j + 1; if (jp1 >= ny) jp1 -= ny;
        jp2 = j + 2; if (jp2 >= ny) jp2 -= ny;
        jp3 = j + 3; if (jp3 >= ny) jp3 -= ny;
        for (auto i = 0; i < nx; ++i)
          ddf[i + j * nx2] =
              one_over_180 *
              (2.0 * f[i + jm3 * nx2] - 27.0 * f[i + jm2 * nx2] +
               270.0 * f[i + jm1 * nx2] - 490.0 * f[i + j * nx2] +
               270.0 * f[i + jp1 * nx2] - 27.0 * f[i + jp2 * nx2] +
               2.0 * f[i + jp3 * nx2]);
      }
    }
    this->r_transpose_jk(gryy, tmp2, tmp3);

    // zz gradient
    for (auto q = 0; q < nq; ++q) {
      int km3, km2, km1, kp1, kp2, kp3;
      const double *f = rho + q * (nx + 2) * nz;
      double *ddf = grzz + q * (nx + 2) * nz;
      for (auto k = 0; k < nz; ++k) {
        km3 = k - 3; if (km3 < 0) km3 += nz;
        km2 = k - 2; if (km2 < 0) km2 += nz;
        km1 = k - 1; if (km1 < 0) km1 += nz;
        kp1 = k + 1; if (kp1 >= nz) kp1 -= nz;
        kp2 = k + 2; if (kp2 >= nz) kp2 -= nz;
        kp3 = k + 3; if (kp3 >= nz) kp3 -= nz;
        for (auto i = 0; i < nx; ++i)
          ddf[i + k * nx2] =
              one_over_180 *
              (2.0 * f[i + km3 * nx2] - 27.0 * f[i + km2 * nx2] +
               270.0 * f[i + km1 * nx2] - 490.0 * f[i + k * nx2] +
               270.0 * f[i + kp1 * nx2] - 27.0 * f[i + kp2 * nx2] +
               2.0 * f[i + kp3 * nx2]);
      }
    }

    this->r_dealloc(a);
    delete[] tmp3;
    delete[] tmp2;
  } else {
    double *a = this->r_alloc();
    double *tmp2 = new (std::nothrow) double[nrft3d]();
    double *tmp3 = new (std::nothrow) double[nrft3d]();

    // xx gradient
    for (auto q = 0; q < nqr1; ++q) {
      int im3, im2, im1, ip1, ip2, ip3;
      const double *f = rho + q * (nx + 2);
      double *ddf = grxx + q * (nx + 2);
      for (auto i = 0; i < nx; ++i) {
        im3 = i - 3; if (im3 < 0) im3 += nx;
        im2 = i - 2; if (im2 < 0) im2 += nx;
        im1 = i - 1; if (im1 < 0) im1 += nx;
        ip1 = i + 1; if (ip1 >= nx) ip1 -= nx;
        ip2 = i + 2; if (ip2 >= nx) ip2 -= nx;
        ip3 = i + 3; if (ip3 >= nx) ip3 -= nx;
        ddf[i] = one_over_180 *
                 (2.0 * f[im3] - 27.0 * f[im2] + 270.0 * f[im1] - 490.0 * f[i] +
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
  this->r_zero_ends(grxx);
  this->r_zero_ends(gryy);
  this->r_zero_ends(grzz);
}

/************************************************
 *           r_tranpose routines block          *
 ************************************************/


/************************************************
 *       Gaussian filters routines block        *
 ************************************************/

 #define PERIODIC_GAUSSIAN_FILTER(nx, nfilter, coeff, a, b) \
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

#define GAUSSIAN_FILTER(nx, nfilter, coeff, a, b) \
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


#define SMOOTH_PERIODIC_2D(nx, ny, nx2, nfilter, coeff, a, b) \
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
                sum += (coeff)[r]*(a)[ii+j*(nx2)]; \
            } \
            (b)[i+j*(nx2)] = sum; \
        } \
    } while(0)

#define SMOOTH_PERIODIC_2D_TRANS(nx, ny, nx2, nfilter, coeff, a, b) \
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
                sum += (coeff)[r]*(a)[i+jj*(nx2)]; \
            } \
            (b)[i+j*(nx2)] = sum; \
        } \
    } while(0)


/**********************************************
 *                                            *
 *     d3db::rr_periodic_gaussian_filter      *
 *                                            *
 **********************************************/
// computes the 3D Gaussian filter function
void d3db::rr_periodic_gaussian_filter(const double sigma, const double *a, double *b)
{
   int nx2 = nx+2;

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
         const double *f  = a + q*(nx2)*nz;
         double       *f1 = b + q*(nx2)*nz;
         SMOOTH_PERIODIC_2D(nx,ny,nx2,nfilter,coeff,f,f1);
      }

      // Apply the filter along the y axis
      for (auto q=0; q<nq; ++q)
      {
         const double *f  = b + q*(nx2)*nz;
         double       *f1 = c + q*(nx2)*nz;
         SMOOTH_PERIODIC_2D_TRANS(nx,ny,nx2,nfilter,coeff,f,f1);
      }

      // transpose y and z
      this->c_transpose_jk(c,tmp2,tmp3);

      // Apply the filter along the z axis
      for (auto q=0; q<nq; ++q)
      {
         const double *f  = c + q*(nx2)*ny;
         double       *f1 = b + q*(nx2)*ny;
         SMOOTH_PERIODIC_2D_TRANS(nx,nz,nx2,nfilter,coeff,f,f1);
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
         const double *f = a + q*(nx2);
         double      *f1 = b + q*(nx2);
         PERIODIC_GAUSSIAN_FILTER(nx,nfilter,coeff,f,f1);
      }

      // tranpose operation #0  (j,k,i) <-- (i,j,k)
      // rotate b(ny,nz,nx) <-- b(nx,ny,nz)
      this->r_transpose_ijk(0,b,tmp2,tmp3);

      // Apply the filter along the y axis
      for (auto q=0; q<nqr2; ++q)
      {
         const double *f = b + q*(ny);
         double      *f1 = c + q*(ny);
         PERIODIC_GAUSSIAN_FILTER(ny,nfilter,coeff,f,f1);
      }

      // tranpose operation #1 (k,i,j) <-- (j,k,i)
      // rotate c(nz,nx,ny) <-- c(ny,nz,nx)
      this->r_transpose_ijk(1,c,tmp2,tmp3);

      // Apply the filter along the z axis
      for (auto q=0; q<nqr3; ++q)
      {
         const double *f = c + q*(nz);
         double      *f1 = b + q*(nz);
         PERIODIC_GAUSSIAN_FILTER(nz,nfilter,coeff,f,f1);
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

/* PGrid.cpp
   Author - Eric Bylaska

        this class is used for defining 3d parallel maps
*/

#include "Control2.hpp"
#include "Lattice.hpp"
#include "fft.h"
//#include "gdevice.hpp"
#include "nwpw_timing.hpp"
#include "util.hpp"
#include <cmath>
#include <cstring> //memset
#include <iostream>

#include "PGrid.hpp"

#include "blas.h"
#define mytaskid 1


// Replace this with the desired alignment value (in bytes)
constexpr std::size_t Alignment = 64;  // Example alignment value

namespace pwdft {

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

// PGrid::PGrid(Parallel *inparall, Lattice *inlattice, Control2& control) :
// d3db(inparall,control.mapping(),control.ngrid(0),control.ngrid(1),control.ngrid(2))

PGrid::PGrid(Parallel *inparall, Lattice *inlattice, int mapping0, int balance0,
             int nx0, int ny0, int nz0, int pfft3_qsize0, bool staged_gpu_fft_pipeline0)
    : d3db(inparall, mapping0, nx0, ny0, nz0) 
{
   int nxh, nyh, nzh, p, q, indx, nb;
   int nwave_in[2], nwave_out[2];
   double *G1, *G2, *G3;
   double ggcut, eps, ggmax, ggmin;
   bool *zero_arow3, *zero_arow2;
   bool yzslab, zrow;
 
   lattice = inlattice;
 
   eps = 1.0e-12;

   // aligned Memory
   std::size_t aligned_size3d = (3 * nfft3d * sizeof(double) + Alignment - 1) & ~(Alignment - 1);
   Garray = new (std::nothrow) double[aligned_size3d]();
   G1 = Garray;
   G2 = Garray + nfft3d;
   G3 = Garray + 2*nfft3d;
   //G2 = (double *)&Garray[nfft3d];
   //G3 = (double *)&Garray[2 * nfft3d];
   nxh = nx/2;
   nyh = ny/2;
   nzh = nz/2;
   ggmin = 9.9e9;
   ggmax = 0.0;
   for (auto k3=(-nzh+1); k3<=nzh; ++k3)
     for (auto k2=(-nyh+1); k2<=nyh; ++k2)
       for (auto k1=0; k1<=nxh; ++k1) 
       {
          auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2);
          auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2);
          auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2);
          auto gg = gx*gx + gy*gy + gz*gz;
          if (gg > ggmax)
            ggmax = gg;
          if ((gg < ggmin) && (gg > 1.0e-6))
            ggmin = gg;
          auto i = k1;
          if (i < 0)
            i = i + nx;
          auto j = k2;
          if (j < 0)
            j = j + ny;
          auto k = k3;
          if (k < 0)
            k = k + nz;
         
          auto indx = ijktoindex(i, j, k);
          auto p = ijktop(i, j, k);
          if (p == parall->taskid_i()) {
            G1[indx] = gx;
            G2[indx] = gy;
            G3[indx] = gz;
          }
       }
   Gmax = sqrt(ggmax);
   Gmin = sqrt(ggmin);
 
 
   // aligned Memory
   std::size_t aligned_size = (nfft3d * sizeof(int) + Alignment - 1) & ~(Alignment - 1);
   masker[0] = new (std::nothrow) int[aligned_size]();
   masker[1] = new (std::nothrow) int[aligned_size]();
   packarray[0] = new (std::nothrow) int[aligned_size]();
   packarray[1] = new (std::nothrow) int[aligned_size]();

   //std::size_t aligned_size = (4 * nfft3d * sizeof(int) + Alignment - 1) & ~(Alignment - 1);
   //masker[1]    = masker[0] + nfft3d;
   //packarray[0] = masker[0] + 2*nfft3d;
   //packarray[1] = masker[0] + 3*nfft3d;
   //masker[1] = (int *)&(masker[0][nfft3d]);
   //packarray[0] = (int *)&(masker[0][nfft3d+nfft3d]);
   //packarray[1] = (int *)&(masker[0][nfft3d+nfft3d+nfft3d]);
 
   parall->Barrier();
   for (int k = 0; k < (nfft3d); ++k) 
   {
      masker[0][k] = 1;
      masker[1][k] = 1;
      packarray[0][k] = 0;
      packarray[1][k] = 0;
   }
 
   for (auto nb = 0; nb <= 1; ++nb) 
   {
      nwave[nb] = 0;
      if (nb == 0)
         ggcut = lattice->eggcut();
      else
         ggcut = lattice->wggcut();
     
      for (auto k3 = (-nzh + 1); k3 < nzh; ++k3)
        for (auto k2 = (-nyh + 1); k2 < nyh; ++k2)
          for (auto k1 = 0; k1 < nxh; ++k1) 
          {
             auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2);
             auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2);
             auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2);
             auto i = k1; if (i < 0) i += nx;
             auto j = k2; if (j < 0) j += ny;
             auto k = k3; if (k < 0) k += nz;
            
             auto indx = ijktoindex(i, j, k);
             auto p = ijktop(i, j, k);
             if (p == parall->taskid_i()) 
             {
                auto gg = gx * gx + gy * gy + gz * gz;
                gg = gg - ggcut;
                if (gg < (-eps)) 
                {
                   masker[nb][indx] = 0;
                   ++nwave[nb];
                }
             }
          }
      nwave_entire[nb] = nwave[nb];
      nwave_entire[nb] = parall->ISumAll(1, nwave_entire[nb]);
   }
 
   /*packarray[0] = new (std::nothrow) int[2 * nfft3d]();
   packarray[1] = (int *)&(packarray[0][nfft3d]);
   parall->Barrier();
   for (int k = 0; k < (nfft3d); ++k) 
   {
      packarray[0][k] = 0;
      packarray[1][k] = 0;
   }
   // packarray[1] = new (std::nothrow) int [nfft3d]();
   */
 
   for (auto nb=0; nb<=1; ++nb) 
   {
      nida[nb] = 0;
      nidb2[nb] = 0;
     
      /* k=(0,0,0)  */
      auto k1 = 0;
      auto k2 = 0;
      auto k3 = 0;
      auto indx = ijktoindex(k1, k2, k3);
      auto p = ijktop(k1, k2, k3);
      if (p == parall->taskid_i())
        if (!masker[nb][indx]) 
        {
           packarray[nb][nida[nb]] = indx;
           ++nida[nb];
        }
     
      // k=(0,0,k3) - neglect (0,0,-k3) points
      for (auto k=1; k<=(nzh-1); ++k) 
      {
         k1 = 0;
         k2 = 0;
         k3 = k;
         indx = ijktoindex(k1, k2, k3);
         p = ijktop(k1, k2, k3);
         if (p == parall->taskid_i())
           if (!masker[nb][indx]) 
           {
              packarray[nb][nida[nb] + nidb2[nb]] = indx;
              ++nidb2[nb];
           }
      }
     
      // k=(0,k2,k3) - neglect (0,-k2, -k3) points
      for (auto k = (-nzh+1); k<=(nzh-1); ++k)
         for (auto j=1; j<=(nyh-1); ++j) 
         {
            k1 = 0;
            k2 = j;
            k3 = k;
            if (k3 < 0) k3 += nz;
            indx = ijktoindex(k1, k2, k3);
            p = ijktop(k1, k2, k3);
            if (p == parall->taskid_i())
              if (!masker[nb][indx]) 
              {
                 packarray[nb][nida[nb] + nidb2[nb]] = indx;
                 ++nidb2[nb];
              }
         }
     
      // k=(k1,k2,k3)
      for (auto k=(-nzh+1); k<=(nzh-1); ++k)
         for (auto j=(-nyh+1); j<=(nyh-1); ++j)
            for (auto i=1; i<=(nxh-1); ++i) 
            {
               k1 = i;
               k2 = j;
               k3 = k;
               if (k2 < 0) k2 += ny;
               if (k3 < 0) k3 += nz;
               indx = ijktoindex(k1, k2, k3);
               p = ijktop(k1, k2, k3);
               if (p == parall->taskid_i())
                  if (!masker[nb][indx]) 
                  {
                     packarray[nb][nida[nb] + nidb2[nb]] = indx;
                     ++nidb2[nb];
                  }
            }
   }
 
   nwave_in[0] = nida[0] + nidb2[0];
   nwave_in[1] = nida[1] + nidb2[1];
 
   // if (control.balance())
   if (balance0) 
   {
      balanced = 1;
      mybalance = new Balance(parall, 2, nwave_in, nwave_out);
   } 
   else 
   {
      balanced = 0;
      nwave_out[0] = nwave_in[0];
      nwave_out[1] = nwave_in[1];
   }
   nidb[0] = nidb2[0] + (nwave_out[0] - nwave_in[0]);
   nidb[1] = nidb2[1] + (nwave_out[1] - nwave_in[1]);
 
   nwave_all[0] = nida[0] + nidb2[0];
   nwave_all[1] = nida[1] + nidb2[1];
   parall->Vector_ISumAll(1, 2, nwave_all);
 
   if (maptype == 1) 
   {
      // Calculate aligned memory size
      //std::size_t aligned_sizeb = ((nxh+1) * (4*nq+2) * sizeof(bool) + Alignment - 1)  & ~(Alignment - 1);
      //zero_row3[0] = new (std::nothrow) bool[aligned_sizeb];
      //zero_row3[1] = zero_row3[0] + (nxh+1)*nq;
      //zero_row2[0] = zero_row3[0] + 2*(nxh+1)*nq;
      //zero_row2[1] = zero_row3[0] + 3*(nxh+1)*nq;
      //zero_slab23[0] = zero_row3[0] + 4*(nxh+1)*nq;
      //zero_slab23[1] = zero_row3[0] + (nxh+1)*(4*nq+1);
 
      zero_row3[0] = new (std::nothrow) bool[((nxh + 1) * nq + Alignment - 1) & ~(Alignment - 1)];
      zero_row3[1] = new (std::nothrow) bool[((nxh + 1) * nq + Alignment - 1) & ~(Alignment - 1)];
      zero_row2[0] = new (std::nothrow) bool[((nxh + 1) * nq + Alignment - 1) & ~(Alignment - 1)];
      zero_row2[1] = new (std::nothrow) bool[((nxh + 1) * nq + Alignment - 1) & ~(Alignment - 1)];
      zero_slab23[0] = new (std::nothrow) bool[nxh + 1];
      zero_slab23[1] = new (std::nothrow) bool[nxh + 1];
     
      zero_arow3 = new bool[((nxh+1)*ny + Alignment - 1) & ~(Alignment - 1)];
      for (auto nb = 0; nb <= 1; ++nb) 
      {
         if (nb == 0)
            ggcut = lattice->eggcut();
         else
            ggcut = lattice->wggcut();
        
         /* find zero_row3 - (i,j,*) rows that are zero */
         for (auto i=0; i<((nxh+1)*nq); ++i)
            zero_row3[nb][i] = true;
 
         for (auto i=0; i<((nxh+1)*ny); ++i)
            zero_arow3[i] = true;
        
         for (auto k2 = (-nyh + 1); k2 < nyh; ++k2)
            for (auto k1 = 0; k1 < nxh; ++k1) 
            {
               auto i = k1;
               auto j = k2;
               if (i < 0) i += nx;
               if (j < 0) j += ny;
               zrow = true;
               for (auto k3 = (-nzh+1); k3<nzh; ++k3) 
               {
                  auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2);
                  auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2);
                  auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2);
                  auto gg = gx*gx + gy*gy + gz*gz;
                  gg = gg - ggcut;
                  if (gg < (-eps))
                     zrow = false;
               }
               if (!zrow) 
               {
                  zero_arow3[i + (nxh + 1) * j] = false;
                  q = ijktoq1(0, j, 0);
                  p = ijktop1(0, j, 0);
                  if (p == parall->taskid_i()) 
                  {
                     zero_row3[nb][i + (nxh + 1) * q] = false;
                  }
               }
            }
         // call D3dB_c_ptranspose_jk_init(nb,log_mb(zero_arow3(1)))
         d3db::c_ptranspose_jk_init(nb, zero_arow3);
        

         /* find zero_slab2 - (i,*,*) slabs that are zero */
         for (auto i=0; i<(nxh+1); ++i)
            zero_slab23[nb][i] = true;
        
         for (auto k1 = 0; k1<nxh; ++k1) 
         {
            auto i = k1; if (i<0) i += nx;
            yzslab = true;
            for (auto k3=(-nzh+1); k3<nzh; ++k3)
               for (auto k2=(-nyh+1); k2<nyh; ++k2) 
               {
                  auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2);
                  auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2);
                  auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2);
                  auto gg = gx*gx + gy*gy + gz*gz;
                  gg = gg - ggcut;
                  if (gg < (-eps))
                     yzslab = false;
               }
            if (!yzslab) zero_slab23[nb][i] = false;
         }

         /* initialize zero_row2 */
         for (auto i=0; i<((nxh+1)*nq); ++i)
            zero_row2[nb][i] = false;
        
         /* find zero_row2 - (i,*,k) rows that are zero after fft of (i,j,*) */
         for (auto k = 0; k<nz; ++k)
            for (auto i=0; i<(nxh+1); ++i) 
            {
               q = ijktoq(i, 0, k);
               p = ijktop(i, 0, k);
               if (p == parall->taskid_i())
                  zero_row2[nb][q] = zero_slab23[nb][i];
            }
      }
     
      delete[] zero_arow3;
   } 
   else 
   {
      zero_row3[0] = new (std::nothrow) bool[(nq3 + Alignment - 1) & ~(Alignment - 1)]();
      zero_row3[1] = new (std::nothrow) bool[(nq3 + Alignment - 1) & ~(Alignment - 1)]();
      zero_row2[0] = new (std::nothrow) bool[(nq2 + Alignment - 1) & ~(Alignment - 1)]();
      zero_row2[1] = new (std::nothrow) bool[(nq2 + Alignment - 1) & ~(Alignment - 1)]();
      zero_slab23[0] = new (std::nothrow) bool[nxh + 1]();
      zero_slab23[1] = new (std::nothrow) bool[nxh + 1]();
     
      zero_arow2 = new (std::nothrow) bool[((nxh+1)*nz + Alignment - 1) & ~(Alignment - 1)]();
      zero_arow3 = new (std::nothrow) bool[((nxh+1)*ny + Alignment - 1) & ~(Alignment - 1)]();
     
      for (auto nb = 0; nb <= 1; ++nb) {
        if (nb == 0)
          ggcut = lattice->eggcut();
        else
          ggcut = lattice->wggcut();
     
        // find zero_row3 - (i,j,*) rows that are zero
        for (auto q = 0; q < nq3; ++q)
          zero_row3[nb][q] = true;
     
        for (auto q = 0; q < (nxh + 1) * ny; ++q)
          zero_arow3[q] = true;
     
        for (auto k2 = (-nyh + 1); k2 < nyh; ++k2)
          for (auto k1 = 0; k1 < nxh; ++k1) {
            auto i = k1;
            auto j = k2;
            if (i < 0) i += nx;
            if (j < 0) j += ny;
            zrow = true;
            for (auto k3 = (-nzh + 1); k3 < nzh; ++k3) 
            {
               auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2);
               auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2);
               auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2);
               auto gg = gx*gx + gy*gy + gz*gz;
               gg = gg - ggcut;
               if (gg < (-eps))
                  zrow = false;
            }
            if (!zrow) {
              // zero_arow3[i-1+(nxh+1)*(j-1)] = 0;
              zero_arow3[i + (nxh + 1) * j] = false;
              q = ijktoq(i, j, 0);
              p = ijktop(i, j, 0);
              if (p == parall->taskid_i()) {
                zero_row3[nb][q] = false;
              }
            }
          }
     
        /* find zero_slab23 - (i,*,*) slabs that are zero */
        for (auto i = 0; i < (nxh + 1); ++i)
          zero_slab23[nb][i] = true;
     
        for (auto k1=0; k1<nxh; ++k1) 
        {
           auto i = k1; if (i < 0) i += nx;
           yzslab = true;
           for (auto k3=(-nzh+1); k3<nzh; ++k3)
              for (auto k2=(-nyh+1); k2<nyh; ++k2) 
              {
                 auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2);
                 auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2);
                 auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2);
                 auto gg = gx*gx + gy*gy + gz*gz;
                 gg = gg - ggcut;
                 if (gg < (-eps))
                   yzslab = false;
              }
           if (!yzslab)
             zero_slab23[nb][i] = false;
        }
     
        // find zero_row2 - (i,*,k) rows that are zero after fft of (i,j,*)
        for (auto k=0; k<nz; ++k)
          for (auto i=0; i<(nxh+1); ++i) 
          {
             q = ijktoq1(i, 0, k);
             p = ijktop1(i, 0, k);
             zero_arow2[i + (nxh + 1) * k] = zero_slab23[nb][i];
             if (p == parall->taskid_i())
                zero_row2[nb][q] = zero_slab23[nb][i];
          }
     
        // call
        // D3dB_c_ptranspose_ijk_init(nb,log_mb(zero_arow2(1)),log_mb(zero_arow3(1)))
        d3db::c_ptranspose_ijk_init(nb, zero_arow2, zero_arow3);
      }
      delete[] zero_arow3;
      delete[] zero_arow2;
   }
 
   //Gpack[0] = new (std::nothrow) double[3*(nida[0] + nidb[0]) + 3*(nida[1] + nidb[1])]();
   //Gpack[1] = (double *)&(Gpack[0][3 * (nida[0] + nidb[0])]);
   Gpack[0] = new (std::nothrow) double[(3*(nida[0] + nidb[0]) + Alignment - 1) & ~(Alignment - 1)]();
   Gpack[1] = new (std::nothrow) double[(3*(nida[1] + nidb[1]) + Alignment - 1) & ~(Alignment - 1)]();
 
   double *Gtmp = new (std::nothrow) double[nfft3d]();
   int one = 1;
   for (auto nb = 0; nb <= 1; ++nb) 
   {
      for (auto i = 0; i < 3; ++i) 
      {
         //DCOPY_PWDFT(nfft3d, &(Garray[i * nfft3d]), one, Gtmp, one);
         //DCOPY_PWDFT(nfft3d, Garray+i*nfft3d, one, Gtmp, one);
         std::memcpy(Gtmp,Garray+i*nfft3d,nfft3d*sizeof(double));
        
         this->t_pack(nb,Gtmp);
         //this->tt_pack_copy(nb, Gtmp, &(Gpack[nb][i * (nida[nb] + nidb[nb])]));
         this->tt_pack_copy(nb,Gtmp,Gpack[nb]+i*(nida[nb]+nidb[nb]));
      }
   }
 
   delete [] Gtmp;
 
   zplane_tmp1 = new (std::nothrow) double[ (2*zplane_size+8 + Alignment - 1) & ~(Alignment - 1)];
   zplane_tmp2 = new (std::nothrow) double[ (2*zplane_size+8 + Alignment - 1) & ~(Alignment - 1)];
 
   /* initialize r_grid */
   has_r_grid = (lattice->aperiodic());
   if (has_r_grid) 
   {
      r_grid = r_nalloc(3);
      r_nzero(3, r_grid);
      double a[9];
      for (auto i = 0; i < 3; ++i) 
      {
         a[i]   = lattice->unita1d(0 + i) / ((double)nx);
         a[3+i] = lattice->unita1d(3 + i) / ((double)ny);
         a[6+i] = lattice->unita1d(6 + i) / ((double)nz);
      }
     
      /* grid points in coordination space */
      for (auto k3 = (-nzh); k3 < nzh; ++k3)
        for (auto k2 = (-nyh); k2 < nyh; ++k2)
          for (auto k1 = (-nxh); k1 < nxh; ++k1) 
          {
             int i = k1 + nxh;
             int j = k2 + nyh;
             int k = k3 + nzh;
             int indx = ijktoindex2(i, j, k);
             int p = ijktop2(i, j, k);
            
             if (p == parall->taskid_i()) {
               r_grid[3*indx]   = a[0]*k1 + a[3]*k2 + a[6]*k3;
               r_grid[3*indx+1] = a[1]*k1 + a[4]*k2 + a[7]*k3;
               r_grid[3*indx+2] = a[2]*k1 + a[5]*k2 + a[8]*k3;
             }
          }
   }
 
   /* initialize pfft3 queues */
   staged_gpu_fft_pipeline = staged_gpu_fft_pipeline0 && d3db::mygdevice.has_gpu();

#ifdef NWPW_SYCL
   staged_gpu_fft_pipeline = false;
#endif

   //aqmax = 5;
   aqmax = pfft3_qsize0;
   if (staged_gpu_fft_pipeline)
   {
      //std::cout << "Using Staged GPU fft pipeline!!!!" << std::endl;
      aqmax += 6;
      d3db::mygdevice.batch_fft_pipeline_mem_init(aqmax,n2ft3d);
   }
 
   aqsize = 0;
   alast_index = aqmax - 1;
   aqindx = new (std::nothrow) int[aqmax]();
   aqstatus = new (std::nothrow) int[aqmax]();
   atmp = new (std::nothrow) double[2*aqmax*n2ft3d]();
 
   bqmax = pfft3_qsize0;
   if (staged_gpu_fft_pipeline) bqmax += 6;
   
   //bqmax = aqmax;
   bqsize = 0;
   blast_index = bqmax - 1;
   bqindx = new (std::nothrow) int[bqmax]();
   bqstatus = new (std::nothrow) int[bqmax]();
   btmp = new (std::nothrow) double[2*bqmax*n2ft3d]();
 
   /* initialize async buffer data for pfft */
   for (auto q=0; q<aqmax; ++q)
      parall->astart(4+q, 2*parall->np_i()+1);

}

PGrid::PGrid(Parallel *inparall, Lattice *inlattice, Control2 &control)
    : PGrid(inparall, inlattice, control.mapping(), control.balance(),
            control.ngrid(0), control.ngrid(1), control.ngrid(2),
            control.pfft3_qsize(), control.staged_gpu_fft()) {}

/********************************
 *                              *
 *       PGrid:c_unpack         *
 *                              *
 ********************************/
void PGrid::c_unpack(const int nb, double *a) 
{
   //int one = 1;
   //double *tmp1, *tmp2;
   //double *tmp = new (std::nothrow) double[2 * nfft3d];

   int nn = 2*(nida[nb]+nidb2[nb]);

   if (balanced)
      mybalance->c_unbalance(nb, a);
 
   // DCOPY_PWDFT(nn,a,one,tmp,one);
   // dcopy_(&n2ft3d,&rzero,&zero,a,&one);
   std::memcpy(d3db::d3db_tmp1, a, nn * sizeof(double));
   std::memset(a, 0, 2*nfft3d * sizeof(double));
   c_bindexcopy(nida[nb]+nidb2[nb],packarray[nb],d3db::d3db_tmp1,a);
   // c_bindexcopy(nida[nb]+nidb[nb],packarray[nb],cpack_tmp,a);
 
   // tmp1 = new (std::nothrow) double[2*zplane_size+1];
   // tmp2 = new (std::nothrow) double[2*zplane_size+1];
   c_timereverse(a, zplane_tmp1, zplane_tmp2);
   // delete [] tmp2;
   // delete [] tmp1;
   //delete[] tmp;
}

/********************************
 *                              *
 *       PGrid:c_pack           *
 *                              *
 ********************************/
void PGrid::c_pack(const int nb, double *a) 
{
   //int one = 1;
   //double *tmp = new (std::nothrow) double[2*nfft3d];
   //DCOPY_PWDFT(n2ft3d,a,one,tmp,one);

   std::memcpy(d3db::d3db_tmp1,a,2*nfft3d*sizeof(double));
   std::memset(a,  0,2*nfft3d*sizeof(double));

   c_aindexcopy(nida[nb]+nidb2[nb],packarray[nb],d3db::d3db_tmp1,a);

   if (balanced)
      mybalance->c_balance(nb, a);

   //delete [] tmp;
   return;
}

/********************************
 *                              *
 *       PGrid:cc_pack_copy     *
 *                              *
 ********************************/
void PGrid::cc_pack_copy(const int nb, const double *a, double *b)
{
   //int one = 1;
   // int ng  = 2*(nida[nb]+nidb[nb]);
   int ng = 2*(nida[nb]+nidb[nb]);

   // DCOPY_PWDFT(ng,a,one,b,one);
   std::memcpy(b,a,ng*sizeof(double));
}

/********************************
 *                              *
 *       PGrid:cc_pack_dot      *
 *                              *
 ********************************/
 /**
 * @brief Computes a modified dot product for complex number arrays in a parallel computation environment.
 *
 * This function calculates a specialized dot product by taking into account the full span of complex number entries
 * (real and imaginary parts) from two input vectors 'a' and 'b'. The computation doubles the dot product of the
 * complete range and then subtracts the dot product of a subset (only the first half of the range), effectively
 * performing an operation that might be part of a larger numerical or analytical method involving complex numbers.
 *
 * The formula for the operation is:
 * \f[
 * \text{tsum} = 2.0 \cdot \text{DDOT_PWDFT}(ng, a, 1, b, 1) - \text{DDOT_PWDFT}(ng0, a, 1, b, 1)
 * \f]
 * where `ng` is the total number of real and imaginary parts combined from both `nida` and `nidb`, and `ng0` is the count
 * from `nida` alone, each doubled to account for complex numbers (real and imaginary parts).
 *
 * @param nb An index used to access specific dimensions or counts from the `nida` and `nidb` arrays, which define the range
 *           of elements involved in the dot product operations.
 * @param a Pointer to the first element of the vector 'a', expected to hold complex numbers represented in interleaved format
 *          (real followed by imaginary).
 * @param b Pointer to the first element of the vector 'b', structured identically to 'a'.
 *
 * @return The global sum of the computed dot product across all parallel processes, ensuring that all parts of the distributed
 *         system contribute to the final result.
 *
 * @note The function uses `DDOT_PWDFT` for dot product calculations, which should be capable of handling distributed arrays
 *       safely in a parallel computation context. The final result is summed across all participating nodes or processes using
 *       `d3db::parall->SumAll`, integral to reducing results in parallel applications.
 */
double PGrid::cc_pack_dot(const int nb, double *a, double *b) 
{
   int one = 1;
   // int ng  = 2*(nida[nb]+nidb[nb]);
   int ng  = 2*(nida[nb] + nidb[nb]);
   int ng0 = 2*nida[nb];
 
   double tsum  = 2.0*DDOT_PWDFT(ng, a,one,b,one);
   tsum -= DDOT_PWDFT(ng0,a,one,b,one);
 
   return d3db::parall->SumAll(1,tsum);
}

/********************************
 *                              *
 *       PGrid:tt_pack_dot      *
 *                              *
 ********************************/
 /**
 * @brief Computes a modified dot product of two vectors suitable for parallel execution environments, specifically designed for handling distributed data arrays.
 *
 * This function performs a specialized dot product operation by initially doubling the result of a dot product over a combined range
 * of elements from two indices arrays (nida and nidb), then subtracting the dot product of a subset of this range. The function is likely
 * used in contexts where adjustments to standard dot product calculations are necessary due to the nature of the data distribution or
 * computational methodology in use (e.g., partial Fourier space representations in electronic structure calculations).
 *
 * The operation is as follows:
 * \f[
 * \text{tsum} = 2 \cdot \text{DDOT_PWDFT}(ng, a, one, b, one) - \text{DDOT_PWDFT}(ng0, a, one, b, one)
 * \f]
 * where `ng` is the sum of elements defined by `nida[nb]` and `nidb[nb]`, and `ng0` is the count of elements defined by `nida[nb]` alone.
 *
 * @param nb An index used to select the particular size settings from the `nida` and `nidb` arrays, which dictate the range of elements
 *           involved in the computation.
 * @param a Pointer to the first element of the vector 'a'.
 * @param b Pointer to the first element of the vector 'b'.
 *
 * @return The global sum of the computed dot product across all parallel processes, ensuring that all parts of the distributed system
 *         contribute to the final result.
 *
 * @note The function relies on `DDOT_PWDFT` for performing the dot product calculation, which is assumed to be a parallel-safe version
 *       typically used in distributed computational frameworks. The function also uses `d3db::parall->SumAll` for summing up the result
 *       across different processors or nodes, which is a part of a parallel reduction operation.
 */
double PGrid::tt_pack_dot(const int nb, double *a, double *b) 
{
   int one = 1;
   int ng = (nida[nb] + nidb[nb]);
   int ng0 = nida[nb];
 
   double tsum = 2.0*DDOT_PWDFT(ng,a,one,b,one);
   tsum -= DDOT_PWDFT(ng0,a,one,b,one);
 
   return d3db::parall->SumAll(1,tsum);
}

/********************************
 *                              *
 *       PGrid:cc_pack_idot     *
 *                              *
 ********************************/
double PGrid::cc_pack_idot(const int nb, double *a, double *b) 
{
   int one = 1;
   // int ng  = 2*(nida[nb]+nidb[nb]);
   int ng = 2 * (nida[nb] + nidb[nb]);
   int ng0 = 2 * nida[nb];
   double tsum = 0.0;
 
   tsum = 2.0 * DDOT_PWDFT(ng, a, one, b, one);
   tsum -= DDOT_PWDFT(ng0, a, one, b, one);
 
   return tsum;
}

/********************************
 *                              *
 *       PGrid:tt_pack_idot     *
 *                              *
 ********************************/
double PGrid::tt_pack_idot(const int nb, double *a, double *b) 
{
   int one = 1;
   // int ng  = 2*(nida[nb]+nidb[nb]);
   int ng = (nida[nb] + nidb[nb]);
   int ng0 = nida[nb];
   double tsum = 0.0;
 
   tsum = 2.0 * DDOT_PWDFT(ng, a, one, b, one);
   tsum -= DDOT_PWDFT(ng0, a, one, b, one);
 
   return tsum;
}

/********************************
 *                              *
 *       PGrid:cc_pack_indot    *
 *                              *
 ********************************/
void PGrid::cc_pack_indot(const int nb, const int nn, double *a, double *b, double *sum) 
{
   int one = 1;
   // int ng  = 2*(nida[nb]+nidb[nb]);
   int ng = 2 * (nida[nb] + nidb[nb]);
   int ng0 = 2 * nida[nb];
 
   for (int i = 0; i < nn; ++i) {
     sum[i] = 2.0 * DDOT_PWDFT(ng, &(a[i * ng]), one, b, one);
     sum[i] -= DDOT_PWDFT(ng0, &(a[i * ng]), one, b, one);
   }
}

/********************************
 *                              *
 *    PGrid:cc_pack_inprjdot    *
 *                              *
 ********************************/
void PGrid::cc_pack_inprjdot(const int nb, int nn, int nprj, double *a,
                             double *b, double *sum) {
  int ng = 2 * (nida[nb] + nidb[nb]);
  int ng0 = 2 * nida[nb];
  int one = 1;
  double rtwo = 2.0;
  double rone = 1.0;
  double rmone = -1.0;
  double rzero = 0.0;

  // DGEMM_PWDFT((char *) "T",(char *) "N",nn,nprj,ng,
  // 	       rtwo,
  // 	       a,ng,
  // 	       b,ng,
  // 	       rzero,
  // 	       sum,nn);
  d3db::mygdevice.TN_dgemm(nn, nprj, ng, rtwo, a, b, rzero, sum);

  if (ng0 > 0) {
    DGEMM_PWDFT((char *)"T", (char *)"N", nn, nprj, ng0, rmone, a, ng, b, ng,
                rone, sum, nn);
  }
}


/********************************
 *                              *
 *   PGrid:ch3t_pack_i3ndot     *
 *                              *
 ********************************/
 //mypneb->ch3t_pack_i3ndot(1,nn,nprj,psi,prj,Gx,Gy,Gz,sum+3*nn(l+nprjall));

/*
void PGrid::ch3t_pack_i3ndot(const int nb, const int nn, const int nprj, 
                             const double *psi,
                             const double *prj,
                             const double Gx, const double Gy, const double Gz,
                             double *sum3)
{
   int ng = 2 * (nida[nb] + nidb[nb]);
   int ng0 = 2 * nida[nb];
   int one = 1;
   double rtwo = 2.0;
   double rone = 1.0;
   double rmone = -1.0;

  d3db::mygdevice.TN_dgemm3(nn,nprj,ng,rtwo,psi,prj,Gx,Gy,Gz,rzero,sum3);

  if (ng0 > 0) 
  {
     DGEMM3_PWDFT((char *)"T", (char *)"N",nn,nprj,ng0,rmone,psi,ng,prj,ng,Gx,Gy,Gz,
                  rone,sum3,nn);
  }
}
*/




/********************************
 *                              *
 *       PGrid:t_unpack         *
 *                              *
 ********************************/
void PGrid::t_unpack(const int nb, double *a) 
{

   int nn = (nida[nb]+nidb2[nb]);
   if (balanced)
      mybalance->t_unbalance(nb, a);
 
   std::memcpy(d3db::d3db_tmp1, a, nn * sizeof(double));
   std::memset(a, 0, nfft3d*sizeof(double));
 
   t_bindexcopy(nida[nb]+nidb2[nb],packarray[nb],d3db::d3db_tmp1,a);
 
   t_timereverse(a, zplane_tmp1, zplane_tmp2);
}

/********************************
 *                              *
 *       PGrid:t_pack           *
 *                              *
 ********************************/
void PGrid::t_pack(const int nb, double *a) 
{
   //int one = 1;
   //int zero = 0;
   //double rzero = 0.0;
   //double *tmp = new (std::nothrow) double[nfft3d];
 
   std::memcpy(d3db::d3db_tmp1,a,nfft3d*sizeof(double));
   std::memset(a,  0,nfft3d*sizeof(double));

   t_aindexcopy(nida[nb]+nidb2[nb],packarray[nb],d3db::d3db_tmp1,a);

   if (balanced)
      mybalance->t_balance(nb, a);
 
   //delete [] tmp;
}

/********************************
 *                              *
 *       PGrid:tt_pack_copy     *
 *                              *
 ********************************/
void PGrid::tt_pack_copy(const int nb, double *a, double *b) 
{
   //int one = 1;
   int ng = nida[nb] + nidb[nb];

   // DCOPY_PWDFT(ng,a,one,b,one);
   std::memcpy(b, a, ng * sizeof(double));
}

/********************************
 *                              *
 *       PGrid:t_pack_nzero     *
 *                              *
 ********************************/
void PGrid::t_pack_nzero(const int nb, const int n, double *a) 
{
   // int one   = 1;
   // int zero = 0;
   // double azero = 0.0;
   int ng = n * (nida[nb] + nidb[nb]);
   std::memset(a, 0, ng * sizeof(double));
}

/********************************
 *                              *
 *       PGrid:i_pack           *
 *                              *
 ********************************/
void PGrid::i_pack(const int nb, int *a) 
{
   int *tmp = new (std::nothrow) int[nfft3d];
 
   for (auto i=0; i<nfft3d; ++i) tmp[i] = a[i];
   for (auto i=0; i<nfft3d; ++i) a[i] = 0;
 
   i_aindexcopy(nida[nb] + nidb2[nb], packarray[nb], tmp, a);
 
   if (balanced)
      mybalance->i_balance(nb, a);
 
   delete[] tmp;
}

/********************************
 *                              *
 *       PGrid:ii_pack_copy     *
 *                              *
 ********************************/
void PGrid::ii_pack_copy(const int nb, int *a, int *b) 
{
   int ng = nida[nb]+nidb[nb];
   for (auto i=0; i<ng; ++i)
      b[i] = a[i];
}

/********************************
 *                              *
 *    PGrid:cr_pfft3b           *
 *                              *
 ********************************/
/* This routine performs the operation of a three dimensional
  complex to complex inverse fft:
   A(nx,ny(nb),nz(nb)) <- FFT3^(-1)[A(kx,ky,kz)]
*/
void PGrid::cr_pfft3b(const int nb, double *a)
{

   nwpw_timing_function ftime(1);
   int  indx0, indx2, nn;
   //double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
   //double *tmp3 = new (std::nothrow) double[2 * nfft3d]();
   double *tmp2 = d3db::d3db_tmp1;
   double *tmp3 = d3db::d3db_tmp2;
 
   auto nxh = nx / 2 + 1;
   auto nxhy = nxh * ny;
   auto nxhz = nxh * nz;
   auto nxh2 = nx + 2;
   auto nxhy2 = nxh2 * ny;
   auto nxhz2 = nxh2 * nz;
 
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      //std::memset(tmp2,0,2*nfft3d*sizeof(double));

     /***************************************************
      ***     do fft along kz dimension               ***
      ***     A(kx,nz,ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
      ***************************************************/
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nxh; ++i) 
         {
            if (!zero_row3[nb][indx2]) 
            {
               auto kk = 0;
               auto indx3 = 2 * i + indx0;
               auto shift = 2 * nz * nn;
               for (auto k=0; k<nz; ++k) 
               {
                  tmp2[kk   + shift] = a[indx3];
                  tmp2[kk+1 + shift] = a[indx3+1];
                  kk += 2;
                  indx3 += nxh2;
               }
               ++nn;
            }
            ++indx2;
         }
         indx0 += nxhz2;
      }
 
      d3db::mygdevice.batch_cfftz_tmpz(d3db::fft_tag,false, nz, nn, n2ft3d, tmp2, d3db::tmpz);
 
      //std::memset(a,0,n2ft3d*sizeof(double));
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nxh; ++i) 
         {
            if (!zero_row3[nb][indx2]) 
            {
               auto kk = 0;
               auto indx3 = 2 * i + indx0;
               auto shift = 2 * nz * nn;
               for (auto k=0; k<nz; ++k) 
               {
                  a[indx3]   = tmp2[kk   + shift];
                  a[indx3+1] = tmp2[kk+1 + shift];
                  kk += 2;
                  indx3 += nxh2;
               }
               ++nn;
            }
            ++indx2;
         }
         indx0 += nxhz2;
      }
      
      
      /***********************************************
       ***         Do a ptranspose of A            ***
       ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
       ************************************************/
      d3db::c_ptranspose1_jk(nb, a, tmp2, tmp3);
      
      /*************************************************
       ***        do fft along ky dimension          ***
       ***    A(kx,ny,nz) <- fft1d^(-1)[A(kx,ky,nz)] ***
       *************************************************/
      //std::memset(tmp2,0,2*nfft3d*sizeof(double));
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nxh; ++i) 
         {
            if (!zero_row2[nb][indx2]) 
            {
               auto jj = 0;
               auto indx3 = 2 * i + indx0;
               auto shift = 2 * ny * nn;
               for (auto j=0; j<ny; ++j) 
               {
                  tmp2[jj   + shift] = a[indx3];
                  tmp2[jj+1 + shift] = a[indx3+1];
                  jj += 2;
                  indx3 += nxh2;
               }
               ++nn;
            }
            ++indx2;
         }
         indx0 += nxhy2;
      }
      
      d3db::mygdevice.batch_cffty_tmpy(d3db::fft_tag,false, ny, nn, n2ft3d, tmp2, d3db::tmpy);
      
      //std::memset(a,0,n2ft3d*sizeof(double));
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nxh; ++i) 
         {
            if (!zero_row2[nb][indx2]) 
            {
               auto jj = 0;
               auto indx3 = 2 * i + indx0;
               auto shift = 2 * ny * nn;
               for (auto j=0; j<ny; ++j) 
               {
                  a[indx3]   = tmp2[jj   + shift];
                  a[indx3+1] = tmp2[jj+1 + shift];
                  jj += 2;
                  indx3 += nxh2;
               }
               ++nn;
            }
            ++indx2;
         }
         indx0 += nxhy2;
      }
      
      
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      /*
       d3db::cshift1_fftb(nx,ny,nq,1,a);
       indx = 0;
       for (q=0; q<nq; ++q)
       for (j=0; j<ny; ++j)
       {
          drfftb_(&nx,a+indx,d3db::tmpx);
          indx += nxh2;
       }
       d3db::zeroend_fftb(nx,ny,nq,1,a);
      */
      d3db::mygdevice.batch_rfftx_tmpx(d3db::fft_tag,false, nx, ny * nq, n2ft3d, a, d3db::tmpx);
      d3db::zeroend_fftb(nx, ny, nq, 1, a);
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
      d3db::mygdevice.batch_cfftz_tmpz_zero(d3db::fft_tag,false, nz, nq3, n2ft3d, a, d3db::tmpz, zero_row3[nb]);
     
      d3db::c_ptranspose_ijk(nb, 2, a, tmp2, tmp3);
      
      /************************************************
       ***     do fft along ky dimension            ***
       ***   A(ny,nz,kx) <- fft1d^(-1)[A(ky,nz,kx)] ***
       ************************************************/
      d3db::mygdevice.batch_cffty_tmpy_zero(d3db::fft_tag,false,ny,nq2,n2ft3d,a,d3db::tmpy,zero_row2[nb]);
      
      d3db::c_ptranspose_ijk(nb, 3, a, tmp2, tmp3);
      
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      d3db::mygdevice.batch_rfftx_tmpx(d3db::fft_tag,false, nx, nq1, n2ft3d, a, d3db::tmpx);
       
      d3db::zeroend_fftb(nx, nq1, 1, 1, a);
      if (n2ft3d_map < n2ft3d)
         std::memset(a + n2ft3d_map, 0, (n2ft3d - n2ft3d_map) * sizeof(double));
   }
 
   //delete[] tmp3;
   //delete[] tmp2;
}

/********************************
 *                              *
 *       PGrid::rc_pfft3f       *
 *                              *
 ********************************/
/*
   This routine performs the operation of a three dimensional
   complex to complex fft:
      A(kx,ky,kz) <- FFT3[A(nx(nb),ny(nb),nz(nb))]
*/
void PGrid::rc_pfft3f(const int nb, double *a)
{
   nwpw_timing_function ftime(1);
   int  indx0, indx2, nn;
   //double *tmp2 = new (std::nothrow) double[2 * nfft3d]();
   //double *tmp3 = new (std::nothrow) double[2 * nfft3d]();
   double *tmp2 = d3db::d3db_tmp1;
   double *tmp3 = d3db::d3db_tmp2;
 
   auto nxh = nx / 2 + 1;
   auto nxhy = nxh * ny;
   auto nxhz = nxh * nz;
   auto nxh2 = nx + 2;
   auto nxhy2 = nxh2 * ny;
   auto nxhz2 = nxh2 * nz;
 
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      /********************************************
       ***     do fft along nx dimension        ***
       ***   A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]  ***
       ********************************************/
      /*
      indx = 0;
      for (q=0; q<nq; ++q)
      for (j=0; j<ny; ++j)
      {
         drfftf_(&nx,a+indx,tmpx);
         indx += nxh2;
      }
      d3db::cshift_fftf(nx,ny,nq,1,a);
      */

      d3db::mygdevice.batch_rfftx_tmpx(d3db::fft_tag,true,nx,ny*nq,n2ft3d,a,d3db::tmpx);

      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
       ********************************************/
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nxh; ++i) 
         {
            if (!zero_row2[nb][indx2]) 
            {
               auto jj = 0;
               auto indx3 = 2 * i + indx0;
               auto shift = 2 * ny * nn;
               for (auto j=0; j<ny; ++j) 
               {
                  tmp2[jj   + shift] = a[indx3];
                  tmp2[jj+1 + shift] = a[indx3+1];
                  jj += 2;
                  indx3 += nxh2;
               }
               ++nn;
            }
            ++indx2;
         }
         indx0 += nxhy2;
      }
     
      d3db::mygdevice.batch_cffty_tmpy(d3db::fft_tag,true, ny, nn, n2ft3d, tmp2, d3db::tmpy);
     
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nxh; ++i) 
         {
            if (!zero_row2[nb][indx2]) 
            {
               auto jj = 0;
               auto indx3 = 2 * i + indx0;
               auto shift = 2 * ny * nn;
               for (auto j=0; j<ny; ++j) 
               {
                  a[indx3]   = tmp2[jj   + shift];
                  a[indx3+1] = tmp2[jj+1 + shift];
                  jj += 2;
                  indx3 += nxh2;
               }
               ++nn;
            }
            ++indx2;
         }
         indx0 += nxhy2;
      }
     
     
      /********************************************
       ***         Do a transpose of A          ***
       ***      A(ky,nz,ky) <- A(kx,ky,nz)      ***
       ********************************************/
      d3db::c_ptranspose2_jk(nb,a,tmp2,tmp3);
      // d3db::c_transpose_jk(a,tmp2,tmp3);

      /********************************************
       ***     do fft along nz dimension        ***
       ***   A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]  ***
       ********************************************/
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nxh; ++i) 
         {
            if (!zero_row3[nb][indx2]) 
            {
               auto kk = 0;
               auto indx3 = 2 * i + indx0;
               auto shift = 2 * nz * nn;
               for (auto k=0; k<nz; ++k) 
               {
                  tmp2[kk   + shift] = a[indx3];
                  tmp2[kk+1 + shift] = a[indx3+1];
                  kk += 2;
                  indx3 += nxh2;
               }
               ++nn;
            }
            ++indx2;
         }
         indx0 += nxhz2;
      }
     
      d3db::mygdevice.batch_cfftz_tmpz(d3db::fft_tag,true, nz, nn, n2ft3d, tmp2, d3db::tmpz);
     
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nxh; ++i) 
         {
            if (!zero_row3[nb][indx2]) 
            {
               auto kk = 0;
               auto indx3 = 2 * i + indx0;
               auto shift = 2 * nz * nn;
               for (auto k=0; k<nz; ++k) 
               {
                  a[indx3]   = tmp2[kk   + shift];
                  a[indx3+1] = tmp2[kk+1 + shift];
                  kk += 2;
                  indx3 += nxh2;
               }
               ++nn;
            }
            ++indx2;
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
      d3db::mygdevice.batch_rfftx_tmpx(d3db::fft_tag,true, nx, nq1, n2ft3d, a, d3db::tmpx);
      
      d3db::c_ptranspose_ijk(nb, 0, a, tmp2, tmp3);
     
      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
       ********************************************/
      d3db::mygdevice.batch_cffty_tmpy_zero(d3db::fft_tag,true,ny,nq2,n2ft3d,a,d3db::tmpy,zero_row2[nb]);
     
      d3db::c_ptranspose_ijk(nb, 1, a, tmp2, tmp3);
     
      /********************************************
       ***     do fft along nz dimension        ***
       ***   A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]  ***
       ********************************************/
      d3db::mygdevice.batch_cfftz_tmpz_zero(d3db::fft_tag,true, nz, nq3, n2ft3d, a, d3db::tmpz, zero_row3[nb]);
   }
 
   //delete[] tmp3;
   //delete[] tmp2;
}

/********************************
 *                              *
 *     PGrid::c_unpack_start    *
 *                              *
 ********************************/
void PGrid::c_unpack_start(const int nb, double *tmp1, double *tmp2,
                           const int request_indx, const int msgtype) {
  if (balanced)
    mybalance->c_unbalance_start(nb, tmp1, request_indx, msgtype);
}

/********************************
 *                              *
 *     PGrid::c_unpack_mid      *
 *                              *
 ********************************/
void PGrid::c_unpack_mid(const int nb, double *tmp1, double *tmp2,
                         const int request_indx, const int msgtype) {
  if (balanced)
    mybalance->c_unbalance_end(nb, tmp1, request_indx);

  std::memcpy(tmp2, tmp1, 2 * (nida[nb] + nidb2[nb]) * sizeof(double));
  std::memset(tmp1, 0, n2ft3d * sizeof(double));

  c_bindexcopy((nida[nb] + nidb2[nb]), packarray[nb], tmp2, tmp1);
  // c_bindexcopy(nida[nb]+nidb[nb],packarray[nb],tmp2,tmp1);

  d3db::c_timereverse_start(tmp1, zplane_tmp1, zplane_tmp2, request_indx,
                            msgtype);
}

/********************************
 *                              *
 *     PGrid::c_unpack_end      *
 *                              *
 ********************************/
void PGrid::c_unpack_end(const int nb, double *tmp1, double *tmp2,
                         const int request_indx) {
  d3db::c_timereverse_end(tmp1, zplane_tmp1, zplane_tmp2, request_indx);
}

/********************************
 *                              *
 *        PGrid::pfftbz         *
 *                              *
 ********************************/
void PGrid::pfftbz(const int nb, double *tmp1, double *tmp2, int request_indx) 
{

   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      auto nxh = nx / 2 + 1;
      auto nxh2 = nx + 2;
      auto nxhz2 = nxh2 * nz;
     
      /***************************************************
       ***     do fft along kz dimension               ***
       ***     A(kx,nz,ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
       ***************************************************/
      int indx0 = 0;
      int indx2 = 0;
      int nn = 0;
      for (auto q = 0; q < nq; ++q) {
        for (auto i = 0; i < nxh; ++i) {
          if (!zero_row3[nb][indx2]) {
            auto kk = 0;
            auto indx3 = 2 * i + indx0;
            auto shift = 2 * nz * nn;
            for (auto k = 0; k < nz; ++k) 
            {
               tmp2[kk   + shift] = tmp1[indx3];
               tmp2[kk+1 + shift] = tmp1[indx3 + 1];
               kk += 2;
               indx3 += nxh2;
            }
            nn += 1;
          }
          ++indx2;
        }
        indx0 += nxhz2;
      }
     
      d3db::mygdevice.batch_cfftz_tmpz(d3db::fft_tag,false, nz, nn, n2ft3d, tmp2, d3db::tmpz);
      // for (auto i=0; i<nn; ++i)
      //    dcfftb_(&nz,tmp2+2*nz*i,d3db::tmpz);
     
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q = 0; q < nq; ++q) {
        for (auto i = 0; i < nxh; ++i) {
          if (!zero_row3[nb][indx2]) {
            auto kk = 0;
            auto indx3 = 2 * i + indx0;
            auto shift = 2 * nz * nn;
            for (auto k = 0; k < nz; ++k) 
            {
               tmp1[indx3]   = tmp2[kk   + shift];
               tmp1[indx3+1] = tmp2[kk+1 + shift];
               kk += 2;
               indx3 += nxh2;
            }
            nn += 1;
          }
          ++indx2;
        }
        indx0 += nxhz2;
      }
     
     
      /***********************************************
       ***         Do a ptranspose of A            ***
       ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
       ************************************************/
      d3db::c_ptranspose1_jk_start(nb, tmp1, tmp2, tmp1, request_indx, 44);
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
      d3db::mygdevice.batch_cfftz_tmpz_zero(d3db::fft_tag,false, nz, nq3, n2ft3d, tmp1, d3db::tmpz,
                                    zero_row3[nb]);
     
      d3db::c_ptranspose_ijk_start(nb, 2, tmp1, tmp2, tmp1, request_indx, 45);
      // d3db::c_ptranspose_ijk(nb,2,tmp1,tmp2,tmp1);
   }
}

/********************************
 *                              *
 *        PGrid::pfftby         *
 *                              *
 ********************************/
void PGrid::pfftby(const int nb, double *tmp1, double *tmp2, int request_indx) 
{
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      auto nxh = nx / 2 + 1;
      auto nxh2 = nx + 2;
      auto nxhy2 = nxh2 * ny;
     
      /***********************************************
       ***         Do a ptranspose of A            ***
       ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
       ************************************************/
      d3db::c_ptranspose1_jk_end(nb, tmp2, tmp1, request_indx);
     
      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
       ********************************************/
      int indx0 = 0;
      int indx2 = 0;
      int nn = 0;
      for (auto q=0; q<nq; ++q) {
        for (auto i=0; i<nxh; ++i) {
          if (!zero_row2[nb][indx2]) 
          {
             auto jj = 0;
             auto indx3 = 2 * i + indx0;
             auto shift = 2 * ny * nn;
             for (auto j=0; j<ny; ++j) 
             {
                tmp1[jj   + shift] = tmp2[indx3];
                tmp1[jj+1 + shift] = tmp2[indx3 + 1];
                jj += 2;
                indx3 += nxh2;
             }
             nn += 1;
          }
          ++indx2;
        }
        indx0 += nxhy2;
      }
     
      d3db::mygdevice.batch_cffty_tmpy(d3db::fft_tag,false, ny, nn, n2ft3d, tmp1, d3db::tmpy);
      // for (auto i=0; i<nn; ++i)
      //    dcfftb_(&ny,tmp1+2*ny*i,d3db::tmpy);
     
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) {
        for (auto i=0; i<nxh; ++i) {
          if (!zero_row2[nb][indx2]) 
          {
             auto jj = 0;
             auto indx3 = 2 * i + indx0;
             auto shift = 2 * ny * nn;
             for (auto j = 0; j < ny; ++j) 
             {
                tmp2[indx3]   = tmp1[jj   + shift];
                tmp2[indx3+1] = tmp1[jj+1 + shift];
                jj += 2;
                indx3 += nxh2;
             }
             nn += 1;
          }
          ++indx2;
        }
        indx0 += nxhy2;
      }
   }
   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      d3db::c_ptranspose_ijk_end(nb, 2, tmp2, tmp1, request_indx);
     
      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
       ********************************************/
      d3db::mygdevice.batch_cffty_tmpy_zero(d3db::fft_tag,false,ny,nq2,n2ft3d,tmp2,d3db::tmpy,zero_row2[nb]);
     
      d3db::c_ptranspose_ijk_start(nb, 3, tmp2, tmp1, tmp2, request_indx, 46);
      // d3db::c_ptranspose_ijk(nb,3,tmp2,tmp1,tmp2);
   }
}

/********************************
 *                              *
 *        PGrid::pfftbx         *
 *                              *
 ********************************/
void PGrid::pfftbx(const int nb, double *tmp1, double *tmp2, int request_indx) 
{
 
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      d3db::mygdevice.batch_rfftx_tmpx(d3db::fft_tag,false, nx, ny * nq, n2ft3d, tmp2, d3db::tmpx);
      d3db::zeroend_fftb(nx, ny, nq, 1, tmp2);
      std::memcpy(tmp1, tmp2, n2ft3d * sizeof(double));
   }
   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      d3db::c_ptranspose_ijk_end(nb, 3, tmp1, tmp2, request_indx);
     
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      d3db::mygdevice.batch_rfftx_tmpx(d3db::fft_tag,false, nx, nq1, n2ft3d, tmp1, d3db::tmpx);
      d3db::zeroend_fftb(nx, nq1, 1, 1, tmp1);
      if (n2ft3d_map < n2ft3d)
         std::memset(tmp1 + n2ft3d_map, 0, (n2ft3d - n2ft3d_map) * sizeof(double));
   }
}

/********************************
 *                              *
 *       PGrid::pfftb_step      *
 *                              *
 ********************************/
void PGrid::pfftb_step(const int step, const int nb, double *a, double *tmp1,
                       double *tmp2, const int request_indx) 
{
   if (step == 0) {
     // parall->astart(request_indx,parall->np_i());
 
     // unpack start, tmp1-->tmp1
     std::memcpy(tmp1, a, 2 * (nida[nb] + nidb[nb]) * sizeof(double));
     this->c_unpack_start(nb, tmp1, tmp2, request_indx, 47);
   } else if (step == 1) {
     // unpack mid
     this->c_unpack_mid(nb, tmp1, tmp2, request_indx, 48);
   } else if (step == 2) {
     // unpack end; mem-->dev,  out=tmp1
     this->c_unpack_end(nb, tmp1, tmp2, request_indx);
   } else if (step == 3) {
     // pfftbz dev-->dev->mem,  tmp1->tmp1
     this->pfftbz(nb, tmp1, tmp2, request_indx);
   } else if (step == 4) {
     // pfftby mem->dev-->dev->mem
     // in=tmp1, tmp2->tmp1, tmp1=in , tmp2=tmp
     pfftby(nb, tmp1, tmp2, request_indx);
   } else if (step == 5) {
     // pfftbx mem->dev->dev->mem
     pfftbx(nb, tmp1, tmp2, request_indx);
     // parall->aend(request_indx);
   }
}

/********************************
 *                              *
 *      PGrid::pfftbz_start     *
 *                              *
 ********************************/
void PGrid::pfftbz_start(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) 
  {
     auto nxh = nx / 2 + 1;
     auto nxh2 = nx + 2;
     auto nxhz2 = nxh2 * nz;
    
     /***************************************************
      ***     do fft along kz dimension               ***
      ***     A(kx,nz,ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
      ***************************************************/
     int indx0 = 0;
     int indx2 = 0;
     int nn = 0;
     for (auto q=0; q<nq; ++q) 
     {
        for (auto i=0; i<nxh; ++i) 
        {
           if (!zero_row3[nb][indx2]) 
           {
              auto kk = 0;
              auto indx3 = 2 * i + indx0;
              auto shift = 2 * nz * nn;
              for (auto k=0; k<nz; ++k) 
              {
                 tmp2[kk   + shift] = tmp1[indx3];
                 tmp2[kk+1 + shift] = tmp1[indx3+1];
                 kk += 2;
                 indx3 += nxh2;
              }
              nn += 1;
           }
           ++indx2;
        }
        indx0 += nxhz2;
     }
    
     d3db::mygdevice.batch_cfftz_stages_tmpz(0,d3db::fft_tag,false, nz, nn, n2ft3d, tmp2, d3db::tmpz,da_indx);
     // for (auto i=0; i<nn; ++i)
     //    dcfftb_(&nz,tmp2+2*nz*i,d3db::tmpz);

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
     d3db::mygdevice.batch_cfftz_stages_tmpz_zero(0,d3db::fft_tag,false, nz, nq3, n2ft3d, tmp1, d3db::tmpz, zero_row3[nb], da_indx);
  }

}

/********************************
 *                              *
 *    PGrid::pfftbz_compute     *
 *                              *
 ********************************/
void PGrid::pfftbz_compute(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) 
  {
     auto nxh = nx / 2 + 1;
     auto nxh2 = nx + 2;
     auto nxhz2 = nxh2 * nz;
    
     /***************************************************
      ***     do fft along kz dimension               ***
      ***     A(kx,nz,ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
      ***************************************************/
     int indx0 = 0;
     int indx2 = 0;
     int nn = 0;
     for (auto q=0; q<nq; ++q) {
       for (auto i = 0; i < nxh; ++i) {
         if (!zero_row3[nb][indx2]) {
           nn += 1;
         }
         ++indx2;
       }
     }
    
     d3db::mygdevice.batch_cfftz_stages_tmpz(1,d3db::fft_tag,false, nz, nn, n2ft3d, tmp2, d3db::tmpz,da_indx);
     // for (auto i=0; i<nn; ++i)
     //    dcfftb_(&nz,tmp2+2*nz*i,d3db::tmpz);
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
      d3db::mygdevice.batch_cfftz_stages_tmpz_zero(1,d3db::fft_tag,false, nz, nq3, n2ft3d, tmp1, d3db::tmpz,zero_row3[nb],da_indx);
  }

}

/********************************
 *                              *
 *      PGrid::pfftbz_end       *
 *                              *
 ********************************/
void PGrid::pfftbz_end(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) 
  {
     auto nxh = nx / 2 + 1;
     auto nxh2 = nx + 2;
     auto nxhz2 = nxh2 * nz;
    
     /***************************************************
      ***     do fft along kz dimension               ***
      ***     A(kx,nz,ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
      ***************************************************/
     int indx0 = 0;
     int indx2 = 0;
     int nn = 0;
     for (auto q=0; q<nq; ++q) 
     {
        for (auto i=0; i<nxh; ++i) 
        {
           if (!zero_row3[nb][indx2]) 
           {
              nn += 1;
           }
           ++indx2;
        }
     }
    
     d3db::mygdevice.batch_cfftz_stages_tmpz(2,d3db::fft_tag,false, nz, nn, n2ft3d, tmp2, d3db::tmpz, da_indx);
     // for (auto i=0; i<nn; ++i)
     //    dcfftb_(&nz,tmp2+2*nz*i,d3db::tmpz);
    
     indx0 = 0;
     indx2 = 0;
     nn = 0;
     for (auto q=0; q<nq; ++q) {
       for (auto i=0; i<nxh; ++i) {
         if (!zero_row3[nb][indx2]) 
         {
            auto kk = 0;
            auto indx3 = 2 * i + indx0;
            auto shift = 2 * nz * nn;
            for (auto k=0; k<nz; ++k) 
            {
               tmp1[indx3]   = tmp2[kk   + shift];
               tmp1[indx3+1] = tmp2[kk+1 + shift];
               kk += 2;
               indx3 += nxh2;
            }
            nn += 1;
         }
         ++indx2;
       }
       indx0 += nxhz2;
     }
    
     /***********************************************
      ***         Do a ptranspose of A            ***
      ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
      ************************************************/
     d3db::c_ptranspose1_jk_start(nb, tmp1, tmp2, tmp1, request_indx, 44);
  }
  /*************************
   **** hilbert mapping ****
   *************************/
  else {

    /************************************************
     ***     do fft along kz dimension            ***
     ***   A(nz,kx,ky) <- fft1d^(-1)[A(kz,kx,ky)] ***
     ************************************************/
    d3db::mygdevice.batch_cfftz_stages_tmpz_zero(2,d3db::fft_tag,false, nz, nq3, n2ft3d, tmp1, d3db::tmpz, zero_row3[nb], da_indx);

    d3db::c_ptranspose_ijk_start(nb, 2, tmp1, tmp2, tmp1, request_indx, 45);
    // d3db::c_ptranspose_ijk(nb,2,tmp1,tmp2,tmp1);
  }

}



/********************************
 *                              *
 *      PGrid::pfftby_start     *
 *                              *
 ********************************/
void PGrid::pfftby_start(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) 
  {
     auto nxh = nx / 2 + 1;
     auto nxh2 = nx + 2;
     auto nxhy2 = nxh2 * ny;
    
     /***********************************************
      ***         Do a ptranspose of A            ***
      ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
      ************************************************/
     d3db::c_ptranspose1_jk_end(nb, tmp2, tmp1, request_indx);
    
     /********************************************
      ***     do fft along ny dimension        ***
      ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
      ********************************************/
     int indx0 = 0;
     int indx2 = 0;
     int nn = 0;
     for (auto q=0; q<nq; ++q) {
       for (auto i=0; i<nxh; ++i) {
         if (!zero_row2[nb][indx2]) 
         {
            auto jj = 0;
            auto indx3 = 2 * i + indx0;
            auto shift = 2 * ny * nn;
            for (auto j = 0; j < ny; ++j) 
            {
               tmp1[jj   + shift] = tmp2[indx3];
               tmp1[jj+1 + shift] = tmp2[indx3 + 1];
               jj += 2;
               indx3 += nxh2;
            }
            nn += 1;
         }
         ++indx2;
       }
       indx0 += nxhy2;
     }
    
     d3db::mygdevice.batch_cffty_stages_tmpy(0,d3db::fft_tag,false,ny,nn,n2ft3d,tmp1,d3db::tmpy,da_indx);
     // for (auto i=0; i<nn; ++i)
     //    dcfftb_(&ny,tmp1+2*ny*i,d3db::tmpy);

  }
  /*************************
   **** hilbert mapping ****
   *************************/
  else 
  {
     d3db::c_ptranspose_ijk_end(nb, 2, tmp2, tmp1, request_indx);
    
     /********************************************
      ***     do fft along ny dimension        ***
      ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
      ********************************************/
     d3db::mygdevice.batch_cffty_stages_tmpy_zero(0,d3db::fft_tag,false,ny,nq2,n2ft3d,tmp2,d3db::tmpy,zero_row2[nb],da_indx);
  }
}


/********************************
 *                              *
 *    PGrid::pfftby_compute     *
 *                              *
 ********************************/
void PGrid::pfftby_compute(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      auto nxh = nx / 2 + 1;
      auto nxh2 = nx + 2;
      auto nxhy2 = nxh2 * ny;
     
      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
       ********************************************/
      int indx0 = 0;
      int indx2 = 0;
      int nn = 0;
      for (auto q=0; q<nq; ++q) {
        for (auto i=0; i<nxh; ++i) {
          if (!zero_row2[nb][indx2]) 
          {
             nn += 1;
          }
          ++indx2;
        }
      }
     
      d3db::mygdevice.batch_cffty_stages_tmpy(1,d3db::fft_tag,false, ny, nn, n2ft3d, tmp1, d3db::tmpy,da_indx);
      // for (auto i=0; i<nn; ++i)
      //    dcfftb_(&ny,tmp1+2*ny*i,d3db::tmpy);
   }
   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
       ********************************************/
      d3db::mygdevice.batch_cffty_stages_tmpy_zero(1,d3db::fft_tag,false,ny,nq2,n2ft3d,tmp2,d3db::tmpy,zero_row2[nb],da_indx);
   }
}


/********************************
 *                              *
 *        PGrid::pfftby_end     *
 *                              *
 ********************************/
void PGrid::pfftby_end(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      auto nxh = nx / 2 + 1;
      auto nxh2 = nx + 2;
      auto nxhy2 = nxh2 * ny;
     
      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
       ********************************************/
      int indx0 = 0;
      int indx2 = 0;
      int nn = 0;
      for (auto q=0; q<nq; ++q) {
        for (auto i=0; i<nxh; ++i) {
          if (!zero_row2[nb][indx2]) 
          {
             nn += 1;
          }
          ++indx2;
        }
      }
     
      d3db::mygdevice.batch_cffty_stages_tmpy(2,d3db::fft_tag,false, ny, nn, n2ft3d, tmp1, d3db::tmpy,da_indx);
      // for (auto i=0; i<nn; ++i)
      //    dcfftb_(&ny,tmp1+2*ny*i,d3db::tmpy);
     
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) {
        for (auto i=0; i<nxh; ++i) {
          if (!zero_row2[nb][indx2]) 
          {
             auto jj = 0;
             auto indx3 = 2 * i + indx0;
             auto shift = 2 * ny * nn;
             for (auto j = 0; j < ny; ++j) 
             {
                tmp2[indx3]   = tmp1[jj   + shift];
                tmp2[indx3+1] = tmp1[jj+1 + shift];
                jj += 2;
                indx3 += nxh2;
             }
             nn += 1;
          }
          ++indx2;
        }
        indx0 += nxhy2;
      }
   }
   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
     /********************************************
      ***     do fft along ny dimension        ***
      ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
      ********************************************/
      d3db::mygdevice.batch_cffty_stages_tmpy_zero(2,d3db::fft_tag,false,ny,nq2,n2ft3d,tmp2,d3db::tmpy,zero_row2[nb],da_indx);
 
      d3db::c_ptranspose_ijk_start(nb, 3, tmp2, tmp1, tmp2, request_indx, 46);
      // d3db::c_ptranspose_ijk(nb,3,tmp2,tmp1,tmp2);
   }
}




/********************************
 *                              *
 *      PGrid::pfftbx_start     *
 *                              *
 ********************************/
void PGrid::pfftbx_start(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      d3db::mygdevice.batch_rfftx_stages_tmpx(0,d3db::fft_tag,false, nx, ny * nq, n2ft3d, tmp2, d3db::tmpx,da_indx);
   }
   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      d3db::c_ptranspose_ijk_end(nb, 3, tmp1, tmp2, request_indx);
     
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      d3db::mygdevice.batch_rfftx_stages_tmpx(0,d3db::fft_tag,false, nx, nq1, n2ft3d, tmp1, d3db::tmpx,da_indx);
   }
}


/********************************
 *                              *
 *    PGrid::pfftbx_compute     *
 *                              *
 ********************************/
void PGrid::pfftbx_compute(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      d3db::mygdevice.batch_rfftx_stages_tmpx(1,d3db::fft_tag,false, nx, ny * nq, n2ft3d, tmp2, d3db::tmpx,da_indx);
   }

   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      d3db::mygdevice.batch_rfftx_stages_tmpx(1,d3db::fft_tag,false, nx, nq1, n2ft3d, tmp1, d3db::tmpx,da_indx);
   }
}



/********************************
 *                              *
 *      PGrid::pfftbx_end       *
 *                              *
 ********************************/
void PGrid::pfftbx_end(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      d3db::mygdevice.batch_rfftx_stages_tmpx(2,d3db::fft_tag,false, nx, ny * nq, n2ft3d, tmp2, d3db::tmpx,da_indx);
      d3db::zeroend_fftb(nx, ny, nq, 1, tmp2);
      std::memcpy(tmp1, tmp2, n2ft3d * sizeof(double));
   }
   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      d3db::mygdevice.batch_rfftx_stages_tmpx(2,d3db::fft_tag,false, nx, nq1, n2ft3d, tmp1, d3db::tmpx,da_indx);
      d3db::zeroend_fftb(nx, nq1, 1, 1, tmp1);
      if (n2ft3d_map < n2ft3d)
         std::memset(tmp1 + n2ft3d_map, 0, (n2ft3d - n2ft3d_map) * sizeof(double));
   }
}




/********************************
 *                              *
 *       PGrid::pfftb_step12    *
 *                              *
 ********************************/
void PGrid::pfftb_step12(const int step, const int nb, double *a, double *tmp1,
                         double *tmp2, const int request_indx, const int indx)
{
   if (step == 0) {
     // parall->astart(request_indx,parall->np_i());

     // unpack start, tmp1-->tmp1
     std::memcpy(tmp1, a, 2 * (nida[nb] + nidb[nb]) * sizeof(double));
     this->c_unpack_start(nb, tmp1, tmp2, request_indx, 47);
   } else if (step == 1) {
     // unpack mid
     this->c_unpack_mid(nb, tmp1, tmp2, request_indx, 48);
   } else if (step == 2) {
     // unpack end; mem-->dev,  out=tmp1
     this->c_unpack_end(nb, tmp1, tmp2, request_indx);

   // pfftbz dev-->dev->mem,  tmp1->tmp1
   } else if (step == 3) {
     this->pfftbz_start(nb, tmp1, tmp2, request_indx,indx);
   } else if (step == 4) {
     this->pfftbz_compute(nb, tmp1, tmp2, request_indx,indx);
   } else if (step == 5) {
     this->pfftbz_end(nb, tmp1, tmp2, request_indx,indx);

   // pfftby mem->dev-->dev->mem
   // in=tmp1, tmp2->tmp1, tmp1=in , tmp2=tmp
   } else if (step == 6) {
     pfftby_start(nb, tmp1, tmp2, request_indx,indx);
   } else if (step == 7) {
     pfftby_compute(nb, tmp1, tmp2, request_indx,indx);
   } else if (step == 8) {
     pfftby_end(nb, tmp1, tmp2, request_indx,indx);

   // pfftbx mem->dev->dev->mem
   } else if (step == 9) {
     pfftbx_start(nb, tmp1, tmp2, request_indx,indx);
   } else if (step == 10) {
     pfftbx_compute(nb, tmp1, tmp2, request_indx,indx);
   } else if (step == 11) {
     pfftbx_end(nb, tmp1, tmp2, request_indx,indx);
     // parall->aend(request_indx);
   }
}

/********************************
 *                              *
 *    PGrid:cr_pfft3b_queuein   *
 *                              *
 ********************************/
void PGrid::cr_pfft3b_queuein(const int nb, double *a) {
  int shift1, shift2;
  int np = parall->np_i();

  for (auto q=0; q<aqsize; ++q) {
    int indx = aqindx[q];
    int status = aqstatus[indx] + 1;
    shift1 = n2ft3d*(2*indx);
    shift2 = n2ft3d*(2*indx + 1);
    if (staged_gpu_fft_pipeline)
       pfftb_step12(status, nb, a, atmp+shift1, atmp+shift2, indx+4,indx);
    else
       pfftb_step(status, nb, a, atmp+shift1, atmp+shift2, indx+4);
    ++aqstatus[indx];
  }

  ++alast_index;
  if (alast_index >= aqmax)
    alast_index = 0;
  ++aqsize;
  aqindx[aqsize - 1] = alast_index;
  aqstatus[alast_index] = 0;

  // status = 0;
  shift1 = n2ft3d*(2*alast_index);
  shift2 = n2ft3d*(2*alast_index+1);

  if (staged_gpu_fft_pipeline)
     pfftb_step12(0,nb,a,atmp+shift1,atmp+shift2, alast_index+4,alast_index);
  else
     pfftb_step(0, nb, a, atmp+shift1, atmp+shift2, alast_index+4);
}

/********************************
 *                              *
 *    PGrid:cr_pfft3b_queueout  *
 *                              *
 ********************************/
void PGrid::cr_pfft3b_queueout(const int nb, double *a) {
  int shift1, shift2;
  int indx1 = aqindx[0];

  //while (aqstatus[indx1] < 5) {
  while (aqstatus[indx1] < aqmax) {

    for (auto q = 0; q < aqsize; ++q) {
      int indx = aqindx[q];
      int status = aqstatus[indx] + 1;
      shift1 = n2ft3d * (2*indx);
      shift2 = n2ft3d * (2*indx+1);
      if (staged_gpu_fft_pipeline)
         pfftb_step12(status,nb,a,atmp+shift1,atmp+shift2,indx+4,indx);
      else
         pfftb_step(status,nb,a,atmp+shift1,atmp+shift2,indx+4);
      ++aqstatus[indx];
    }
  }
  double scal1 = 1.0 / ((double)((nx) * (ny) * (nz)));
  double enrr0 = scal1 * d3db::rr_dot(atmp, atmp);

  shift1 = n2ft3d * (2 * indx1);
  std::memcpy(a, atmp+shift1, n2ft3d*sizeof(double));
  --aqsize;
  for (auto q = 0; q < aqsize; ++q)
    aqindx[q] = aqindx[q+1];
}

/********************************
 *                              *
 * PGrid:cr_pfft3b_queuefilled  *
 *                              *
 ********************************/
int PGrid::cr_pfft3b_queuefilled() { return (aqsize >= aqmax); }

/********************************
 *                              *
 *        PGrid::pfftfx         *
 *                              *
 ********************************/
void PGrid::pfftfx(const int nb, double *a, double *tmp1, double *tmp2, int request_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      // do fft along nx dimension
      d3db::mygdevice.batch_rfftx_tmpx(d3db::fft_tag,true, nx, ny*nq, n2ft3d, a, d3db::tmpx);
      std::memcpy(tmp1, a, n2ft3d * sizeof(double));
   }
   /**** hilbert mapping ****/
   else 
   {
      // do fft along nx dimension
      // A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]
      d3db::mygdevice.batch_rfftx_tmpx(d3db::fft_tag,true, nx, nq1, n2ft3d, a, d3db::tmpx);
      d3db::c_ptranspose_ijk_start(nb, 0, a, tmp1, tmp2, request_indx, 40);
   }
}

/********************************
 *                              *
 *        PGrid::pfftfy         *
 *                              *
 ********************************/
void PGrid::pfftfy(const int nb, double *tmp1, double *tmp2, int request_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      auto nxh = nx / 2 + 1;
      auto nxh2 = nx + 2;
      auto nxhy2 = nxh2 * ny;
     
      // do fft along ny dimension
      // A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]
      int indx0 = 0;
      int indx2 = 0;
      int nn = 0;
      for (auto q = 0; q < nq; ++q) {
        for (auto i = 0; i < nxh; ++i) {
          if (!zero_row2[nb][indx2]) 
          {
             auto jj = 0;
             auto indx3 = 2 * i + indx0;
             auto shift = 2 * ny * nn;
             for (auto j = 0; j < ny; ++j) 
             {
                tmp2[jj   + shift] = tmp1[indx3];
                tmp2[jj+1 + shift] = tmp1[indx3 + 1];
                jj += 2;
                indx3 += nxh2;
             }
             ++nn;
          }
          ++indx2;
        }
        indx0 += nxhy2;
      }
     
      d3db::mygdevice.batch_cffty_tmpy(d3db::fft_tag,true, ny, nn, n2ft3d, tmp2, d3db::tmpy);
     
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) {
        for (auto i=0; i<nxh; ++i) {
          if (!zero_row2[nb][indx2]) 
          {
             auto jj = 0;
             auto indx3 = 2 * i + indx0;
             auto shift = 2 * ny * nn;
             for (auto j = 0; j < ny; ++j) 
             {
                tmp1[indx3] = tmp2[jj + shift];
                tmp1[indx3 + 1] = tmp2[jj + 1 + shift];
                jj += 2;
                indx3 += nxh2;
             }
             ++nn;
          }
          ++indx2;
        }
        indx0 += nxhy2;
      }
     
      // Do a transpose of A
      // A(ky,nz,ky) <- A(kx,ky,nz)
      d3db::c_ptranspose2_jk_start(nb, tmp1, tmp2, tmp1, request_indx, 41);
   }
 
   /**** hilbert mapping ****/
   else 
   {
      // in=tmp1, out=tmp2
      d3db::c_ptranspose_ijk_end(nb, 0, tmp1, tmp2, request_indx);
     
      // do fft along ny dimension
      // A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]
      d3db::mygdevice.batch_cffty_tmpy_zero(d3db::fft_tag,true,ny,nq2,n2ft3d,tmp1,d3db::tmpy,zero_row2[nb]);
     
      // in=tmp2, out=tmp2
      d3db::c_ptranspose_ijk_start(nb, 1, tmp1, tmp2, tmp1, request_indx, 42);
   }
}

/********************************
 *                              *
 *        PGrid::pfftfz         *
 *                              *
 ********************************/
void PGrid::pfftfz(const int nb, double *tmp1, double *tmp2, int request_indx) 
{

   /**** slab mapping ****/
   if (maptype == 1) 
   {
      d3db::c_ptranspose2_jk_end(nb, tmp2, tmp1, request_indx);
     
      auto nxh = nx / 2 + 1;
      auto nxh2 = nx + 2;
      auto nxhz2 = nxh2 * nz;
     
      // do fft along nz dimension
      // A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]
      int indx0 = 0;
      int indx2 = 0;
      int nn = 0;
      for (auto q=0; q<nq; ++q) {
        for (auto i=0; i<nxh; ++i) {
          if (!zero_row3[nb][indx2]) 
          {
             auto kk = 0;
             auto indx3 = 2 * i + indx0;
             auto shift = 2 * nz * nn;
             for (auto k=0; k<nz; ++k) 
             {
                tmp1[kk   + shift] = tmp2[indx3];
                tmp1[kk+1 + shift] = tmp2[indx3 + 1];
                kk += 2;
                indx3 += nxh2;
             }
             ++nn;
          }
          ++indx2;
        }
        indx0 += nxhz2;
      }
     
      d3db::mygdevice.batch_cfftz_tmpz(d3db::fft_tag,true, nz, nn, n2ft3d, tmp1, d3db::tmpz);
     
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) {
        for (auto i=0; i<nxh; ++i) {
          if (!zero_row3[nb][indx2]) 
          {
             auto kk = 0;
             auto indx3 = 2 * i + indx0;
             auto shift = 2 * nz * nn;
             for (auto k = 0; k < nz; ++k) {
               tmp2[indx3]   = tmp1[kk   + shift];
               tmp2[indx3+1] = tmp1[kk+1 + shift];
               kk += 2;
               indx3 += nxh2;
             }
             ++nn;
          }
          ++indx2;
        }
        indx0 += nxhz2;
      }
   }
 
   /**** hilbert mapping ****/
   else 
   {
      // in=tmp1, out=tmp2
      d3db::c_ptranspose_ijk_end(nb, 1, tmp2, tmp1, request_indx);
     
      // do fft along nz dimension
      // A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]
      d3db::mygdevice.batch_cfftz_tmpz_zero(d3db::fft_tag,true, nz, nq3, n2ft3d, tmp2, d3db::tmpz,zero_row3[nb]);
   }
}

/********************************
 *                              *
 *       PGrid::pfftf_step      *
 *                              *
 ********************************/
void PGrid::pfftf_step(const int step, const int nb, double *a, double *tmp1, double *tmp2, int request_indx)
{
   if (step==0)
   {
      // pfftfx mem-->device, in=a out=tmp2
      pfftfx(nb, a, tmp1, tmp2, request_indx);
   }
   else if (step==1)
   {
      // pfftfy device, in=tmp1
      pfftfy(nb, tmp1, tmp2, request_indx);
   }
   else if (step==2)
   {
      // pfftfz device-->mem
      pfftfz(nb, tmp1, tmp2, request_indx);
      this->c_pack_start(nb, tmp2, tmp1, request_indx, 47);
   }
   else if (step==3)
   {
      // pfftf final
      this->c_pack_end(nb, tmp2, request_indx);
   }
}

/********************************
 *                              *
 *    PGrid::pfftfx_start       *
 *                              *
 ********************************/
void PGrid::pfftfx_start(const int nb, double *a, double *tmp1, double *tmp2, int request_indx, int da_indx)
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      // do fft along nx dimension
      std::memcpy(tmp2, a, n2ft3d * sizeof(double));
      d3db::mygdevice.batch_rfftx_stages_tmpx(0,d3db::fft_tag,true, nx, ny*nq, n2ft3d, tmp2, d3db::tmpx,da_indx);
   }
   /**** hilbert mapping ****/
   else 
   {
      // do fft along nx dimension
      // A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]
      std::memcpy(tmp2, a, n2ft3d * sizeof(double));
      d3db::mygdevice.batch_rfftx_stages_tmpx(0,d3db::fft_tag,true, nx, nq1, n2ft3d, tmp2, d3db::tmpx,da_indx);
   }
}

/********************************
 *                              *
 *    PGrid::pfftfx_compute     *
 *                              *
 ********************************/
void PGrid::pfftfx_compute(const int nb, double *a, double *tmp1, double *tmp2, int request_indx, int da_indx)
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      // do fft along nx dimension
      d3db::mygdevice.batch_rfftx_stages_tmpx(1,d3db::fft_tag,true, nx, ny*nq, n2ft3d, tmp2, d3db::tmpx,da_indx);
   }
   /**** hilbert mapping ****/
   else 
   {
      // do fft along nx dimension
      // A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]
      d3db::mygdevice.batch_rfftx_stages_tmpx(1,d3db::fft_tag,true, nx, nq1, n2ft3d, tmp2, d3db::tmpx,da_indx);
   }
}


/********************************
 *                              *
 *      PGrid::pfftfx_end       *
 *                              *
 ********************************/
void PGrid::pfftfx_end(const int nb, double *a, double *tmp1, double *tmp2, int request_indx, int da_indx)
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      // do fft along nx dimension
      d3db::mygdevice.batch_rfftx_stages_tmpx(2,d3db::fft_tag,true, nx, ny*nq, n2ft3d, tmp2, d3db::tmpx,da_indx);
      std::memcpy(tmp1, tmp2, n2ft3d * sizeof(double));
   }
   /**** hilbert mapping ****/
   else 
   {
      // do fft along nx dimension
      // A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]
      d3db::mygdevice.batch_rfftx_stages_tmpx(2,d3db::fft_tag,true, nx, nq1, n2ft3d, tmp2, d3db::tmpx,da_indx);
 
      d3db::c_ptranspose_ijk_start(nb, 0, tmp2, tmp1, tmp2, request_indx, 40);
 
   }
}



/********************************
 *                              *
 *      PGrid::pfftfy_start     *
 *                              *
 ********************************/
void PGrid::pfftfy_start(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      auto nxh = nx / 2 + 1;
      auto nxh2 = nx + 2;
      auto nxhy2 = nxh2 * ny;
     
      // do fft along ny dimension
      // A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]
      int indx0 = 0;
      int indx2 = 0;
      int nn = 0;
      for (auto q=0; q<nq; ++q) {
        for (auto i=0; i<nxh; ++i) {
          if (!zero_row2[nb][indx2]) {
            auto jj = 0;
            auto indx3 = 2 * i + indx0;
            auto shift = 2 * ny * nn;
            for (auto j=0; j<ny; ++j) 
            {
               tmp2[jj   + shift] = tmp1[indx3];
               tmp2[jj+1 + shift] = tmp1[indx3 + 1];
               jj += 2;
               indx3 += nxh2;
            }
            ++nn;
          }
          ++indx2;
        }
        indx0 += nxhy2;
      }
     
      d3db::mygdevice.batch_cffty_stages_tmpy(0,d3db::fft_tag,true, ny, nn, n2ft3d, tmp2, d3db::tmpy,da_indx);
   }
 
   /**** hilbert mapping ****/
   else 
   {
      // in=tmp1, out=tmp2
      d3db::c_ptranspose_ijk_end(nb, 0, tmp1, tmp2, request_indx);
     
      // do fft along ny dimension
      // A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]
      d3db::mygdevice.batch_cffty_stages_tmpy_zero(0,d3db::fft_tag,true,ny,nq2,n2ft3d,tmp1,d3db::tmpy,zero_row2[nb],da_indx);
   }
}

/********************************
 *                              *
 *      PGrid::pfftfy_compute   *
 *                              *
 ********************************/
void PGrid::pfftfy_compute(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      auto nxh = nx / 2 + 1;
      auto nxh2 = nx + 2;
      auto nxhy2 = nxh2 * ny;
     
      // do fft along ny dimension
      // A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]
      int indx0 = 0;
      int indx2 = 0;
      int nn = 0;
      for (auto q=0; q<nq; ++q) {
        for (auto i=0; i<nxh; ++i) {
          if (!zero_row2[nb][indx2]) 
          {
             ++nn;
          }
          ++indx2;
        }
      }
     
      d3db::mygdevice.batch_cffty_stages_tmpy(1,d3db::fft_tag,true, ny, nn, n2ft3d, tmp2, d3db::tmpy,da_indx);
   }
 
   /**** hilbert mapping ****/
   else 
   {
      // in=tmp1, out=tmp2
      // do fft along ny dimension
      // A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]
      d3db::mygdevice.batch_cffty_stages_tmpy_zero(1,d3db::fft_tag,true,ny,nq2,n2ft3d,tmp1,d3db::tmpy,zero_row2[nb],da_indx);
   }
}




/********************************
 *                              *
 *      PGrid::pfftfy_end       *
 *                              *
 ********************************/
void PGrid::pfftfy_end(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      auto nxh = nx / 2 + 1;
      auto nxh2 = nx + 2;
      auto nxhy2 = nxh2 * ny;
     
      // do fft along ny dimension
      // A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]
      int indx0 = 0;
      int indx2 = 0;
      int nn = 0;
      for (auto q=0; q<nq; ++q) {
        for (auto i=0; i<nxh; ++i) {
          if (!zero_row2[nb][indx2]) 
          {
             ++nn;
          }
          ++indx2;
        }
      }
     
      d3db::mygdevice.batch_cffty_stages_tmpy(2,d3db::fft_tag,true, ny, nn, n2ft3d, tmp2, d3db::tmpy,da_indx);
     
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) {
        for (auto i=0; i<nxh; ++i) {
          if (!zero_row2[nb][indx2]) 
          {
             auto jj = 0;
             auto indx3 = 2 * i + indx0;
             auto shift = 2 * ny * nn;
             for (auto j=0; j<ny; ++j) 
             {
               tmp1[indx3]   = tmp2[jj   + shift];
               tmp1[indx3+1] = tmp2[jj+1 + shift];
               jj += 2;
               indx3 += nxh2;
             }
             ++nn;
          }
          ++indx2;
        }
        indx0 += nxhy2;
      }
     
      // Do a transpose of A
      // A(ky,nz,ky) <- A(kx,ky,nz)
      d3db::c_ptranspose2_jk_start(nb, tmp1, tmp2, tmp1, request_indx, 41);
   }
   /**** hilbert mapping ****/
   else 
   {
      // do fft along ny dimension
      // A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]
      d3db::mygdevice.batch_cffty_stages_tmpy_zero(2,d3db::fft_tag,true,ny,nq2,n2ft3d,tmp1,d3db::tmpy,zero_row2[nb],da_indx);
     
      // in=tmp2, out=tmp2
      d3db::c_ptranspose_ijk_start(nb, 1, tmp1, tmp2, tmp1, request_indx, 42);
   }
}



/********************************
 *                              *
 *      PGrid::pfftfz_start     *
 *                              *
 ********************************/
void PGrid::pfftfz_start(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      d3db::c_ptranspose2_jk_end(nb, tmp2, tmp1, request_indx);
     
      auto nxh = nx / 2 + 1;
      auto nxh2 = nx + 2;
      auto nxhz2 = nxh2 * nz;
     
      // do fft along nz dimension
      // A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]
      int indx0 = 0;
      int indx2 = 0;
      int nn = 0;
      for (auto q=0; q<nq; ++q) {
        for (auto i=0; i<nxh; ++i) {
          if (!zero_row3[nb][indx2]) 
          {
            auto kk = 0;
            auto indx3 = 2 * i + indx0;
            auto shift = 2 * nz * nn;
            for (auto k=0; k<nz; ++k) 
            {
              tmp1[kk   + shift] = tmp2[indx3];
              tmp1[kk+1 + shift] = tmp2[indx3 + 1];
              kk += 2;
              indx3 += nxh2;
            }
            ++nn;
          }
          ++indx2;
        }
        indx0 += nxhz2;
      }
     
      d3db::mygdevice.batch_cfftz_stages_tmpz(0,d3db::fft_tag,true, nz, nn, n2ft3d, tmp1, d3db::tmpz,da_indx);
   }
   /**** hilbert mapping ****/
   else 
   {
      // in=tmp1, out=tmp2
      d3db::c_ptranspose_ijk_end(nb, 1, tmp2, tmp1, request_indx);
     
      // do fft along nz dimension
      // A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]
      d3db::mygdevice.batch_cfftz_stages_tmpz_zero(0,d3db::fft_tag,true, nz, nq3, n2ft3d, tmp2, d3db::tmpz, zero_row3[nb],da_indx);
   }
}


/********************************
 *                              *
 *      PGrid::pfftfz_compute   *
 *                              *
 ********************************/
void PGrid::pfftfz_compute(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      auto nxh = nx / 2 + 1;
      auto nxh2 = nx + 2;
      auto nxhz2 = nxh2 * nz;
     
      // do fft along nz dimension
      // A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]
      int indx0 = 0;
      int indx2 = 0;
      int nn = 0;
      for (auto q=0; q<nq; ++q) {
        for (auto i=0; i<nxh; ++i) {
          if (!zero_row3[nb][indx2]) 
          {
             ++nn;
          }
          ++indx2;
        }
      }
      d3db::mygdevice.batch_cfftz_stages_tmpz(1,d3db::fft_tag,true, nz, nn, n2ft3d, tmp1, d3db::tmpz,da_indx);
   }
   /**** hilbert mapping ****/
   else 
   {
      // do fft along nz dimension
      // A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]
      d3db::mygdevice.batch_cfftz_stages_tmpz_zero(1,d3db::fft_tag,true, nz, nq3, n2ft3d, tmp2, d3db::tmpz, zero_row3[nb],da_indx);
   }
}

/********************************
 *                              *
 *      PGrid::pfftfz_end       *
 *                              *
 ********************************/
void PGrid::pfftfz_end(const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
     auto nxh = nx / 2 + 1;
     auto nxh2 = nx + 2;
     auto nxhz2 = nxh2 * nz;
 
     // do fft along nz dimension
     // A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]
     int indx0 = 0;
     int indx2 = 0;
     int nn = 0;
     for (auto q=0; q<nq; ++q) {
       for (auto i=0; i<nxh; ++i) {
         if (!zero_row3[nb][indx2]) 
         {
            ++nn;
         }
         ++indx2;
       }
     }
 
     d3db::mygdevice.batch_cfftz_stages_tmpz(2,d3db::fft_tag,true, nz, nn, n2ft3d, tmp1, d3db::tmpz,da_indx);
 
     indx0 = 0;
     indx2 = 0;
     nn = 0;
     for (auto q=0; q<nq; ++q) {
       for (auto i=0; i<nxh; ++i) {
         if (!zero_row3[nb][indx2]) 
         {
            auto kk = 0;
            auto indx3 = 2 * i + indx0;
            auto shift = 2 * nz * nn;
            for (auto k = 0; k < nz; ++k) 
            {
               tmp2[indx3]   = tmp1[kk   + shift];
               tmp2[indx3+1] = tmp1[kk+1 + shift];
               kk += 2;
               indx3 += nxh2;
            }
            ++nn;
         }
         ++indx2;
       }
       indx0 += nxhz2;
     }
 
   }
 
   /**** hilbert mapping ****/
   else 
   {
      // do fft along nz dimension
      // A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]
      d3db::mygdevice.batch_cfftz_stages_tmpz_zero(2,d3db::fft_tag,true, nz, nq3, n2ft3d, tmp2, d3db::tmpz, zero_row3[nb],da_indx);
   }
}


/********************************
 *                              *
 *       PGrid::pfftf_step10    *
 *                              *
 ********************************/
void PGrid::pfftf_step10(const int step, const int nb, double *a, double *tmp1,
                         double *tmp2, int request_indx, int da_indx)
{
   // pfftfx mem-->device, in=a out=tmp2
   if (step==0)
   {
     pfftfx_start(nb, a, tmp1, tmp2, request_indx,da_indx);
   }
   else if (step==1)
   {
     pfftfx_compute(nb, a, tmp1, tmp2, request_indx,da_indx);
   }
   else if (step==2)
   {
     pfftfx_end(nb, a, tmp1, tmp2, request_indx,da_indx);
   }


   // pfftfy device, in=tmp1
   else if (step==3)
   {
      pfftfy_start(nb, tmp1, tmp2, request_indx,da_indx);
   }
   else if (step==4)
   {
      pfftfy_compute(nb, tmp1, tmp2, request_indx,da_indx);
   }
   else if (step==5)
   {
      pfftfy_end(nb, tmp1, tmp2, request_indx,da_indx);
   }


   // pfftfz device-->mem
   else if (step==6)
   {
      pfftfz_start(nb, tmp1, tmp2, request_indx,da_indx);
   }
   else if (step==7)
   {
      pfftfz_compute(nb, tmp1, tmp2, request_indx,da_indx);
   }
   else if (step==8)
   {
      pfftfz_end(nb, tmp1, tmp2, request_indx,da_indx);
      this->c_pack_start(nb, tmp2, tmp1, request_indx, 47);
   }


   else if (step==9)
   {
      // pfftf final
      this->c_pack_end(nb, tmp2, request_indx);
   }
}


/********************************
 *                              *
 *       PGrid:c_pack_start     *
 *                              *
 ********************************/
void PGrid::c_pack_start(const int nb, double *a, double *tmp1,
                         const int request_indx, const int msgtype) {
  // int one=1;

  // DCOPY_PWDFT(n2ft3d,a,one,tmp,one);
  std::memcpy(tmp1, a, n2ft3d * sizeof(double));
  std::memset(a, 0, n2ft3d * sizeof(double));

  c_aindexcopy(nida[nb] + nidb2[nb], packarray[nb], tmp1, a);

  if (balanced)
    mybalance->c_balance_start(nb, a, request_indx, msgtype);

  return;
}

/********************************
 *                              *
 *       PGrid:c_pack_end       *
 *                              *
 ********************************/
void PGrid::c_pack_end(const int nb, double *tmp1, const int request_indx) {

  if (balanced)
    mybalance->c_balance_end(nb, tmp1, request_indx);

  return;
}

/********************************
 *                              *
 *    PGrid:rc_pfft3f_queuein   *
 *                              *
 ********************************/
void PGrid::rc_pfft3f_queuein(const int nb, double *b) {
  int shift1, shift2;
  int np = parall->np_i();

  for (auto q = 0; q < bqsize; ++q) {
     int indx = bqindx[q];
     int status = bqstatus[indx] + 1;
     shift1 = n2ft3d * (2*indx);
     shift2 = n2ft3d * (2*indx + 1);
     if (staged_gpu_fft_pipeline)
        pfftf_step10(status, nb, b, btmp + shift1, btmp + shift2, indx+4,indx);
     else
        pfftf_step(status, nb, b, btmp + shift1, btmp + shift2, indx+4);
     ++bqstatus[indx];
  }

  ++blast_index;
  if (blast_index >= bqmax)
    blast_index = 0;
  ++bqsize;
  bqindx[bqsize - 1] = blast_index;
  bqstatus[blast_index] = 0;

  // status = 0;
  shift1 = n2ft3d * (2*blast_index);
  shift2 = n2ft3d * (2*blast_index + 1);

  if (staged_gpu_fft_pipeline)
     pfftf_step10(0, nb, b, btmp + shift1, btmp + shift2, blast_index+4,blast_index);
  else
     pfftf_step(0, nb, b, btmp + shift1, btmp + shift2, blast_index+4);
}

/********************************
 *                              *
 *    PGrid:rc_pfft3f_queueout  *
 *                              *
 ********************************/
void PGrid::rc_pfft3f_queueout(const int nb, double *b) {
  int shift1, shift2;
  int indx1 = bqindx[0];

  //while (bqstatus[indx1] < 5) {
  while (bqstatus[indx1] < bqmax) {

    for (auto q = 0; q < bqsize; ++q) {
      int indx = bqindx[q];
      int status = bqstatus[indx] + 1;
      shift1 = n2ft3d * (2*indx);
      shift2 = n2ft3d * (2*indx + 1);
      if (staged_gpu_fft_pipeline)
         pfftf_step10(status, nb, b, btmp + shift1, btmp + shift2, indx+4,indx);
      else
         pfftf_step(status, nb, b, btmp + shift1, btmp + shift2, indx+4);
      ++bqstatus[indx];
    }
  }
  double scal1 = 1.0 / ((double)((nx) * (ny) * (nz)));
  double enrr0 = scal1 * d3db::rr_dot(btmp, btmp);

  shift2 = n2ft3d * (2 * indx1 + 1);
  std::memcpy(b, btmp + shift2, n2ft3d * sizeof(double));
  --bqsize;
  for (auto q = 0; q < bqsize; ++q)
    bqindx[q] = bqindx[q + 1];
}

/********************************
 *                              *
 * PGrid:rc_pfft3f_queuefilled  *
 *                              *
 ********************************/
int PGrid::rc_pfft3f_queuefilled() { return (bqsize >= bqmax); }

/********************************
 *                              *
 *     PGrid:tc_pack_copy       *
 *                              *
 ********************************/
void PGrid::tc_pack_copy(const int nb, double *a, double *b) {
  int i, ii;
  int ng = nida[nb] + nidb[nb];

  ii = 0;
  for (i = 0; i < ng; ++i) {
    b[ii] = a[i];
    b[ii + 1] = 0.0;
    ii += 2;
  }
}

/********************************
 *                              *
 *      PGrid:tcc_pack_Mul      *
 *                              *
 ********************************/
void PGrid::tcc_pack_Mul(const int nb, const double *a, const double *b, double *c)
{
   int i, ii;
   int ng = nida[nb]+nidb[nb];

   ii = 0;
   for (i=0; i<ng; ++i)
   {
      c[ii]   = b[ii]*  a[i];
      c[ii+1] = b[ii+1]*a[i];
      ii += 2;
   }
}

/********************************
 *                              *
 *      PGrid:tcc_pack_aMul     *
 *                              *
 ********************************/
void PGrid::tcc_pack_aMul(const int nb, const double alpha, const double *a, const double *b, double *c)
{
   int i, ii;
   int ng = nida[nb]+nidb[nb];

   ii = 0;
   for (i=0; i<ng; ++i)
   {
      c[ii]   = alpha*b[ii]*  a[i];
      c[ii+1] = alpha*b[ii+1]*a[i];
      ii += 2;
   }
}

/********************************
 *                              *
 *      PGrid:tc_pack_Mul       *
 *                              *
 ********************************/
void PGrid::tc_pack_Mul(const int nb, const double *a, double *c) {
  int i, ii;
  int ng = nida[nb] + nidb[nb];

  ii = 0;
  for (i = 0; i < ng; ++i) {
    c[ii] = c[ii] * a[i];
    c[ii + 1] = c[ii + 1] * a[i];
    ii += 2;
  }
}

/********************************
 *                              *
 *    PGrid:tcc_pack_aMulAdd    *
 *                              *
 ********************************/
void PGrid::tcc_pack_aMulAdd(const int nb, const double alpha, const double *a, const double *b, double *c)
{
   int i, ii;
   int ng = nida[nb] + nidb[nb];

   ii=0;
   for (i=0; i<ng; ++i)
   {
      c[ii]   += alpha*b[ii] * a[i];
      c[ii+1] += alpha*b[ii+1]*a[i];
      ii += 2;
   }
}

/********************************
 *                              *
 *      PGrid:tcc_pack_iMul     *
 *                              *
 ********************************/
void PGrid::tcc_pack_iMul(const int nb, const double *a, const double *b,
                          double *c) {
  int i, ii;
  int ng = nida[nb] + nidb[nb];

  ii = 0;
  for (i = 0; i < ng; ++i) {
    c[ii] = -b[ii + 1] * a[i];
    c[ii + 1] = b[ii] * a[i];
    ii += 2;
  }
}

/*******************************************
 *                                         *
 *     PGrid:tcr_pack_iMul_unpack_fft      *
 *                                         *
 *******************************************/
void PGrid::tcr_pack_iMul_unpack_fft(const int nb, const double *a, const double *b, double *c)
{
   int i, ii;
   int ng = nida[nb] + nidb[nb];

   ii = 0;
   for (i=0; i<ng; ++i)
   {
      c[ii]   = -b[ii+1]* a[i];
      c[ii+1] = b[ii]   * a[i];
      ii += 2;
   }
   this->c_unpack(nb,c);
   this->cr_pfft3b(nb,c);
   this->d3db::r_zero_ends(c);
}

/********************************
 *                              *
 *     PGrid:tc_pack_iMul       *
 *                              *
 ********************************/
void PGrid::tc_pack_iMul(const int nb, const double *a, double *c) 
{
   double x, y;
   int ng = nida[nb] + nidb[nb];
   int ii = 0;
   for (auto i=0; i<ng; ++i) 
   {
      x = c[ii];
      y = c[ii+1];
     
      c[ii]   = -y*a[i];
      c[ii+1] =  x*a[i];
      ii += 2;
   }
}

/********************************
 *                              *
 *    PGrid:ttc_pack_MulSum2    *
 *                              *
 ********************************/
void PGrid::tcc_pack_MulSum2(const int nb, const double *a, const double *b, double *c) 
{
   int ng = nida[nb] + nidb[nb];
   int ii = 0;
   for (auto i=0; i<ng; ++i) 
   {
      c[ii]   += b[ii]   * a[i];
      c[ii+1] += b[ii+1] * a[i];
      ii += 2;
   }
}

/********************************
 *                              *
 *      PGrid:cc_pack_Sum2      *
 *                              *
 ********************************/
void PGrid::cc_pack_Sum2(const int nb, const double *a, double *b) 
{
   int ng = 2*(nida[nb] + nidb[nb]);

   for (auto i=0; i<ng; ++i)
      b[i] += a[i];
}

/********************************
 *                              *
 *     PGrid:cccc_pack_Sum      *
 *                              *
 ********************************/
void PGrid::cccc_pack_Sum(const int nb, const double *a, const double *b, const double *c, double *d) 
{
   int ng = 2*(nida[nb] + nidb[nb]);

   for (auto i=0; i<ng; ++i)
      d[i] = (a[i] + b[i] + c[i]);
}

/********************************
 *                              *
 *     PGrid:c_pack_addzero     *
 *                              *
 ********************************/
void PGrid::c_pack_addzero(const int nb, const double vzero, double *a) {
  int pzero = ijktop(0, 0, 0);
  if (pzero == parall->taskid_i())
    a[0] += vzero;
}

/********************************
 *                              *
 *   PGrid:c_pack_noimagzero    *
 *                              *
 ********************************/

void PGrid::c_pack_noimagzero(const int nb, double *a) 
{
   int pzero = ijktop(0, 0, 0);
   if (pzero == parall->taskid_i())
      a[1] = 0.0;
}

/********************************
 *                              *
 *     PGrid:c_pack_zero        *
 *                              *
 ********************************/
void PGrid::c_pack_zero(const int nb, double *b) 
{
   int ng = 2*(nida[nb] + nidb[nb]);

   for (auto i = 0; i < ng; ++i)
      b[i] = 0.0;
}

/********************************
 *                              *
 *       PGrid:c_pack_SMul      *
 *                              *
 ********************************/
void PGrid::c_pack_SMul(const int nb, const double alpha, double *b) 
{
   int ng = 2 * (nida[nb] + nidb[nb]);

   for (auto i = 0; i < ng; ++i)
      b[i] *= alpha;
}

/********************************
 *                              *
 *     PGrid:cc_pack_SMul       *
 *                              *
 ********************************/
void PGrid::cc_pack_SMul(const int nb, const double alpha, const double *a, double *b) 
{
   int ng = 2*(nida[nb] + nidb[nb]);
 
   for (auto i=0; i<ng; ++i)
     b[i] = alpha*a[i];
}

/********************************
 *                              *
 *      PGrid:cc_pack_daxpy     *
 *                              *
 ********************************/
/**
 * @brief Performs the DAXPY operation on complex number arrays, optimized for parallel computation.
 *
 * This function scales the complex number vector 'a' by a scalar 'alpha' and adds the result to the complex number vector 'b'.
 * It operates on interleaved complex numbers, where each complex number is represented by two consecutive double values (the real part followed by the imaginary part).
 * The operation is vectorized across the total number of real and imaginary parts indicated by 'nida' and 'nidb' for the given index 'nb'.
 *
 * The mathematical operation performed is:
 * \f[
 * b[i] := b[i] + \alpha \times a[i]
 * \f]
 * for each element \( i \) in the vector, where \( i \) ranges from 0 to \( 2 \times (\text{nida}[nb] + \text{nidb}[nb]) - 1 \), covering all parts of the complex numbers in 'a' and 'b'.
 *
 * @param nb An index used to determine the number of complex elements in 'a' and 'b', with the count retrieved from the arrays `nida` and `nidb`.
 * @param alpha A double precision scalar by which each element in the vector 'a' is scaled.
 * @param a Pointer to the first element of the input vector 'a', where complex numbers are stored in an interleaved format.
 * @param b Pointer to the first element of the output vector 'b', structured identically to 'a', to which the scaled elements of 'a' are added.
 *
 * @note Assumes that arrays 'a' and 'b' are pre-allocated with sufficient size to hold all necessary elements, which should include space for both the real and imaginary parts of the complex numbers.
 */
void PGrid::cc_pack_daxpy(const int nb, const double alpha, const double *a, double *b) 
{
   int ng = 2 * (nida[nb] + nidb[nb]);

   for (auto i=0; i<ng; ++i)
      b[i] += alpha*a[i];
}

/********************************
 *                              *
 *   PGrid:cct_pack_iconjgMul   *
 *                              *
 ********************************/
 /**
 * @brief Performs element-wise multiplication of complex numbers from array 'a' with the conjugates of complex numbers from array 'b', storing the imaginary part of each product in array 'c'.
 *
 * This function iterates over pairs of complex numbers taken from the arrays 'a' and 'b'. For each pair, it computes the imaginary part of the product of the complex number from 'a' and the conjugate of the complex number from 'b'. The result of this operation for each pair is stored in 'c'. The operation performed is:
 * \f[
 * C[i] = A[i] \times \overline{B[i]} = (a_{re} + i \times a_{im}) \times (b_{re} - i \times b_{im}) = (a_{re} \times b_{re} + a_{im} \times b_{im}) + i \times (a_{im} \times b_{re} - a_{re} \times b_{im})
 * \f]
 * Here, only the imaginary part \f$ (a_{im} \times b_{re} - a_{re} \times b_{im}) \f$ is stored in array 'c'.
 *
 * @param nb An index used to determine the number of complex pairs to process, based on the sum of entries from arrays `nida` and `nidb`.
 * @param a Pointer to the first element of the input array 'a', where each complex number is represented by consecutive elements (first the real part, then the imaginary part).
 * @param b Pointer to the first element of the input array 'b', structured identically to 'a'.
 * @param c Pointer to the first element of the output array where the results (imaginary parts of the products) are stored.
 *
 * @note Assumes that arrays 'a', 'b', and 'c' are allocated with sufficient size to hold the necessary number of elements, and that 'nida[nb]' and 'nidb[nb]' correctly reflect the counts of complex elements to be processed. This function calculates only the imaginary parts of the complex products and stores them in 'c'.
 */
void PGrid::cct_pack_iconjgMul(const int nb, const double *a, const double *b, double *c)
{
   for (auto i=0; i<(nida[nb]+nidb[nb]); ++i)
      c[i] = a[2*i]*b[2*i+1] - a[2*i+1]*b[2*i];
}

/********************************
 *                              *
 *  PGrid:cct_pack_iconjgMulb   *
 *                              *
 ********************************/
/**
 * @brief Performs element-wise multiplication of complex numbers from array 'a' with the conjugates of complex numbers from array 'b', specifically calculating a component typically associated with the imaginary part of the resultant complex product.
 *
 * This function iterates through pairs of complex numbers taken from the arrays 'a' and 'b'. For each pair, it computes:
 * \f[
 * C[i] = a_{im} \times b_{re} - a_{re} \times b_{im}
 * \f]
 * which is the negative of the imaginary component of the complex product if the second complex number were conjugated. The computed values are stored in the output array 'c'. This component is crucial in many applications involving complex arithmetic where only specific components of the product are required.
 *
 * @param nb An index used to determine the number of complex pairs to process, based on the sum of entries from arrays `nida` and `nidb`.
 * @param a Pointer to the first element of the input array 'a', where each complex number is represented by consecutive elements (first the real part, then the imaginary part).
 * @param b Pointer to the first element of the input array 'b', structured identically to 'a'.
 * @param c Pointer to the first element of the output array where the results are stored.
 *
 * @note Assumes that arrays 'a', 'b', and 'c' are pre-allocated with sufficient size to hold the necessary number of elements, and that 'nida[nb]' and 'nidb[nb]' correctly reflect the counts of complex elements to be processed. This function is optimized for scenarios where only a specific component of the complex product is needed, enhancing computational efficiency in such cases.
 */
void PGrid::cct_pack_iconjgMulb(const int nb, const double *a, const double *b, double *c)
{
   for (auto i=0; i<(nida[nb]+nidb[nb]); ++i)
      c[i] = a[2*i+1]*b[2*i] - a[2*i]*b[2*i+1];
}

/**********************************
 *                                *
 *    PGrid::regenerate_r_grid    *
 *                                *
 **********************************/
void PGrid::regenerate_r_grid() {
  int nxh = nx / 2;
  int nyh = ny / 2;
  int nzh = nz / 2;
  double a[9];
  for (auto i = 0; i < 3; ++i) {
    a[i] = lattice->unita1d(0 + i) / ((double)nx);
    a[3 + i] = lattice->unita1d(3 + i) / ((double)ny);
    a[6 + i] = lattice->unita1d(6 + i) / ((double)nz);
  }

  r_nzero(3, r_grid);

  /* grid points in coordination space */
  for (auto k3 = (-nzh); k3 < nzh; ++k3)
    for (auto k2 = (-nyh); k2 < nyh; ++k2)
      for (auto k1 = (-nxh); k1 < nxh; ++k1) {
        int i = k1 + nxh;
        int j = k2 + nyh;
        int k = k3 + nzh;
        int indx = ijktoindex2(i, j, k);
        int p = ijktop2(i, j, k);

        if (p == parall->taskid_i()) {
          r_grid[3 * indx] = a[0] * k1 + a[3] * k2 + a[6] * k3;
          r_grid[3 * indx + 1] = a[1] * k1 + a[4] * k2 + a[7] * k3;
          r_grid[3 * indx + 2] = a[2] * k1 + a[5] * k2 + a[8] * k3;
        }
      }
}

/************************************
 *                                  *
 *    PGrid::generate_r_sym_grid    *
 *                                  *
 ************************************/
void PGrid::generate_r_sym_grid(double *r_sym_grid) {
  int nxh = nx / 2;
  int nyh = ny / 2;
  int nzh = nz / 2;
  double a[9];
  for (auto i = 0; i < 3; ++i) {
    a[i] = lattice->unita1d(0 + i) / ((double)nx);
    a[3 + i] = lattice->unita1d(3 + i) / ((double)ny);
    a[6 + i] = lattice->unita1d(6 + i) / ((double)nz);
  }

  r_nzero(3, r_sym_grid);

  /* grid points in coordination space */
  for (auto k3 = (-nzh + 1); k3 < nzh; ++k3)
    for (auto k2 = (-nyh + 1); k2 < nyh; ++k2)
      for (auto k1 = (-nxh + 1); k1 < nxh; ++k1) {
        int i = k1 + nxh;
        int j = k2 + nyh;
        int k = k3 + nzh;
        int indx = ijktoindex2(i, j, k);
        int p = ijktop2(i, j, k);

        if (p == parall->taskid_i()) {
          r_sym_grid[3 * indx] = a[0] * k1 + a[3] * k2 + a[6] * k3;
          r_sym_grid[3 * indx + 1] = a[1] * k1 + a[4] * k2 + a[7] * k3;
          r_sym_grid[3 * indx + 2] = a[2] * k1 + a[5] * k2 + a[8] * k3;
        }
      }
}

/************************************
 *                                  *
 *    PGrid::generate_r_sym_mask    *
 *                                  *
 ************************************/
void PGrid::generate_r_sym_mask(double *rmask) {
  int nxh = nx / 2;
  int nyh = ny / 2;
  int nzh = nz / 2;
  r_zero(rmask);

  /* grid points in coordination space */
  for (auto k3 = (-nzh); k3 < nzh; ++k3)
    for (auto k2 = (-nyh); k2 < nyh; ++k2)
      for (auto k1 = (-nxh); k1 < nxh; ++k1) {
        int i = k1 + nxh;
        int j = k2 + nyh;
        int k = k3 + nzh;
        int indx = ijktoindex2(i, j, k);
        int p = ijktop2(i, j, k);

        if (p == parall->taskid_i())
          rmask[indx] = 1.0;
      }
}

/************************************
 *                                  *
 *       PGrid::c_Laplacian         *
 *                                  *
 ************************************/
void PGrid::c_Laplacian(const int nb, double *w) {
  int npack0 = this->npack(nb);
  double *Gx = this->Gpackxyz(nb, 0);
  double *Gy = this->Gpackxyz(nb, 1);
  double *Gz = this->Gpackxyz(nb, 2);
  int kk = 0;
  for (auto k = 0; k < npack0; ++k) {
    auto gg = Gx[k] * Gx[k] + Gy[k] * Gy[k] + Gz[k] * Gz[k];
    w[kk] *= (-gg);
    w[kk + 1] *= (-gg);
    kk += 2;
  }
}

/************************************
 *                                  *
 *       PGrid::cc_Laplacian        *
 *                                  *
 ************************************/
void PGrid::cc_Laplacian(const int nb, const double *w0, double *w) {
  int npack0 = this->npack(nb);
  double *Gx = this->Gpackxyz(nb, 0);
  double *Gy = this->Gpackxyz(nb, 1);
  double *Gz = this->Gpackxyz(nb, 2);
  int kk = 0;
  for (auto k = 0; k < npack0; ++k) {
    auto gg = Gx[k] * Gx[k] + Gy[k] * Gy[k] + Gz[k] * Gz[k];
    w[kk] = -w0[kk] * gg;
    w[kk + 1] = -w0[kk + 1] * gg;
    kk += 2;
  }
}

/************************************
 *                                  *
 *      PGrid::rr_Laplacian         *
 *                                  *
 ************************************/
void PGrid::rr_Laplacian(const int nb, const double *w0, double *w) {
  double scal1 = 1.0 / ((double)((nx) * (ny) * (nz)));

  rr_SMul(scal1, w0, w);
  this->rc_pfft3f(nb, w);
  this->c_pack(nb, w);

  this->c_Laplacian(nb, w);

  this->c_unpack(nb, w);
  this->cr_pfft3b(nb, w);
}

/************************************
 *                                  *
 *       PGrid::rr_Helmholtz        *
 *                                  *
 ************************************/
void PGrid::rr_Helmholtz(const int nb, const double *k2, const double *w,
                         double *Hw) {
  this->rr_Laplacian(nb, w, Hw);
  this->rrr_Mul2Add(k2, w, Hw);
}

/************************************
 *                                  *
 *   PGrid::rrr_solve_Helmholtz     *
 *                                  *
 ************************************/
/* The routine solves the inhomogeneous Helmholtz equation

 ∇^2 w(r) + k2(r)w(r) = b(r)  using CG.

 Entry - nb: 0-density grid, 1-wvfnc grid
         k2(r): wavenumber function
         b(r): source term

 Entry/Exit - w0(r) : initial guess / output solution

*/
void PGrid::rrr_solve_Helmholtz(const int nb, const double *k2, const double *b,
                                double *w) {
  double alpha, alpha0;
  double delta = 1.0;
  double *Aw = r_alloc();
  double *R = r_alloc();
  double *HR = r_alloc();

  double omega = this->lattice->omega();
  double scal1 = 1.0 / ((double)((nx) * (ny) * (nz)));
  double dv = omega * scal1;

  this->r_zero(w);

  int it = 0;
  while ((it < 10000) && (delta > 0.01)) {
    // compute the Aw and the residual
    this->rr_Helmholtz(nb, k2, w, Aw);
    this->rrr_Minus(b, Aw, R);

    this->rr_Helmholtz(nb, k2, R, HR);
    delta = rr_dot(R, R) * dv;
    // alpha = -1.0e-6;
    alpha0 = delta / rr_dot(R, HR);
    alpha = 1.0e-5 * alpha0;

    std::cout << "it=" << it << " delta=" << delta << " alpha0=" << alpha0
              << " alpha=" << alpha << std::endl;

    rr_daxpy(alpha, R, w);
    ++it;
  }

  r_dealloc(HR);
  r_dealloc(R);
  r_dealloc(Aw);
}

/************************************
 *                                  *
 *   PGrid::rrrr_FD_gradient        *
 *                                  *
 ************************************/
void PGrid::rrrr_FD_gradient(const double *rho, double *rhox, double *rhoy,
                             double *rhoz) {
  double ua[9];
  for (auto i = 0; i < 3; ++i) {
    ua[i] = lattice->unita1d(0 + i) / ((double)nx);
    ua[3 + i] = lattice->unita1d(3 + i) / ((double)ny);
    ua[6 + i] = lattice->unita1d(6 + i) / ((double)nz);
  }
  double dx = std::sqrt(ua[0] * ua[0] + ua[1] * ua[1] + ua[2] * ua[2]);
  double dy = std::sqrt(ua[3] * ua[3] + ua[4] * ua[4] + ua[5] * ua[5]);
  double dz = std::sqrt(ua[6] * ua[6] + ua[7] * ua[7] + ua[8] * ua[8]);

  this->rrrr_periodic_gradient(rho, rhox, rhoy, rhoz);
  for (auto i = 0; i < n2ft3d; ++i) {
    rhox[i] /= dx;
    rhoy[i] /= dy;
    rhoz[i] /= dz;
  }
}

/************************************
 *                                  *
 *   PGrid::rrrr_FD_laplacian       *
 *                                  *
 ************************************/
void PGrid::rrrr_FD_laplacian(const double *rho, double *rhoxx, double *rhoyy,
                              double *rhozz) {
  double ua[9];
  for (auto i = 0; i < 3; ++i) {
    ua[i] = lattice->unita1d(0 + i) / ((double)nx);
    ua[3 + i] = lattice->unita1d(3 + i) / ((double)ny);
    ua[6 + i] = lattice->unita1d(6 + i) / ((double)nz);
  }
  double dxx = (ua[0] * ua[0] + ua[1] * ua[1] + ua[2] * ua[2]);
  double dyy = (ua[3] * ua[3] + ua[4] * ua[4] + ua[5] * ua[5]);
  double dzz = (ua[6] * ua[6] + ua[7] * ua[7] + ua[8] * ua[8]);

  this->rrrr_periodic_laplacian(rho, rhoxx, rhoyy, rhozz);
  for (auto i = 0; i < n2ft3d; ++i) {
    rhoxx[i] /= (dxx);
    rhoyy[i] /= (dyy);
    rhozz[i] /= (dzz);
  }
}

} // namespace pwdft

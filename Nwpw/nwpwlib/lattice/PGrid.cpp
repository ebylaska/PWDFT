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

void print_tmp(int taskid, double *tmp1)
{     
   
   std::cout << "taskid = " << taskid;
   std::cout << " TMP=";
   for (auto k=0; k<32; ++k)
      std::cout << tmp1[k] << " " ;
   std::cout << std::endl;
}  

/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

// PGrid::PGrid(Parallel *inparall, Lattice *inlattice, Control2& control) :
// d3db(inparall,control.mapping(),control.ngrid(0),control.ngrid(1),control.ngrid(2))

PGrid::PGrid(Parallel *inparall, Lattice *inlattice, int mapping0, int balance0,
             int nx0, int ny0, int nz0, int pfft3_qsize0, bool staged_gpu_fft_pipeline0, int fft_container_size0)
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

   nffts_max = fft_container_size0;
 
   zplane_tmp1 = new (std::nothrow) double[ (nffts_max*2*zplane_size+8 + Alignment - 1) & ~(Alignment - 1)];
   zplane_tmp2 = new (std::nothrow) double[ (nffts_max*2*zplane_size+8 + Alignment - 1) & ~(Alignment - 1)];
 
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
   aqindx   = new (std::nothrow) int[aqmax]();
   aqstatus = new (std::nothrow) int[aqmax]();
   aqnffts  = new (std::nothrow) int[aqmax]();
   //atmp = new (std::nothrow) double[2*aqmax*n2ft3d]();
   atmp = new (std::nothrow) double[2*aqmax*n2ft3d*nffts_max]();
   //std::vector<double> atmp_vector(2*aqmax*n2ft3d*nffts_max);
   //atmp = atmp_vector.data();



   bqmax = pfft3_qsize0;
   if (staged_gpu_fft_pipeline) bqmax += 6;
   
   //bqmax = aqmax;
   bqsize = 0;
   blast_index = bqmax - 1;
   bqindx   = new (std::nothrow) int[bqmax]();
   bqstatus = new (std::nothrow) int[bqmax]();
   bqnffts  = new (std::nothrow) int[bqmax]();
   //btmp = new (std::nothrow) double[2*bqmax*n2ft3d]();
   btmp = new (std::nothrow) double[2*bqmax*n2ft3d*nffts_max]();
   //std::vector<double> btmp_vector(2*bqmax*n2ft3d*nffts_max);
   //btmp = btmp_vector.data();

 
   /* initialize async buffer data for pfft */
   for (auto q=0; q<aqmax; ++q)
      parall->astart(4+q, 2*parall->np_i()+1);

}

PGrid::PGrid(Parallel *inparall, Lattice *inlattice, Control2 &control)
    : PGrid(inparall, inlattice, control.mapping(), control.balance(),
            control.ngrid(0), control.ngrid(1), control.ngrid(2),
            control.pfft3_qsize(), control.staged_gpu_fft(), control.fft_container_size()) {}

/********************************
 *                              *
 *       PGrid:c_unpack         *
 *                              *
 ********************************/
/**
 * @brief Unpacks FFT data, applies balancing, and performs post-processing transformations.
 *
 * This function unpacks FFT data from a compact format into a full array with optional data redistribution
 * if balancing is required. It clears the input array and refills it using indexed data from a global temporary array.
 * Post unpacking, it applies a time reversal transformation for data normalization.
 *
 * @param nb Index used for accessing specific configuration arrays and determining the scope of the operation.
 * @param a Pointer to the data array containing packed FFT data.
 */
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
/**
 * @brief Packs FFT data into a compact format and applies optional load balancing.
 *
 * This function copies FFT data into temporary storage, clears the original data array, and reorganizes it
 * according to predefined indices. It is typically used to prepare data for efficient FFT processing or to
 * reorganize data post-FFT. Load balancing may be applied conditionally to optimize data distribution in
 * parallel computing environments.
 *
 * @param nb The index used to access specific configuration and indexing arrays.
 * @param a Pointer to the data array to be packed.
 */
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
/**
 * @brief Copies a specified number of elements from one array to another based on indexed configuration.
 *
 * This function performs a direct memory copy of elements from source array `a` to destination array `b`.
 * The number of elements copied is determined by configuration indices specific to the provided index `nb`.
 * This is typically used for handling subsets of data in large-scale numerical computations.
 *
 * @param nb Index used to determine the number of elements to copy based on application-specific configuration arrays.
 * @param a Pointer to the source data array.
 * @param b Pointer to the destination data array.
 */
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
/**
 * @brief Computes a modified dot product of two data arrays based on configuration indices.
 *
 * This function calculates a weighted sum of the dot products of two arrays, specifically designed
 * to handle complex number arrays in grid-based computations. The function computes the dot product 
 * of the entire configured range, scales it, and subtracts the dot product of a subset of that range.
 *
 * @param nb Index used to access configuration-specific indices that determine the range of data processed.
 * @param a Pointer to the first array.
 * @param b Pointer to the second array.
 * @return double The calculated modified dot product.
 */
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
/**
 * @brief Packs and computes the dot product of projections using computeTrans_Mult.
 *
 * This function utilizes the `computeTrans_Mult` function to perform transformations
 * on projection data. It calculates the dot product of projections (`a`, `b`) and 
 * stores the results in `sum`. After the initial transformation, it aggregates the
 * results across all processes using the `Vector_SumAll` function. This function is used 
 * by non-local pseudopotential routines.
 *
 * @param nb   Index of the data partition.
 * @param nn   Number of functions (e.g., `psi` functions) to process.
 * @param nprj Number of projections to process.
 * @param a    Pointer to the first projection data array.
 * @param b    Pointer to the second projection data array.
 * @param sum  Array to store the computed dot product results.
 *
 * Temporary Variables:
 * - `ng`: Total grid size, calculated as the sum of `nida[nb]` and `nidb[nb]`.
 * - `ng0`: Size of the first sub-grid, based on `nida[nb]`.
 * - `rtwo`: Scalar value set to 2.0, used as a scaling factor in matrix multiplication.
 * - `rone`: Scalar value set to 1.0, typically used to preserve computed results.
 * - `rmone`: Scalar value set to -1.0, used to negate values during matrix operations.
 * - `rzero`: Scalar value set to 0.0, used to initialize or reset matrices or results.
 *
 */
/********************************
 *                              *
 *    PGrid:cc_pack_inprjdot    *
 *                              *
 ********************************/
void PGrid::cc_pack_inprjdot(const int nb, int nn, int nprj, double *a, double *b, double *sum) 
{
   int ng = 2 * (nida[nb] + nidb[nb]);
   int ng0 = 2 * nida[nb];
 
   int one = 1;
   double rtwo = 2.0;
   double rone = 1.0;
   double rmone = -1.0;
   double rzero = 0.0;
 
   d3db::mygdevice.TN_dgemm(nn, nprj, ng, rtwo, a, b, rzero, sum);
   if (ng0 > 0) 
   {
      DGEMM_PWDFT((char *)"T", (char *)"N", nn, nprj, ng0, rmone, a, ng, b, ng, rone, sum, nn);
   }

   // Aggregate the results across all processes
   parall->Vector_SumAll(1, nn*nprj, sum);
}


/********************************
 *                              *
 *   PGrid:n2ccttt_pack_i3ndot  *
 *                              *
 ********************************/
 /**
 * @brief Computes a series of dot products between transformed vectors with parallelization.
 *
 * This function performs a set of dot product calculations on complex numbers stored in the
 * arrays `psi` and `prj`, with results stored in `sum3`. The input arrays `Gx`, `Gy`, and `Gz`
 * represent the data against which the transformed vectors are compared. The function is optimized
 * for parallel execution, allowing the workload to be distributed across multiple cores.
 *
 * @param nb An integer index used to determine the range of elements from the global arrays, influencing the size of the operations.
 * @param nn The number of unique psi arrays to be processed.
 * @param nprj The number of projection matrices or elements to be considered.
 * @param psi Pointer to the first element of the psi array, which stores data for multiplication.
 * @param prj Pointer to the first element of the prj array, which contains projection matrices or data.
 * @param Gx Pointer to the first element of the Gx array, representing transformation data.
 * @param Gy Pointer to the first element of the Gy array, representing transformation data.
 * @param Gz Pointer to the first element of the Gz array, representing transformation data.
 * @param sum3 Pointer to the first element of the output array where the results are accumulated.
 *
 * @note This function can be optimized for AVX2 instruction sets, ARM NEON intrinsics, OpenMP, and GPUs for more efficient computation.
 */
void PGrid::n2ccttt_pack_i3ndot(const int nb, const int nn, const int nprj, 
                                const double *psi,
                                const double *prj,
                                double *Gx, double *Gy, double *Gz,
                                double *sum3)
{
   double *xtmp1 = d3db::d3db_tmp1;

   int one = 1;
   int ng = (nida[nb] + nidb[nb]);
   int ng0 = nida[nb];
   int nshift = 2*ng;
   int count3 = 0;


   //d3db::mygdevice.computeTrans3_Mult(ne,nprj,psi,prj,ng,ng0,Gx,Gy,Gz,xtmp1,sum3)

   for (auto l=0; l<nprj; ++l)
   for (auto n=0; n<nn; ++n)
   {
      //cct_pack_iconjgMul(1, prj+l*ng, psi + n*ng, xtmp1);
      //sum3[3*n + (3*nn*nprj)*l]     = tt_pack_idot(1,Gx, xtmp1);
      //sum3[3*n + (3*nn*nprj)*l + 1] = tt_pack_idot(1,Gy, xtmp1);
      //sum3[3*n + (3*nn*nprj)*l + 2] = tt_pack_idot(1,Gz, xtmp1);

      // Perform cct_pack_iconjgMul
      const double *a = prj + l*nshift;
      const double *b = psi + n*nshift;
      for (int i = 0; i < ng; ++i) 
         xtmp1[i] = a[2*i]*b[2*i+1] - a[2*i+1]*b[2*i];
      
      double tsumx = 2.0*DDOT_PWDFT(ng,Gx,one,xtmp1,one);
      double tsumy = 2.0*DDOT_PWDFT(ng,Gy,one,xtmp1,one);
      double tsumz = 2.0*DDOT_PWDFT(ng,Gz,one,xtmp1,one);
      tsumx -= DDOT_PWDFT(ng0,Gx,one,xtmp1,one);
      tsumy -= DDOT_PWDFT(ng0,Gy,one,xtmp1,one);
      tsumz -= DDOT_PWDFT(ng0,Gz,one,xtmp1,one);

      sum3[count3]   = tsumx;
      sum3[count3+1] = tsumy;
      sum3[count3+2] = tsumz;
      count3 += 3;
   }

   parall->Vector_SumAll(1,3*nn*nprj,sum3);
}


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
/**
 * @brief Pack a real-space integer array `a` into a compact, 
 * load-balanced layout for FFT or matvec operations.
 *
 * This function extracts and rearranges a subset of elements from
 * the input array `a` based on a precomputed index mapping stored 
 * in `packarray[nb]`. The selected elements are first gathered into 
 * a temporary buffer (`tmp`) and then packed into the correct positions 
 * in `a` for subsequent parallel operations.
 *
 * If the `balanced` flag is set, the packed array is further permuted 
 * by `mybalance->i_balance()` to ensure a load-balanced distribution 
 * across processors or GPUs.
 *
 * @param[in]  nb   Index of the packing block (e.g., orbital or FFT slice).
 * @param[in,out] a Integer array of length `nfft3d` that will be overwritten 
 *                  with its packed and optionally balanced version.
 *
 * @note 
 *  - Original data in `a` is overwritten.
 *  - Packing layout is determined by `packarray[nb]` and potentially 
 *    modified by the balancing map `mybalance`.
 *  - Used in routines that expect compact, balanced G-vector layouts
 *    (e.g., 3D FFTs, matrix-vector products, and data exchange).
 */
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
/**
 * @brief Copy a packed segment of integer data from array `a` to array `b`.
 *
 * This routine copies `ng = nida[nb] + nidb[nb]` integers from the beginning 
 * of the packed input array `a` into output array `b`. The packed form is 
 * typically created by first zero-initializing a full-size array (length `nfft3d`), 
 * then applying `PGrid::i_pack()` to retain only data relevant to the index block `nb`.
 *
 * @param[in]  nb   Index block ID indicating which packed segment to copy.
 * @param[in]  a    Input array of size `nfft3d`, assumed to have been packed beforehand.
 * @param[out] b    Output array of size `ng = nida[nb] + nidb[nb]`, receiving packed values.
 *
 * @note
 * - This is a direct memory copy with no additional index mapping or balancing.
 * - Especially useful for truncating and transferring packed integer index arrays
 *   after sparse masking or FFT data preparation.
 * - Enables efficient redistribution of index maps across processor grids or FFT tasks.
 */

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
      d3db::mygdevice.batch_cfftz_tmpz_zero(d3db::fft_tag,false, nz, nq3, 1, n2ft3d, a, d3db::tmpz, zero_row3[nb]);
      //d3db::mygdevice.batch_cfftz_tmpz(d3db::fft_tag,false, nz, nq3, n2ft3d, a, d3db::tmpz );
     
      d3db::c_ptranspose_ijk(nb, 2, a, tmp2, tmp3);
      
      /************************************************
       ***     do fft along ky dimension            ***
       ***   A(ny,nz,kx) <- fft1d^(-1)[A(ky,nz,kx)] ***
       ************************************************/
      d3db::mygdevice.batch_cffty_tmpy_zero(d3db::fft_tag,false,ny,nq2,1,n2ft3d,a,d3db::tmpy,zero_row2[nb]);
      
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
      d3db::mygdevice.batch_cffty_tmpy_zero(d3db::fft_tag,true,ny,nq2,1,n2ft3d,a,d3db::tmpy,zero_row2[nb]);
     
      d3db::c_ptranspose_ijk(nb, 1, a, tmp2, tmp3);
     
      /********************************************
       ***     do fft along nz dimension        ***
       ***   A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]  ***
       ********************************************/
      d3db::mygdevice.batch_cfftz_tmpz_zero(d3db::fft_tag,true, nz, nq3,1,n2ft3d, a, d3db::tmpz, zero_row3[nb]);
      //d3db::mygdevice.batch_cfftz_tmpz_zero(d3db::fft_tag,true, nz, nq3, n2ft3d, a, d3db::tmpz );
   }
 
   //delete[] tmp3;
   //delete[] tmp2;
}

/********************************
 *                              *
 *     PGrid::c_unpack_start    *
 *                              *
 ********************************/
void PGrid::c_unpack_start(const int nffts, const int nb, double *tmp1, double *tmp2,
                           const int request_indx, const int msgtype) 
{
   if (balanced)
   {
      for (auto s=0; s<nffts; ++s)
         mybalance->c_unbalance_start(1, 1, tmp1+s*n2ft3d, request_indx, msgtype);
   }
}

/********************************
 *                              *
 *     PGrid::c_unpack_mid      *
 *                              *
 ********************************/
void PGrid::c_unpack_mid(const int nffts, const int nb, double *tmp1, double *tmp2,
                         const int request_indx, const int msgtype) 
{
   if (balanced)
      mybalance->c_unbalance_end(nffts, nb, tmp1, request_indx);
 
   for (auto s=0; s<nffts; ++s)
      std::memcpy(tmp2 + s*n2ft3d, tmp1 + s*n2ft3d, 2*(nida[nb]+nidb2[nb])*sizeof(double));
 
   std::memset(tmp1, 0, nffts*n2ft3d*sizeof(double));
   for (auto s=0; s<nffts; ++s)
      c_bindexcopy((nida[nb]+nidb2[nb]),packarray[nb], tmp2 + s*n2ft3d, tmp1 + s*n2ft3d);
   // c_bindexcopy(nida[nb]+nidb[nb],packarray[nb],tmp2,tmp1);
   
 
   d3db::c_timereverse_start(nffts, tmp1, zplane_tmp1, zplane_tmp2, request_indx, msgtype);
   //for (auto s=0; s<nffts; ++s)
  // {
  //    d3db::c_timereverse_start(1, tmp1+s*n2ft3d, zplane_tmp1, zplane_tmp2, request_indx, msgtype);
  //    d3db::c_timereverse_end(1, tmp1+s*n2ft3d, zplane_tmp1, zplane_tmp2, request_indx);
  // }
}

/********************************
 *                              *
 *     PGrid::c_unpack_end      *
 *                              *
 ********************************/
void PGrid::c_unpack_end(const int nffts, const int nb, double *tmp1, double *tmp2,
                         const int request_indx) 
{
   d3db::c_timereverse_end(nffts, tmp1, zplane_tmp1, zplane_tmp2, request_indx);

}

/********************************
 *                              *
 *        PGrid::pfftbz         *
 *                              *
 ********************************/
void PGrid::pfftbz(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx) 
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
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nxh; ++i) 
            {
               if (!zero_row3[nb][indx2]) 
               {
                  auto kk = 0;
                  auto indx3 = 2*i + indx0;
                  auto shift = 2*nz*nn;
                  for (auto k=0; k<nz; ++k) 
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
      }
     
      d3db::mygdevice.batch_cfftz_tmpz(d3db::fft_tag,false, nz, nn, n2ft3d, tmp2, d3db::tmpz);
      // for (auto i=0; i<nn; ++i)
      //    dcfftb_(&nz,tmp2+2*nz*i,d3db::tmpz);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nxh; ++i) 
            {
               if (!zero_row3[nb][indx2]) 
               {
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
      }
     
      /***********************************************
       ***         Do a ptranspose of A            ***
       ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
       ************************************************/
      d3db::c_ptranspose1_jk_start(nffts, nb, tmp1, tmp2, tmp1, request_indx, 44);
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
      d3db::mygdevice.batch_cfftz_tmpz_zero(d3db::fft_tag,false, nz, nq3,nffts, n2ft3d, tmp1, d3db::tmpz, zero_row3[nb]);
      //d3db::mygdevice.batch_cfftz_tmpz(d3db::fft_tag, false, nz, nffts*nq3, n2ft3d, tmp1, d3db::tmpz);

      // GOOD
      d3db::c_ptranspose_ijk_start(nffts, nb, 2, tmp1, tmp2, tmp1, request_indx, 45);
      // BAD 
      // d3db::c_ptranspose_ijk(nb,2,tmp1,tmp2,tmp1);
   }
}

/********************************
 *                              *
 *        PGrid::pfftby         *
 *                              *
 ********************************/
void PGrid::pfftby(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx) 
{
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      auto nxh = nx/2 + 1;
      auto nxh2 = nx + 2;
      auto nxhy2 = nxh2 * ny;
     
      /***********************************************
       ***         Do a ptranspose of A            ***
       ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
       ************************************************/
      d3db::c_ptranspose1_jk_end(nffts, nb, tmp2, tmp1, request_indx);
     
      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
       ********************************************/
      int indx0 = 0;
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
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
      }
     
      d3db::mygdevice.batch_cffty_tmpy(d3db::fft_tag,false, ny, nn, n2ft3d, tmp1, d3db::tmpy);
      // for (auto i=0; i<nn; ++i)
      //    dcfftb_(&ny,tmp1+2*ny*i,d3db::tmpy);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nxh; ++i) 
            {
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
   }

   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      d3db::c_ptranspose_ijk_end(nffts, nb, 2, tmp2, tmp1, request_indx);
     
      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
       ********************************************/
      d3db::mygdevice.batch_cffty_tmpy_zero(d3db::fft_tag,false,ny,nq2,nffts,n2ft3d,tmp2,d3db::tmpy,zero_row2[nb]);
      //d3db::mygdevice.batch_cffty_tmpy(d3db::fft_tag,false,ny,nffts*nq2,n2ft3d,tmp2,d3db::tmpy);

     
      d3db::c_ptranspose_ijk_start(nffts, nb, 3, tmp2, tmp1, tmp2, request_indx, 46);
      //for (s=0; s<nffts; ++s)
      //   d3db::c_ptranspose_ijk(nb,3,tmp2+s*n2ft3d,tmp1+s*n2ft3d,tmp2+s*n2ft3d);
   }
}

/********************************
 *                              *
 *        PGrid::pfftbx         *
 *                              *
 ********************************/
void PGrid::pfftbx(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx) 
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
      d3db::mygdevice.batch_rfftx_tmpx(d3db::fft_tag,false,nx,nffts*ny*nq,n2ft3d,tmp2,d3db::tmpx);

      d3db::zeroend_fftb(nx,ny,nq,nffts,tmp2);
      std::memcpy(tmp1,tmp2,nffts*n2ft3d*sizeof(double));
   }

   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      d3db::c_ptranspose_ijk_end(nffts,nb, 3, tmp1, tmp2, request_indx);

     
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      d3db::mygdevice.batch_rfftx_tmpx(d3db::fft_tag,false, nx, nffts*nq1, n2ft3d, tmp1, d3db::tmpx);
      //for (auto s=0; s<nffts; ++s)
      //   d3db::mygdevice.batch_rfftx_tmpx(d3db::fft_tag,false, nx, nq1, n2ft3d, tmp1+s*n2ft3d, d3db::tmpx);


      for (auto s=0; s<nffts; ++s)
      d3db::zeroend_fftb(nx, nq1, 1, 1, tmp1+s*n2ft3d);

      if (n2ft3d_map < n2ft3d)
         for (auto s=0; s<nffts; ++s)
            std::memset(tmp1+s*n2ft3d + n2ft3d_map,0,(n2ft3d-n2ft3d_map)*sizeof(double));
   }
}

/********************************
 *                              *
 *       PGrid::pfftb_step      *
 *                              *
 ********************************/
void PGrid::pfftb_step(const int step, const int nffts, const int nb, double *a, double *tmp1,
                       double *tmp2, const int request_indx) 
{
   if (step == 0) 
   {
      // parall->astart(request_indx,parall->np_i());
      std::memset(tmp1,0,nffts*n2ft3d*sizeof(double));
      // unpack start, tmp1-->tmp1
      for (auto s=0; s<nffts; ++s)
         std::memcpy(tmp1 + s*n2ft3d, a + s*2*(nida[nb]+nidb[nb]), 2*(nida[nb]+nidb[nb])*sizeof(double));

      d3db:parall->comm_Barrier(1);
      this->c_unpack_start(nffts, nb, tmp1, tmp2, request_indx, 47);

   } 
   else if (step == 1) 
   {
      parall->comm_Barrier(1);
      // unpack mid
      this->c_unpack_mid(nffts, nb, tmp1, tmp2, request_indx, 48);

   } 
   else if (step == 2) 
   {
      parall->comm_Barrier(1);
      // unpack end; mem-->dev,  out=tmp1
      this->c_unpack_end(nffts, nb, tmp1, tmp2, request_indx);
   } 
   else if (step == 3) 
   {
      parall->comm_Barrier(1);
      // pfftbz dev-->dev->mem,  tmp1->tmp1
      this->pfftbz(nffts, nb, tmp1, tmp2, request_indx);
   } 
   else if (step == 4) 
   {
      parall->comm_Barrier(1);
      // pfftby mem->dev-->dev->mem
      // in=tmp1, tmp2->tmp1, tmp1=in , tmp2=tmp
      pfftby(nffts, nb, tmp1, tmp2, request_indx);
   } 
   else if (step == 5) 
   {
      parall->comm_Barrier(1);
      // pfftbx mem->dev->dev->mem
      pfftbx(nffts, nb, tmp1, tmp2, request_indx);
      // parall->aend(request_indx);

   }
}

/********************************
 *                              *
 *      PGrid::pfftbz_start     *
 *                              *
 ********************************/
void PGrid::pfftbz_start(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
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
     int nn = 0;
     for (auto s=0; s<nffts; ++s)
     {
        int indx2 = 0;
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
     }
    
     d3db::mygdevice.batch_cfftz_stages_tmpz(0,d3db::fft_tag,false, nz, nffts*nn, n2ft3d, tmp2, d3db::tmpz,da_indx);
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
     d3db::mygdevice.batch_cfftz_stages_tmpz_zero(0,d3db::fft_tag,false, nz, nq3,nffts, n2ft3d, tmp1, d3db::tmpz, zero_row3[nb], da_indx);
  }

}

/********************************
 *                              *
 *    PGrid::pfftbz_compute     *
 *                              *
 ********************************/
void PGrid::pfftbz_compute(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
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
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i = 0; i < nxh; ++i) 
            {
               if (!zero_row3[nb][indx2]) 
               {
                  nn += 1;
               }
               ++indx2;
            }
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
      d3db::mygdevice.batch_cfftz_stages_tmpz_zero(1,d3db::fft_tag,false, nz, nq3,nffts, n2ft3d, tmp1, d3db::tmpz,zero_row3[nb],da_indx);
   }
}

/********************************
 *                              *
 *      PGrid::pfftbz_end       *
 *                              *
 ********************************/
void PGrid::pfftbz_end(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
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
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
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
      }
     
      d3db::mygdevice.batch_cfftz_stages_tmpz(2,d3db::fft_tag,false, nz, nn, n2ft3d, tmp2, d3db::tmpz, da_indx);
      // for (auto i=0; i<nn; ++i)
      //    dcfftb_(&nz,tmp2+2*nz*i,d3db::tmpz);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
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
      }
     
      /***********************************************
       ***         Do a ptranspose of A            ***
       ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
       ************************************************/
      d3db::c_ptranspose1_jk_start(nffts,nb,tmp1,tmp2,tmp1,request_indx, 44);
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
      d3db::mygdevice.batch_cfftz_stages_tmpz_zero(2,d3db::fft_tag,false, nz, nq3,nffts, n2ft3d, tmp1, d3db::tmpz, zero_row3[nb], da_indx);
     
      d3db::c_ptranspose_ijk_start(nffts,nb,2,tmp1,tmp2,tmp1,request_indx, 45);
      // d3db::c_ptranspose_ijk(nb,2,tmp1,tmp2,tmp1);
   }
}



/********************************
 *                              *
 *      PGrid::pfftby_start     *
 *                              *
 ********************************/
void PGrid::pfftby_start(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      auto nxh = nx/2 + 1;
      auto nxh2 = nx + 2;
      auto nxhy2 = nxh2 * ny;
     
      /***********************************************
       ***         Do a ptranspose of A            ***
       ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
       ************************************************/
      d3db::c_ptranspose1_jk_end(nffts, nb, tmp2, tmp1, request_indx);
     
      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
       ********************************************/
      int indx0 = 0;
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nxh; ++i) 
            {
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
      d3db::c_ptranspose_ijk_end(nffts,nb, 2, tmp2, tmp1, request_indx);
     
      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
       ********************************************/
      d3db::mygdevice.batch_cffty_stages_tmpy_zero(0,d3db::fft_tag,false,ny,nq2,nffts,n2ft3d,tmp2,d3db::tmpy,zero_row2[nb],da_indx);
   }
}


/********************************
 *                              *
 *    PGrid::pfftby_compute     *
 *                              *
 ********************************/
void PGrid::pfftby_compute(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      auto nxh = nx/2 + 1;
      auto nxh2 = nx + 2;
      auto nxhy2 = nxh2 * ny;
     
      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
       ********************************************/
      int indx0 = 0;
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nxh; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  nn += 1;
               }
               ++indx2;
            }
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
      d3db::mygdevice.batch_cffty_stages_tmpy_zero(1,d3db::fft_tag,false,ny,nq2,nffts,n2ft3d,tmp2,d3db::tmpy,zero_row2[nb],da_indx);
   }
}


/********************************
 *                              *
 *        PGrid::pfftby_end     *
 *                              *
 ********************************/
void PGrid::pfftby_end(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
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
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nxh; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  nn += 1;
               }
               ++indx2;
            }
         }
      }
     
      d3db::mygdevice.batch_cffty_stages_tmpy(2,d3db::fft_tag,false, ny, nn, n2ft3d, tmp1, d3db::tmpy,da_indx);
      // for (auto i=0; i<nn; ++i)
      //    dcfftb_(&ny,tmp1+2*ny*i,d3db::tmpy);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nxh; ++i) 
            {
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
      d3db::mygdevice.batch_cffty_stages_tmpy_zero(2,d3db::fft_tag,false,ny,nq2,nffts,n2ft3d,tmp2,d3db::tmpy,zero_row2[nb],da_indx);
 
      d3db::c_ptranspose_ijk_start(nffts, nb, 3, tmp2, tmp1, tmp2, request_indx, 46);
      // d3db::c_ptranspose_ijk(nb,3,tmp2,tmp1,tmp2);
   }
}




/********************************
 *                              *
 *      PGrid::pfftbx_start     *
 *                              *
 ********************************/
void PGrid::pfftbx_start(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
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
      d3db::mygdevice.batch_rfftx_stages_tmpx(0,d3db::fft_tag,false, nx, nffts*ny*nq, n2ft3d, tmp2, d3db::tmpx,da_indx);
   }
   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      d3db::c_ptranspose_ijk_end(nffts, nb, 3, tmp1, tmp2, request_indx);
     
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      d3db::mygdevice.batch_rfftx_stages_tmpx(0,d3db::fft_tag,false, nx, nffts*nq1, n2ft3d, tmp1, d3db::tmpx,da_indx);
   }
}


/********************************
 *                              *
 *    PGrid::pfftbx_compute     *
 *                              *
 ********************************/
void PGrid::pfftbx_compute(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
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
      d3db::mygdevice.batch_rfftx_stages_tmpx(1,d3db::fft_tag,false, nx, nffts*ny*nq, n2ft3d, tmp2, d3db::tmpx,da_indx);
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
      d3db::mygdevice.batch_rfftx_stages_tmpx(1,d3db::fft_tag,false, nx, nffts*nq1, n2ft3d, tmp1, d3db::tmpx,da_indx);
   }
}



/********************************
 *                              *
 *      PGrid::pfftbx_end       *
 *                              *
 ********************************/
void PGrid::pfftbx_end(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
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
      d3db::mygdevice.batch_rfftx_stages_tmpx(2,d3db::fft_tag,false, nx, nffts*ny*nq, n2ft3d, tmp2, d3db::tmpx,da_indx);
      d3db::zeroend_fftb(nx, ny, nq, nffts, tmp2);
      std::memcpy(tmp1,tmp2,nffts*n2ft3d*sizeof(double));
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
      d3db::mygdevice.batch_rfftx_stages_tmpx(2,d3db::fft_tag,false, nx, nffts*nq1, n2ft3d, tmp1, d3db::tmpx,da_indx);
      d3db::zeroend_fftb(nx, nq1, 1, nffts, tmp1);
      if (n2ft3d_map < n2ft3d)
         for (auto s=0; s<nffts; ++s)
            std::memset(tmp1+s*n2ft3d + n2ft3d_map, 0, (n2ft3d-n2ft3d_map)*sizeof(double));
   }
}




/********************************
 *                              *
 *       PGrid::pfftb_step12    *
 *                              *
 ********************************/
void PGrid::pfftb_step12(const int step, const int nffts, const int nb, double *a, double *tmp1,
                         double *tmp2, const int request_indx, const int indx)
{
   if (step == 0) 
   {
      // parall->astart(request_indx,parall->np_i());
      // unpack start, tmp1-->tmp1
      for (auto s=0; s<nffts; ++s)
         std::memcpy(tmp1+s*n2ft3d, a+s*n2ft3d, 2*(nida[nb]+nidb[nb])*sizeof(double));
      this->c_unpack_start(nffts, nb, tmp1, tmp2, request_indx, 47);
   } 
   else if (step == 1) 
   {
      // unpack mid
      this->c_unpack_mid(nffts, nb, tmp1, tmp2, request_indx, 48);
   } 
   else if (step == 2) 
   {
      // unpack end; mem-->dev,  out=tmp1
      this->c_unpack_end(nffts, nb, tmp1, tmp2, request_indx);
   } 
   // pfftbz dev-->dev->mem,  tmp1->tmp1
   else if (step == 3) 
   {
      this->pfftbz_start(nffts, nb, tmp1, tmp2, request_indx,indx);
   } 
   else if (step == 4) 
   {
      this->pfftbz_compute(nffts, nb, tmp1, tmp2, request_indx,indx);
   } 
   else if (step == 5) 
   {
      this->pfftbz_end(nffts, nb, tmp1, tmp2, request_indx,indx);
   } 
   // pfftby mem->dev-->dev->mem
   // in=tmp1, tmp2->tmp1, tmp1=in , tmp2=tmp
   else if (step == 6) 
   {
      pfftby_start(nffts, nb, tmp1, tmp2, request_indx,indx);
   } 
   else if (step == 7) 
   {
      pfftby_compute(nffts, nb, tmp1, tmp2, request_indx,indx);
   } 
   else if (step == 8) 
   {
      pfftby_end(nffts, nb, tmp1, tmp2, request_indx,indx);
   } 
   // pfftbx mem->dev->dev->mem
   else if (step == 9) 
   {
      pfftbx_start(nffts, nb, tmp1, tmp2, request_indx,indx);
   } 
   else if (step == 10) 
   {
      pfftbx_compute(nffts, nb, tmp1, tmp2, request_indx,indx);
   } 
   else if (step == 11) 
   {
      pfftbx_end(nffts, nb, tmp1, tmp2, request_indx,indx);
      // parall->aend(request_indx);
   }
}

/********************************
 *                              *
 *    PGrid:cr_pfft3b_queuein   *
 *                              *
 ********************************/
 /**
 * @brief Queues complex-to-real FFT operations, handling both GPU-accelerated and non-accelerated pipelines.
 *
 * This method processes a batch of FFT operations by updating operation statuses and managing a queue
 * of FFT tasks. Depending on the configuration, it either uses a GPU-accelerated pipeline or a standard
 * FFT process. It supports wrapping of queue indices and manages temporary storage for intermediate results.
 *
 * @param nb The fermi sphere packing used.
 * @param a Pointer to the data array on which FFT operations are performed.
 *
 * @note This function is designed for high-throughput FFT processing in parallel computing environments,
 *       particularly using GPUs. It is part of the PGrid class which encapsulates grid and FFT operation management
 *       in a parallelized scientific computation context.
 */
void PGrid::cr_pfft3b_queuein(const int nb, const int nffts_in, double *a) 
{
   //int nffts_in = 1;
   int shift1, shift2;
   int np = parall->np_i();
   //std::cout << "cr queuein HERA" << std::endl;
   //std::cout << "cr_queuein nb=" << nb << " nffts_in=" << nffts_in << std::endl;
   //std::cout << "       aqsize=" << aqsize << " nffts_max=" << nffts_max << std::endl;
   //std::cout << "       alast_index=" << alast_index << std::endl;
   //std::cout << "src ptr:" << a << std::endl;

   for (auto q=0; q<aqsize; ++q) 
   {
      int indx   = aqindx[q];
      int status = aqstatus[indx] + 1;
      int nffts  = aqnffts[indx];
      shift1 = nffts_max*n2ft3d*(2*indx);
      shift2 = nffts_max*n2ft3d*(2*indx+1);
      if (staged_gpu_fft_pipeline)
         pfftb_step12(status,nffts,nb,a,atmp+shift1,atmp+shift2,indx+4,indx);
      else
         pfftb_step(status,nffts,nb,a,atmp+shift1,atmp+shift2,indx+4);
      ++aqstatus[indx];
   }
 
   ++alast_index;
   if (alast_index >= aqmax)
      alast_index = 0;
   ++aqsize;
   aqindx[aqsize - 1] = alast_index;
   aqstatus[alast_index] = 0;
   aqnffts[alast_index]  = nffts_in;
 
   // status = 0;
   shift1 = nffts_max*n2ft3d*(2*alast_index);
   shift2 = nffts_max*n2ft3d*(2*alast_index+1);
 
   if (staged_gpu_fft_pipeline)
      pfftb_step12(0,nffts_in,nb,a,atmp+shift1,atmp+shift2, alast_index+4,alast_index);
   else
      pfftb_step(0,nffts_in,nb,a,atmp+shift1,atmp+shift2,alast_index+4);
}


/********************************
 *                              *
 *    PGrid:cr_pfft3b_queueout  *
 *                              *
 ********************************/
void PGrid::cr_pfft3b_queueout(const int nb, const int nffts_out, double *a) 
{
   int shift1, shift2;
   int indx1 = aqindx[0];

   //while (aqstatus[indx1] < 5) {
   while (aqstatus[indx1] < aqmax) 
   {
      for (auto q = 0; q < aqsize; ++q) 
      {
         int indx   = aqindx[q];
         int status = aqstatus[indx] + 1;
         int nffts  = aqnffts[indx];
         shift1 = nffts_max*n2ft3d*(2*indx);
         shift2 = nffts_max*n2ft3d*(2*indx+1);
         if (staged_gpu_fft_pipeline)
            pfftb_step12(status,nffts,nb,a,atmp+shift1,atmp+shift2,indx+4,indx);
         else
            pfftb_step(status,nffts,nb,a,atmp+shift1,atmp+shift2,indx+4);
         ++aqstatus[indx];
      }
   }
   double scal1 = 1.0 / ((double)((nx) * (ny) * (nz)));
   double enrr0 = scal1 * d3db::rr_dot(atmp, atmp);
 
   shift1 = nffts_max*n2ft3d * (2 * indx1);
   std::memcpy(a, atmp+shift1, nffts_out*n2ft3d*sizeof(double));
   --aqsize;
   for (auto q=0; q<aqsize; ++q)
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
void PGrid::pfftfx(const int nffts, const int nb, double *a, double *tmp1, double *tmp2, int request_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      // do fft along nx dimension
      d3db::mygdevice.batch_rfftx_tmpx(d3db::fft_tag,true, nx, nffts*ny*nq, n2ft3d, a, d3db::tmpx);
      std::memcpy(tmp1, a, nffts*n2ft3d*sizeof(double));
   }
   /**** hilbert mapping ****/
   else 
   {
      // do fft along nx dimension
      // A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]
      d3db::mygdevice.batch_rfftx_tmpx(d3db::fft_tag,true, nx, nffts*nq1, n2ft3d, a, d3db::tmpx);
      d3db::c_ptranspose_ijk_start(nffts, nb, 0, a, tmp1, tmp2, request_indx, 40);
   }
}

/********************************
 *                              *
 *        PGrid::pfftfy         *
 *                              *
 ********************************/
/**
 * Performs the 3D Fast Fourier Transform (FFT) pipeline along the ny dimension of the dataset as part of
 * a multi-dimensional FFT processing sequence. This function handles FFT computations differently
 * based on the specified mapping strategy (slab or Hilbert mapping).
 *
 * @param nffts        Number of 3D FFTs to perform in this batch, affecting the scale of transformation.
 * @param nb           Fermi sphere packing number, indicating energy cutoffs:
 *                     - '1' for lower cutoff energy used in wavefunctions,
 *                     - '0' for higher cutoff energy used in density grids.
 * @param tmp1         Temporary storage array used as input data for the FFT or for data transposition.
 * @param tmp2         Temporary storage array where the FFT output or transposed data is stored.
 * @param request_indx Identifier for this specific FFT operation sequence, used to manage and track
 *                     FFT operations across different stages and configurations.
 *
 * Processing Details:
 * - Slab Mapping (maptype == 1):
 *   - Transposes data for alignment along the ny dimension before performing FFT.
 *   - Data is reorganized based on zero row filtering, followed by FFT along ny.
 *   - The function uses optimized data handling to efficiently perform transformations and memory operations.
 * - Hilbert Mapping (maptype != 1):
 *   - Data is transposed from tmp1 to tmp2 for subsequent FFT along nz, using configurations that adapt
 *     the process to Hilbert mapping strategies.
 *
 * Note:
 * - 'maptype' determines the FFT configuration and influences how the function processes the data.
 * - This function is critical in a high-performance computing context, designed to handle complex,
 *   large-scale data arrays across multiple dimensions.
 */
void PGrid::pfftfy(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      auto nxh = nx/2 + 1;
      auto nxh2 = nx + 2;
      auto nxhy2 = nxh2 * ny;
     
      // do fft along ny dimension
      // A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]
      int indx0 = 0;
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i = 0; i < nxh; ++i) 
            {
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
      }
     
      d3db::mygdevice.batch_cffty_tmpy(d3db::fft_tag,true, ny, nn, n2ft3d, tmp2, d3db::tmpy);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nxh; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  auto jj = 0;
                  auto indx3 = 2*i + indx0;
                  auto shift = 2*ny*nn;
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
      }
     
      // Do a transpose of A
      // A(ky,nz,ky) <- A(kx,ky,nz)
      d3db::c_ptranspose2_jk_start(nffts, nb, tmp1, tmp2, tmp1, request_indx, 41);
   }
 
   /**** hilbert mapping ****/
   else 
   {
      // in=tmp1, out=tmp2
      d3db::c_ptranspose_ijk_end(nffts, nb, 0, tmp1, tmp2, request_indx);
     
      // do fft along ny dimension
      // A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]
      d3db::mygdevice.batch_cffty_tmpy_zero(d3db::fft_tag,true,ny,nq2,nffts,n2ft3d,tmp1,d3db::tmpy,zero_row2[nb]);
     
      // in=tmp2, out=tmp2
      d3db::c_ptranspose_ijk_start(nffts, nb, 1, tmp1, tmp2, tmp1, request_indx, 42);
   }
}

/********************************
 *                              *
 *        PGrid::pfftfz         *
 *                              *
 ********************************/
/**
 * Performs the 3D Fast Fourier Transform (FFT) pipeline along the nz dimension of the dataset as part of
 * a multi-dimensional FFT processing sequence. This function handles FFT computations differently
 * based on the specified mapping strategy (slab or Hilbert mapping).
 *
 * @param nffts        Number of 3D FFTs to perform in this batch, affecting the scale of transformation.
 * @param nb           Fermi sphere packing number, indicating energy cutoffs:
 *                     - '1' for lower cutoff energy used in wavefunctions,
 *                     - '0' for higher cutoff energy used in density grids.
 * @param tmp1         Temporary storage array used as input data for the FFT or for data transposition.
 * @param tmp2         Temporary storage array where the FFT output or transposed data is stored.
 * @param request_indx Identifier for this specific FFT operation sequence, used to manage and track
 *                     FFT operations across different stages and configurations.
 *
 * Processing Details:
 * - Slab Mapping (maptype == 1):
 *   - Transposes data for alignment along the nz dimension before performing FFT.
 *   - Data is reorganized based on zero row filtering, followed by FFT along nz.
 *   - The function uses optimized data handling to efficiently perform transformations and memory operations.
 * - Hilbert Mapping (maptype != 1):
 *   - Data is transposed from tmp1 to tmp2 for subsequent FFT along nz, using configurations that adapt
 *     the process to Hilbert mapping strategies.
 *
 * Note:
 * - 'maptype' determines the FFT configuration and influences how the function processes the data.
 * - This function is critical in a high-performance computing context, designed to handle complex,
 *   large-scale data arrays across multiple dimensions.
 */
void PGrid::pfftfz(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx) 
{

   /**** slab mapping ****/
   if (maptype == 1) 
   {
      d3db::c_ptranspose2_jk_end(nffts, nb, tmp2, tmp1, request_indx);
     
      auto nxh = nx / 2 + 1;
      auto nxh2 = nx + 2;
      auto nxhz2 = nxh2 * nz;
     
      // do fft along nz dimension
      // A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]
      int indx0 = 0;
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
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
      }
     
      d3db::mygdevice.batch_cfftz_tmpz(d3db::fft_tag,true, nz, nn, n2ft3d, tmp1, d3db::tmpz);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nxh; ++i) 
            {
               if (!zero_row3[nb][indx2]) 
               {
                  auto kk = 0;
                  auto indx3 = 2*i + indx0;
                  auto shift = 2*nz*nn;
                  for (auto k=0; k<nz; ++k) 
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
   }
 
   /**** hilbert mapping ****/
   else 
   {
      // in=tmp1, out=tmp2
      d3db::c_ptranspose_ijk_end(nffts, nb, 1, tmp2, tmp1, request_indx);
     
      // do fft along nz dimension
      // A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]
      d3db::mygdevice.batch_cfftz_tmpz_zero(d3db::fft_tag,true,nz,nq3,nffts,n2ft3d,tmp2,d3db::tmpz,zero_row3[nb]);
   }
}

/********************************
 *                              *
 *       PGrid::pfftf_step      *
 *                              *
 ********************************/
/**
 * Controls the execution of multi-dimensional Fast Fourier Transform (FFT) pipeline operations 
 * on a dataset by orchestrating steps across different dimensions (x, y, z). This function
 * sequentially triggers FFT operations or packing steps based on the 'step' parameter.
 *
 * @param step         The specific stage of the FFT process to execute, dictating which FFT
 *                     or packing function is called:
 *                     - 0: Initial transfer of data for FFT along x-axis.
 *                     - 1: Perform FFT along y-axis.
 *                     - 2: Finalize FFT along z-axis and start packing data.
 *                     - 3: Complete packing of data.
 * @param nffts        Number of 3D FFTs to perform, affects the scope of transformation.
 * @param nb           Fermi sphere packing number, influencing configurations like:
 *                     - '1' for lower energy cutoff used in wavefunctions.
 *                     - '0' for higher energy cutoff used in density grids.
 * @param a            Input data array for the initial FFT step.
 * @param tmp1         Temporary storage array used as an intermediary in FFT transformations.
 * @param tmp2         Temporary storage array used to hold the output of FFT operations.
 * @param request_indx Identifier for this specific FFT operation sequence, used to track
 *                     and manage the request through various transformation stages.
 *
 * Processing Flow:
 * - step 0: Moves memory to device and prepares for FFT along the x-dimension.
 * - step 1: Processes FFT along the y-dimension.
 * - step 2: Completes FFT along the z-dimension and begins data packing.
 * - step 3: Finalizes data packing and prepares it for further processing or output.
 *
 * Note:
 * - Each step corresponds to a specific operation and must be executed in sequence for correct
 *   processing. Misalignment in step ordering can lead to incorrect results or data corruption.
 * - Designed for high-performance computing contexts dealing with large-scale data sets in multi-dimensional spaces.
 */
void PGrid::pfftf_step(const int step, const int nffts,  const int nb, double *a, double *tmp1, double *tmp2, int request_indx)
{
   if (step==0)
   {
      parall->comm_Barrier(1);
      // pfftfx mem-->device, in=a out=tmp2
      pfftfx(nffts, nb, a, tmp1, tmp2, request_indx);
   }
   else if (step==1)
   {
      parall->comm_Barrier(1);
      // pfftfy device, in=tmp1
      pfftfy(nffts, nb, tmp1, tmp2, request_indx);
   }
   else if (step==2)
   {
      parall->comm_Barrier(1);
      // pfftfz device-->mem
      pfftfz(nffts, nb, tmp1, tmp2, request_indx);
      this->c_pack_start(nffts, nb, tmp2, tmp1, request_indx, 47);
   }
   else if (step==3)
   {
      parall->comm_Barrier(1);
      // pfftf final
      this->c_pack_end(nffts, nb, tmp2, request_indx);
   }
}

/********************************
 *                              *
 *    PGrid::pfftfx_start       *
 *                              *
 ********************************/
/**
 * Initiates batch processing of 3D Fast Fourier Transforms (FFT) pipeline  on a given dataset, setting up data and
 * performing the initial stages of FFT. This function supports two types of FFT mappings, each affecting
 * the execution along the nx dimension.
 *
 * @param nffts        Number of 3D FFTs to perform in this batch, affecting the transformation scale.
 * @param nb           Fermi sphere packing number, with possible values:
 *                     - '1' for lower cutoff energy used in wavefunctions,
 *                     - '0' for higher cutoff energy used in density grids.
 * @param a            Pointer to the input data array.
 * @param tmp1         Pointer to a temporary storage array (unused in this function).
 * @param tmp2         Pointer to a temporary storage array, used for storing data before FFT processing.
 * @param request_indx Index or identifier for this specific FFT request (unused directly in this function).
 * @param da_indx      Dataset identifier used in GPU operations to specify a subset of data.
 *
 * Processing Details:
 * - Slab Mapping (maptype == 1):
 *   - Copies data from 'a' to 'tmp2'.
 *   - Performs FFT along the nx dimension using the copied data.
 * - Hilbert Mapping (maptype != 1):
 *   - Copies data from 'a' to 'tmp2'.
 *   - Executes FFT along the nx dimension, tailored for the Hilbert mapping strategy.
 *
 * Note:
 * - 'maptype' needs to be set externally and determines the FFT configuration.
 * - Designed for high-performance computing contexts for processing extensive multidimensional datasets.
 */
void PGrid::pfftfx_start(const int nffts, const int nb, double *a, double *tmp1, double *tmp2, int request_indx, int da_indx)
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      // do fft along nx dimension
      std::memcpy(tmp2, a, nffts*n2ft3d*sizeof(double));
      d3db::mygdevice.batch_rfftx_stages_tmpx(0,d3db::fft_tag,true,nx,nffts*ny*nq, n2ft3d,tmp2,d3db::tmpx,da_indx);
   }
   /**** hilbert mapping ****/
   else 
   {
      // do fft along nx dimension
      // A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]
      std::memcpy(tmp2, a, nffts*n2ft3d*sizeof(double));
      d3db::mygdevice.batch_rfftx_stages_tmpx(0,d3db::fft_tag,true,nx,nffts*nq1,n2ft3d,tmp2,d3db::tmpx,da_indx);
   }
}

/********************************
 *                              *
 *    PGrid::pfftfx_compute     *
 *                              *
 ********************************/
/**
 * Executes the main computation step of the 3D Fast Fourier Transforms (FFT) pipeline  for a given dataset,
 * which is part of a batch processing sequence initiated by `pfftfx_start`. This function specifically
 * handles the computation stages of FFT based on the selected mapping strategy.
 *
 * @param nffts        Number of 3D FFTs to perform in this batch, determines the scale of transformation.
 * @param nb           Fermi sphere packing number, indicating energy cutoffs:
 *                     - '1' for lower cutoff energy used in wavefunctions,
 *                     - '0' for higher cutoff energy used in density grids.
 * @param a            Pointer to the input data array (not directly used in this function).
 * @param tmp1         Pointer to a temporary storage array (unused in this function).
 * @param tmp2         Pointer to the temporary storage array, used for holding data during FFT processing.
 * @param request_indx Index or identifier for this specific FFT request (not directly used in this function).
 * @param da_indx      Dataset identifier used in GPU operations to specify the data subset being processed.
 *
 * Processing Details:
 * - Slab Mapping (maptype == 1):
 *   - Performs FFT computation along the nx dimension using the data in 'tmp2'.
 * - Hilbert Mapping (maptype != 1):
 *   - Executes FFT computation along the nx dimension, tailored for the Hilbert mapping strategy.
 *
 * Note:
 * - 'maptype' is externally defined and selects the FFT configuration.
 * - This function is crucial within a high-performance computing framework intended for extensive multidimensional data arrays.
 */
void PGrid::pfftfx_compute(const int nffts, const int nb, double *a, double *tmp1, double *tmp2, int request_indx, int da_indx)
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      // do fft along nx dimension
      d3db::mygdevice.batch_rfftx_stages_tmpx(1,d3db::fft_tag,true,nx,nffts*ny*nq,n2ft3d,tmp2,d3db::tmpx,da_indx);
   }
   /**** hilbert mapping ****/
   else 
   {
      // do fft along nx dimension
      // A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]
      d3db::mygdevice.batch_rfftx_stages_tmpx(1,d3db::fft_tag,true,nx,nffts*nq1,n2ft3d,tmp2,d3db::tmpx,da_indx);
   }
}

/********************************
 *                              *
 *      PGrid::pfftfx_end       *
 *                              *
 ********************************/
/**
 * Finalizes the batch of 3D Fast Fourier Transform (FFT) pipeline operations initiated by `pfftfx_start`.
 * This function performs the last stages of FFT and handles the necessary post-processing for the data.
 * It supports different processing strategies based on the 'maptype' setting.
 *
 * @param nffts        Number of 3D FFTs to be processed in this batch.
 * @param nb           Fermi sphere packing number, indicating energy cutoffs:
 *                     - '1' for lower cutoff energy used in wavefunctions.
 *                     - '0' for higher cutoff energy used in density grids.
 * @param a            Pointer to the initial data array (unused directly in this function).
 * @param tmp1         Pointer to the temporary array where the final processed data is stored.
 * @param tmp2         Pointer to the temporary array for intermediate data processing.
 * @param request_indx Identifier for this FFT request, used for data transposition in Hilbert mapping.
 * @param da_indx      Dataset identifier used in GPU operations.
 *
 * Processing Details:
 * - Slab Mapping (maptype == 1):
 *   - Completes FFT along the nx dimension.
 *   - Copies processed data from 'tmp2' to 'tmp1'.
 * - Hilbert Mapping (maptype != 1):
 *   - Completes FFT along the nx dimension.
 *   - Invokes 'c_ptranspose_ijk_start' for complex data transposition.
 *
 * Note:
 * - Assumes external setting of 'maptype' to distinguish between Slab and Hilbert mappings.
 * - Integral to a high-performance computing framework for multidimensional FFT operations.
 */
void PGrid::pfftfx_end(const int nffts, const int nb, double *a, double *tmp1, double *tmp2, int request_indx, int da_indx)
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      // do fft along nx dimension
      d3db::mygdevice.batch_rfftx_stages_tmpx(2,d3db::fft_tag,true,nx,nffts*ny*nq,n2ft3d,tmp2,d3db::tmpx,da_indx);
      std::memcpy(tmp1,tmp2,nffts*n2ft3d*sizeof(double));
   }
   /**** hilbert mapping ****/
   else 
   {
      // do fft along nx dimension
      // A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]
      d3db::mygdevice.batch_rfftx_stages_tmpx(2,d3db::fft_tag,true, nx, nffts*nq1, n2ft3d, tmp2, d3db::tmpx,da_indx);
 
      d3db::c_ptranspose_ijk_start(nffts,nb,0,tmp2,tmp1,tmp2,request_indx,40);
 
   }
}

/********************************
 *                              *
 *      PGrid::pfftfy_start     *
 *                              *
 ********************************/
/**
 * Initiates the batch processing of 3D Fast Fourier Transforms (FFT) pipeline  along the ny dimension,
 * catering to different FFT mapping strategies: slab and Hilbert. This function sets up and
 * executes initial stages of FFT tailored for large-scale data arrays.
 *
 * @param nffts        Number of 3D FFTs to perform in this batch, affecting the transformation scale.
 * @param nb           Fermi sphere packing number, with possible values:
 *                     - '1' for lower cutoff energy used in wavefunctions,
 *                     - '0' for higher cutoff energy used in density grids.
 * @param tmp1         Pointer to the input array of data.
 * @param tmp2         Pointer to the output array where transformed data is stored.
 * @param request_indx Index or identifier for this specific FFT request, used in transposition in Hilbert mapping.
 * @param da_indx      Dataset identifier used in GPU operations to specify a subset of data.
 *
 * Processing Details:
 * - Slab Mapping (maptype == 1):
 *   - Arranges and copies data based on zero row filtering before performing FFT along ny.
 *   - Data preparation involves optimized memory handling to set up for batched FFT operations.
 * - Hilbert Mapping (maptype != 1):
 *   - Data transposition is performed from tmp1 to tmp2 using pre-defined configurations.
 *   - Post-transposition, executes FFT along ny dimension, designed for the Hilbert mapping approach.
 *
 * Note:
 * - 'maptype' is externally defined, guiding the FFT configuration and execution approach.
 * - This function is integral to a high-performance computing framework focused on multidimensional FFT operations.
 */
void PGrid::pfftfy_start(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      auto nxh = nx/2 + 1;
      auto nxh2 = nx + 2;
      auto nxhy2 = nxh2 * ny;
     
      // do fft along ny dimension
      // A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]
      int indx0 = 0;
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<(nffts*nq); ++q) 
         {
            for (auto i=0; i<nxh; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  auto jj = 0;
                  auto indx3 = 2*i + indx0;
                  auto shift = 2*ny*nn;
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
      }
     
      d3db::mygdevice.batch_cffty_stages_tmpy(0,d3db::fft_tag,true,ny,nn,n2ft3d,tmp2,d3db::tmpy,da_indx);
   }
 
   /**** hilbert mapping ****/
   else 
   {
      // in=tmp1, out=tmp2
      d3db::c_ptranspose_ijk_end(nffts,nb,0,tmp1,tmp2,request_indx);
     
      // do fft along ny dimension
      // A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]
      d3db::mygdevice.batch_cffty_stages_tmpy_zero(0,d3db::fft_tag,true,ny,nq2,nffts,n2ft3d,tmp1,d3db::tmpy,zero_row2[nb],da_indx);
   }
}

/********************************
 *                              *
 *      PGrid::pfftfy_compute   *
 *                              *
 ********************************/
/**
 * Executes the main computation step of the 3D Fast Fourier Transforms (FFT) pipeline along the ny dimension
 * for a given dataset. This function is part of a batch processing sequence that deals with large-scale
 * data arrays, applying FFT based on selected mapping strategies.
 *
 * @param nffts        Number of 3D FFTs to perform in this batch, determines the scale of transformation.
 * @param nb           Fermi sphere packing number, with possible values:
 *                     - '1' for lower cutoff energy used in wavefunctions,
 *                     - '0' for higher cutoff energy used in density grids.
 * @param tmp1         Pointer to the input array of data used in Hilbert mapping.
 * @param tmp2         Pointer to the output array where transformed data is stored.
 * @param request_indx Index or identifier for this specific FFT request, used in configurations for Hilbert mapping.
 * @param da_indx      Dataset identifier used in GPU operations to specify the subset of data.
 *
 * Processing Details:
 * - Slab Mapping (maptype == 1):
 *   - Arranges data based on zero row filtering before performing FFT along ny.
 *   - Efficiently handles memory and data flow to optimize FFT computations.
 * - Hilbert Mapping (maptype != 1):
 *   - Directly uses pre-configured data from tmp1 for FFT computations along ny, adapted for Hilbert mapping strategy.
 *
 * Note:
 * - 'maptype' must be externally set and determines the FFT configuration and execution approach.
 * - This function is essential in high-performance computing frameworks focused on extensive multidimensional FFT operations.
 */
void PGrid::pfftfy_compute(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
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
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q)
         {
            for (auto i=0; i<nxh; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  ++nn;
               }
               ++indx2;
            }
         }
      }
     
      d3db::mygdevice.batch_cffty_stages_tmpy(1,d3db::fft_tag,true,ny,nn,n2ft3d,tmp2,d3db::tmpy,da_indx);
   }
 
   /**** hilbert mapping ****/
   else 
   {
      // in=tmp1, out=tmp2
      // do fft along ny dimension
      // A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]
      d3db::mygdevice.batch_cffty_stages_tmpy_zero(1,d3db::fft_tag,true,ny,nq2,nffts,n2ft3d,tmp1,d3db::tmpy,zero_row2[nb],da_indx);
   }
}




/********************************
 *                              *
 *      PGrid::pfftfy_end       *
 *                              *
 ********************************/
void PGrid::pfftfy_end(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
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
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nxh; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  ++nn;
               }
               ++indx2;
            }
         }
      }
     
      d3db::mygdevice.batch_cffty_stages_tmpy(2,d3db::fft_tag,true,ny,nn,n2ft3d,tmp2,d3db::tmpy,da_indx);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
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
      }
     
      // Do a transpose of A
      // A(ky,nz,ky) <- A(kx,ky,nz)
      d3db::c_ptranspose2_jk_start(nffts,nb,tmp1,tmp2,tmp1,request_indx,41);
   }
   /**** hilbert mapping ****/
   else 
   {
      // do fft along ny dimension
      // A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]
      d3db::mygdevice.batch_cffty_stages_tmpy_zero(2,d3db::fft_tag,true,ny,nq2,nffts,n2ft3d,tmp1,d3db::tmpy,zero_row2[nb],da_indx);
     
      // in=tmp2, out=tmp2
      d3db::c_ptranspose_ijk_start(nffts, nb, 1, tmp1, tmp2, tmp1, request_indx, 42);
   }
}



/********************************
 *                              *
 *      PGrid::pfftfz_start     *
 *                              *
 ********************************/
void PGrid::pfftfz_start(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      d3db::c_ptranspose2_jk_end(nffts, nb, tmp2, tmp1, request_indx);
     
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
      d3db::c_ptranspose_ijk_end(nffts, nb, 1, tmp2, tmp1, request_indx);
     
      // do fft along nz dimension
      // A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]
      d3db::mygdevice.batch_cfftz_stages_tmpz_zero(0,d3db::fft_tag,true, nz, nq3,nffts, n2ft3d, tmp2, d3db::tmpz, zero_row3[nb],da_indx);
   }
}


/********************************
 *                              *
 *      PGrid::pfftfz_compute   *
 *                              *
 ********************************/
void PGrid::pfftfz_compute(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
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
      d3db::mygdevice.batch_cfftz_stages_tmpz_zero(1,d3db::fft_tag,true, nz, nq3,nffts, n2ft3d, tmp2, d3db::tmpz, zero_row3[nb],da_indx);
   }
}

/********************************
 *                              *
 *      PGrid::pfftfz_end       *
 *                              *
 ********************************/
void PGrid::pfftfz_end(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
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
      d3db::mygdevice.batch_cfftz_stages_tmpz_zero(2,d3db::fft_tag,true, nz, nq3,nffts, n2ft3d, tmp2, d3db::tmpz, zero_row3[nb],da_indx);
   }
}


/********************************
 *                              *
 *       PGrid::pfftf_step10    *
 *                              *
 ********************************/
void PGrid::pfftf_step10(const int step, const int nffts, const int nb, double *a, double *tmp1,
                         double *tmp2, int request_indx, int da_indx)
{
   // pfftfx mem-->device, in=a out=tmp2
   if (step==0)
   {
     pfftfx_start(nffts, nb, a, tmp1, tmp2, request_indx,da_indx);
   }
   else if (step==1)
   {
     pfftfx_compute(nffts, nb, a, tmp1, tmp2, request_indx,da_indx);
   }
   else if (step==2)
   {
     pfftfx_end(nffts, nb, a, tmp1, tmp2, request_indx,da_indx);
   }


   // pfftfy device, in=tmp1
   else if (step==3)
   {
      pfftfy_start(nffts, nb, tmp1, tmp2, request_indx,da_indx);
   }
   else if (step==4)
   {
      pfftfy_compute(nffts, nb, tmp1, tmp2, request_indx,da_indx);
   }
   else if (step==5)
   {
      pfftfy_end(nffts, nb, tmp1, tmp2, request_indx,da_indx);
   }


   // pfftfz device-->mem
   else if (step==6)
   {
      pfftfz_start(nffts, nb, tmp1, tmp2, request_indx,da_indx);
   }
   else if (step==7)
   {
      pfftfz_compute(nffts, nb, tmp1, tmp2, request_indx,da_indx);
   }
   else if (step==8)
   {
      pfftfz_end(nffts, nb, tmp1, tmp2, request_indx,da_indx);
      this->c_pack_start(nffts, nb, tmp2, tmp1, request_indx, 47);
   }


   else if (step==9)
   {
      // pfftf final
      this->c_pack_end(nffts, nb, tmp2, request_indx);
   }
}


/********************************
 *                              *
 *       PGrid:c_pack_start     *
 *                              *
 ********************************/
void PGrid::c_pack_start(const int nffts, const int nb, double *a, double *tmp1,
                         const int request_indx, const int msgtype) 
{
   // int one=1;
 
   // DCOPY_PWDFT(n2ft3d,a,one,tmp,one);
   std::memcpy(tmp1, a, nffts*n2ft3d * sizeof(double));
   std::memset(a, 0, nffts*n2ft3d * sizeof(double));
 
   for (auto s=0; s<nffts; ++s)
      c_aindexcopy((nida[nb] + nidb2[nb]), packarray[nb], tmp1+s*n2ft3d, a+s*n2ft3d);
 
   if (balanced)
      for (auto s=0; s<nffts; ++s)
         mybalance->c_balance_start(1, nb, a+s*n2ft3d, request_indx, msgtype);
 
   return;
}

/********************************
 *                              *
 *       PGrid:c_pack_end       *
 *                              *
 ********************************/
void PGrid::c_pack_end(const int nffts, const int nb, double *tmp1, const int request_indx) 
{
   if (balanced)
      for (auto s=0; s<nffts; ++s)
         mybalance->c_balance_end(1, nb, tmp1+s*n2ft3d, request_indx);

   return;
}

/********************************
 *                              *
 *    PGrid:rc_pfft3f_queuein   *
 *                              *
 ********************************/
void PGrid::rc_pfft3f_queuein(const int nb, const int nffts_in,  double *b) 
{
   //int nffts_in = 1;
   int shift1, shift2;
   int np = parall->np_i();
 
   for (auto q = 0; q < bqsize; ++q) 
   {
      int indx   = bqindx[q];
      int status = bqstatus[indx] + 1;
      int nffts  = bqnffts[indx];
      shift1 = nffts_max*n2ft3d * (2*indx);
      shift2 = nffts_max*n2ft3d * (2*indx + 1);
      if (staged_gpu_fft_pipeline)
         pfftf_step10(status,nffts,nb,b,btmp+shift1,btmp+shift2,indx+4,indx);
      else
         pfftf_step(status,nffts,nb,b,btmp+shift1,btmp+shift2,indx+4);
      ++bqstatus[indx];
   }
 
   ++blast_index;
   if (blast_index >= bqmax)
      blast_index = 0;
   ++bqsize;
   bqindx[bqsize - 1] = blast_index;
   bqstatus[blast_index] = 0;
   bqnffts[blast_index]  = nffts_in;
 
   // status = 0;
   shift1 = nffts_max*n2ft3d * (2*blast_index);
   shift2 = nffts_max*n2ft3d * (2*blast_index + 1);
 
   if (staged_gpu_fft_pipeline)
      pfftf_step10(0,nffts_in,nb,b,btmp+shift1,btmp+shift2,blast_index+4,blast_index);
   else
      pfftf_step(0,nffts_in,nb,b,btmp+shift1,btmp+shift2,blast_index+4);
}

/********************************
 *                              *
 *    PGrid:rc_pfft3f_queueout  *
 *                              *
 ********************************/
void PGrid::rc_pfft3f_queueout(const int nb, const int nffts_out, double *b) 
{
   int shift1, shift2;
   int indx1 = bqindx[0];
 
   //while (bqstatus[indx1] < 5) {
   while (bqstatus[indx1] < bqmax) 
   {
      for (auto q = 0; q < bqsize; ++q) 
      {
         int indx   = bqindx[q];
         int status = bqstatus[indx] + 1;
         int nffts  = bqnffts[indx];
         shift1 = nffts_max*n2ft3d * (2*indx);
         shift2 = nffts_max*n2ft3d * (2*indx + 1);
         if (staged_gpu_fft_pipeline)
            pfftf_step10(status,nffts,nb,b,btmp+shift1,btmp+shift2,indx+4,indx);
         else
            pfftf_step(status,nffts,nb,b,btmp+shift1,btmp+shift2,indx+4);
         ++bqstatus[indx];
      }
   }
   double scal1 = 1.0 / ((double)((nx) * (ny) * (nz)));
   double enrr0 = scal1 * d3db::rr_dot(btmp, btmp);
 
   shift2 = nffts_max*n2ft3d * (2 * indx1 + 1);
   std::memcpy(b, btmp + shift2, nffts_out*n2ft3d*sizeof(double));
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
/**
 * @brief Copies real numbers from one array to another, converting them into complex format.
 *
 * This function transforms a real-valued array into a complex array where each real number is followed 
 * by an imaginary part set to zero. It is typically used to prepare real data for computations that 
 * require complex number inputs, such as FFTs.
 *
 * @param nb An index that determines the number of elements to process, based on the configuration arrays
 *           `nida` and `nidb`, which provide bounds for the operations.
 * @param a Pointer to the source array containing real values.
 * @param b Pointer to the destination array where real values from `a` are stored as the real parts of complex numbers,
 *          with each subsequent imaginary part initialized to zero.
 *
 * @note The destination array `b` must be at least twice the size of the number of elements specified by `nb` to accommodate
 *       both real and imaginary parts. It is assumed that appropriate memory allocation for `b` has been managed externally.
 */
void PGrid::tc_pack_copy(const int nb, double *a, double *b) 
{
   int i, ii;
   int ng = nida[nb] + nidb[nb];
 
   ii = 0;
   for (i = 0; i < ng; ++i) 
   {
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
/**
 * @brief Multiplies a scalar array with a complex array element-wise.
 *
 * This function performs an element-wise multiplication of a scalar array `a` with a complex array `b`,
 * storing the results in a complex array `c`. Each scalar in `a` multiplies both the real and imaginary 
 * components of the corresponding complex number in `b`. The complex numbers are expected to be stored 
 * in a contiguous block of memory with real and imaginary parts interleaved.
 *
 * @param nb An index that determines the number of elements to process, based on the configuration arrays
 *           `nida` and `nidb`, which provide bounds for the operations.
 * @param a Pointer to the scalar array `a` of size `nida[nb] + nidb[nb]`.
 * @param b Pointer to the complex array `b` which stores complex numbers as interleaved doubles.
 *           Thus, for a complex number at index `i`, `b[2*i]` is the real part, and `b[2*i + 1]` is the imaginary part.
 * @param c Pointer to the complex array `c` where the results of the multiplication are stored in the same
 *           interleaved format as `b`.
 *
 * @note It's assumed that the memory allocation for `a`, `b`, and `c` has been appropriately managed externally,
 *       with `c` having at least as many elements as `b`. The function does not handle cases where `a`, `b`, or `c`
 *       are `nullptr`. If such error checking is desired, it should be added externally before calling this function.
 *
 * @warning This function does not perform bounds checking on the array accesses. Ensure the size specifications
 *          via `nb` correctly map to the sizes of `a`, `b`, and `c` to prevent out-of-bounds errors.
 */
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
 *     PGrid:c_pack_addzeros    *
 *                              *
 ********************************/
void PGrid::c_pack_addzeros(const int nb, const double vzero, double *a) 
{
   int pzero0  = ijktop(0, 0, 0);
   int pzero1  = ijktop(1, 0, 0);
   int pzero2  = ijktop(0, 1, 0);
   int pzero3  = ijktop(0, 0, 1);
   int index0 = ijktoindex(0, 0, 0);
   int index1 = ijktoindex(1, 0, 0);
   int index2 = ijktoindex(0, 1, 0);
   int index3 = ijktoindex(0, 0, 1);
   if (pzero0 == parall->taskid_i()) a[2*index0] += 0.25*vzero;
   if (pzero1 == parall->taskid_i()) a[2*index1] += 0.25*vzero;
   if (pzero2 == parall->taskid_i()) a[2*index2] += 0.25*vzero;
   if (pzero3 == parall->taskid_i()) a[2*index3] += 0.25*vzero;
}

/********************************
 *                              *
 *     PGrid:c_pack_addzero     *
 *                              *
 ********************************/
void PGrid::c_pack_addzero(const int nb, const double vzero, double *a) 
{ 
   int pzero0  = ijktop(0, 0, 0);
   int index0 = ijktoindex(0, 0, 0);
   if (pzero0 == parall->taskid_i()) a[2*index0] += 0.25*vzero;
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

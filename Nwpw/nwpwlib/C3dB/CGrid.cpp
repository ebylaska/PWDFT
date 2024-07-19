/* CGrid.cpp
   Author - Eric Bylaska

        this class is used for defining 3d c3db::parallel maps
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

#include "CGrid.hpp"

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

// CGrid::CGrid(Parallel *inparall, Lattice *inlattice, Control2& control) :
// c3db(inparall,control.mapping(),control.ngrid(0),control.ngrid(1),control.ngrid(2))

CGrid::CGrid(Parallel *inparall, Lattice *inlattice, int mapping0, int balance0,
             int nx0, int ny0, int nz0, int pfft3_qsize0, bool staged_gpu_fft_pipeline0, Brillouin *inbrillouin)
    : k1db(inparall, mapping0, inbrillouin->nbrillouin),
      c3db(inparall, mapping0, nx0, ny0, nz0, nbrillq) 
{
   int nxh, nyh, nzh, p, q, indx, nb;
   double *G1, *G2, *G3;
   double ggcut, eps, ggmax, ggmin;
   bool *zero_arow3, *zero_arow2;
   bool yzslab, zrow;


   int nwave_in[nbrillq+1], nwave_out[nbrillq+1];
 
   lattice = inlattice;
   mybrillouin = inbrillouin;
 
   eps = 1.0e-12;

   // aligned Memory
   //std::size_t aligned_size3d = (3 * nfft3d * sizeof(double) + Alignment - 1) & ~(Alignment - 1);
   std::size_t aligned_size3d = (3 * nfft3d);
   Garray = new (std::nothrow) double[aligned_size3d]();
   G1 = Garray;
   G2 = Garray + nfft3d;
   G3 = Garray + 2*nfft3d;
   nxh = nx/2;
   nyh = ny/2;
   nzh = nz/2;
   ggmin = 9.9e9;
   ggmax = 0.0;
   for (auto k3=(-nzh+1); k3<=nzh; ++k3)
   for (auto k2=(-nyh+1); k2<=nyh; ++k2)
   for (auto k1=(-nxh+1); k1<=nxh; ++k1)
   {
      auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2);
      auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2);
      auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2);
      auto gg = gx*gx + gy*gy + gz*gz;
      if (gg > ggmax)
        ggmax = gg;
      if ((gg < ggmin) && (gg > 1.0e-6))
        ggmin = gg;
      auto i = k1; if (i<0) i = i + nx;
      auto j = k2; if (j<0) j = j + ny;
      auto k = k3; if (k<0) k = k + nz;
     
      auto indx = cijktoindex(i, j, k);
      auto p = cijktop(i, j, k);
      if (p == c3db::parall->taskid_i()) 
      {
         G1[indx] = gx;
         G2[indx] = gy;
         G3[indx] = gz;
      }
   }
   Gmax = sqrt(ggmax);
   Gmin = sqrt(ggmin);

   // for printing grid values
   //nidb_print      = new (std::nothrow) int[nbrillouin]();
   //nwave_all_print = new (std::nothrow) int[nbrillouin]();
   //nidb_print.resize(nbrillouin,0);
   //nwave_all_print.resize(nbrillouin,0);
 
   // aligned Memory
   //nidb         = new (std::nothrow) int[nbrillq+1]();
   //nidb2        = new (std::nothrow) int[nbrillq+1]();
   //nwave        = new (std::nothrow) int[nbrillq+1]();
   //nwave_all    = new (std::nothrow) int[nbrillq+1]();
   //nwave_entire = new (std::nothrow) int[nbrillq+1]();

 
   // aligned Memory
   //nidb         = new (std::nothrow) int[nbrillq+1]();
   //nidb2        = new (std::nothrow) int[nbrillq+1]();
   //nwave        = new (std::nothrow) int[nbrillq+1]();
   //nwave_all    = new (std::nothrow) int[nbrillq+1]();
   //nwave_entire = new (std::nothrow) int[nbrillq+1]();

   nidb.resize(nbrillq + 1);
   nidb2.resize(nbrillq + 1);
   nwave.resize(nbrillq + 1);
   nwave_all.resize(nbrillq + 1);
   nwave_entire.resize(nbrillq + 1);
   Gpack.resize(nbrillq + 1);

   p_weight  = new (std::nothrow) double[nbrillq];
   p_kvector = new (std::nothrow) double[3*nbrillq];
   for (auto nbq=0; nbq<nbrillq; ++nbq)
   {
      auto nb = k1db::ktoptok(nbq);
      p_weight[nbq]      = mybrillouin->weight[nb];
      p_kvector[3*nbq]   = mybrillouin->kvector[3*nb];
      p_kvector[3*nbq+1] = mybrillouin->kvector[3*nb+1];
      p_kvector[3*nbq+2] = mybrillouin->kvector[3*nb+2];
   }

   //std::size_t aligned_size = (nfft3d * sizeof(int) + Alignment - 1) & ~(Alignment - 1);
   std::size_t aligned_size = (nfft3d);
   //masker    = new (std::nothrow) int*[nbrillq+1]();
   //packarray = new (std::nothrow) int*[nbrillq+1]();
   masker.resize(nbrillq+1);
   packarray.resize(nbrillq+1);
   for (auto nbq=0; nbq<=nbrillq; ++nbq)
   {
      //masker[nbq]    = new (std::nothrow) int[aligned_size]();
      //packarray[nbq] = new (std::nothrow) int[aligned_size]();
      masker[nbq].resize(aligned_size, 0.0);
      packarray[nbq].resize(aligned_size, 0.0);
      for (auto k=0; k<(nfft3d); ++k) 
      {
         masker[nbq][k] = 1;
         packarray[nbq][k] = 0;
      }
   }
 
   c3db::parall->Barrier();

   double kv[3];
   for (auto nb=0; nb<=nbrillq; ++nb) 
   {
      nwave[nb] = 0;
      if (nb == 0)
      {
         ggcut = lattice->eggcut();
         kv[0] = 0.0;
         kv[1] = 0.0;
         kv[2] = 0.0;
      }
      else
      {
         ggcut = lattice->wggcut();
         kv[0] = p_kvector[3*(nb-1)];
         kv[1] = p_kvector[3*(nb-1)+1];
         kv[2] = p_kvector[3*(nb-1)+2];
      }
     
      for (auto k3 = (-nzh+1); k3<nzh; ++k3)
      for (auto k2 = (-nyh+1); k2<nyh; ++k2)
      for (auto k1 = (-nxh+1); k1<nxh; ++k1) 
      {
         auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2) + kv[0];
         auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2) + kv[1];
         auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2) + kv[2];
         auto i = k1; if (i < 0) i += nx;
         auto j = k2; if (j < 0) j += ny;
         auto k = k3; if (k < 0) k += nz;
        
         auto indx = cijktoindex(i, j, k);
         auto p = cijktop(i, j, k);
         if (p == c3db::parall->taskid_i()) 
         {
            auto gg = gx*gx + gy*gy + gz*gz;
            gg = gg - ggcut;
            if (gg < (-eps)) 
            {
               masker[nb][indx] = 0;
               ++nwave[nb];
            }
         }
      }
      nwave_entire[nb] = nwave[nb];
      nwave_entire[nb] = c3db::parall->ISumAll(1, nwave_entire[nb]);
   }
 
   for (auto nb=0; nb<=nbrillq; ++nb) 
   {
      nidb2[nb] = 0;
     
      // k=(k1,k2,k3)
      for (auto k=(-nzh+1); k<=(nzh-1); ++k)
      for (auto j=(-nyh+1); j<=(nyh-1); ++j)
      for (auto i=(-nxh+1); i<=(nxh-1); ++i) 
      {
         int k1 = i;
         int k2 = j;
         int k3 = k;
         if (k1 < 0) k1 += nx;
         if (k2 < 0) k2 += ny;
         if (k3 < 0) k3 += nz;
         indx = cijktoindex(k1, k2, k3);
         p = cijktop(k1, k2, k3);
         if (p == c3db::parall->taskid_i())
         {
            if (!masker[nb][indx]) 
            {
               packarray[nb][nidb2[nb]] = indx;
               ++nidb2[nb];
            }
         }
      }
      nwave_in[nb] = nidb2[nb];
   }
 
   // if (control.balance())
   if (balance0) 
   {
      balanced = 1;
      mybalance = new CBalance(c3db::parall, nbrillq+1, nwave_in, nwave_out);
   } 
   else 
   {
      balanced = 0;
      for (auto nbq=0; nbq<=nbrillq; ++nbq)
          nwave_out[nbq] = nwave_in[0];
   }

   for (auto nbq=0; nbq<=nbrillq; ++nbq)
   {
      nidb[nbq] = nidb2[nbq] + (nwave_out[nbq] - nwave_in[nbq]);
      nwave_all[nbq] = nidb2[nbq];
   }
   c3db::parall->Vector_ISumAll(1, nbrillq+1, nwave_all.data());

   nidb1_max = 0; 
   for (auto nbq=0; nbq<nbrillq; ++nbq)
      if (nidb[nbq+1] > nidb1_max) nidb1_max = nidb[nbq+1];


   zero_row3 = new (std::nothrow) bool*[nbrillq+1];
   zero_row2 = new (std::nothrow) bool*[nbrillq+1];
   zero_slab23 = new (std::nothrow) bool*[nbrillq+1];
   //zero_row3.resize(nbrillq+1);
   //zero_row2.resize(nbrillq+1);
   //zero_slab23.resize(nbrillq+1);

   if (maptype == 1) 
   {
      // Calculate aligned memory size
      for (auto nbq=0; nbq<=nbrillq; ++nbq)
      {
         //zero_row3[nbq] = new (std::nothrow) bool[(nx * nq + Alignment - 1) & ~(Alignment - 1)];
         //zero_row2[nbq] = new (std::nothrow) bool[(nx * nq + Alignment - 1) & ~(Alignment - 1)];
         zero_row3[nbq] = new (std::nothrow) bool[(nx * nq)];
         zero_row2[nbq] = new (std::nothrow) bool[(nx * nq)];
         zero_slab23[nbq] = new (std::nothrow) bool[nx];
         //zero_row3[nbq].resize(nx*nq,false);
         //zero_row2[nbq].resize(nx*nq,false);
         //zero_slab23[nbq].resize(nx,false);
      }
     
      //zero_arow3 = new bool[(nx*ny + Alignment - 1) & ~(Alignment - 1)];
      zero_arow3 = new bool[(nx*ny)];
      for (auto nb=0; nb<=nbrillq; ++nb) 
      {
         if (nb == 0)
         {
            ggcut = lattice->eggcut();
            kv[0] = 0.0;
            kv[1] = 0.0;
            kv[2] = 0.0;
         }
         else
         {
            ggcut = lattice->wggcut();
            kv[0] = p_kvector[3*(nb-1)];
            kv[1] = p_kvector[3*(nb-1)+1];
            kv[2] = p_kvector[3*(nb-1)+2];
         }
      }
     
      //zero_arow3 = new bool[(nx*ny + Alignment - 1) & ~(Alignment - 1)];
      zero_arow3 = new bool[(nx*ny)];
      for (auto nb=0; nb<=nbrillq; ++nb) 
      {
         if (nb == 0)
         {
            ggcut = lattice->eggcut();
            kv[0] = 0.0;
            kv[1] = 0.0;
            kv[2] = 0.0;
         }
         else
         {
            ggcut = lattice->wggcut();
            kv[0] = p_kvector[3*(nb-1)];
            kv[1] = p_kvector[3*(nb-1)+1];
            kv[2] = p_kvector[3*(nb-1)+2];
         }
        
         /* find zero_row3 - (i,j,*) rows that are zero */
         for (auto i=0; i<((nxh+1)*nq); ++i)
            zero_row3[nb][i] = true;
 
         for (auto i=0; i<((nxh+1)*ny); ++i)
            zero_arow3[i] = true;
        
         for (auto k2 = (-nyh+1); k2<nyh; ++k2)
         for (auto k1 = (-nxh+1); k1<nxh; ++k1) 
            {
               auto i = k1;
               auto j = k2;
               if (i < 0) i += nx;
               if (j < 0) j += ny;
               zrow = true;
               for (auto k3 = (-nzh+1); k3<nzh; ++k3) 
               {
                  auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2) + kv[0];
                  auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2) + kv[1];
                  auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2) + kv[2];
                  auto gg = gx*gx + gy*gy + gz*gz;
                  gg = gg - ggcut;
                  if (gg < (-eps))
                     zrow = false;
               }
               if (!zrow) 
               {
                  zero_arow3[i + nx*j] = false;
                  q = cijktoq1(0, j, 0);
                  p = cijktop1(0, j, 0);
                  if (p == c3db::parall->taskid_i()) 
                  {
                     zero_row3[nb][i + nx*q] = false;
                  }
               }
            }
         // call c3db_c_ptranspose_jk_init(nb,log_mb(zero_arow3(1)))
         c3db::c_pctranspose_jk_init(nb, zero_arow3);
        

         /* find zero_slab2 - (i,*,*) slabs that are zero */
         for (auto i=0; i<nx; ++i)
            zero_slab23[nb][i] = true;
        
         for (auto k1=(-nxh+1); k1<nxh; ++k1)
         {
            auto i = k1; if (i<0) i += nx;
            yzslab = true;
            for (auto k3=(-nzh+1); k3<nzh; ++k3)
            for (auto k2=(-nyh+1); k2<nyh; ++k2) 
            {
               auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2) + kv[0];
               auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2) + kv[1];
               auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2) + kv[2];
               auto gg = gx*gx + gy*gy + gz*gz;
               gg = gg - ggcut;
               if (gg < (-eps))
                  yzslab = false;
            }
            if (!yzslab) zero_slab23[nb][i] = false;
         }

         /* initialize zero_row2 */
         for (auto i=0; i<(nx*nq); ++i)
            zero_row2[nb][i] = false;
        
         /* find zero_row2 - (i,*,k) rows that are zero after fft of (i,j,*) */
         for (auto k=0; k<nz; ++k)
         for (auto i=0; i<nx; ++i) 
         {
            q = cijktoq(i, 0, k);
            p = cijktop(i, 0, k);
            if (p == c3db::parall->taskid_i())
               zero_row2[nb][q] = zero_slab23[nb][i];
         }
      }
     
      delete[] zero_arow3;
   } 
   else 
   {
      for (auto nbq=0; nbq<=nbrillq; ++nbq)
      {
         //zero_row3[nbq]   = new (std::nothrow) bool[(nq3 + Alignment - 1) & ~(Alignment - 1)]();
         //zero_row2[nbq]   = new (std::nothrow) bool[(nq2 + Alignment - 1) & ~(Alignment - 1)]();
         zero_row3[nbq]   = new (std::nothrow) bool[(nq3)]();
         zero_row2[nbq]   = new (std::nothrow) bool[(nq2)]();
         zero_slab23[nbq] = new (std::nothrow) bool[nx]();
         //zero_row3[nbq].resize(nq3,false);
         //zero_row2[nbq].resize(nq2,false);
         //zero_slab23[nbq].resize(nx,false);
      }
     
      //zero_arow2 = new (std::nothrow) bool[(nx*nz + Alignment - 1) & ~(Alignment - 1)]();
      //zero_arow3 = new (std::nothrow) bool[(nx*ny + Alignment - 1) & ~(Alignment - 1)]();
      zero_arow2 = new (std::nothrow) bool[(nx*nz)]();
      zero_arow3 = new (std::nothrow) bool[(nx*ny)]();
     
      for (auto nb=0; nb<=nbrillq; ++nb) 
      {
         if (nb==0)
         {
            ggcut = lattice->eggcut();
            kv[0] = 0.0;
            kv[1] = 0.0;
            kv[2] = 0.0;
         }
         else
         {
            ggcut = lattice->wggcut();
            kv[0] = p_kvector[3*(nb-1)];
            kv[1] = p_kvector[3*(nb-1)+1];
            kv[2] = p_kvector[3*(nb-1)+2];
         }

         // find zero_row3 - (i,j,*) rows that are zero
         for (auto q=0; q<nq3; ++q)
            zero_row3[nb][q] = true;

         for (auto q=0; q < nx*ny; ++q)
            zero_arow3[q] = true;

         for (auto k2 = (-nyh+1); k2<nyh; ++k2)
         for (auto k1 = (-nxh+1); k1<nxh; ++k1)
         {
            auto i = k1;
            auto j = k2;
            if (i < 0) i += nx;
            if (j < 0) j += ny;
            zrow = true;
            for (auto k3=(-nzh+1); k3<nzh; ++k3) 
            {
               auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2) + kv[0];
               auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2) + kv[1];
               auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2) + kv[2];
               auto gg = gx*gx + gy*gy + gz*gz;
               gg = gg - ggcut;
               if (gg < (-eps))
                  zrow = false;
            }

            if (!zrow) 
            {
              // zero_arow3[i-1+(nxh+1)*(j-1)] = 0;
              zero_arow3[i + nx*j] = false;
              q = cijktoq(i, j, 0);
              p = cijktop(i, j, 0);
              if (p == c3db::parall->taskid_i()) 
              {
                 zero_row3[nb][q] = false;
              }
            }
         }

         /* find zero_slab23 - (i,*,*) slabs that are zero */
         for (auto i=0; i<nx; ++i)
            zero_slab23[nb][i] = true;
        
         for (auto k1=(-nxh+1); k1<nxh; ++k1)
         {
            auto i = k1; if (i < 0) i += nx;
            yzslab = true;
            for (auto k3=(-nzh+1); k3<nzh; ++k3)
            for (auto k2=(-nyh+1); k2<nyh; ++k2) 
            {
               auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2) + kv[0];
               auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2) + kv[1];
               auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2) + kv[2];
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
         for (auto i=0; i<nx; ++i) 
         {
            q = cijktoq1(i, 0, k);
            p = cijktop1(i, 0, k);
            zero_arow2[i + nx*k] = zero_slab23[nb][i];
            if (p == c3db::parall->taskid_i())
               zero_row2[nb][q] = zero_slab23[nb][i];
         }
        
         c3db::c_pctranspose_ijk_init(nb, zero_arow2, zero_arow3);
      }
      delete[] zero_arow3;
      delete[] zero_arow2;
   }

      //Gpack[nbq] = new (std::nothrow) double[(3*(nidb[nbq]) + Alignment - 1) & ~(Alignment - 1)]();
      //Gpack[nbq] = (std::nothrow) double[(3*(nidb[nbq]))]();
      //Gpack[nbq] = new (std::nothrow) double[(3*(nidb[nbq]))]();
   for (auto nbq=0; nbq<=nbrillq; ++nbq)
      Gpack[nbq].resize(3*nidb[nbq], 0.0);

   double *Gtmp = new (std::nothrow) double[nfft3d]();

   int one = 1;
   for (auto nb=0; nb<=nbrillq; ++nb) 
   {
      for (auto i=0; i<3; ++i) 
      {
         std::memcpy(Gtmp,Garray+i*nfft3d,nfft3d*sizeof(double));
        
         this->r_pack(nb,Gtmp);
         this->rr_pack_copy(nb,Gtmp,Gpack[nb].data()+i*(nidb[nb]));
      }
   }
 
   delete [] Gtmp;
 
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
         int indx = cijktoindex2(i, j, k);
         int p = cijktop2(i, j, k);
        
         if (p == c3db::parall->taskid_i()) 
         {
            r_grid[3*indx]   = a[0]*k1 + a[3]*k2 + a[6]*k3;
            r_grid[3*indx+1] = a[1]*k1 + a[4]*k2 + a[7]*k3;
            r_grid[3*indx+2] = a[2]*k1 + a[5]*k2 + a[8]*k3;
         }
      }
   }

   /* initialize print_grid */
   //std::memset(nidb_print.data(),0,nbrillouin*sizeof(int));
   //std::memset(nwave_all_print.data(),0,nbrillouin*sizeof(int));
   nidb_print.resize(nbrillouin,0);
   nwave_all_print.resize(nbrillouin,0);
   for (auto nbq=0; nbq<nbrillq; ++nbq)
   {
      auto nb = k1db::ktoptok(nbq);
      nidb_print[nb]      = nidb[1+nbq];
      nwave_all_print[nb] = nwave_all[1+nbq];
   }
   c3db::parall->Vector_ISumAll(3,nbrillouin,nidb_print.data());
   c3db::parall->Vector_ISumAll(3,nbrillouin,nwave_all_print.data());

 
   /* initialize pfft3 queues */
   nffts_max = 1;
   staged_gpu_fft_pipeline = staged_gpu_fft_pipeline0 && c3db::mygdevice.has_gpu();

#ifdef NWPW_SYCL
   staged_gpu_fft_pipeline = false;
#endif

   //aqmax = 5;
   aqmax = pfft3_qsize0;
   if (staged_gpu_fft_pipeline)
   {
      //std::cout << "Using Staged GPU fft pipeline!!!!" << std::endl;
      aqmax += 6;
      c3db::mygdevice.batch_fft_pipeline_mem_init(aqmax,2*nfft3d);
   }
 
   aqsize = 0;
   alast_index = aqmax - 1;
   aqindx   = new (std::nothrow) int[aqmax]();
   aqstatus = new (std::nothrow) int[aqmax]();
   aqnffts  = new (std::nothrow) int[aqmax]();
   aqnbb    = new (std::nothrow) int[aqmax]();
   atmp = new (std::nothrow) double[2*aqmax*2*nfft3d*nffts_max]();
 
   bqmax = pfft3_qsize0;
   if (staged_gpu_fft_pipeline) bqmax += 6;
   
   //bqmax = aqmax;
   bqsize = 0;
   blast_index = bqmax - 1;
   bqindx   = new (std::nothrow) int[bqmax]();
   bqstatus = new (std::nothrow) int[bqmax]();
   bqnffts  = new (std::nothrow) int[bqmax]();
   bqnbb    = new (std::nothrow) int[bqmax]();
   btmp = new (std::nothrow) double[2*bqmax*2*nfft3d*nffts_max]();
 
   /* initialize async buffer data for pfft */
   for (auto q=0; q<aqmax; ++q)
      c3db::parall->astart(4+q, 2*c3db::parall->np_i()+1);

}

CGrid::CGrid(Parallel *inparall, Lattice *inlattice, Control2 &control, Brillouin *inbrillouin)
    : CGrid(inparall, inlattice, control.mapping(), control.balance(),
            control.ngrid(0), control.ngrid(1), control.ngrid(2),
            control.pfft3_qsize(), control.staged_gpu_fft(), inbrillouin) {}

/********************************
 *                              *
 *       CGrid:c_unpack         *
 *                              *
 ********************************/
void CGrid::c_unpack(const int nb, double *a) 
{

   int nn = 2*(nidb2[nb]);

   if (balanced)
      mybalance->c_unbalance(nb, a);
 
   std::memcpy(c3db::c3db_tmp1,a,nn*sizeof(double));
   std::memset(a, 0, 2*nfft3d*sizeof(double));
   c_bindexcopy(nidb2[nb],packarray[nb].data(),c3db::c3db_tmp1,a);
 
}

/********************************
 *                              *
 *       CGrid:c_pack           *
 *                              *
 ********************************/
void CGrid::c_pack(const int nb, double *a) 
{

   std::memcpy(c3db::c3db_tmp1,a,2*nfft3d*sizeof(double));
   std::memset(a,  0,2*nfft3d*sizeof(double));

   c_aindexcopy(nidb2[nb],packarray[nb].data(),c3db::c3db_tmp1,a);

   if (balanced)
      mybalance->c_balance(nb, a);

   //delete [] tmp;
   return;
}

/********************************
 *                              *
 *       CGrid:cc_pack_copy     *
 *                              *
 ********************************/
void CGrid::cc_pack_copy(const int nb, const double *a, double *b)
{
   //int one = 1;
   // int ng  = 2*(nidb[nb]);
   int ng = 2*(nidb[nb]);

   // DCOPY_PWDFT(ng,a,one,b,one);
   std::memcpy(b,a,ng*sizeof(double));
}

/********************************
 *                              *
 *       CGrid:cc_pack_dot      *
 *                              *
 ********************************/
double CGrid::cc_pack_dot(const int nb, double *a, double *b) 
{
   int one = 1;
   // int ng  = 2*(nidb[nb]);
   int ng  = 2*(nidb[nb]);
 
   double tsum  = DDOT_PWDFT(ng, a,one,b,one);
 
   return c3db::parall->SumAll(1,tsum);
}


/********************************
 *                              *
 *       CGrid:cc_pack_idot     *
 *                              *
 ********************************/
double CGrid::cc_pack_idot(const int nb, double *a, double *b) 
{
   int one = 1;
   int ng = 2*(nidb[nb]);
   double tsum = 0.0;
 
   tsum = DDOT_PWDFT(ng, a, one, b, one);
 
   return tsum;
}

/********************************
 *                              *
 *       CGrid:cc_pack_zdot     *
 *                              *
 ********************************/
std::complex<double> CGrid::cc_pack_zdot(const int nb, double *a, double *b) 
{
   int one = 1;
   // int ng  = 2*(nidb[nb]);
   int ng  = (nidb[nb]);
 
   std::complex<double> tsum  = ZDOTC_PWDFT(ng, a,one,b,one);
   c3db::parall->Vector_SumAll(1,2,reinterpret_cast<double*>(&tsum));

   return tsum;
}

/********************************
 *                              *
 *      CGrid:cc_pack_izdot     *
 *                              *
 ********************************/
std::complex<double> CGrid::cc_pack_izdot(const int nb, double *a, double *b) 
{
   int one = 1;
   // int ng  = 2*(nidb[nb]);
   int ng  = (nidb[nb]);
 
   std::complex<double> tsum  = ZDOTC_PWDFT(ng, a,one,b,one);

   return tsum;
}


/********************************
 *                              *
 *      CGrid:cc_pack_inzdot    *
 *                              *
 ********************************/
void CGrid::cc_pack_inzdot(const int nb, const int nn, double *a, double *b, double *sum)
{  
   int one = 1;
   int ng  = (nidb[nb]);

   for (int i=0; i<nn; ++i) 
   {
      std::complex<double> tsum = ZDOTC_PWDFT(ng, a+i*2*ng,one,b,one);
      sum[2*i]   = tsum.real();
      sum[2*i+1] = tsum.imag();
   }
}


/********************************
 *                              * 
 *   CGrid:cc_pack_inprjzdot    *
 *                              *
 ********************************/
void CGrid::cc_pack_inprjzdot(const int nb, int nn, int nprj, double *a,
                              double *b, double *sum) 
{  
   int npack1_max = (nidb1_max);
   int ng = (nidb[nb]);
   double rone[2] = {1.0,0.0};
   double rzero[2] = {0.0,0.0}; 
   
   c3db::mygdevice.CN2_zgemm(nn, nprj, ng, npack1_max, rone, a, b, rzero, sum);
   //c3db::mygdevice.CN2_stride_zgemm(nn, nprj, npack1, ng, rone, a, b, rzero, sum);
}  



/********************************
 *                              *
 *       CGrid:rr_pack_copy     *
 *                              *
 ********************************/
void CGrid::rr_pack_copy(const int nb, const double *a, double *b)
{
   //int one = 1;
   // int ng  = 2*(nidb[nb]);
   int ng = (nidb[nb]);

   // DCOPY_PWDFT(ng,a,one,b,one);
   std::memcpy(b,a,ng*sizeof(double));
}

/********************************
 *                              *
 *       CGrid:tt_pack_copy     *
 *                              *
 ********************************/
/*
void CGrid::tt_pack_copy(const int nb, const double *a, double *b)

   //int one = 1;
   // int ng  = 2*(nidb[nb]);
   int ng = (nidb[nb]);

   // DCOPY_PWDFT(ng,a,one,b,one);
   std::memcpy(b,a,ng*sizeof(double));
}
*/

/********************************
 *                              *
 *       CGrid::tt_pack_copy    * 
 *                              *
 ********************************/
void CGrid::tt_pack_copy(const int nb, const double *a, double *b) {
    // Ensure nb is within bounds
    if (nb < 0 ) {
        std::cerr << "Error: 'nb' index out of bounds in tt_pack_copy" << std::endl;
        return;
    }

    // Ensure a and b are valid pointers
    if (a == nullptr || b == nullptr) {
        std::cerr << "Error: Null pointer passed to tt_pack_copy" << std::endl;
        return;
    }

    // Calculate the number of elements to copy
    int ng = nidb[nb];

    // Ensure the size of the source and destination buffers is sufficient
    // This requires knowing the size of the allocated buffers a and b, assuming they are valid and properly sized
    if (ng <= 0) {
        std::cerr << "Error: 'ng' must be positive in tt_pack_copy" << std::endl;
        return;
    }


    // Perform the copy safely
    std::memcpy(b, a, ng * sizeof(double));
}




/********************************
 *                              *
 *       CGrid:tt_pack_dot      *
 *                              *
 ********************************/
double CGrid::tt_pack_dot(const int nb, double *a, double *b)
{
   int one = 1;
   // int ng  = 2*(nidb[nb]);
   int ng  = (nidb[nb]);

   double tsum  = DDOT_PWDFT(ng, a,one,b,one);

   return c3db::parall->SumAll(1,tsum);
}


/********************************
 *                              *
 *       CGrid:tt_pack_idot     *
 *                              *
 ********************************/
double CGrid::tt_pack_idot(const int nb, double *a, double *b)
{
   int one = 1;
   int ng = (nidb[nb]);
   double tsum = 0.0;

   tsum = DDOT_PWDFT(ng, a, one, b, one);

   return tsum;
}


/********************************
 *                              *
 *       CGrid:rr_pack_idot     *
 *                              *
 ********************************/
double CGrid::rr_pack_idot(const int nb, double *a, double *b) 
{
   int one = 1;
   int ng = (nidb[nb]);
   double tsum = 0.0;
 
   tsum = DDOT_PWDFT(ng, a, one, b, one);
 
   return tsum;
}

/********************************
 *                              *
 *       CGrid:cc_pack_indot    *
 *                              *
 ********************************/
void CGrid::cc_pack_indot(const int nb, const int nn, double *a, double *b, double *sum) 
{
   int one = 1;
   int ng = 2 * (nidb[nb]);
 
   for (int i = 0; i < nn; ++i) {
     sum[i] = DDOT_PWDFT(ng, a+i*ng, one, b, one);
   }
}

/********************************
 *                              *
 *    CGrid:cc_pack_inprjdot    *
 *                              *
 ********************************/
void CGrid::cc_pack_inprjdot(const int nb, int nn, int nprj, double *a,
                             double *b, double *sum) 
{
   int ng = 2 * (nidb[nb]);
   int one = 1;
   double rone = 1.0;
   double rzero = 0.0;

   c3db::mygdevice.TN_dgemm(nn, nprj, ng, rone, a, b, rzero, sum);
}

/********************************
 *                              *
 *       CGrid:r_unpack         *
 *                              *
 ********************************/
void CGrid::r_unpack(const int nb, double *a) 
{
   int nn = (nidb2[nb]);
   if (balanced)
      mybalance->r_unbalance(nb, a);
 
   std::memcpy(c3db::c3db_tmp1, a, nn * sizeof(double));
   std::memset(a, 0, nfft3d*sizeof(double));
 
   t_bindexcopy(nidb2[nb],packarray[nb].data(),c3db::c3db_tmp1,a);
}

/********************************
 *                              *
 *       CGrid:t_unpack         *
 *                              *
 ********************************/
void CGrid::t_unpack(const int nb, double *a)
{
   int nn = (nidb2[nb]);
   if (balanced)
      mybalance->t_unbalance(nb, a);

   std::memcpy(c3db::c3db_tmp1, a, nn * sizeof(double));
   std::memset(a, 0, nfft3d*sizeof(double));

   t_bindexcopy(nidb2[nb],packarray[nb].data(),c3db::c3db_tmp1,a);
}



/********************************
 *                              *
 *       CGrid:r_pack           *
 *                              *
 ********************************/
void CGrid::r_pack(const int nb, double *a) 
{
   std::memcpy(c3db::c3db_tmp1,a,nfft3d*sizeof(double));
   std::memset(a,0,nfft3d*sizeof(double));

   t_aindexcopy(nidb2[nb],packarray[nb].data(),c3db::c3db_tmp1,a);

   if (balanced)
      mybalance->r_balance(nb, a);
}

/********************************
 *                              *
 *       CGrid:t_pack           *
 *                              *
 ********************************/
void CGrid::t_pack(const int nb, double *a)
{
   std::memcpy(c3db::c3db_tmp1,a,nfft3d*sizeof(double));
   std::memset(a,0,nfft3d*sizeof(double));

   t_aindexcopy(nidb2[nb],packarray[nb].data(),c3db::c3db_tmp1,a);

   if (balanced)
      mybalance->t_balance(nb, a);
}



/********************************
 *                              *
 *       CGrid:r_pack_nzero     *
 *                              *
 ********************************/
void CGrid::r_pack_nzero(const int nb, const int n, double *a) 
{
   int ng = n * (nidb[nb]);
   std::memset(a, 0, ng * sizeof(double));
}

/********************************
 *                              *
 *       CGrid:t_pack_nzero     *
 *                              *
 ********************************/
void CGrid::t_pack_nzero(const int nb, const int n, double *a)
{
   int ng = n * (nidb[nb]);
   std::memset(a, 0, ng * sizeof(double));
}

/********************************
 *                              *
 *     CGrid:t_pack_max_nzero   *
 *                              *
 ********************************/
void CGrid::t_pack_max_nzero(const int n, double *a)
{
   int ng = n*(nidb1_max);
   std::memset(a, 0, ng*sizeof(double));
}




/********************************
 *                              *
 *       CGrid:i_pack           *
 *                              *
 ********************************/
void CGrid::i_pack(const int nb, int *a) 
{
   int *tmp = new (std::nothrow) int[nfft3d];
 
   for (auto i=0; i<nfft3d; ++i) tmp[i] = a[i];
   for (auto i=0; i<nfft3d; ++i) a[i] = 0;
 
   i_aindexcopy(nidb2[nb], packarray[nb].data(), tmp, a);
 
   if (balanced)
      mybalance->i_balance(nb, a);
 
   delete[] tmp;
}

/********************************
 *                              *
 *       CGrid:ii_pack_copy     *
 *                              *
 ********************************/
void CGrid::ii_pack_copy(const int nb, int *a, int *b) 
{
   int ng = nidb[nb];
   for (auto i=0; i<ng; ++i)
      b[i] = a[i];
}

/********************************
 *                              *
 *    CGrid:cr_pfft3b           *
 *                              *
 ********************************/
/* This routine performs the operation of a three dimensional
  complex to complex inverse fft:
   A(nx,ny(nb),nz(nb)) <- FFT3^(-1)[A(kx,ky,kz)]
*/
void CGrid::cr_pfft3b(const int nb, double *a)
{

   nwpw_timing_function ftime(1);
   int  indx0, indx2, nn;
   double *tmp2 = c3db::c3db_tmp1;
   double *tmp3 = c3db::c3db_tmp2;
 
   auto nxy = nx * ny;
   auto nxz = nx * nz;
 
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {

     /***************************************************
      ***     do fft along kz dimension               ***
      ***     A(kx,nz,ky) <- fft1d^(-1)[A(kx,kz,ky)]  ***
      ***************************************************/
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nx; ++i) 
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
                  indx3 += nx;
               }
               ++nn;
            }
            ++indx2;
         }
         indx0 += nxz;
      }
 
      c3db::mygdevice.batch_cfftz_tmpz(c3db::fft_tag,false, nz, nn, 2*nfft3d, tmp2, c3db::tmpz);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,false,nz,nn,nz,tmp2,c3db::backward_z,c3db::tmpz,2);
 
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nx; ++i) 
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
                  indx3 += nx;
               }
               ++nn;
            }
            ++indx2;
         }
         indx0 += nxz;
      }
      
      
      /***********************************************
       ***         Do a ptranspose of A            ***
       ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
       ************************************************/
      c3db::c_pctranspose1_jk(nb, a, tmp2, tmp3);
      
      /*************************************************
       ***        do fft along ky dimension          ***
       ***    A(kx,ny,nz) <- fft1d^(-1)[A(kx,ky,nz)] ***
       *************************************************/
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nx; ++i) 
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
                  indx3 += nx;
               }
               ++nn;
            }
            ++indx2;
         }
         indx0 += nxy;
      }
      
      c3db::mygdevice.batch_cffty_tmpy(c3db::fft_tag,false, ny, nn, 2*nfft3d, tmp2, c3db::tmpy);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,false,ny,nn,ny,tmp2,c3db::backward_y,c3db::tmpy,1);
      
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nx; ++i) 
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
                  indx3 += nx;
               }
               ++nn;
            }
            ++indx2;
         }
         indx0 += nxy;
      }
      
      
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      c3db::mygdevice.batch_cfftx_tmpx(c3db::fft_tag,false, nx, ny * nq, 2*nfft3d, a, c3db::tmpx);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,false,nx,ny*nq,nx,a,backward_x,c3db::tmpx,0);
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
      c3db::mygdevice.batch_cfftz_tmpz_zero(c3db::fft_tag,false, nz, nq3, 1,2*nfft3d, a, c3db::tmpz, zero_row3[nb]);
      //c3db::mygdevice.batch_cfft_zero(c3db::fft_tag,false,nz,nq3,1,2*nfft3d,a,c3db::backward_z,c3db::tmpz,zero_row3[nb],2);
     
      c3db::c_pctranspose_ijk(nb, 2, a, tmp2, tmp3);
      
      /************************************************
       ***     do fft along ky dimension            ***
       ***   A(ny,nz,kx) <- fft1d^(-1)[A(ky,nz,kx)] ***
       ************************************************/
      c3db::mygdevice.batch_cffty_tmpy_zero(c3db::fft_tag,false,ny,nq2,1,2*nfft3d,a,c3db::tmpy,zero_row2[nb]);
      //c3db::mygdevice.batch_cfft_zero(c3db::fft_tag,false,ny,nq2,1,2*nfft3d,a,c3db::backward_y,c3db::tmpy,zero_row2[nb],1);
      
      c3db::c_pctranspose_ijk(nb, 3, a, tmp2, tmp3);
      
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/

      c3db::mygdevice.batch_cfftx_tmpx(c3db::fft_tag,false, nx, nq1, 2*nfft3d, a, c3db::tmpx);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,false,nx,nq1,nx,a,c3db::backward_x,c3db::tmpx,0);
       
      if (nfft3d_map < nfft3d)
         std::memset(a + 2*nfft3d_map, 0, 2*(nfft3d - nfft3d_map) * sizeof(double));
   }
 
}

/********************************
 *                              *
 *       CGrid::rc_pfft3f       *
 *                              *
 ********************************/
/*
   This routine performs the operation of a three dimensional
   complex to complex fft:
      A(kx,ky,kz) <- FFT3[A(nx(nb),ny(nb),nz(nb))]
*/
void CGrid::rc_pfft3f(const int nb, double *a)
{
   nwpw_timing_function ftime(1);
   int  indx0, indx2, nn;
   double *tmp2 = c3db::c3db_tmp1;
   double *tmp3 = c3db::c3db_tmp2;
 
   auto nxy = nx * ny;
   auto nxz = nx * nz;
 
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      /********************************************
       ***     do fft along nx dimension        ***
       ***   A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]  ***
       ********************************************/

      c3db::mygdevice.batch_cfftx_tmpx(c3db::fft_tag,true,nx,ny*nq,2*nfft3d,a,c3db::tmpx);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,true,nx,ny*nq,nx,a,forward_x,c3db::tmpx,0);

      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]  ***
       ********************************************/
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nx; ++i) 
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
                  indx3 += nx;
               }
               ++nn;
            }
            ++indx2;
         }
         indx0 += nxy;
      }
     
      c3db::mygdevice.batch_cffty_tmpy(c3db::fft_tag,true, ny, nn, 2*nfft3d, tmp2, c3db::tmpy);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,true,ny,nn,ny,tmp2,forward_y,c3db::tmpy,1);
     
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nx; ++i) 
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
                  indx3 += nx;
               }
               ++nn;
            }
            ++indx2;
         }
         indx0 += nxy;
      }
     
     
      /********************************************
       ***         Do a transpose of A          ***
       ***      A(ky,nz,ky) <- A(kx,ky,nz)      ***
       ********************************************/
      c3db::c_pctranspose2_jk(nb,a,tmp2,tmp3);

      /********************************************
       ***     do fft along nz dimension        ***
       ***   A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]  ***
       ********************************************/
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nx; ++i) 
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
                  indx3 += nx;
               }
               ++nn;
            }
            ++indx2;
         }
         indx0 += nxz;
      }
     
      c3db::mygdevice.batch_cfftz_tmpz(c3db::fft_tag,true, nz, nn, 2*nfft3d, tmp2, c3db::tmpz);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,true,nz,nn,nz,tmp2,c3db::forward_z,c3db::tmpz,2);
     
      indx0 = 0;
      indx2 = 0;
      nn = 0;
      for (auto q=0; q<nq; ++q) 
      {
         for (auto i=0; i<nx; ++i) 
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
                  indx3 += nx;
               }
               ++nn;
            }
            ++indx2;
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
      c3db::mygdevice.batch_cfftx_tmpx(c3db::fft_tag,true, nx, nq1, 2*nfft3d, a, c3db::tmpx);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,true,nx,nq1,nx,a,c3db::forward_x,c3db::tmpx,0);

      c3db::c_pctranspose_ijk(nb, 0, a, tmp2, tmp3);

      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
       ********************************************/
      c3db::mygdevice.batch_cffty_tmpy_zero(c3db::fft_tag,true,ny,nq2,1,2*nfft3d,a,c3db::tmpy,zero_row2[nb]);
      //c3db::mygdevice.batch_cfft_zero(c3db::fft_tag,true,ny,nq2,1,2*nfft3d,a,c3db::forward_y,c3db::tmpy,zero_row2[nb],1);

      c3db::c_pctranspose_ijk(nb, 1, a, tmp2, tmp3);

      /********************************************
       ***     do fft along nz dimension        ***
       ***   A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]  ***
       ********************************************/
      c3db::mygdevice.batch_cfftz_tmpz_zero(c3db::fft_tag,true, nz, nq3, 1,2*nfft3d, a, c3db::tmpz, zero_row3[nb]);
      //c3db::mygdevice.batch_cfft_zero(c3db::fft_tag,true,nz,nq3,1,2*nfft3d,a,c3db::forward_z,c3db::tmpz,zero_row3[nb],2);
   }
 
}

/********************************
 *                              *
 *     CGrid::c_unpack_start    *
 *                              *
 ********************************/
void CGrid::c_unpack_start(const int nffts, const int nb, double *tmp1, double *tmp2,
                           const int request_indx, const int msgtype) 
{
   if (balanced)
      for (auto s=0; s<nffts; ++s)
         mybalance->c_unbalance_start(nffts, nb, tmp1+s*n2ft3d, request_indx, msgtype);
}

/********************************
 *                              *
 *     CGrid::c_unpack_mid      *
 *                              *
 ********************************/
void CGrid::c_unpack_mid(const int nffts, const int nb, double *tmp1, double *tmp2,
                         const int request_indx, const int msgtype) 
{
   if (balanced)
     mybalance->c_unbalance_end(nffts, nb, tmp1, request_indx);
 
   std::memcpy(tmp2, tmp1, nffts*2*(nidb2[nb])*sizeof(double));
   std::memset(tmp1, 0, nffts*2*nfft3d*sizeof(double));
 
   for (auto s=0; s<nffts; ++s)
      c_bindexcopy((nidb2[nb]), packarray[nb].data(), tmp2+s*2*nfft3d, tmp1+s*2*nfft3d);
}

/********************************
 *                              *
 *     CGrid::c_unpack_end      *
 *                              *
 ********************************/
void CGrid::c_unpack_end(const int nffts, const int nb, double *tmp1, double *tmp2,
                         const int request_indx) {
}

/********************************
 *                              *
 *        CGrid::pfftbz         *
 *                              *
 ********************************/
void CGrid::pfftbz(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx) 
{

   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      auto nxz = nx*nz;
     
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
            for (auto i=0; i<nx; ++i) 
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
                     indx3 += nx;
                  }
                  nn += 1;
               }
               ++indx2;
            }
            indx0 += nxz;
         }
      }
     
      c3db::mygdevice.batch_cfftz_tmpz(c3db::fft_tag,false, nz, nn, 2*nfft3d, tmp2, c3db::tmpz);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,false,nz,nn,nz,tmp2,c3db::backward_z,c3db::tmpz,2);
      // for (auto i=0; i<nn; ++i)
      //    dcfftb_(&nz,tmp2+2*nz*i,c3db::tmpz);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
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
                     indx3 += nx;
                  }
                  nn += 1;
               }
               ++indx2;
            }
            indx0 += nxz;
         }
      }
     
     
      /***********************************************
       ***         Do a ptranspose of A            ***
       ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
       ************************************************/
      c3db::c_pctranspose1_jk_start(nffts, nb, tmp1, tmp2, tmp1, request_indx, 44);
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
      c3db::mygdevice.batch_cfftz_tmpz_zero(c3db::fft_tag,false,nz,nq3,1,2*nfft3d,tmp1,c3db::tmpz,zero_row3[nb]);
     // std::cout << "nz=" << nz << " nq3=" << nq3 << " nffts=" << nffts << " n2ft3d=" << 2*nfft3d 
     //           << " nfft3d=" << nfft3d << std::endl;
     
      //c3db::mygdevice.batch_cfft_zero(c3db::fft_tag,false,nz,nq3,nffts,2*nfft3d,tmp1,c3db::backward_z,c3db::tmpz,zero_row3[nb],2);
     
      c3db::c_pctranspose_ijk_start(nffts,nb,2,tmp1,tmp2,tmp1,request_indx, 45);
      // c3db::c_ptranspose_ijk(nb,2,tmp1,tmp2,tmp1);
   }
}

/********************************
 *                              *
 *        CGrid::pfftby         *
 *                              *
 ********************************/
void CGrid::pfftby(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx) 
{
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      auto nxy = nx*ny;
     
      /***********************************************
       ***         Do a ptranspose of A            ***
       ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
       ************************************************/
      c3db::c_pctranspose1_jk_end(nffts, nb, tmp2, tmp1, request_indx);
     
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
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  auto jj = 0;
                  auto indx3 = 2*i + indx0;
                  auto shift = 2*ny*nn;
                  for (auto j=0; j<ny; ++j) 
                  {
                     tmp1[jj   + shift] = tmp2[indx3];
                     tmp1[jj+1 + shift] = tmp2[indx3 + 1];
                     jj += 2;
                     indx3 += nx;
                  }
                  nn += 1;
               }
               ++indx2;
            }
            indx0 += nxy;
         }
      }
     
      c3db::mygdevice.batch_cffty_tmpy(c3db::fft_tag,false, ny, nn, 2*nfft3d, tmp1, c3db::tmpy);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,false,ny,nn,ny,tmp1,c3db::backward_y,c3db::tmpy,1);
      // for (auto i=0; i<nn; ++i)
      //    dcfftb_(&ny,tmp1+2*ny*i,c3db::tmpy);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  auto jj = 0;
                  auto indx3 = 2*i + indx0;
                  auto shift = 2*ny*nn;
                  for (auto j = 0; j < ny; ++j) 
                  {
                     tmp2[indx3]   = tmp1[jj   + shift];
                     tmp2[indx3+1] = tmp1[jj+1 + shift];
                     jj += 2;
                     indx3 += nx;
                  }
                  nn += 1;
               }
               ++indx2;
            }
            indx0 += nxy;
         }
      }

   }

   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      c3db::c_pctranspose_ijk_end(nffts, nb, 2, tmp2, tmp1, request_indx);
     
      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
       ********************************************/
      c3db::mygdevice.batch_cffty_tmpy_zero(c3db::fft_tag,false,ny,nq2,1,2*nfft3d,tmp2,c3db::tmpy,zero_row2[nb]);
      //c3db::mygdevice.batch_cfft_zero(c3db::fft_tag,false,ny,nq2,nffts,1,2*nfft3d,tmp2,c3db::backward_y,c3db::tmpy,zero_row2[nb],1);
     
      c3db::c_pctranspose_ijk_start(nffts,nb,3,tmp2,tmp1,tmp2,request_indx, 46);
      // c3db::c_ptranspose_ijk(nb,3,tmp2,tmp1,tmp2);
   }
}

/********************************
 *                              *
 *        CGrid::pfftbx         *
 *                              *
 ********************************/
void CGrid::pfftbx(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx) 
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
      c3db::mygdevice.batch_cfftx_tmpx(c3db::fft_tag,false, nx, ny * nq, 2*nfft3d, tmp2, c3db::tmpx);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,false,nx,nffts*ny*nq,nx,tmp2,c3db::backward_x,c3db::tmpx,0);
      std::memcpy(tmp1,tmp2,nffts*2*nfft3d*sizeof(double));
   }
   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      c3db::c_pctranspose_ijk_end(nffts,nb,3,tmp1,tmp2,request_indx);
     
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      c3db::mygdevice.batch_cfftx_tmpx(c3db::fft_tag,false, nx, nq1, 2*nfft3d, tmp1, c3db::tmpx);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,false,nx,nffts*nq1,nx,tmp1,c3db::backward_x,c3db::tmpx,0);
      if (2*nfft3d_map < 2*nfft3d)
         for (auto s=0; s<nffts; ++s)
            std::memset(tmp1+s*2*nfft3d + 2*nfft3d_map, 0, (2*nfft3d-2*nfft3d_map)*sizeof(double));
   }
}

/********************************
 *                              *
 *       CGrid::pfftb_step      *
 *                              *
 ********************************/
void CGrid::pfftb_step(const int step, const int nffts, const int nb, double *a, double *tmp1,
                       double *tmp2, const int request_indx) 
{
   if (step == 0) 
   {
      // c3db::parall->astart(request_indx,parall->np_i());
      // unpack start, tmp1-->tmp1
      for (auto s=0; s<nffts; ++s)
      std::memcpy(tmp1, a, nffts*2*(nidb[nb])*sizeof(double));
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
   else if (step == 3) 
   {
      // pfftbz dev-->dev->mem,  tmp1->tmp1
      this->pfftbz(nffts, nb, tmp1, tmp2, request_indx);
   } 
   else if (step == 4) 
   {
      // pfftby mem->dev-->dev->mem
      // in=tmp1, tmp2->tmp1, tmp1=in , tmp2=tmp
      pfftby(nffts, nb, tmp1, tmp2, request_indx);
   } 
   else if (step == 5) 
   {
      // pfftbx mem->dev->dev->mem
      pfftbx(nffts, nb, tmp1, tmp2, request_indx);
      // c3db::parall->aend(request_indx);
   }
}

/********************************
 *                              *
 *      CGrid::pfftbz_start     *
 *                              *
 ********************************/
void CGrid::pfftbz_start(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{

  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) 
  {
     auto nxz = nx*nz;
    
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
           for (auto i=0; i<nx; ++i) 
           {
              if (!zero_row3[nb][indx2]) 
              {
                 auto kk = 0;
                 auto indx3 = 2*i + indx0;
                 auto shift = 2*nz*nn;
                 for (auto k=0; k<nz; ++k) 
                 {
                    tmp2[kk   + shift] = tmp1[indx3];
                    tmp2[kk+1 + shift] = tmp1[indx3+1];
                    kk += 2;
                    indx3 += nx;
                 }
                 nn += 1;
              }
              ++indx2;
           }
           indx0 += nxz;
        }
     }
    
     c3db::mygdevice.batch_cfftz_stages_tmpz(0,c3db::fft_tag,false, nz, nn, 2*nfft3d, tmp2, c3db::tmpz,da_indx);
     // for (auto i=0; i<nn; ++i)
     //    dcfftb_(&nz,tmp2+2*nz*i,c3db::tmpz);
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
     c3db::mygdevice.batch_cfftz_stages_tmpz_zero(0,c3db::fft_tag,false, nz, nq3, nffts, 2*nfft3d, tmp1, c3db::tmpz, zero_row3[nb], da_indx);
  }

}

/********************************
 *                              *
 *    CGrid::pfftbz_compute     *
 *                              *
 ********************************/
void CGrid::pfftbz_compute(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
  /**********************
   **** slab mapping ****
   **********************/
  if (maptype == 1) 
  {
     auto nxz = nx*nz;
    
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
           for (auto i=0; i<nx; ++i) 
           {
              if (!zero_row3[nb][indx2]) 
              {
                 nn += 1;
              }
              ++indx2;
           }
        }
     }
    
     c3db::mygdevice.batch_cfftz_stages_tmpz(1,c3db::fft_tag,false, nz, nn, 2*nfft3d, tmp2, c3db::tmpz,da_indx);
     // for (auto i=0; i<nn; ++i)
     //    dcfftb_(&nz,tmp2+2*nz*i,c3db::tmpz);
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
      c3db::mygdevice.batch_cfftz_stages_tmpz_zero(1,c3db::fft_tag,false, nz, nq3, nffts,2*nfft3d, tmp1, c3db::tmpz,zero_row3[nb],da_indx);
  }

}

/********************************
 *                              *
 *      CGrid::pfftbz_end       *
 *                              *
 ********************************/
void CGrid::pfftbz_end(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      auto nxz = nx*nz;
     
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
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row3[nb][indx2]) 
               {
                  nn += 1;
               }
               ++indx2;
            }
         }
      }
     
      c3db::mygdevice.batch_cfftz_stages_tmpz(2,c3db::fft_tag,false, nz, nn, 2*nfft3d, tmp2, c3db::tmpz, da_indx);
      // for (auto i=0; i<nn; ++i)
      //    dcfftb_(&nz,tmp2+2*nz*i,c3db::tmpz);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
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
                     indx3 += nx;
                  }
                  nn += 1;
               }
               ++indx2;
            }
            indx0 += nxz;
         }
      }
     
      /***********************************************
       ***         Do a pctranspose of A            ***
       ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
       ************************************************/
      c3db::c_pctranspose1_jk_start(nffts, nb, tmp1, tmp2, tmp1, request_indx, 44);
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
     c3db::mygdevice.batch_cfftz_stages_tmpz_zero(2,c3db::fft_tag,false,nz,nq3,nffts,2*nfft3d,tmp1,c3db::tmpz,zero_row3[nb],da_indx);
 
     c3db::c_pctranspose_ijk_start(nffts,nb,2,tmp1,tmp2,tmp1,request_indx, 45);
     // c3db::c_ptranspose_ijk(nb,2,tmp1,tmp2,tmp1);
   }
}

/********************************
 *                              *
 *      CGrid::pfftby_start     *
 *                              *
 ********************************/
void CGrid::pfftby_start(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      auto nxy = nx*ny;
     
      /***********************************************
       ***         Do a ptranspose of A            ***
       ***       A(kx,ky,nz) <- A(kx,nz,ky)        ***
       ************************************************/
      c3db::c_pctranspose1_jk_end(nffts, nb, tmp2, tmp1, request_indx);
     
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
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  auto jj = 0;
                  auto indx3 = 2*i + indx0;
                  auto shift = 2*ny*nn;
                  for (auto j = 0; j < ny; ++j) 
                  {
                     tmp1[jj   + shift] = tmp2[indx3];
                     tmp1[jj+1 + shift] = tmp2[indx3 + 1];
                     jj += 2;
                     indx3 += nx;
                  }
                  nn += 1;
               }
               ++indx2;
            }
            indx0 += nxy;
         }
      }
     
      c3db::mygdevice.batch_cffty_stages_tmpy(0,c3db::fft_tag,false,ny,nn,2*nfft3d,tmp1,c3db::tmpy,da_indx);
      // for (auto i=0; i<nn; ++i)
      //    dcfftb_(&ny,tmp1+2*ny*i,c3db::tmpy);
   }

   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      c3db::c_pctranspose_ijk_end(nffts,nb,2,tmp2,tmp1,request_indx);
     
      /********************************************
       ***     do fft along ny dimension        ***
       ***   A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]  ***
       ********************************************/
      c3db::mygdevice.batch_cffty_stages_tmpy_zero(0,c3db::fft_tag,false,ny,nq2,nffts,2*nfft3d,tmp2,c3db::tmpy,zero_row2[nb],da_indx);
   }
}


/********************************
 *                              *
 *    CGrid::pfftby_compute     *
 *                              *
 ********************************/
void CGrid::pfftby_compute(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      auto nxy = nx*ny;
     
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
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  nn += 1;
               }
               ++indx2;
            }
         }
      }
     
      c3db::mygdevice.batch_cffty_stages_tmpy(1,c3db::fft_tag,false, ny, nn,  2*nfft3d, tmp1, c3db::tmpy,da_indx);
      // for (auto i=0; i<nn; ++i)
      //    dcfftb_(&ny,tmp1+2*ny*i,c3db::tmpy);
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
      c3db::mygdevice.batch_cffty_stages_tmpy_zero(1,c3db::fft_tag,false,ny,nq2,1,2*nfft3d,tmp2,c3db::tmpy,zero_row2[nb],da_indx);
   }
}


/********************************
 *                              *
 *        CGrid::pfftby_end     *
 *                              *
 ********************************/
void CGrid::pfftby_end(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**********************
    **** slab mapping ****
    **********************/
   if (maptype == 1) 
   {
      auto nxy = nx*ny;
     
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
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  nn += 1;
               }
               ++indx2;
            }
         }
      }
     
      c3db::mygdevice.batch_cffty_stages_tmpy(2,c3db::fft_tag,false, ny, nn, 2*nfft3d, tmp1, c3db::tmpy,da_indx);
      // for (auto i=0; i<nn; ++i)
      //    dcfftb_(&ny,tmp1+2*ny*i,c3db::tmpy);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  auto jj = 0;
                  auto indx3 = 2*i + indx0;
                  auto shift = 2*ny*nn;
                  for (auto j=0; j<ny; ++j) 
                  {
                     tmp2[indx3]   = tmp1[jj   + shift];
                     tmp2[indx3+1] = tmp1[jj+1 + shift];
                     jj += 2;
                     indx3 += nx;
                  }
                  nn += 1;
               }
               ++indx2;
            }
            indx0 += nxy;
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
      c3db::mygdevice.batch_cffty_stages_tmpy_zero(2,c3db::fft_tag,false,ny,nq2,nffts,2*nfft3d,tmp2,c3db::tmpy,zero_row2[nb],da_indx);
 
      c3db::c_pctranspose_ijk_start(nffts, nb, 3, tmp2, tmp1, tmp2, request_indx, 46);
      // c3db::c_ptranspose_ijk(nb,3,tmp2,tmp1,tmp2);
   }
}




/********************************
 *                              *
 *      CGrid::pfftbx_start     *
 *                              *
 ********************************/
void CGrid::pfftbx_start(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
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
      c3db::mygdevice.batch_cfftx_stages_tmpx(0,c3db::fft_tag,false,nx,nffts*ny*nq,2*nfft3d,tmp2,c3db::tmpx,da_indx);
   }
   /*************************
    **** hilbert mapping ****
    *************************/
   else 
   {
      c3db::c_pctranspose_ijk_end(nffts, nb, 3, tmp1, tmp2, request_indx);
     
      /************************************************
       ***     do fft along kx dimension            ***
       ***   A(nx,ny,nz) <- fft1d^(-1)[A(kx,ny,nz)] ***
       ************************************************/
      c3db::mygdevice.batch_cfftx_stages_tmpx(0,c3db::fft_tag,false,nx,nffts*nq1,2*nfft3d,tmp1,c3db::tmpx,da_indx);
   }
}


/********************************
 *                              *
 *    CGrid::pfftbx_compute     *
 *                              *
 ********************************/
void CGrid::pfftbx_compute(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
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
      c3db::mygdevice.batch_cfftx_stages_tmpx(1,c3db::fft_tag,false,nx,nffts*ny*nq,2*nfft3d,tmp2,c3db::tmpx,da_indx);
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
      c3db::mygdevice.batch_cfftx_stages_tmpx(1,c3db::fft_tag,false,nx,nffts*nq1,2*nfft3d,tmp1,c3db::tmpx,da_indx);
   }
}



/********************************
 *                              *
 *      CGrid::pfftbx_end       *
 *                              *
 ********************************/
void CGrid::pfftbx_end(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
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
      c3db::mygdevice.batch_cfftx_stages_tmpx(2,c3db::fft_tag,false, nx, nffts*ny*nq,2*nfft3d,tmp2,c3db::tmpx,da_indx);
      std::memcpy(tmp1,tmp2,nffts*2*nfft3d*sizeof(double));
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
      c3db::mygdevice.batch_cfftx_stages_tmpx(2,c3db::fft_tag,false,nx,nffts*nq1,2*nfft3d,tmp1,c3db::tmpx,da_indx);
      if (2*nfft3d_map < 2*nfft3d)
         std::memset(tmp1 + 2*nfft3d_map, 0, (2*nfft3d - 2*nfft3d_map) * sizeof(double));
   }
}




/********************************
 *                              *
 *       CGrid::pfftb_step12    *
 *                              *
 ********************************/
void CGrid::pfftb_step12(const int step, const int nffts, const int nb, double *a, double *tmp1,
                         double *tmp2, const int request_indx, const int indx)
{
   if (step == 0) 
   {
      // c3db::parall->astart(request_indx,parall->np_i());
      // unpack start, tmp1-->tmp1
      for (auto s=0; s<nffts; ++s)
         std::memcpy(tmp1 + s*n2ft3d, a + s*2*nidb[nb], 2*(nidb[nb])*sizeof(double));
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
      // c3db::parall->aend(request_indx);
   }
}

/********************************
 *                              *
 *    CGrid:cr_pfft3b_queuein   *
 *                              *
 ********************************/
void CGrid::cr_pfft3b_queuein(const int nb, const int nffts_in, double *a) 
{
   int shift1, shift2;
   int np = c3db::parall->np_i();
 
   for (auto q=0; q<aqsize; ++q) 
   {
      int indx = aqindx[q];
      int status = aqstatus[indx] + 1;
      int nffts  = aqnffts[indx];
      shift1 = nffts_max*2*nfft3d*(2*indx);
      shift2 = nffts_max*2*nfft3d*(2*indx + 1);
      if (staged_gpu_fft_pipeline)
         pfftb_step12(status, nffts, aqnbb[indx], a, atmp+shift1, atmp+shift2, indx+4,indx);
      else
         pfftb_step(status, nffts, aqnbb[indx], a, atmp+shift1, atmp+shift2, indx+4);
      ++aqstatus[indx];
   }
 
   ++alast_index;
   if (alast_index >= aqmax)
      alast_index = 0;
   ++aqsize;
   aqindx[aqsize - 1] = alast_index;
   aqstatus[alast_index] = 0;
   aqnffts[alast_index]  = nffts_in;
   aqnbb[alast_index] = nb;
 
   // status = 0;
   shift1 = nffts_max*2*nfft3d*(2*alast_index);
   shift2 = nffts_max*2*nfft3d*(2*alast_index+1);
 
   if (staged_gpu_fft_pipeline)
      pfftb_step12(0,nffts_in, nb,a,atmp+shift1,atmp+shift2, alast_index+4,alast_index);
   else
      pfftb_step(0, nffts_in, nb, a, atmp+shift1, atmp+shift2, alast_index+4);
}

/********************************
 *                              *
 *    CGrid:cr_pfft3b_queueout  *
 *                              *
 ********************************/
void CGrid::cr_pfft3b_queueout(const int nb, const int nffts_out, double *a) 
{
   int shift1, shift2;
   int indx1 = aqindx[0];
 
   while (aqstatus[indx1] < aqmax) 
   {
      for (auto q=0; q<aqsize; ++q) 
      {
         int indx = aqindx[q];
         int status = aqstatus[indx] + 1;
         int nffts  = aqnffts[indx];
         shift1 = nffts_max*2*nfft3d * (2*indx);
         shift2 = nffts_max*2*nfft3d * (2*indx+1);
         if (staged_gpu_fft_pipeline)
            pfftb_step12(status,nffts,aqnbb[indx],a,atmp+shift1,atmp+shift2,indx+4,indx);
         else
            pfftb_step(status,nffts,aqnbb[indx],a,atmp+shift1,atmp+shift2,indx+4);
         ++aqstatus[indx];
      }
   }
   double scal1 = 1.0 / ((double)((nx) * (ny) * (nz)));
   double enrr0 = scal1 * c3db::rr_dot(atmp, atmp);
 
   shift1 = nffts_max*2*nfft3d * (2 * indx1);
   std::memcpy(a, atmp+shift1, nffts_out*2*nfft3d*sizeof(double));
   --aqsize;
   for (auto q = 0; q < aqsize; ++q)
     aqindx[q] = aqindx[q+1];
}

/********************************
 *                              *
 * CGrid:cr_pfft3b_queuefilled  *
 *                              *
 ********************************/
int CGrid::cr_pfft3b_queuefilled() { return (aqsize >= aqmax); }

/********************************
 *                              *
 *        CGrid::pfftfx         *
 *                              *
 ********************************/
void CGrid::pfftfx(const int nffts, const int nb, double *a, double *tmp1, double *tmp2, int request_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      // do fft along nx dimension
      c3db::mygdevice.batch_cfftx_tmpx(c3db::fft_tag,true, nx, ny*nq, 2*nfft3d, a, c3db::tmpx);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,true,nx,nffts*ny*nq,nx,a,c3db::forward_x,c3db::tmpx,0);
      std::memcpy(tmp1, a, nffts*2*nfft3d * sizeof(double));
   }
   /**** hilbert mapping ****/
   else 
   {
      // do fft along nx dimension
      // A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]
      c3db::mygdevice.batch_cfftx_tmpx(c3db::fft_tag,true, nx, nq1, 2*nfft3d, a, c3db::tmpx);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,true,nx,nffts*nq1,nx,a,c3db::forward_x,c3db::tmpx,0);
      c3db::c_pctranspose_ijk_start(nffts,nb,0,a,tmp1,tmp2,request_indx, 40);
   }
}

/********************************
 *                              *
 *        CGrid::pfftfy         *
 *                              *
 ********************************/
void CGrid::pfftfy(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      auto nxy = nx*ny;
     
      // do fft along ny dimension
      // A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]
      int indx0 = 0;
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
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
                     indx3 += nx;
                  }
                  ++nn;
               }
               ++indx2;
            }
            indx0 += nxy;
         }
      }
     
      c3db::mygdevice.batch_cffty_tmpy(c3db::fft_tag,true, ny, nn, 2*nfft3d, tmp2, c3db::tmpy);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,true,ny,nn,ny,tmp2,c3db::forward_y,c3db::tmpy,1);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  auto jj = 0;
                  auto indx3 = 2*i + indx0;
                  auto shift = 2*ny*nn;
                  for (auto j=0; j<ny; ++j) 
                  {
                     tmp1[indx3] = tmp2[jj + shift];
                     tmp1[indx3 + 1] = tmp2[jj + 1 + shift];
                     jj += 2;
                     indx3 += nx;
                  }
                  ++nn;
               }
               ++indx2;
            }
            indx0 += nxy;
         }
      }
     
      // Do a transpose of A
      // A(ky,nz,ky) <- A(kx,ky,nz)
      c3db::c_pctranspose2_jk_start(nffts, nb, tmp1, tmp2, tmp1, request_indx, 41);
   }
 
   /**** hilbert mapping ****/
   else 
   {
      // in=tmp1, out=tmp2
      c3db::c_pctranspose_ijk_end(nffts, nb, 0, tmp1, tmp2, request_indx);
     
      // do fft along ny dimension
      // A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]
      c3db::mygdevice.batch_cffty_tmpy_zero(c3db::fft_tag,true,ny,nq2,1,2*nfft3d,tmp1,c3db::tmpy,zero_row2[nb]);
      //c3db::mygdevice.batch_cfft_zero(c3db::fft_tag,true,ny,nq2,nffts,ny,tmp1,c3db::forward_y,c3db::tmpy,zero_row2[nb],1);
     
      // in=tmp2, out=tmp2
      c3db::c_pctranspose_ijk_start(nffts,nb,1,tmp1,tmp2,tmp1,request_indx, 42);
   }
}

/********************************
 *                              *
 *        CGrid::pfftfz         *
 *                              *
 ********************************/
void CGrid::pfftfz(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx) 
{

   /**** slab mapping ****/
   if (maptype == 1) 
   {
      c3db::c_pctranspose2_jk_end(nffts, nb, tmp2, tmp1, request_indx);
     
      auto nxz = nx*nz;
     
      // do fft along nz dimension
      // A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]
      int indx0 = 0;
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row3[nb][indx2]) 
               {
                  auto kk = 0;
                  auto indx3 = 2*i + indx0;
                  auto shift = 2*nz*nn;
                  for (auto k=0; k<nz; ++k) 
                  {
                     tmp1[kk   + shift] = tmp2[indx3];
                     tmp1[kk+1 + shift] = tmp2[indx3 + 1];
                     kk += 2;
                     indx3 += nx;
                  }
                  ++nn;
               }
               ++indx2;
            }
            indx0 += nxz;
         }
      }
     
      c3db::mygdevice.batch_cfftz_tmpz(c3db::fft_tag,true, nz, nn, 2*nfft3d, tmp1, c3db::tmpz);
      //c3db::mygdevice.batch_cfft(c3db::fft_tag,true,nz,nn,nz,tmp1,c3db::forward_z,c3db::tmpz,2);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row3[nb][indx2]) 
               {
                  auto kk = 0;
                  auto indx3 = 2 * i + indx0;
                  auto shift = 2 * nz * nn;
                  for (auto k=0; k<nz; ++k) 
                  {
                     tmp2[indx3]   = tmp1[kk   + shift];
                     tmp2[indx3+1] = tmp1[kk+1 + shift];
                     kk += 2;
                     indx3 += nx;
                  }
                  ++nn;
               }
               ++indx2;
            }
            indx0 += nxz;
         }
      }

   }
 
   /**** hilbert mapping ****/
   else 
   {
      // in=tmp1, out=tmp2
      c3db::c_pctranspose_ijk_end(nffts, nb, 1, tmp2, tmp1, request_indx);
     
      // do fft along nz dimension
      // A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]
      c3db::mygdevice.batch_cfftz_tmpz_zero(c3db::fft_tag,true, nz, nq3, 1,2*nfft3d, tmp2, c3db::tmpz,zero_row3[nb]);
      //c3db::mygdevice.batch_cfft_zero(c3db::fft_tag,true,nz,nq3,nffts,2*nfft3d,tmp2,c3db::forward_z,c3db::tmpz,zero_row3[nb],2);
   }
}

/********************************
 *                              *
 *       CGrid::pfftf_step      *
 *                              *
 ********************************/
void CGrid::pfftf_step(const int step, const int nffts, const int nb, double *a, double *tmp1, double *tmp2, int request_indx)
{
   if (step==0)
   {
      // pfftfx mem-->device, in=a out=tmp2
      pfftfx(nffts, nb, a, tmp1, tmp2, request_indx);
   }
   else if (step==1)
   {
      // pfftfy device, in=tmp1
      pfftfy(nffts, nb, tmp1, tmp2, request_indx);
   }
   else if (step==2)
   {
      // pfftfz device-->mem
      pfftfz(nffts, nb, tmp1, tmp2, request_indx);
      this->c_pack_start(nffts, nb, tmp2, tmp1, request_indx, 47);
   }
   else if (step==3)
   {
      // pfftf final
      this->c_pack_end(nffts, nb, tmp2, request_indx);
   }
}

/********************************
 *                              *
 *    CGrid::pfftfx_start       *
 *                              *
 ********************************/
void CGrid::pfftfx_start(const int nffts, const int nb, double *a, double *tmp1, double *tmp2, int request_indx, int da_indx)
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      // do fft along nx dimension
      std::memcpy(tmp2,a,nffts*2*nfft3d*sizeof(double));
      c3db::mygdevice.batch_cfftx_stages_tmpx(0,c3db::fft_tag,true,nx,nffts*ny*nq,2*nfft3d,tmp2,c3db::tmpx,da_indx);
   }
   /**** hilbert mapping ****/
   else 
   {
      // do fft along nx dimension
      // A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]
      std::memcpy(tmp2,a,nffts*2*nfft3d*sizeof(double));
      c3db::mygdevice.batch_cfftx_stages_tmpx(0,c3db::fft_tag,true,nx,nffts*nq1,2*nfft3d,tmp2,c3db::tmpx,da_indx);
   }
}

/********************************
 *                              *
 *    CGrid::pfftfx_compute     *
 *                              *
 ********************************/
void CGrid::pfftfx_compute(const int nffts, const int nb, double *a, double *tmp1, double *tmp2, int request_indx, int da_indx)
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      // do fft along nx dimension
      c3db::mygdevice.batch_cfftx_stages_tmpx(1,c3db::fft_tag,true,nx,nffts*ny*nq,2*nfft3d,tmp2,c3db::tmpx,da_indx);
   }

   /**** hilbert mapping ****/
   else 
   {
      // do fft along nx dimension
      // A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]
      c3db::mygdevice.batch_cfftx_stages_tmpx(1,c3db::fft_tag,true,nx,nffts*nq1,2*nfft3d,tmp2,c3db::tmpx,da_indx);
   }
}


/********************************
 *                              *
 *      CGrid::pfftfx_end       *
 *                              *
 ********************************/
void CGrid::pfftfx_end(const int nffts, const int nb, double *a, double *tmp1, double *tmp2, int request_indx, int da_indx)
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      // do fft along nx dimension
      c3db::mygdevice.batch_cfftx_stages_tmpx(2,c3db::fft_tag,true,nx,nffts*ny*nq,2*nfft3d,tmp2,c3db::tmpx,da_indx);
      std::memcpy(tmp1, tmp2, nffts*2*nfft3d*sizeof(double));
   }
   /**** hilbert mapping ****/
   else 
   {
      // do fft along nx dimension
      // A(kx,ny,nz) <- fft1d[A(nx,ny,nz)]
      c3db::mygdevice.batch_cfftx_stages_tmpx(2,c3db::fft_tag,true,nx,nffts*nq1,2*nfft3d,tmp2,c3db::tmpx,da_indx);
 
      c3db::c_pctranspose_ijk_start(nffts,nb,0,tmp2,tmp1,tmp2,request_indx, 40);
 
   }
}

/********************************
 *                              *
 *      CGrid::pfftfy_start     *
 *                              *
 ********************************/
void CGrid::pfftfy_start(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      auto nxy = nx*ny;
     
      // do fft along ny dimension
      // A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]
      int indx0 = 0;
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
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
                     indx3 += nx;
                  }
                  ++nn;
               }
               ++indx2;
            }
            indx0 += nxy;
         }
      }
     
      c3db::mygdevice.batch_cffty_stages_tmpy(0,c3db::fft_tag,true, ny, nn, 2*nfft3d, tmp2, c3db::tmpy,da_indx);
   }
 
   /**** hilbert mapping ****/
   else 
   {
      // in=tmp1, out=tmp2
      c3db::c_pctranspose_ijk_end(nffts,nb,0,tmp1,tmp2,request_indx);
     
      // do fft along ny dimension
      // A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]
      c3db::mygdevice.batch_cffty_stages_tmpy_zero(0,c3db::fft_tag,true,ny,nq2,nffts,2*nfft3d,tmp1,c3db::tmpy,zero_row2[nb],da_indx);
   }
}

/********************************
 *                              *
 *      CGrid::pfftfy_compute   *
 *                              *
 ********************************/
void CGrid::pfftfy_compute(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      auto nxy = nx*ny;
     
      // do fft along ny dimension
      // A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]
      int indx0 = 0;
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  ++nn;
               }
               ++indx2;
            }
         }
      }
     
      c3db::mygdevice.batch_cffty_stages_tmpy(1,c3db::fft_tag,true, ny, nn, 2*nfft3d, tmp2, c3db::tmpy,da_indx);
   }
 
   /**** hilbert mapping ****/
   else 
   {
      // in=tmp1, out=tmp2
      // do fft along ny dimension
      // A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]
      c3db::mygdevice.batch_cffty_stages_tmpy_zero(1,c3db::fft_tag,true,ny,nq2,nffts,2*nfft3d,tmp1,c3db::tmpy,zero_row2[nb],da_indx);
   }
}

/********************************
 *                              *
 *      CGrid::pfftfy_end       *
 *                              *
 ********************************/
void CGrid::pfftfy_end(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      auto nxy = nx*ny;
     
      // do fft along ny dimension
      // A(kx,ky,nz) <- fft1d[A(kx,ny,nz)]
      int indx0 = 0;
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  ++nn;
               }
               ++indx2;
            }
         }
      }
     
      c3db::mygdevice.batch_cffty_stages_tmpy(2,c3db::fft_tag,true, ny, nn, 2*nfft3d, tmp2, c3db::tmpy,da_indx);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row2[nb][indx2]) 
               {
                  auto jj = 0;
                  auto indx3 = 2*i + indx0;
                  auto shift = 2*ny*nn;
                  for (auto j=0; j<ny; ++j) 
                  {
                    tmp1[indx3]   = tmp2[jj   + shift];
                    tmp1[indx3+1] = tmp2[jj+1 + shift];
                    jj += 2;
                    indx3 += nx;
                  }
                  ++nn;
               }
               ++indx2;
            }
            indx0 += nxy;
         }
      }
     
      // Do a transpose of A
      // A(ky,nz,ky) <- A(kx,ky,nz)
      c3db::c_pctranspose2_jk_start(nffts, nb, tmp1, tmp2, tmp1, request_indx, 41);
   }

   /**** hilbert mapping ****/
   else 
   {
      // do fft along ny dimension
      // A(ky,nz,kx) <- fft1d[A(ny,nz,kx)]
      c3db::mygdevice.batch_cffty_stages_tmpy_zero(2,c3db::fft_tag,true,ny,nq2,nffts,2*nfft3d,tmp1,c3db::tmpy,zero_row2[nb],da_indx);
     
      // in=tmp2, out=tmp2
      c3db::c_pctranspose_ijk_start(nffts,nb,1,tmp1,tmp2,tmp1,request_indx, 42);
   }
}



/********************************
 *                              *
 *      CGrid::pfftfz_start     *
 *                              *
 ********************************/
void CGrid::pfftfz_start(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      c3db::c_pctranspose2_jk_end(nffts, nb, tmp2, tmp1, request_indx);
     
      auto nxz = nx*nz;
     
      // do fft along nz dimension
      // A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]
      int indx0 = 0;
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row3[nb][indx2]) 
               {
                  auto kk = 0;
                  auto indx3 = 2*i + indx0;
                  auto shift = 2*nz*nn;
                  for (auto k=0; k<nz; ++k) 
                  {
                     tmp1[kk   + shift] = tmp2[indx3];
                     tmp1[kk+1 + shift] = tmp2[indx3 + 1];
                     kk += 2;
                     indx3 += nx;
                  }
                  ++nn;
               }
               ++indx2;
            }
            indx0 += nxz;
         }
      }
     
      c3db::mygdevice.batch_cfftz_stages_tmpz(0,c3db::fft_tag,true,nz,nn,2*nfft3d,tmp1,c3db::tmpz,da_indx);
   }

   /**** hilbert mapping ****/
   else 
   {
      // in=tmp1, out=tmp2
      c3db::c_pctranspose_ijk_end(nffts, nb, 1, tmp2, tmp1, request_indx);
     
      // do fft along nz dimension
      // A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]
      c3db::mygdevice.batch_cfftz_stages_tmpz_zero(0,c3db::fft_tag,true,nz,nq3,nffts,2*nfft3d,tmp2,c3db::tmpz,zero_row3[nb],da_indx);
   }
}


/********************************
 *                              *
 *      CGrid::pfftfz_compute   *
 *                              *
 ********************************/
void CGrid::pfftfz_compute(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      auto nxz = nx*nz;
     
      // do fft along nz dimension
      // A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]
      int indx0 = 0;
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row3[nb][indx2]) 
               {
                  ++nn;
               }
               ++indx2;
            }
         }
      }

      c3db::mygdevice.batch_cfftz_stages_tmpz(1,c3db::fft_tag,true,nz,nn,2*nfft3d,tmp1,c3db::tmpz,da_indx);
   }

   /**** hilbert mapping ****/
   else 
   {
      // do fft along nz dimension
      // A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]
      c3db::mygdevice.batch_cfftz_stages_tmpz_zero(1,c3db::fft_tag,true,nz,nq3,nffts,2*nfft3d,tmp2,c3db::tmpz,zero_row3[nb],da_indx);
   }
}

/********************************
 *                              *
 *      CGrid::pfftfz_end       *
 *                              *
 ********************************/
void CGrid::pfftfz_end(const int nffts, const int nb, double *tmp1, double *tmp2, int request_indx, int da_indx) 
{
   /**** slab mapping ****/
   if (maptype == 1) 
   {
      auto nxz = nx*nz;
     
      // do fft along nz dimension
      // A(kx,kz,ky) <- fft1d[A(kx,nz,ky)]
      int indx0 = 0;
      int nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
            {
               if (!zero_row3[nb][indx2]) 
               {
                  ++nn;
               }
               ++indx2;
            }
         }
      }
     
      c3db::mygdevice.batch_cfftz_stages_tmpz(2,c3db::fft_tag,true, nz, nn, 2*nfft3d, tmp1, c3db::tmpz,da_indx);
     
      indx0 = 0;
      nn = 0;
      for (auto s=0; s<nffts; ++s)
      {
         int indx2 = 0;
         for (auto q=0; q<nq; ++q) 
         {
            for (auto i=0; i<nx; ++i) 
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
                     indx3 += nx;
                  }
                  ++nn;
               }
               ++indx2;
            }
            indx0 += nxz;
         }
      }

   }
 
   /**** hilbert mapping ****/
   else 
   {
      // do fft along nz dimension
      // A(kz,kx,ky) <- fft1d[A(nz,kx,ky)]
      c3db::mygdevice.batch_cfftz_stages_tmpz_zero(2,c3db::fft_tag,true, nz, nq3, 1,2*nfft3d, tmp2, c3db::tmpz, zero_row3[nb],da_indx);
   }
}


/********************************
 *                              *
 *       CGrid::pfftf_step10    *
 *                              *
 ********************************/
void CGrid::pfftf_step10(const int step, const int nffts, const int nb, double *a, double *tmp1,
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
 *       CGrid:c_pack_start     *
 *                              *
 ********************************/
void CGrid::c_pack_start(const int nffts, const int nb, double *a, double *tmp1,
                         const int request_indx, const int msgtype) {
  // int one=1;

  // DCOPY_PWDFT(2*nfft3d,a,one,tmp,one);
  std::memcpy(tmp1, a, nffts*2*nfft3d * sizeof(double));
  std::memset(a, 0, nffts*2*nfft3d * sizeof(double));

  for (auto s=0; s<nffts; ++s)
     c_aindexcopy(nidb2[nb], packarray[nb].data(), tmp1 + s*n2ft3d, a + s*n2ft3d);

  if (balanced)
     mybalance->c_balance_start(nffts, nb, a, request_indx, msgtype);

  return;
}

/********************************
 *                              *
 *       CGrid:c_pack_end       *
 *                              *
 ********************************/
void CGrid::c_pack_end(const int nffts, const int nb, double *tmp1, const int request_indx) {

  if (balanced)
    mybalance->c_balance_end(nffts, nb, tmp1, request_indx);

  return;
}

/********************************
 *                              *
 *    CGrid:rc_pfft3f_queuein   *
 *                              *
 ********************************/
void CGrid::rc_pfft3f_queuein(const int nb, const int nffts_in, double *b) 
{
   int shift1, shift2;
   int np = c3db::parall->np_i();
 
   for (auto q = 0; q < bqsize; ++q) 
   {
      int indx = bqindx[q];
      int status = bqstatus[indx] + 1;
      int nffts  = bqnffts[indx];
      shift1 = nffts_max*2*nfft3d * (2*indx);
      shift2 = nffts_max*2*nfft3d * (2*indx + 1);
      if (staged_gpu_fft_pipeline)
         pfftf_step10(status, nffts, bqnbb[indx], b, btmp + shift1, btmp + shift2, indx+4,indx);
      else
         pfftf_step(status, nffts, bqnbb[indx], b, btmp + shift1, btmp + shift2, indx+4);
      ++bqstatus[indx];
   }
 
   ++blast_index;
   if (blast_index >= bqmax)
      blast_index = 0;
   ++bqsize;
   bqindx[bqsize - 1] = blast_index;
   bqstatus[blast_index] = 0;
   bqnffts[blast_index]  = nffts_in;
   bqnbb[blast_index] = nb;
 
   // status = 0;
   shift1 = nffts_max*2*nfft3d * (2*blast_index);
   shift2 = nffts_max*2*nfft3d * (2*blast_index + 1);
 
   if (staged_gpu_fft_pipeline)
      pfftf_step10(0, nffts_in, nb, b, btmp + shift1, btmp + shift2, blast_index+4,blast_index);
   else
      pfftf_step(0, nffts_in, nb, b, btmp + shift1, btmp + shift2, blast_index+4);
}


/********************************
 *                              *
 *    CGrid:rc_pfft3f_queueout  *
 *                              *
 ********************************/
void CGrid::rc_pfft3f_queueout(const int nb, const int nffts_out, double *b) 
{
   int shift1, shift2;
   int indx1 = bqindx[0];
 
   while (bqstatus[indx1] < bqmax) 
   {
      for (auto q = 0; q < bqsize; ++q) 
      {
         int indx = bqindx[q];
         int status = bqstatus[indx] + 1;
         int nffts  = bqnffts[indx];
         shift1 = nffts_max*2*nfft3d * (2*indx);
         shift2 = nffts_max*2*nfft3d * (2*indx + 1);
         if (staged_gpu_fft_pipeline)
            pfftf_step10(status, nffts, bqnbb[indx], b, btmp + shift1, btmp + shift2, indx+4,indx);
         else
            pfftf_step(status, nffts, bqnbb[indx], b, btmp + shift1, btmp + shift2, indx+4);
         ++bqstatus[indx];
      }
   }
   double scal1 = 1.0 / ((double)((nx) * (ny) * (nz)));
   double enrr0 = scal1 * c3db::rr_dot(btmp, btmp);
 
   shift2 = nffts_max*2*nfft3d * (2 * indx1 + 1);
   std::memcpy(b, btmp + shift2, nffts_out*2*nfft3d * sizeof(double));
   --bqsize;
   for (auto q = 0; q < bqsize; ++q)
     bqindx[q] = bqindx[q + 1];
}


/********************************
 *                              *
 * CGrid:rc_pfft3f_queuefilled  *
 *                              *
 ********************************/
int CGrid::rc_pfft3f_queuefilled() { return (bqsize >= bqmax); }

/********************************
 *                              *
 *     CGrid:rc_pack_copy       *
 *                              *
 ********************************/
void CGrid::rc_pack_copy(const int nb, double *a, double *b) {
  int i, ii;
  int ng = nidb[nb];

  ii = 0;
  for (i = 0; i < ng; ++i) {
    b[ii] = a[i];
    b[ii + 1] = 0.0;
    ii += 2;
  }
}


/********************************
 *                              *
 *     CGrid:tc_pack_copy       *
 *                              *
 ********************************/
void CGrid::tc_pack_copy(const int nb, double *a, double *b) 
{
   int i, ii;
   int ng = nidb[nb];

   ii = 0;
   for (i = 0; i < ng; ++i) 
   {
      b[ii]   = a[i];
      b[ii+1] = 0.0;
      ii += 2;
   }
}


/********************************
 *                              *
 *      CGrid:rcc_pack_Mul      *
 *                              *
 ********************************/
void CGrid::rcc_pack_Mul(const int nb, const double *a, const double *b, double *c)
{
   int i, ii;
   int ng = nidb[nb];

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
 *      CGrid:tcc_pack_Mul      *
 *                              *
 ********************************/
void CGrid::tcc_pack_Mul(const int nb, const double *a, const double *b, double *c)
{
   int ng = nidb[nb];

   int ii = 0;
   for (auto i=0; i<ng; ++i)
   {
      c[ii]   = b[ii]*  a[i];
      c[ii+1] = b[ii+1]*a[i];
      ii += 2;
   }
}

/********************************
 *                              *
 *      CGrid:rcc_pack_aMul     *
 *                              *
 ********************************/
void CGrid::rcc_pack_aMul(const int nb, const double alpha, const double *a, const double *b, double *c)
{
   int i, ii;
   int ng = nidb[nb];

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
 *      CGrid:rc_pack_Mul       *
 *                              *
 ********************************/
void CGrid::rc_pack_Mul(const int nb, const double *a, double *c) {
  int i, ii;
  int ng = nidb[nb];

  ii = 0;
  for (i = 0; i < ng; ++i) {
    c[ii] = c[ii] * a[i];
    c[ii + 1] = c[ii + 1] * a[i];
    ii += 2;
  }
}

/********************************
 *                              *
 *    CGrid:rcc_pack_aMulAdd    *
 *                              *
 ********************************/
void CGrid::rcc_pack_aMulAdd(const int nb, const double alpha, const double *a, const double *b, double *c)
{
   int i, ii;
   int ng = nidb[nb];

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
 *      CGrid:rcc_pack_iMul     *
 *                              *
 ********************************/
void CGrid::rcc_pack_iMul(const int nb, const double *a, const double *b, double *c) {
  int i, ii;
  int ng = nidb[nb];

  ii = 0;
  for (i = 0; i < ng; ++i) 
  {
     c[ii]   = -b[ii+1]*a[i];
     c[ii+1] =  b[ii]  *a[i];
     ii += 2;
  }
}


/********************************
 *                              *
 *      CGrid:tcc_pack_iMul     *
 *                              *
 ********************************/
void CGrid::tcc_pack_iMul(const int nb, const double *a, const double *b, double *c) 
{
   int i, ii;
   int ng = nidb[nb];

   ii = 0;
   for (i=0; i<ng; ++i) 
   {
      c[ii]   = -b[ii+1]*a[i];
      c[ii+1] = b[ii]   *a[i];
      ii += 2;
   }
}


/*******************************************
 *                                         *
 *     CGrid:rcr_pack_iMul_unpack_fft      *
 *                                         *
 *******************************************/
void CGrid::rcr_pack_iMul_unpack_fft(const int nb, const double *a, const double *b, double *c)
{
   int i, ii;
   int ng = nidb[nb];

   ii = 0;
   for (i=0; i<ng; ++i)
   {
      c[ii]   = -b[ii+1]* a[i];
      c[ii+1] = b[ii]   * a[i];
      ii += 2;
   }
   this->c_unpack(nb,c);
   this->cr_pfft3b(nb,c);
}

/********************************
 *                              *
 *     CGrid:rc_pack_iMul       *
 *                              *
 ********************************/
void CGrid::rc_pack_iMul(const int nb, const double *a, double *c) 
{
   double x, y;
   int ng =  nidb[nb];
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
 *    CGrid:rcc_pack_MulSum2    *
 *                              *
 ********************************/
void CGrid::rcc_pack_MulSum2(const int nb, const double *a, const double *b, double *c) 
{
   int ng = nidb[nb];
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
 *      CGrid:cc_pack_Sum2      *
 *                              *
 ********************************/
void CGrid::cc_pack_Sum2(const int nb, const double *a, double *b) 
{
   int ng = 2*(nidb[nb]);

   for (auto i=0; i<ng; ++i)
      b[i] += a[i];
}

/********************************
 *                              *
 *     CGrid:cccc_pack_Sum      *
 *                              *
 ********************************/
void CGrid::cccc_pack_Sum(const int nb, const double *a, const double *b, const double *c, double *d) 
{
   int ng = 2*(nidb[nb]);

   for (auto i=0; i<ng; ++i)
      d[i] = (a[i] + b[i] + c[i]);
}

/********************************
 *                              *
 *     CGrid:c_pack_addzero     *
 *                              *
 ********************************/
void CGrid::c_pack_addzero(const int nb, const double vzero, double *a) {
  int pzero = cijktop(0, 0, 0);
  if (pzero == c3db::parall->taskid_i())
    a[0] += vzero;
}

/********************************
 *                              *
 *   CGrid:c_pack_noimagzero    *
 *                              *
 ********************************/

void CGrid::c_pack_noimagzero(const int nb, double *a) 
{
   int pzero = cijktop(0, 0, 0);
   if (pzero == c3db::parall->taskid_i())
      a[1] = 0.0;
}

/********************************
 *                              *
 *     CGrid:c_pack_zero        *
 *                              *
 ********************************/
void CGrid::c_pack_zero(const int nb, double *b) 
{
   int ng = 2*(nidb[nb]);

   for (auto i=0; i<ng; ++i)
      b[i] = 0.0;
}

/********************************
 *                              *
 *       CGrid:c_pack_SMul      *
 *                              *
 ********************************/
void CGrid::c_pack_SMul(const int nb, const double alpha, double *b) 
{
   int ng = 2 * (nidb[nb]);

   for (auto i = 0; i < ng; ++i)
      b[i] *= alpha;
}

/********************************
 *                              *
 *     CGrid:cc_pack_SMul       *
 *                              *
 ********************************/
void CGrid::cc_pack_SMul(const int nb, const double alpha, const double *a, double *b) 
{
   int ng = 2*(nidb[nb]);
 
   for (auto i=0; i<ng; ++i)
     b[i] = alpha*a[i];
}

/********************************
 *                              *
 *      CGrid:cc_pack_daxpy     *
 *                              *
 ********************************/
void CGrid::cc_pack_daxpy(const int nb, const double alpha, const double *a, double *b) 
{
   int ng = 2 * (nidb[nb]);

   for (auto i=0; i<ng; ++i)
      b[i] += alpha*a[i];
}

/********************************
 *                              *
 *      CGrid:cc_pack_zaxpy     *
 *                              *
 ********************************/
/**
 * @brief Performs a complex scaled vector addition (ZAXPY operation) on interleaved complex arrays.
 *
 * This function multiplies each element of an array 'a' (representing complex numbers in interleaved format)
 * by a complex scalar 'alpha' and then adds the result to an array 'b'. The operation is defined as:
 * b[i] = alpha * a[i] + b[i], where 'a' and 'b' are arrays of complex numbers represented in interleaved format.
 *
 * @param nb An integer index used to determine the number of complex elements in 'a' and 'b'.
 *           The actual number of elements processed is derived from the 'nidb' array of the CGrid class.
 * @param alpha A std::complex<double> scalar by which each element of 'a' is scaled.
 * @param a Pointer to the first element of the input array, representing complex numbers in interleaved format.
 * @param b Pointer to the first element of the output array, also in interleaved format.
 *
 * @note This function assumes that 'a' and 'b' are arrays of double representing complex numbers,
 * where each complex number is stored as consecutive real and imaginary parts.
 *
 */
void CGrid::cc_pack_zaxpy(const int nb, const std::complex<double> alpha, const double *a, double *b) 
{
   double realPart = alpha.real(); // Gets the real part of alpha
   double imagPart = alpha.imag(); // Gets the imaginary part of alpha
   int ng = (nidb[nb]);

   int i1 = 0; int i2 = 1;
   for (auto i=0; i<ng; ++i)
   {
      b[i1] += realPart*a[i1] - imagPart*a[i2];
      b[i2] += imagPart*a[i1] + realPart*a[i2];
      i1 += 2;
      i2 += 2;
   }
}


/********************************
 *                              *
 *   CGrid:ccr_pack_iconjgMul   *
 *                              *
 ********************************/
 /**
 * @brief Computes the real part of the product of the imaginary unit i and the product
 *        of a complex number A and the conjugate of another complex number B.
 *
 * This function performs a complex operation where it computes the expression:
 *
 * Re(i*A*conjg(B)) = Re(i*(a + ib)*(c - id)) = Re(i*(ac + bd) + (bc - ad)i) = bc - ad.
 *
 * Specifically, the function calculates C = Re(A) × Im(B) − Im(A) × Re(B) for each
 * pair of complex numbers in the provided arrays 'a' and 'b', storing the result in 'c'.
 *
 * @param nb An integer index used to determine the number of complex elements to process.
 *           The actual number of elements processed is derived from the 'nidb' array of the CGrid class.
 * @param a Pointer to the first element of the input array, representing complex numbers
 *          in interleaved format, where each complex number consists of a real part followed by an imaginary part.
 * @param b Pointer to the first element of the input array, representing complex numbers
 *          in interleaved format, structured similarly to 'a'.
 * @param c Pointer to the first element of the output array where the results are stored.
 *
 * @note This function assumes that 'a', 'b', and 'c' are arrays of double representing complex numbers,
 *       where each complex number is stored as consecutive real and imaginary parts. This operation is
 *       part of complex conjugate multiplication but focuses on computing one specific component of the full result.
 *
 * Usage Example:
 *   double a[] = {1.0, 2.0, 3.0, 4.0};  // Two complex numbers: 1+2i, 3+4i
 *   double b[] = {5.0, 6.0, 7.0, 8.0};  // Two complex numbers: 5+6i, 7+8i
 *   double c[2];                        // Result array
 *   int nb = 0;                         // Assuming nidb[0] is set appropriately
 *   CGrid grid;
 *   grid.cct_pack_iconjgMul(nb, a, b, c);
 *   // Now c contains the results: [12 - 16 for each complex number pair]
 */
void CGrid::ccr_pack_iconjgMul(const int nb, const double *a, const double *b, double *c)
{
   for (auto i=0; i<nidb[nb]; ++i)
      c[i] = a[2*i]*b[2*i+1] - a[2*i+1]*b[2*i];
}


/********************************
 *                              *
 *   CGrid:cct_pack_iconjgMul   *
 *                              *
 ********************************/
/**
 same as above
 */
void CGrid::cct_pack_iconjgMul(const int nb, const double *a, const double *b, double *c)
{
   for (auto i=0; i<nidb[nb]; ++i)
      c[i] = a[2*i]*b[2*i+1] - a[2*i+1]*b[2*i];
}



/********************************
 *                              *
 *  CGrid:ccr_pack_iconjgMulb   *
 *                              *
 ********************************/
 /**
 * @brief Performs a specific complex arithmetic operation by computing the expression
 * derived from the multiplication of a complex number with the conjugate of another,
 * followed by multiplying the result by the imaginary unit and extracting the negative
 * of the real part.
 *
 * The expression computed is:
 * \f[
 * -\text{Re}(i \cdot A \cdot \text{conjg}(B)) = -\text{Re}(i \cdot (a + i \cdot b) \cdot (c - i \cdot d))
 * = -\text{Re}((-b + i \cdot a) \cdot (c - i \cdot d)) = -\text{Re}((-b \cdot c + a \cdot d) + i \cdot (a \cdot c + b \cdot d))
 * = b \cdot c - a \cdot d
 * \f]
 * This computes \f$ b \cdot c - a \cdot d \f$ for each pair of complex elements from the arrays.
 *
 * @param nb An integer index used to access the number of elements to process in the 'nidb' array,
 *           which dictates the loop's iteration limit.
 * @param a Pointer to the first element of the input array 'a', where each complex number is represented
 *          by consecutive elements (real part followed by imaginary part).
 * @param b Pointer to the first element of the input array 'b', structured identically to 'a'.
 * @param c Pointer to the first element of the output array 'c', where the results are stored.
 *
 * @note The arrays 'a', 'b', and 'c' must be pre-allocated with sufficient size to hold the necessary
 *       number of elements, with 'a' and 'b' containing at least 2*nidb[nb] elements, and 'c' having
 *       at least nidb[nb] elements. The function assumes that 'nidb[nb]' is correctly set to reflect
 *       the number of complex elements to be processed. This function performs a reduction from two
 *       multiplications and one subtraction per element pair.
 */
void CGrid::ccr_pack_iconjgMulb(const int nb, const double *a, const double *b, double *c)
{
   for (auto i=0; i<(nidb[nb]); ++i)
      c[i] = a[2*i+1]*b[2*i] - a[2*i]*b[2*i+1];
}

/********************************
 *                              *
 *  CGrid:cct_pack_iconjgMulb   *
 *                              *
 ********************************/
/**
 same as above
 */
void CGrid::cct_pack_iconjgMulb(const int nb, const double *a, const double *b, double *c)
{
   for (auto i=0; i<(nidb[nb]); ++i)
      c[i] = a[2*i+1]*b[2*i] - a[2*i]*b[2*i+1];
}


/********************************
 *                              *
 *   CGrid:zccr_pack_iconjgMul  *
 *                              *
 ********************************/
/**

    C(i) = Imag (alpha*A(i)*dconjg(B(i)))

 * @brief Computes a scaled complex number multiplication with the conjugate of another complex number,
 * extracting specific components from the resulting complex product.
 *
 * This function calculates the component C(i) as:
 * \f[
 * C(i) = \text{Imag}(\alpha \cdot A(i) \cdot \text{dconjg}(B(i))) =
 * \alpha[0] \cdot (a[2i+1] \cdot b[2i] - a[2i] \cdot b[2i+1]) +
 * \alpha[1] \cdot (a[2i] \cdot b[2i] + a[2i+1] \cdot b[2i+1])
 * \f]
 * Here, \f$\alpha\f$ is a complex scalar, and \f$A(i)\f$ and \f$B(i)\f$ are the i-th complex numbers from arrays 'a' and 'b'.
 * The function performs this computation for each complex element indexed by 'i', dictated by the size stored in 'nidb[nb]'.
 *
 * @param nb An integer index used to access the number of elements to process in the 'nidb' array.
 * @param alpha Pointer to a two-element array representing the real and imaginary parts of the complex scalar alpha.
 * @param a Pointer to the first element of the input array 'a', where each complex number is represented by consecutive
 *          elements (real part followed by an imaginary part).
 * @param b Pointer to the first element of the input array 'b', structured identically to 'a'.
 * @param c Pointer to the first element of the output array where the computed results are stored.
 *
 * @note The arrays 'a', 'b', and 'c' must be pre-allocated with sufficient size to hold the necessary number of elements,
 *       with 'a' and 'b' containing at least 2*nidb[nb] elements, and 'c' having at least nidb[nb] elements.
 *       Ensure that 'nidb[nb]' is correctly set to reflect the number of complex elements to be processed.
 *       The function assumes each complex number involves two double elements (one for the real part and one for the imaginary part).
 */
void CGrid::zccr_pack_iconjgMul(const int nb, const double *alpha, const double *a, const double *b, double *c)
{
   for (auto i=0; i<nidb[nb]; ++i)
      c[i] = alpha[0]*(a[2*i+1]*b[2*i] - a[2*i]*b[2*i+1]) 
           + alpha[1]*(a[2*i]*b[2*i] + a[2*i+1]*b[2*i+1]);
}


/**********************************
 *                                *
 *    CGrid::regenerate_r_grid    *
 *                                *
 **********************************/
void CGrid::regenerate_r_grid() 
{
   int nxh = nx / 2;
   int nyh = ny / 2;
   int nzh = nz / 2;
   double a[9];
   for (auto i = 0; i < 3; ++i) 
   {
      a[i]   = lattice->unita1d(0+i) / ((double)nx);
      a[3+i] = lattice->unita1d(3+i) / ((double)ny);
      a[6+i] = lattice->unita1d(6+i) / ((double)nz);
   }
 
   r_nzero(3, r_grid);
 
   /* grid points in coordination space */
   for (auto k3 = (-nzh); k3 < nzh; ++k3)
   for (auto k2 = (-nyh); k2 < nyh; ++k2)
   for (auto k1 = (-nxh); k1 < nxh; ++k1) 
   {
      int i = k1 + nxh;
      int j = k2 + nyh;
      int k = k3 + nzh;
      int indx = cijktoindex2(i, j, k);
      int p = cijktop2(i, j, k);
 
      if (p == c3db::parall->taskid_i()) 
      {
         r_grid[3 * indx] = a[0] * k1 + a[3] * k2 + a[6] * k3;
         r_grid[3 * indx + 1] = a[1] * k1 + a[4] * k2 + a[7] * k3;
         r_grid[3 * indx + 2] = a[2] * k1 + a[5] * k2 + a[8] * k3;
      }
   }
}

/************************************
 *                                  *
 *    CGrid::generate_r_sym_grid    *
 *                                  *
 ************************************/
void CGrid::generate_r_sym_grid(double *r_sym_grid) 
{
   int nxh = nx / 2;
   int nyh = ny / 2;
   int nzh = nz / 2;
   double a[9];
   for (auto i = 0; i < 3; ++i) 
   {
      a[i] = lattice->unita1d(0 + i) / ((double)nx);
      a[3 + i] = lattice->unita1d(3 + i) / ((double)ny);
      a[6 + i] = lattice->unita1d(6 + i) / ((double)nz);
   }
 
   r_nzero(3, r_sym_grid);
 
   /* grid points in coordination space */
   for (auto k3 = (-nzh + 1); k3 < nzh; ++k3)
   for (auto k2 = (-nyh + 1); k2 < nyh; ++k2)
   for (auto k1 = (-nxh + 1); k1 < nxh; ++k1) 
   {
      int i = k1 + nxh;
      int j = k2 + nyh;
      int k = k3 + nzh;
      int indx = cijktoindex2(i, j, k);
      int p = cijktop2(i, j, k);
 
      if (p == c3db::parall->taskid_i()) 
      {
         r_sym_grid[3 * indx] = a[0] * k1 + a[3] * k2 + a[6] * k3;
         r_sym_grid[3 * indx + 1] = a[1] * k1 + a[4] * k2 + a[7] * k3;
         r_sym_grid[3 * indx + 2] = a[2] * k1 + a[5] * k2 + a[8] * k3;
      }
   }
}

/************************************
 *                                  *
 *    CGrid::generate_r_sym_mask    *
 *                                  *
 ************************************/
void CGrid::generate_r_sym_mask(double *rmask) 
{
   int nxh = nx / 2;
   int nyh = ny / 2;
   int nzh = nz / 2;
   r_zero(rmask);
 
   /* grid points in coordination space */
   for (auto k3 = (-nzh); k3 < nzh; ++k3)
   for (auto k2 = (-nyh); k2 < nyh; ++k2)
   for (auto k1 = (-nxh); k1 < nxh; ++k1) 
   {
      int i = k1 + nxh;
      int j = k2 + nyh;
      int k = k3 + nzh;
      int indx = cijktoindex2(i, j, k);
      int p = cijktop2(i, j, k);
 
      if (p == c3db::parall->taskid_i())
         rmask[indx] = 1.0;
   }
}

/************************************
 *                                  *
 *       CGrid::c_Laplacian         *
 *                                  *
 ************************************/
void CGrid::c_Laplacian(const int nb, double *w) 
{
   int npack0 = this->npack(nb);
   double *Gx = this->Gpackxyz(nb, 0);
   double *Gy = this->Gpackxyz(nb, 1);
   double *Gz = this->Gpackxyz(nb, 2);
   int kk = 0;
   for (auto k = 0; k < npack0; ++k) 
   {
      auto gg = Gx[k] * Gx[k] + Gy[k] * Gy[k] + Gz[k] * Gz[k];
      w[kk] *= (-gg);
      w[kk + 1] *= (-gg);
      kk += 2;
   }
}

/************************************
 *                                  *
 *       CGrid::cc_Laplacian        *
 *                                  *
 ************************************/
void CGrid::cc_Laplacian(const int nb, const double *w0, double *w) 
{
   int npack0 = this->npack(nb);
   double *Gx = this->Gpackxyz(nb, 0);
   double *Gy = this->Gpackxyz(nb, 1);
   double *Gz = this->Gpackxyz(nb, 2);
   int kk = 0;
   for (auto k = 0; k < npack0; ++k) 
   {
      auto gg = Gx[k] * Gx[k] + Gy[k] * Gy[k] + Gz[k] * Gz[k];
      w[kk] = -w0[kk] * gg;
      w[kk + 1] = -w0[kk + 1] * gg;
      kk += 2;
   }
}

/************************************
 *                                  *
 *      CGrid::rr_Laplacian         *
 *                                  *
 ************************************/
void CGrid::rr_Laplacian(const int nb, const double *w0, double *w) 
{
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
 *       CGrid::rr_Helmholtz        *
 *                                  *
 ************************************/
void CGrid::rr_Helmholtz(const int nb, const double *k2, const double *w, double *Hw) 
{
   this->rr_Laplacian(nb, w, Hw);
   this->rrr_Mul2Add(k2, w, Hw);
}

/************************************
 *                                  *
 *   CGrid::rrr_solve_Helmholtz     *
 *                                  *
 ************************************/
/* The routine solves the inhomogeneous Helmholtz equation

 ∇^2 w(r) + k2(r)w(r) = b(r)  using CG.

 Entry - nb: 0-density grid, 1-wvfnc grid
         k2(r): wavenumber function
         b(r): source term

 Entry/Exit - w0(r) : initial guess / output solution

*/
void CGrid::rrr_solve_Helmholtz(const int nb, const double *k2, const double *b, double *w) 
{
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
   while ((it < 10000) && (delta > 0.01)) 
   {
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
 *   CGrid::rrrr_FD_gradient        *
 *                                  *
 ************************************/
void CGrid::rrrr_FD_gradient(const double *rho, double *rhox, double *rhoy, double *rhoz) 
{
   double ua[9];
   for (auto i = 0; i < 3; ++i) 
   {
      ua[i] = lattice->unita1d(0 + i) / ((double)nx);
      ua[3 + i] = lattice->unita1d(3 + i) / ((double)ny);
      ua[6 + i] = lattice->unita1d(6 + i) / ((double)nz);
   }
   double dx = std::sqrt(ua[0] * ua[0] + ua[1] * ua[1] + ua[2] * ua[2]);
   double dy = std::sqrt(ua[3] * ua[3] + ua[4] * ua[4] + ua[5] * ua[5]);
   double dz = std::sqrt(ua[6] * ua[6] + ua[7] * ua[7] + ua[8] * ua[8]);
 
   this->rrrr_periodic_gradient(rho, rhox, rhoy, rhoz);
   for (auto i = 0; i < nfft3d; ++i) 
   {
      rhox[i] /= dx;
      rhoy[i] /= dy;
      rhoz[i] /= dz;
   }
}

/************************************
 *                                  *
 *   CGrid::rrrr_FD_laplacian       *
 *                                  *
 ************************************/
void CGrid::rrrr_FD_laplacian(const double *rho, double *rhoxx, double *rhoyy, double *rhozz) 
{
   double ua[9];
   for (auto i = 0; i < 3; ++i) 
   {
      ua[i] = lattice->unita1d(0 + i) / ((double)nx);
      ua[3 + i] = lattice->unita1d(3 + i) / ((double)ny);
      ua[6 + i] = lattice->unita1d(6 + i) / ((double)nz);
   }
   double dxx = (ua[0] * ua[0] + ua[1] * ua[1] + ua[2] * ua[2]);
   double dyy = (ua[3] * ua[3] + ua[4] * ua[4] + ua[5] * ua[5]);
   double dzz = (ua[6] * ua[6] + ua[7] * ua[7] + ua[8] * ua[8]);
 
   this->rrrr_periodic_laplacian(rho, rhoxx, rhoyy, rhozz);
   for (auto i = 0; i < nfft3d; ++i) 
   {
      rhoxx[i] /= (dxx);
      rhoyy[i] /= (dyy);
      rhozz[i] /= (dzz);
   }
}

} // namespace pwdft

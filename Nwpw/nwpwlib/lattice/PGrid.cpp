/* PGrid.cpp
   Author - Eric Bylaska

	this class is used for defining 3d parallel maps
*/



#include        <cmath>
#include 	<cstring> //memset
#include        <iostream>
#include	"Control2.hpp"
#include	"Lattice.hpp"
#include	"util.hpp"
#include	"nwpw_timing.hpp"
#include	"gdevice.hpp"

#include	"PGrid.hpp"

#include "blas.h"

namespace pwdft {


/********************************
 *                              *
 *         Constructors         *
 *                              *
 ********************************/

//PGrid::PGrid(Parallel *inparall, Lattice *inlattice, Control2& control) : d3db(inparall,control.mapping(),control.ngrid(0),control.ngrid(1),control.ngrid(2))

PGrid::PGrid(Parallel *inparall, Lattice *inlattice, int mapping0, int balance0, int nx0, int ny0, int nz0) : d3db(inparall,mapping0,nx0,ny0,nz0)
{
   int nxh,nyh,nzh,p,q,indx,nb;
   int nwave_in[2],nwave_out[2];
   double *G1, *G2,*G3;
   double ggcut,eps,ggmax,ggmin;
   bool *zero_arow3,*zero_arow2;
   bool yzslab,zrow;

   lattice = inlattice;

   eps = 1.0e-12;
   Garray = new (std::nothrow) double [3*nfft3d]();
   G1 = Garray;
   G2 = (double *) &Garray[nfft3d];
   G3 = (double *) &Garray[2*nfft3d];
   nxh = nx/2;
   nyh = ny/2;
   nzh = nz/2;
   ggmin = 9.9e9;
   ggmax = 0.0;
   for (auto k3 = (-nzh+1); k3<= nzh; ++k3)
   for (auto k2 = (-nyh+1); k2<= nyh; ++k2)
   for (auto k1 = 0;        k1<= nxh; ++k1)
   {
      auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2);
      auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2);
      auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2);
      auto gg = gx*gx + gy*gy + gz*gz;
      if (gg>ggmax) ggmax = gg;
      if ((gg<ggmin)&&(gg>1.0e-6)) ggmin = gg;
      auto i=k1; if (i < 0) i = i + nx;
      auto j=k2; if (j < 0) j = j + ny;
      auto k=k3; if (k < 0) k = k + nz;

      auto indx = ijktoindex(i,j,k);
      auto p    = ijktop(i,j,k);
      if (p==parall->taskid_i())
      {
         G1[indx] = gx;
         G2[indx] = gy;
         G3[indx] = gz;
      }

   }
   Gmax = sqrt(ggmax);
   Gmin = sqrt(ggmin);
   masker[0] = new (std::nothrow) int [2*nfft3d]();
   masker[1] = (int *) &(masker[0][nfft3d]);
   parall->Barrier();
   for (int k=0; k<(nfft3d); ++k)
   {
      masker[0][k] = 1;
      masker[1][k] = 1;
   }

   for (auto nb=0; nb<=1; ++nb)
   {
      nwave[nb] = 0;
      if (nb==0)
         ggcut = lattice->eggcut();
      else
         ggcut = lattice->wggcut();

      for (auto k3 = (-nzh+1); k3<  nzh; ++k3)
      for (auto k2 = (-nyh+1); k2<  nyh; ++k2)
      for (auto k1 = 0;        k1<  nxh; ++k1)
      {
         auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2);
         auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2);
         auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2);
         auto i=k1; if (i < 0) i = i + nx;
         auto j=k2; if (j < 0) j = j + ny;
         auto k=k3; if (k < 0) k = k + nz;

         auto indx = ijktoindex(i,j,k);
         auto p    = ijktop(i,j,k);
         if (p==parall->taskid_i())
         {
            auto gg = gx*gx + gy*gy + gz*gz;
            gg = gg-ggcut;
            if (gg < (-eps))
            {
               masker[nb][indx] = 0;
               ++nwave[nb] ;
            }
         }
      }
      nwave_entire[nb] = nwave[nb];
      nwave_entire[nb] = parall->ISumAll(1,nwave_entire[nb]);
   }


   packarray[0] = new (std::nothrow) int [2*nfft3d]();
   packarray[1] = (int *) &(packarray[0][nfft3d]);

   for (auto nb=0; nb<=1; ++nb)
   {
      nida[nb]  = 0;
      nidb2[nb] = 0;

      /* k=(0,0,0)  */
      auto k1=0;
      auto k2=0;
      auto k3=0;
      auto indx = ijktoindex(k1,k2,k3);
      auto p    = ijktop(k1,k2,k3);
      if (p==parall->taskid_i())
         if (!masker[nb][indx])
         {
            packarray[nb][nida[nb]] = indx;
            ++nida[nb];
         }

      // k=(0,0,k3) - neglect (0,0,-k3) points
      for (auto k=1; k<=(nzh-1); ++k)
      {
         k1=0;
         k2=0;
         k3=k;
         indx = ijktoindex(k1,k2,k3);
         p    = ijktop(k1,k2,k3);
         if (p==parall->taskid_i())
            if (!masker[nb][indx])
            {
               packarray[nb][nida[nb]+nidb2[nb]] = indx;
               ++nidb2[nb];
            }
      }

      // k=(0,k2,k3) - neglect (0,-k2, -k3) points 
      for (auto k=(-nzh+1); k<=(nzh-1); ++k)
      for (auto j=1;        j<=(nyh-1); ++j)
      {
         k1=0;
         k2=j;
         k3=k;
         if (k3 < 0) k3 = k3 + nz;
         indx = ijktoindex(k1,k2,k3);
         p    = ijktop(k1,k2,k3);
         if (p==parall->taskid_i())
            if (!masker[nb][indx])
            {
               packarray[nb][nida[nb]+nidb2[nb]] = indx;
               ++nidb2[nb];
            }
      }


      // k=(k1,k2,k3)
      for (auto k=(-nzh+1); k<=(nzh-1); ++k)
      for (auto j=(-nyh+1); j<=(nyh-1); ++j)
      for (auto i=1;        i<=(nxh-1); ++i)
      {
         k1=i;
         k2=j;
         k3=k;
         if (k2 < 0) k2 = k2 + ny;
         if (k3 < 0) k3 = k3 + nz;
         indx = ijktoindex(k1,k2,k3);
         p    = ijktop(k1,k2,k3);
         if (p==parall->taskid_i())
            if (!masker[nb][indx])
            {
               packarray[nb][nida[nb]+nidb2[nb]] = indx;
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
      mybalance = new Balance(parall,2,nwave_in,nwave_out);
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
   parall->Vector_ISumAll(1,2,nwave_all);


   if (maptype==1)
   {

      zero_row3[0]   = new (std::nothrow) bool[(nxh+1)*nq];
      zero_row3[1]   = new (std::nothrow) bool[(nxh+1)*nq];
      zero_row2[0]   = new (std::nothrow) bool[(nxh+1)*nq];
      zero_row2[1]   = new (std::nothrow) bool[(nxh+1)*nq];
      zero_slab23[0] = new (std::nothrow) bool[nxh+1];
      zero_slab23[1] = new (std::nothrow) bool[nxh+1];

      zero_arow3 = new bool [(nxh+1)*ny];
      for (auto nb=0; nb<=1; ++nb)
      {
         if (nb==0)
            ggcut = lattice->eggcut();
         else
            ggcut = lattice->wggcut();

         /* find zero_row3 - (i,j,*) rows that are zero */
         for (auto i=0; i<((nxh+1)*nq); ++i) zero_row3[nb][i] = true;
         for (auto i=0; i<((nxh+1)*ny); ++i) zero_arow3[i] = true;

         for (auto k2=(-nyh+1); k2<nyh; ++k2)
         for (auto k1=0;        k1<nxh; ++k1)
         {
            auto i = k1; auto j = k2;
            if (i<0) i = i + nx;
            if (j<0) j = j + ny;
            zrow = true;
            for (auto k3=(-nzh+1); k3<nzh; ++k3)
            {
               auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2);
               auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2);
               auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2);
               auto gg = gx*gx + gy*gy + gz*gz;
                    gg = gg-ggcut;
               if (gg<(-eps)) zrow = false;
            }
            if (!zrow)
            {
               zero_arow3[i+(nxh+1)*j] = false;
               q = ijktoq1(0,j,0);
               p = ijktop1(0,j,0);
               if (p==parall->taskid_i())
               {
                 zero_row3[nb][i+(nxh+1)*q] = false;
               }
             }
         }
         // call D3dB_c_ptranspose_jk_init(nb,log_mb(zero_arow3(1)))
         c_ptranspose_jk_init(nb,zero_arow3);


         /* find zero_slab23 - (i,*,*) slabs that are zero */
         for (auto i=0;  i<nxh; ++i)
            zero_slab23[nb][i] = 1;

         for (auto k1=0;  k1<nxh; ++k1)
         {
            auto i = k1;
            if (i<0) i = i + nx;
            yzslab = true;
            for (auto k3=(-nzh+1); k3<nzh; ++k3)
            for (auto k2=(-nyh+1); k2<nyh; ++k2)
            {
               auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2);
               auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2);
               auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2);
               auto gg = gx*gx + gy*gy + gz*gz;
                    gg = gg-ggcut;
               if (gg<(-eps)) yzslab = false;
            }
            if (!yzslab)
               zero_slab23[nb][i] = false;
         }

         /* find zero_row2 - (i,*,k) rows that are zero after fft of (i,j,*) */
         for (auto k=0; k<nz;      ++k)
         for (auto i=0; i<(nxh+1); ++i)
         {
            q = ijktoq(i,0,k);
            p = ijktop(i,0,k);
            if (p==parall->taskid_i())
               zero_row2[nb][q] = zero_slab23[nb][i];
         }



      }

      delete [] zero_arow3;
   }
   else
   {
      zero_row3[0] = new (std::nothrow) bool[nq3]();
      zero_row3[1] = new (std::nothrow) bool[nq3]();
      zero_row2[0] = new (std::nothrow) bool[nq2]();
      zero_row2[1] = new (std::nothrow) bool[nq2]();
      zero_slab23[0] = new (std::nothrow) bool[nxh+1]();
      zero_slab23[1] = new (std::nothrow) bool[nxh+1]();

      zero_arow2 = new (std::nothrow) bool [(nxh+1)*nz]();
      zero_arow3 = new (std::nothrow) bool [(nxh+1)*ny]();

      for (auto nb=0; nb<=1; ++nb)
      {
         if (nb==0)
            ggcut = lattice->eggcut();
         else
            ggcut = lattice->wggcut();

         // find zero_row3 - (i,j,*) rows that are zero 
         for (auto q=0; q<nq3; ++q)
            zero_row3[nb][q] = true;

         for (auto q=0; q<(nxh+1)*ny; ++q)
            zero_arow3[q] = true;


         for (auto k2=(-nyh+1); k2<nyh; ++k2)
         for (auto k1=0;        k1<nxh; ++k1)
         {
            auto i = k1; auto j = k2;
            if (i<0) i = i + nx;
            if (j<0) j = j + ny;
            zrow = true;
            for (auto k3=(-nzh+1); k3<nzh; ++k3)
            {
               auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2);
               auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2);
               auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2);
               auto gg = gx*gx + gy*gy + gz*gz;
               gg= gg-ggcut;
               if (gg<(-eps)) zrow = false;
            }
            if (!zrow)
            {
               //zero_arow3[i-1+(nxh+1)*(j-1)] = 0;
               zero_arow3[i+(nxh+1)*j] = false;
               q = ijktoq(i,j,0);
               p = ijktop(i,j,0);
               if (p==parall->taskid_i())
               {
                 zero_row3[nb][q] = false;
               }
             }
         }

         /* find zero_slab23 - (i,*,*) slabs that are zero */
         for (auto i=0;  i<nxh; ++i)
            zero_slab23[nb][i] = 1;

         for (auto k1=0;  k1<nxh; ++k1)
         {
            auto i = k1;
            if (i<0) i = i + nx;
            yzslab = true;
            for (auto k3=(-nzh+1); k3<nzh; ++k3)
            for (auto k2=(-nyh+1); k2<nyh; ++k2)
            {
               auto gx = k1*lattice->unitg(0,0) + k2*lattice->unitg(0,1) + k3*lattice->unitg(0,2);
               auto gy = k1*lattice->unitg(1,0) + k2*lattice->unitg(1,1) + k3*lattice->unitg(1,2);
               auto gz = k1*lattice->unitg(2,0) + k2*lattice->unitg(2,1) + k3*lattice->unitg(2,2);
               auto gg = gx*gx + gy*gy + gz*gz;
               gg= gg-ggcut;
               if (gg<(-eps)) yzslab = false;
            }
            if (!yzslab)
               zero_slab23[nb][i] = false;
         }

         // find zero_row2 - (i,*,k) rows that are zero after fft of (i,j,*)
         for (auto k=0; k<nz;      ++k)
         for (auto i=0; i<(nxh+1); ++i)
         {
            q = ijktoq1(i,0,k);
            p = ijktop1(i,0,k);
            if (p==parall->taskid_i())
               zero_row2[nb][q] = zero_slab23[nb][i];
         }


         //call D3dB_c_ptranspose_ijk_init(nb,log_mb(zero_arow2(1)),log_mb(zero_arow3(1)))
         c_ptranspose_ijk_init(nb,zero_arow2,zero_arow3);

      }


      delete [] zero_arow3;
      delete [] zero_arow2;

   }

   Gpack[0] = new (std::nothrow) double [3*(nida[0]+nidb[0]) + 3*(nida[1]+nidb[1])]();
   Gpack[1] = (double *) &(Gpack[0][3*(nida[0]+nidb[0])]);

   double *Gtmp = new (std::nothrow) double [nfft3d]();
   int one      = 1;
   for (auto nb=0; nb<=1; ++nb)
   {
      for (auto i=0; i<3; ++i)
      {
         DCOPY_PWDFT(nfft3d,&(Garray[i*nfft3d]),one,Gtmp,one);

         this->t_pack(nb,Gtmp);
         this->tt_pack_copy(nb,Gtmp,&(Gpack[nb][i*(nida[nb]+nidb[nb])]));
      }
   }

   delete [] Gtmp;

   zplane_tmp1 = new (std::nothrow) double[2*zplane_size+8];
   zplane_tmp2 = new (std::nothrow) double[2*zplane_size+8];


   /* initialize r_grid */
   has_r_grid = (lattice->aperiodic());
   if (has_r_grid)
   {
      r_grid = r_nalloc(3);
      r_nzero(3,r_grid);
      double a[9];
      for (auto i=0; i<3; ++i) {
         a[i]   = lattice->unita1d(0+i)/((double) nx);
         a[3+i] = lattice->unita1d(3+i)/((double) ny);
         a[6+i] = lattice->unita1d(6+i)/((double) nz);
      }

      /* grid points in coordination space */
      for (auto k3=(-nzh); k3<nzh; ++k3)
      for (auto k2=(-nyh); k2<nyh; ++k2)
      for (auto k1=(-nxh); k1<nxh; ++k1)
      {
         int i = k1 + nxh;
         int j = k2 + nyh;
         int k = k3 + nzh;
         int indx = ijktoindex2(i,j,k);
         int p    = ijktop2(i,j,k);
   
         if (p==parall->taskid_i())
         {
            r_grid[3*indx]   = a[0]*k1 + a[3]*k2 + a[6]*k3;
            r_grid[3*indx+1] = a[1]*k1 + a[4]*k2 + a[7]*k3;
            r_grid[3*indx+2] = a[2]*k1 + a[5]*k2 + a[8]*k3;
         }
      }
   }

}

PGrid::PGrid(Parallel *inparall, Lattice *inlattice, Control2& control) : PGrid(inparall,inlattice,control.mapping(),control.balance(),control.ngrid(0),control.ngrid(1),control.ngrid(2)) {}


/*
void c_indexcopy(const int n, const int *indx, double *A, double *B)
{
   int ii,jj;
   ii = 0;
   for (int i=0; i<n; ++i)
   {
      jj      = 2*indx[i];
      B[ii]   = A[jj];
      B[ii+1] = A[jj+1];
      ii +=2;
   }
}
void t_indexcopy(const int n, const int *indx, double *A, double *B)
{
   for (int i=0; i<n; ++i)
      B[i] = A[indx[i]];
}
*/

/********************************
 *                              *
 *       PGrid:c_unpack         *
 *                              *
 ********************************/
void PGrid::c_unpack(const int nb, double *a)
{
   int one=1;
   int nn = 2*(nida[nb]+nidb2[nb]);
   double *tmp1,*tmp2;
   double *tmp = new (std::nothrow) double [2*nfft3d];
   if (balanced)
      mybalance->c_unbalance(nb,a);

   DCOPY_PWDFT(nn,a,one,tmp,one);

   //dcopy_(&n2ft3d,&rzero,&zero,a,&one);
   std::memset(a, 0, n2ft3d * sizeof(double));

   c_bindexcopy(nida[nb]+nidb2[nb],packarray[nb],tmp,a);

   //tmp1 = new (std::nothrow) double[2*zplane_size+1];
   //tmp2 = new (std::nothrow) double[2*zplane_size+1];
   c_timereverse(a,zplane_tmp1,zplane_tmp2);
   //delete [] tmp2;
   //delete [] tmp1;
   delete [] tmp;
}

/********************************
 *                              *
 *       PGrid:c_pack           *
 *                              *
 ********************************/
void PGrid::c_pack(const int nb, double *a)
{
   int one=1;
   double *tmp = new (std::nothrow) double [n2ft3d];

   DCOPY_PWDFT(n2ft3d,a,one,tmp,one);

   std::memset(a, 0, n2ft3d * sizeof(double));

   c_aindexcopy(nida[nb]+nidb2[nb],packarray[nb],tmp,a);

   if (balanced)
      mybalance->c_balance(nb,a);

   delete [] tmp;
   return;
}

/********************************
 *                              *
 *       PGrid:cc_pack_copy     *
 *                              *
 ********************************/
void PGrid::cc_pack_copy(const int nb, double *a, double *b)
{
   int one = 1;
   //int ng  = 2*(nida[nb]+nidb[nb]);
   int ng  = 2*(nida[nb]+nidb[nb]);

   DCOPY_PWDFT(ng,a,one,b,one);
}

/********************************
 *                              *
 *       PGrid:cc_pack_dot      *
 *                              *
 ********************************/
double PGrid::cc_pack_dot(const int nb, double *a, double *b)
{
   int one = 1;
   //int ng  = 2*(nida[nb]+nidb[nb]);
   int ng  = 2*(nida[nb]+nidb[nb]);
   int ng0 = 2*nida[nb];
   double tsum;

   tsum = 2.0*DDOT_PWDFT(ng,a,one,b,one);
   tsum -= DDOT_PWDFT(ng0,a,one,b,one);

   return d3db::parall->SumAll(1,tsum);
}

/********************************
 *                              *
 *       PGrid:tt_pack_dot      *
 *                              *
 ********************************/
double PGrid::tt_pack_dot(const int nb, double *a, double *b)
{
   int one = 1;
   int ng  = (nida[nb]+nidb[nb]);
   int ng0 = nida[nb];
   double tsum;

   tsum = 2.0*DDOT_PWDFT(ng,a,one,b,one);
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
   //int ng  = 2*(nida[nb]+nidb[nb]);
   int ng  = 2*(nida[nb]+nidb[nb]);
   int ng0 = 2*nida[nb];
   double tsum=0.0;

   tsum = 2.0*DDOT_PWDFT(ng,a,one,b,one);
   tsum -= DDOT_PWDFT(ng0,a,one,b,one);

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
   //int ng  = 2*(nida[nb]+nidb[nb]);
   int ng  = (nida[nb]+nidb[nb]);
   int ng0 = nida[nb];
   double tsum=0.0;

   tsum = 2.0*DDOT_PWDFT(ng,a,one,b,one);
   tsum -= DDOT_PWDFT(ng0,a,one,b,one);

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
   //int ng  = 2*(nida[nb]+nidb[nb]);
   int ng  = 2*(nida[nb]+nidb[nb]);
   int ng0 = 2*nida[nb];

   for (int i=0; i<nn; ++i)
   {
      sum[i] = 2.0 * DDOT_PWDFT(ng,&(a[i*ng]),one,b,one);
      sum[i] -= DDOT_PWDFT(ng0,&(a[i*ng]),one,b,one);

   }

}


/********************************
 *                              *
 *    PGrid:cc_pack_inprjdot    *
 *                              *
 ********************************/
void PGrid::cc_pack_inprjdot(const int nb, int nn, int nprj, double *a, double *b, double *sum)
{
   int ng  = 2*(nida[nb]+nidb[nb]);
   int ng0 = 2*nida[nb];
   int one = 1;
   double rtwo  = 2.0;
   double rone  = 1.0;
   double rmone = -1.0;
   double rzero = 0.0;

   // DGEMM_PWDFT((char *) "T",(char *) "N",nn,nprj,ng,
   // 	       rtwo,
   // 	       a,ng,
   // 	       b,ng,
   // 	       rzero,
   // 	       sum,nn);
   gdevice_TN_dgemm(nn,nprj,ng,rtwo,a,b,rzero,sum);

   if (ng0>0)
   {
      DGEMM_PWDFT((char *) "T",(char *) "N",nn,nprj,ng0,
                 rmone,
                 a,ng,
                 b,ng,
                 rone,
                 sum,nn);
   }
}

/********************************
 *                              *
 *       PGrid:t_unpack         *
 *                              *
 ********************************/
void PGrid::t_unpack(const int nb, double *a)
{
   int one=1;
   int nn = (nida[nb]+nidb2[nb]);
   double *tmp1,*tmp2;
   double *tmp = new (std::nothrow) double [nfft3d];
   if (balanced)
      mybalance->t_unbalance(nb,a);

   DCOPY_PWDFT(nn,a,one,tmp,one);

   //dcopy_(&n2ft3d,&rzero,&zero,a,&one);
   std::memset(a, 0, nfft3d * sizeof(double));

   t_bindexcopy(nida[nb]+nidb2[nb],packarray[nb],tmp,a);

   //tmp1 = new (std::nothrow) double[2*zplane_size+1];
   //tmp2 = new (std::nothrow) double[2*zplane_size+1];
   t_timereverse(a,zplane_tmp1,zplane_tmp2);
   //delete [] tmp2;
   //delete [] tmp1;
   delete [] tmp;
}



/********************************
 *                              *
 *       PGrid:t_pack           *
 *                              *
 ********************************/
void PGrid::t_pack(const int nb, double *a)
{
   int one      = 1;
   int zero     = 0;
   double rzero = 0.0;
   double *tmp  = new (std::nothrow) double [nfft3d];

   DCOPY_PWDFT(nfft3d,a,one,tmp,one);

   //dcopy_(&nfft3d,&rzero,&zero,a,&one);
   std::memset(a, 0, nfft3d * sizeof(double));

   t_aindexcopy(nida[nb]+nidb2[nb],packarray[nb],tmp,a);

   if (balanced)
      mybalance->t_balance(nb,a);

   delete [] tmp;
}

/********************************
 *                              *
 *       PGrid:tt_pack_copy     *
 *                              *
 ********************************/
void PGrid::tt_pack_copy(const int nb, double *a, double *b)
{
   int one = 1;
   int ng  = nida[nb]+nidb[nb];

   DCOPY_PWDFT(ng,a,one,b,one);
}

/********************************
 *                              *
 *       PGrid:t_pack_nzero     *
 *                              *
 ********************************/
void PGrid::t_pack_nzero(const int nb, const int n, double *a)
{
   //int one   = 1;
  // int zero = 0;
  // double azero = 0.0;
   int ng  = n*(nida[nb]+nidb[nb]);
   std::memset(a, 0, ng * sizeof(double));
}



/********************************
 *                              *
 *       PGrid:i_pack           *
 *                              *
 ********************************/
void PGrid::i_pack(const int nb, int *a)
{
   int i;
   int *tmp  = new (std::nothrow) int [nfft3d];

   for (i=0; i<nfft3d; ++i) tmp[i] = a[i];
   for (i=0; i<nfft3d; ++i) a[i]    = 0;

   i_aindexcopy(nida[nb]+nidb2[nb],packarray[nb],tmp,a);

   if (balanced)
      mybalance->i_balance(nb,a);

   delete [] tmp;
}

/********************************
 *                              *
 *       PGrid:ii_pack_copy     *
 *                              *
 ********************************/
void PGrid::ii_pack_copy(const int nb, int *a, int *b)
{
   int i;
   int ng  = nida[nb]+nidb[nb];
   for (i=0; i<ng; ++i) b[i] = a[i];
}




/********************************
 *                              *
 *    PGrid:cr_pfft3b_queuein   *
 *                              *
 ********************************/
void PGrid::cr_pfft3b_queuein(const int nb, double *a)
{
   int np = parall->np_i();
}

/********************************
 *                              *
 *    PGrid:cr_pfft3b_queueout  *
 *                              *
 ********************************/
void PGrid::cr_pfft3b_queueout(const int nb, double *a)
{
}


/********************************
 *                              *
 * PGrid:cr_pfft3b_queuefilled  *
 *                              *
 ********************************/
int  PGrid::cr_pfft3b_queuefilled()
{
   //return (qsize>=qmax);
   return true;
}



/********************************
 *                              *
 *     PGrid:tc_pack_copy       *
 *                              *
 ********************************/
void PGrid::tc_pack_copy(const int nb, double *a, double *b)
{
   int i,ii;
   int ng  = nida[nb]+nidb[nb];

   ii = 0;
   for (i=0; i<ng; ++i)
   {
      b[ii]   = a[i];
      b[ii+1] = 0.0;
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
   int i,ii;
   int ng  = nida[nb]+nidb[nb];

   ii = 0;
   for (i=0; i<ng; ++i)
   {
      c[ii]   = b[ii]  *a[i];
      c[ii+1] = b[ii+1]*a[i];
      ii += 2;
   }
}



/********************************
 *                              *
 *      PGrid:tc_pack_Mul       *
 *                              *
 ********************************/
void PGrid::tc_pack_Mul(const int nb, const double *a, double *c)
{
   int i,ii;
   int ng  = nida[nb]+nidb[nb];

   ii = 0;
   for (i=0; i<ng; ++i)
   {
      c[ii]   = c[ii]  *a[i];
      c[ii+1] = c[ii+1]*a[i];
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
   int i,ii;
   int ng  = nida[nb]+nidb[nb];

   ii = 0;
   for (i=0; i<ng; ++i)
   {
      c[ii]   += alpha*b[ii]  *a[i];
      c[ii+1] += alpha*b[ii+1]*a[i];
      ii += 2;
   }
}



/********************************
 *                              *
 *      PGrid:ttc_pack_iMul     *
 *                              *
 ********************************/
void PGrid::tcc_pack_iMul(const int nb, const double *a, const double *b, double *c)
{
   int i,ii;
   int ng  = nida[nb]+nidb[nb];

   ii = 0;
   for (i=0; i<ng; ++i)
   {
      c[ii]   = -b[ii+1]*a[i];
      c[ii+1] =  b[ii]  *a[i];
      ii += 2;
   }
}

/********************************
 *                              *
 *     PGrid:tc_pack_iMul       *
 *                              *
 ********************************/
void PGrid::tc_pack_iMul(const int nb, const double *a, double *c)
{
   int i,ii;
   int ng  = nida[nb]+nidb[nb];
   double x,y;

   ii = 0;
   for (i=0; i<ng; ++i)
   {
      x = c[ii]; y = c[ii+1];

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
   int i,ii;
   int ng  = nida[nb]+nidb[nb];

   ii = 0;
   for (i=0; i<ng; ++i)
   {
      c[ii]   += b[ii]  *a[i];
      c[ii+1] += b[ii+1]*a[i];
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
   int i;
   int ng  = 2*(nida[nb]+nidb[nb]);

   for (i=0; i<ng; ++i) b[i] += a[i];
}

/********************************
 *                              *
 *     PGrid:cccc_pack_Sum      *
 *                              *
 ********************************/
void PGrid::cccc_pack_Sum(const int nb, const double *a, const double *b, const double *c, double *d)
{
   int i;
   int ng  = 2*(nida[nb]+nidb[nb]);

   for (i=0; i<ng; ++i) d[i] = (a[i]+b[i]+c[i]);
}

/********************************
 *                              *
 *     PGrid:c_pack_addzero     *
 *                              *
 ********************************/
void PGrid::c_pack_addzero(const int nb, const double vzero, double *a)
{
   int pzero = ijktop(0,0,0);
   if (pzero==parall->taskid_i()) a[0] += vzero;
}


/********************************
 *                              *
 *     PGrid:c_pack_zero        *
 *                              *
 ********************************/
void PGrid::c_pack_zero(const int nb, double *b)
{
   int i;
   int ng  = 2*(nida[nb]+nidb[nb]);

   for (i=0; i<ng; ++i) b[i] = 0.0;
}

/********************************
 *                              *
 *       PGrid:c_pack_SMul      *
 *                              *
 ********************************/
void PGrid::c_pack_SMul(const int nb, const double alpha, double *b)
{
   int i;
   int ng  = 2*(nida[nb]+nidb[nb]);

   for (i=0; i<ng; ++i) b[i] *= alpha;
}

/********************************
 *                              *
 *     PGrid:cc_pack_SMul       *
 *                              *
 ********************************/
void PGrid::cc_pack_SMul(const int nb, const double alpha, const double *a, double *b)
{
   int i;
   int ng  = 2*(nida[nb]+nidb[nb]);

   for (i=0; i<ng; ++i) b[i] = alpha*a[i];
}

/********************************
 *                              *
 *      PGrid:cc_pack_daxpy     *
 *                              *
 ********************************/
void PGrid::cc_pack_daxpy(const int nb, const double alpha, const double *a, double *b)
{
   int i;
   int ng  = 2*(nida[nb]+nidb[nb]);

   for (i=0; i<ng; ++i) b[i] += alpha*a[i];
}

/********************************
 *                              *
 *   PGrid:cct_pack_iconjgMul   *
 *                              *
 ********************************/
void PGrid::cct_pack_iconjgMul(const int nb, const double *a, const double *b, double *c)
{
   for (int i=0; i<(nida[nb]+nidb[nb]); ++i)
      c[i] = a[2*i]*b[2*i+1] - a[2*i+1]*b[2*i];
}


/********************************
 *                              *
 *  PGrid:cct_pack_iconjgMulb   *
 *                              *
 ********************************/
void PGrid::cct_pack_iconjgMulb(const int nb, const double *a, const double *b, double *c)
{
   for (int i=0; i<(nida[nb]+nidb[nb]); ++i)
      c[i] = a[2*i+1]*b[2*i] - a[2*i]*b[2*i+1];
}


/**********************************
 *                                *
 *    PGrid::regenerate_r_grid    *
 *                                *
 **********************************/
void PGrid::regenerate_r_grid() 
{
   int nxh = nx/2;
   int nyh = ny/2;
   int nzh = nz/2;
   double a[9];
   for (auto i=0; i<3; ++i) {
      a[i]   = lattice->unita1d(0+i)/((double) nx);
      a[3+i] = lattice->unita1d(3+i)/((double) ny);
      a[6+i] = lattice->unita1d(6+i)/((double) nz);
   }

   r_nzero(3,r_grid);

   /* grid points in coordination space */
   for (auto k3=(-nzh); k3<nzh; ++k3)
   for (auto k2=(-nyh); k2<nyh; ++k2)
   for (auto k1=(-nxh); k1<nxh; ++k1)
   { 
      int i = k1 + nxh;
      int j = k2 + nyh;
      int k = k3 + nzh;
      int indx = ijktoindex2(i,j,k);
      int p    = ijktop2(i,j,k);

      if (p==parall->taskid_i())
      { 
         r_grid[3*indx]   = a[0]*k1 + a[3]*k2 + a[6]*k3;
         r_grid[3*indx+1] = a[1]*k1 + a[4]*k2 + a[7]*k3;
         r_grid[3*indx+2] = a[2]*k1 + a[5]*k2 + a[8]*k3;
      }
   }
}


/************************************
 *                                  *
 *    PGrid::generate_r_sym_grid    *
 *                                  *
 ************************************/
void PGrid::generate_r_sym_grid(double *r_sym_grid) 
{
   int nxh = nx/2;
   int nyh = ny/2;
   int nzh = nz/2;
   double a[9];
   for (auto i=0; i<3; ++i) {
      a[i]   = lattice->unita1d(0+i)/((double) nx);
      a[3+i] = lattice->unita1d(3+i)/((double) ny);
      a[6+i] = lattice->unita1d(6+i)/((double) nz);
   }

   r_nzero(3,r_sym_grid);

   /* grid points in coordination space */
   for (auto k3=(-nzh+1); k3<nzh; ++k3)
   for (auto k2=(-nyh+1); k2<nyh; ++k2)
   for (auto k1=(-nxh+1); k1<nxh; ++k1)
   {
      int i = k1 + nxh;
      int j = k2 + nyh;
      int k = k3 + nzh;
      int indx = ijktoindex2(i,j,k);
      int p    = ijktop2(i,j,k);

      if (p==parall->taskid_i())
      {
         r_sym_grid[3*indx]   = a[0]*k1 + a[3]*k2 + a[6]*k3;
         r_sym_grid[3*indx+1] = a[1]*k1 + a[4]*k2 + a[7]*k3;
         r_sym_grid[3*indx+2] = a[2]*k1 + a[5]*k2 + a[8]*k3;
      }
   }
}

/************************************
 *                                  *
 *    PGrid::generate_r_sym_mask    *
 *                                  *
 ************************************/
void PGrid::generate_r_sym_mask(double *rmask) 
{
   int nxh = nx/2;
   int nyh = ny/2;
   int nzh = nz/2;
   r_zero(rmask);

   /* grid points in coordination space */
   for (auto k3=(-nzh); k3<nzh; ++k3)
   for (auto k2=(-nyh); k2<nyh; ++k2)
   for (auto k1=(-nxh); k1<nxh; ++k1)
   {
      int i = k1 + nxh;
      int j = k2 + nyh;
      int k = k3 + nzh;
      int indx = ijktoindex2(i,j,k);
      int p    = ijktop2(i,j,k);

      if (p==parall->taskid_i()) rmask[indx] = 1.0;
   }
}



}


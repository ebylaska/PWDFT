#ifndef _CGRID_HPP_
#define _CGRID_HPP_

#pragma once

/* CGrid.hpp
   Author - Eric Bylaska

*/

#pragma once

#include "CBalance.hpp"
#include "Lattice.hpp"
#include "Parallel.hpp"
#include "Brillouin.hpp"
#include "c3db.hpp"
#include "k1db.hpp"
#include <cmath>

namespace pwdft {

class CGrid : public k1db, public c3db {

   CBalance  *mybalance;
   int balanced;
 
   /* Brillouin */
   Brillouin *mybrillouin;
 
   /* G grid data */
   double *Garray, **Gpack;
   double Gmax, Gmin;
   int **masker, **packarray;
   int *nwave, *nwave_entire, *nwave_all,  *nidb, *nidb2, nidb1_max;

   double *p_kvector, *p_weight;
 
   /* pfft data */
   bool **zero_row2, **zero_row3, **zero_slab23;
 
   /* pfft_queue data */
   int aqmax, aqsize, alast_index;
   int *aqindx, *aqstatus;
   double *atmp;
 
   int bqmax, bqsize, blast_index;
   int *bqindx, *bqstatus;
   double *btmp;
 
   /* zplane data */
   double *zplane_tmp1, *zplane_tmp2;

public:
   bool staged_gpu_fft_pipeline;
 
   /* lattice pointer */
   Lattice *lattice;
 
   /* r_grid data */
   bool has_r_grid = false;
   double *r_grid;
 
   /* constructor */
   CGrid(Parallel *, Lattice *, int, int, int, int, int, int, bool, Brillouin *);
   CGrid(Parallel *, Lattice *, Control2 &, Brillouin *);
 
   /* destructor */
   ~CGrid() {
      delete [] p_weight;
      delete [] p_kvector;
      delete [] Garray;
      delete [] nwave;
      delete [] nwave_entire;
      delete [] nwave_all;
      delete [] nidb;
      delete [] nidb2;
      for (auto nb=0; nb<=nbrillq; ++nb)
      {
         delete [] Gpack[nb];
         delete [] masker[nb];
         delete [] packarray[nb];
         delete [] zero_row3[nb];
         delete [] zero_row2[nb];
         delete [] zero_slab23[nb];
      }
      
      delete [] Gpack;
      delete [] masker;
      delete [] packarray;
      delete [] zero_row3;
      delete [] zero_row2;
      delete [] zero_slab23;

      if (balanced)
        delete mybalance;

      delete [] zplane_tmp1;
      delete [] zplane_tmp2;
      if (has_r_grid)
        delete [] r_grid;
      delete [] atmp;
      delete [] aqindx;
      delete [] aqstatus;
      delete [] btmp;
      delete [] bqindx;
      delete [] bqstatus;
     
      // deallocate async buffer data
      for (auto q=0; q<aqmax; ++q)
        c3db::parall->aend(3+q);
   }
 
   double *Gxyz(const int i) { return Garray + i*nfft3d; }
   double *Gpackxyz(const int nb, const int i) { return  Gpack[nb] + i*(nidb[nb]); }
   double Gmax_ray() { return Gmax; }
   double Gmin_ray() { return Gmin; }
   double dGmin_ray() { return 0.01 * Gmin; }
   int n_ray() {
     int nray0 = (int)ceil(100 * (Gmax / Gmin) + 1.0);
     if (nray0 < 10)
       nray0 = 10;
     return nray0;
   }
   double *generate_G_ray() {
     int nray0 = (int)ceil(100 * (Gmax / Gmin) + 1.0);
     if (nray0 < 10)
       nray0 = 10;
     double *g_ray = new (std::nothrow) double[nray0]();
     double dGmin = 0.01 * Gmin;
     for (auto i = 0; i < nray0; ++i)
       g_ray[i] = dGmin * i;
     return g_ray;
   }

   double *pbrill_kvector(const int i) { return p_kvector + 3*i; }
   double pbrill_weight(const int i)   { return p_weight[i]; }
 
   int npack1_max() { return (nidb1_max); }
   int npack(const int nb) { return (nidb[nb]); }
   int npack_all(const int nb) { return nwave_all[nb]; }
   int isbalanced() { return balanced; }
 
   double *c_pack_allocate(const int nb) {
     double *ptr;
     ptr = new (std::nothrow) double[2*npack(nb)]();
     return ptr;
   }
   void c_pack_deallocate(double *ptr) { delete[] ptr; }
 
   double *r_pack_allocate(const int nb) {
     double *ptr;
     ptr = new (std::nothrow) double[npack(nb)]();
     return ptr;
   }
   void r_pack_deallocate(double *ptr) { delete[] ptr; }
 
   void c_unpack(const int, double *);
   void c_pack(const int, double *);
   void cc_pack_copy(const int, const double *, double *);
   double cc_pack_dot(const int, double *, double *);
   double cc_pack_idot(const int, double *, double *);
   void cc_pack_indot(const int, const int, double *, double *, double *);
   double rr_pack_dot(const int, double *, double *);
   double rr_pack_idot(const int, double *, double *);
 
   void cc_pack_inprjdot(const int, int, int, double *, double *, double *);
 
   void r_unpack(const int, double *);
   void r_pack(const int, double *);
   void rr_pack_copy(const int, const double *, double *);
   void r_pack_nzero(const int, const int, double *);

   double tt_pack_dot(const int, double *, double *);
   double tt_pack_idot(const int, double *, double *);
   void t_unpack(const int, double *);
   void t_pack(const int, double *);
   void tt_pack_copy(const int, const double *, double *);
   void t_pack_nzero(const int, const int, double *);

   void tcc_pack_Mul(const int, const double *, const double *, double *);
   void tcc_pack_iMul(const int, const double *, const double *, double *);
 
   void rc_pack_copy(const int, double *, double *);
   void tc_pack_copy(const int, double *, double *);
 
   void rcc_pack_Mul(const int, const double *, const double *, double *);
   void rcc_pack_aMul(const int, const double, const double *, const double *, double *);
   void rc_pack_Mul(const int, const double *, double *);
   void rcc_pack_iMul(const int, const double *, const double *, double *);
   void rc_pack_iMul(const int, const double *, double *);
 
   void rcc_pack_MulSum2(const int, const double *, const double *, double *);
   void cc_pack_Sum2(const int, const double *, double *);
   void cccc_pack_Sum(const int, const double *, const double *, const double *,
                      double *);
   void rcc_pack_aMulAdd(const int, const double, const double *, const double *,
                         double *);
 
   void c_pack_addzero(const int, const double, double *);
   void c_pack_noimagzero(const int, double *);
 
   void c_pack_zero(const int, double *);
   void c_pack_SMul(const int, const double, double *);
   void cc_pack_SMul(const int, const double, const double *, double *);
   void cc_pack_daxpy(const int, const double, const double *, double *);
   void ccr_pack_iconjgMul(const int, const double *, const double *, double *);
   void ccr_pack_iconjgMulb(const int, const double *, const double *, double *);
 
   void cct_pack_iconjgMul(const int, const double *, const double *, double *);
   void cct_pack_iconjgMulb(const int, const double *, const double *, double *);

   void i_pack(const int, int *);
   void ii_pack_copy(const int, int *, int *);
 
   void cr_pfft3b_queuein(const int, double *);
   void cr_pfft3b_queueout(const int, double *);
   int cr_pfft3b_queuefilled();
   void cr_pfft3b(const int, double *);
   void pfftb_step(const int, const int, double *, double *, double *, const int);
   void pfftb_step12(const int, const int, double *, double *, double *, const int,const int);
 
   void c_unpack_start(const int, double *, double *, const int, const int);
   void c_unpack_mid(const int, double *, double *, const int, const int);
   void c_unpack_end(const int, double *, double *, const int);
   void pfftbz(const int, double *, double *, int);
   void pfftby(const int, double *, double *, int);
   void pfftbx(const int, double *, double *, int);
 
   void rc_pfft3f_queuein(const int, double *);
   void rc_pfft3f_queueout(const int, double *);
   int rc_pfft3f_queuefilled();
   void rc_pfft3f(const int, double *);
   void pfftf_step(const int, const int, double *, double *, double *, int);
   void pfftf_step10(const int, const int, double *, double *, double *, int,int);
   void c_pack_start(const int, double *, double *, const int, const int);
   void c_pack_end(const int, double *, const int);
   void pfftfx(const int, double *, double *, double *, int);
   void pfftfy(const int, double *, double *, int);
   void pfftfz(const int, double *, double *, int);
   void pfftf_final(const int, double *, double *, int);
 
   void pfftfx_start(const int, double *, double *, double *, int,int);
   void pfftfx_compute(const int, double *, double *, double *, int,int);
   void pfftfx_end(const int, double *, double *, double *, int,int);
 
   void pfftfy_start(const int, double *, double *, int,int);
   void pfftfy_compute(const int, double *, double *, int,int);
   void pfftfy_end(const int, double *, double *, int,int);
 
   void pfftfz_start(const int, double *, double *, int,int);
   void pfftfz_compute(const int, double *, double *, int,int);
   void pfftfz_end(const int, double *, double *, int,int);
 
 
   void pfftbx_start(const int, double *, double *, int,int);
   void pfftbx_compute(const int, double *, double *, int,int);
   void pfftbx_end(const int, double *, double *, int,int);
 
   void pfftby_start(const int, double *, double *, int, int);
   void pfftby_compute(const int, double *, double *, int, int);
   void pfftby_end(const int, double *, double *, int, int);
 
   void pfftbz_start(const int, double *, double *, int, int);
   void pfftbz_compute(const int, double *, double *, int, int);
   void pfftbz_end(const int, double *, double *, int, int);
 
   void rcr_pack_iMul_unpack_fft(const int, const double *, const double *,
                                 double *);
 
   void regenerate_r_grid();
   void initialize_r_grid() {
     if (!has_r_grid) {
       has_r_grid = true;
       r_grid = r_nalloc(3);
       this->regenerate_r_grid();
     }
   }
   void generate_r_sym_grid(double *);
   void generate_r_sym_mask(double *);
 
   void c_Laplacian(const int, double *);
   void cc_Laplacian(const int, const double *, double *);
   void rr_Laplacian(const int, const double *, double *);
   void rr_Helmholtz(const int, const double *, const double *, double *);
   void rrr_solve_Helmholtz(const int, const double *, const double *, double *);
 
   void rrrr_FD_gradient(const double *, double *, double *, double *);
   void rrrr_FD_laplacian(const double *, double *, double *, double *);
 
   void r_pack_gaussian(const int nb, const double rcut, double *gauss) 
   {
      double w = rcut * rcut / 4.0;
      int npack0 = nidb[nb];
      double *gx = Gpackxyz(nb, 0);
      double *gy = Gpackxyz(nb, 1);
      double *gz = Gpackxyz(nb, 2);
      for (auto k = 0; k < npack0; ++k) {
        auto gg = gx[k] * gx[k] + gy[k] * gy[k] + gz[k] * gz[k];
        gauss[k] = std::exp(-w * gg);
      }
   }
};

} // namespace pwdft

#endif

#ifndef _C3dB_HPP_
#define _C3dB_HPP_

/* c3db.h
   Author - Eric Bylaska

  this class is container of operations for distributed 3d blocks
*/

#pragma once

#include "Mapping3c.hpp"
#include "Parallel.hpp"
#include "gdevice2.hpp"

#include <sstream>
#include <stdexcept>

namespace pwdft {


/**
 * @class c3db
 * @brief Container class for distributed 3D blocks operations.
 */
class c3db : public Mapping3c {

   /* transpose indexings */
   int **iq_to_i1, **iq_to_i2;
   int **i1_start, **i2_start;
 
   bool initialized_r_transpose = false;
   int **iq_to_ir1, **iq_to_ir2;
   int **ir1_start, **ir2_start;
 
 
   /* ptrans indexings */
   //int **p_iq_to_i1[2], **p_iq_to_i2[2], **p_iz_to_i2[2], p_iz_to_i2_count[2][6];
   //int **p_i1_start[2], **p_i2_start[2];
   //int **p_jq_to_i1[2], **p_jq_to_i2[2], **p_jz_to_i2[2];
   //int **p_j1_start[2], **p_j2_start[2];

   int p_nbrillq0=1;
   int ***p_iq_to_i1, ***p_iq_to_i2, ***p_iz_to_i2, **p_iz_to_i2_count;
   int ***p_i1_start, ***p_i2_start;

   int ***p_jq_to_i1, ***p_jq_to_i2, ***p_jz_to_i2;
   int ***p_j1_start, ***p_j2_start;


public:
   gdevice2 mygdevice;
 
   /* fft tabulations of of trigonometry functions */
   int fft_tag=-1;
   double *tmpx, *tmpy, *tmpz;
   double *forward_x,  *forward_y,  *forward_z;
   double *backward_x, *backward_y, *backward_z;
 
   Parallel *parall;
   int zplane_size;
 
   /* c3db_tmp data */
   double *c3db_tmp1,*c3db_tmp2;
 
   /* constructor */
   c3db(Parallel *, const int, const int, const int, const int, const int);
 
   /* destructor */

    /**
     * @brief Destructor for the d3db class.
     */
   ~c3db();
 
   /* c array operators */
   double *c_alloc();
   double *c_nalloc(const int);
   void c_dealloc(double *);
   void c_zero(double *);
   void c_nzero(int, double *);
   void cc_copy(const double *, double *);

   double cc_dot(const double *, const double *);

   /* r array operators */
   double *r_alloc();
   double *r_nalloc(const int);
   void r_dealloc(double *);
   void r_zero(double *);
   void r_nzero(int, double *);
   void rr_copy(const double *, double *);
   void rr_SMul(const double, const double *, double *);
   void r_SMul(const double, double *);
   void rrr_SMulAdd(const double, const double *, const double *, double *);
   double r_dsum(const double *);
   double rr_dot(const double *, const double *);
 
   void nrr_vdot(const int, const double *, const double *, double *);
   void r_zero_ends(double *);
   void r_zero_mends(double *);
   void r_abs(double *);
 
   void r_sqr(double *);
   void rr_sqr(const double *, double *);
   void rr_addsqr(const double *, double *);
   void r_sqrt(double *);
 
   void rrr_Sum2Add(const double *, const double *, double *);
   void rrrr_Sum(const double *, const double *, const double *, double *);
   void rrr_Sum(const double *, const double *, double *);
   void rr_Sum(const double *, double *);
   void rrr_Mul(const double *, const double *, double *);
   void rrr_Mul2Add(const double *, const double *, double *);
   void rr_Mul(const double *, double *);
   void arrr_Minus(const double, const double *, const double *, double *);
   void rrr_Minus(const double *, const double *, double *);
   void rr_Minus(const double *, double *);
 
   void rr_Divide(const double *, double *);
   void rrr_Divide(const double *, const double *, double *);
   void rr_screen0(const double *, double *);
 
   void rr_daxpy(const double, const double *, double *);
 
   void rrrrr_SumMulAdd(const double *, const double *, const double *,
                        const double *, double *);

   void rc_Mul(const double *, double *);
   void rrc_Sum(const double *, const double *, double *);
   void c_SMul(const double, double *);
   void cc_SMul(const double, const double *, double *);
   void cc_Mul(const double *, double *);

   void rcc_Sum(const double *, const double *, double *);
   void bb_Mul(const double *, double *);

   void c_ZMul(const double, const double, double *);


   void t_read(const int, double *, const int, const int);
   void t_write(const int, double *, const int, const int);
   void t_write_buffer(const int, double *, const int, const int);

   //void t_write_buffer_max(const int, double *, const int, const int, const int, int &, double *);
   //void t_write_buffer_max_final(const int, int &, double *);
 
   // void  	 r_read(const int, const int, double *);
   void c_read(const int, double *, const int, const int);
   void c_write(const int, double *, const int, const int);
   void c_write_buffer(const int, double *, const int, const int);
   void c_write_buffer_max(const int, double *, const int, const int, const int, int &, double *);
   void c_write_buffer_max_final(const int, int &, double *);

   void r_read(const int, double *, const int, const int, const bool);
   void r_write(const int, double *, const int, const int, const bool);
 
   void cr_fft3d(double *);
   void rc_fft3d(double *);
   void cshift1_fftb(const int, const int, const int, const int, double *);
   void cshift_fftf(const int, const int, const int, const int, double *);
   void zeroend_fftb(const int, const int, const int, const int, double *);
 
   void c_transpose_jk(double *, double *, double *);
   void c_transpose_ijk(const int, double *, double *, double *);
 
 
   void c_setpw(const int *, const double *, double *);
   void c_addrandom(double *);
   void r_setrandom(double *);


   /* expand and contract operators */
   void hr2r_expand(const double *, double *);
   void r2hr_contract(const double *, double *);

 
 
   /* ptranspose operators */
   void c_ptranspose_jk_init(const int, bool *);
   void c_ptranspose_ijk_init(const int, bool *, bool *);
   void c_ptranspose1_jk(const int, double *, double *, double *);
   void c_ptranspose2_jk(const int, double *, double *, double *);
   void c_ptranspose_ijk(const int, const int, double *, double *, double *);
 
   void c_ptranspose1_jk_start(const int, double *, double *, double *,
                               const int, const int);
   void c_ptranspose2_jk_start(const int, double *, double *, double *,
                               const int, const int);
   void c_ptranspose_ijk_start(const int, const int, double *, double *,
                               double *, const int, const int);
   void c_ptranspose1_jk_end(const int, double *, double *, const int);
   void c_ptranspose2_jk_end(const int, double *, double *, const int);
   void c_ptranspose_ijk_end(const int, const int, double *, double *,
                             const int);
 
   /* gcube io */
   std::string r_formatwrite_reverse(double *);
   std::string r_formatwrite(double *);
 
   /* real-space transposes and gradients */
   void r_transpose_ijk_init();
   void r_transpose_ijk_end();
   void r_transpose_jk(double *, double *, double *);
   void r_transpose_ijk(const int, double *, double *, double *);
   void rrrr_periodic_gradient(const double *, double *, double *, double *);
   void rrrr_periodic_laplacian(const double *, double *, double *, double *);
 
   /* real-space special gradients */
   void rrr_SqrMulAdd(const double *, const double *, double *);
   void rrrrrrr_Sqr3MulPlusMul2(const double *, const double *, const double *,
                                const double *, const double *, const double *,
                                double *);
 
   /* real-space Gaussian filters */
   void rr_periodic_gaussian_filter(const double, const double *, double *);

};
} // namespace pwdft

#endif

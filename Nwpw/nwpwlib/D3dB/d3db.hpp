#ifndef _D3dB_H_
#define _D3dB_H_
/* d3db.h
   Author - Eric Bylaska

  this class is container of operations for distributed 3d blocks
*/

#pragma once

#include	"Parallel.hpp"
#include	"Mapping3.hpp"

#ifdef NWPW_SYCL
#include "gdevice.hpp"
#include "oneapi/mkl.hpp"

typedef oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::SINGLE,
				     oneapi::mkl::dft::domain::REAL> desc_real_t;
typedef oneapi::mkl::dft::descriptor<oneapi::mkl::dft::precision::SINGLE,
				     oneapi::mkl::dft::domain::COMPLEX> desc_cmplx_t;
#endif

class d3db : public Mapping3 {

   /* transpose indexings */
   int **iq_to_i1,**iq_to_i2;
   int **i1_start,**i2_start;

   /* timereverse indexings */
   int *t_iq_to_i1,*t_iq_to_i2;
   int *t_i1_start,*t_i2_start;

   /* ptrans indexings */
   int **p_iq_to_i1[2],**p_iq_to_i2[2],**p_iz_to_i2[2], p_iz_to_i2_count[2][6];
   int **p_i1_start[2],**p_i2_start[2];

   int **p_jq_to_i1[2],**p_jq_to_i2[2],**p_jz_to_i2[2];
   int **p_j1_start[2],**p_j2_start[2];

   float *tmpx,*tmpy,*tmpz;

#ifdef NWPW_SYCL
  float *a_dev;
  desc_real_t *desc_x;
  desc_cmplx_t *desc_y, *desc_z;
#endif

public:
        Parallel  *parall;
        int zplane_size;

        /* constructor */
	d3db(Parallel *, const int,const int,const int,const int);

        /* destructor */
        ~d3db();

        /* r array operators */
        float * r_alloc();
        float * r_nalloc(const int);
        void     r_dealloc(float *);
        void     r_zero(float *);
        void     r_nzero(int, float *);
        void     rr_copy(const float *,float *);
        void     tt_copy(const float *,float *);
        void     rr_SMul(const float, const float *, float *);
        void	 r_SMul(const float, float *);
        void     rrr_SMulAdd(const float, const float *, const float *, float *);
        float	 r_dsum(const float *);
        float	 rr_dot(const float *, const float *);
        void  	 r_zero_ends(float *);
        void  	 r_abs(float *);
        void     r_sqr(float *);
        void  	 rrr_Sum(const float *, const float *, float *);
        void  	 rr_Sum(const float *, float *);
        void  	 rrr_Mul(const float *, const float *, float *);
        void  	 rrr_Minus(const float *, const float *, float *);
        void  	 rrr_Divide(const float *, const float *, float *);
        void  	 rr_daxpy(const float, const float *, float *);
        //void  	 r_read(const int, const int, float *);
	void     c_read(const int, float *, const int);
	void     c_write(const int, float *, const int);

	void     cr_fft3d(float *);
	void     rc_fft3d(float *);

	void     c_transpose_jk(float *, float *, float *);
	void     t_transpose_jk(float *, float *, float *);

#ifdef NWPW_SYCL
  void     c_transpose_ijk_sycl(const int, const int*, const int*, float *, float *, float *);
#endif
	void     c_transpose_ijk(const int, float *, float *, float *);

	void     t_transpose_ijk(const int, float *, float *, float *);

	void     t_timereverse(float *, float *, float *);
	void     c_timereverse(float *, float *, float *);
        int	 timereverse_size();

        void     c_setpw(const int *, const float *, float *);
        void     c_addrandom(float *);



        /* t array operators */
        float * t_alloc();
        void     t_dealloc(float *);
	void     t_read(const int, float *, const int);
	void     t_write(const int, float *, const int);
        void     t_nzero(int, float *);

        /* ptranspose operators */
 	void     c_ptranspose_jk_init(const int, int  *);
	void     c_ptranspose_ijk_init(const int, int *, int *);
};

#endif

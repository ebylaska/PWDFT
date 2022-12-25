/* d3db.h
   Author - Eric Bylaska

  this class is container of operations for distributed 3d blocks
*/

#pragma once

#include	"Parallel.hpp"
#include	"Mapping3.hpp"

#include <stdexcept>
#include <sstream>

namespace pwdft {


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

public:

    /* fft tabulations of of trigonometry functions */
    double *tmpx,*tmpy,*tmpz;

    Parallel  *parall;
    int zplane_size;

    /* constructor */
    d3db(Parallel *, const int,const int,const int,const int);

    /* destructor */
    ~d3db();

    /* r array operators */
    double * r_alloc();
    double * r_nalloc(const int);
    void     r_dealloc(double *);
    void     r_zero(double *);
    void     r_nzero(int, double *);
    void     rr_copy(const double *,double *);
    void     tt_copy(const double *,double *);
    void     rr_SMul(const double, const double *, double *);
    void      r_SMul(const double, double *);
    void     rrr_SMulAdd(const double, const double *, const double *, double *);
    double   r_dsum(const double *);
    double  rr_dot(const double *, const double *);
    void     nrr_vdot(const int, const double *, const double *, double *);
    void     r_zero_ends(double *);
    void     r_abs(double *);

    void     r_sqr(double *);
    void     rr_sqr(const double *, double *);
    void     rr_addsqr(const double *, double *);
    void     r_sqrt(double *);

    void  	 rrr_Sum2Add(const double *, const double *, double *);
    void  	 rrr_Sum(const double *, const double *, double *);
    void  	 rr_Sum(const double *, double *);
    void  	 rrr_Mul(const double *, const double *, double *);
    void  	 rr_Mul(const double *, double *);
    void  	 arrr_Minus(const double, const double *, const double *, double *);
    void  	 rrr_Minus(const double *, const double *, double *);
    void  	 rr_Minus(const double *, double *);

    void  	 rr_Divide(const double *, double *);
    void  	 rrr_Divide(const double *, const double *, double *);

    void  	 rr_daxpy(const double, const double *, double *);

    void     rrrrr_SumMulAdd(const double *, const double *, const double *, const double *, double *);

    //void  	 r_read(const int, const int, double *);
	 void     c_read(const int, double *, const int);
    void     c_write(const int, double *, const int);

    void     cr_fft3d(double *);
    void     rc_fft3d(double *);
    void     cshift1_fftb(const int,const int, const int, const int, double *);
    void     cshift_fftf(const int,const int, const int, const int, double *);
    void     zeroend_fftb(const int,const int, const int, const int, double *);

    void     c_transpose_jk(double *, double *, double *);
    void     t_transpose_jk(double *, double *, double *);

    void     c_transpose_ijk(const int, double *, double *, double *);

    void     t_transpose_ijk(const int, double *, double *, double *);

    void     t_timereverse(double *, double *, double *);
    void     c_timereverse(double *, double *, double *);
    void     c_timereverse_start(double *, double *, double *, const int, const int);
    void     c_timereverse_end(double *, double *, double *, const int);
    int      timereverse_size();

    void     c_setpw(const int *, const double *, double *);
    void     c_addrandom(double *);
    void     r_setrandom(double *);

    /* expand and contract operators */
    void     hr2r_expand(const double *, double *);
    void     r2hr_contract(const double *, double *);

    /* t array operators */
    double * t_alloc();
    void     t_dealloc(double *);
    void     t_read(const int, double *, const int);
    void     t_write(const int, double *, const int);
    void     t_nzero(int, double *);
    void     tc_Mul(const double *, double *);

    /* ptranspose operators */
    void     c_ptranspose_jk_init(const int, bool  *);
    void     c_ptranspose_ijk_init(const int, bool *, bool *);
    void     c_ptranspose1_jk(const int, double *, double *, double *);
    void     c_ptranspose2_jk(const int, double *, double *, double *);
    void     c_ptranspose_ijk(const int, const int, double *, double *, double *);

    void     c_ptranspose1_jk_start(const int, double *, double *, double *, const int, const int);
    void     c_ptranspose_ijk_start(const int, const int, double *, double *, double *, const int, const int);
    void     c_ptranspose1_jk_end(const int, double *, double *, const int);
    void     c_ptranspose_ijk_end(const int, const int, double *, double *, const int);

};
}

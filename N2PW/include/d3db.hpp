#ifndef _D3dB_H_
#define _D3dB_H_
/* d3db.h
   Author - Eric Bylaska

  this class is container of operations for distributed 3d blocks
*/

#include	"Parallel.hpp"
#include	"Mapping3.hpp"

class d3db : public Mapping3 {

   /* transpose indexings */
   int **iq_to_i1,**iq_to_i2;
   int **i1_start,**i2_start;

   /* timereverse indexings */
   int *t_iq_to_i1,*t_iq_to_i2;
   int *t_i1_start,*t_i2_start;

   /* ptrans indexings */
   int ***p_iq_to_i1,***p_iq_to_i2,***p_iz_to_i2;
   int ***p_i1_start,***p_i2_start;

   float *tmpx,*tmpy,*tmpz;

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
        void     rr_copy(const float *,float *);
        void     rr_SMul(const float, const float *, float *);
        void	 r_SMul(const float, float *);
        float	 r_dsum(const float *);
        float	 rr_dot(const float *, const float *);
        void  	 r_zero_ends(float *);
        void  	 r_abs(float *);
        void  	 rrr_Sum(const float *, const float *, float *);
        void  	 rrr_Mul(const float *, const float *, float *);
        void  	 rrr_Minus(const float *, const float *, float *);
        void  	 rrr_Divide(const float *, const float *, float *);
        void  	 rr_daxpy(const float, const float *, float *);
        //void  	 r_read(const int, const int, float *);
	void     c_read(const int, float *, const int);

	void     cr_fft3d(float *);
	void     rc_fft3d(float *);

	void     c_transpose_jk(float *, float *, float *);
	void     t_transpose_jk(float *, float *, float *);

	void     c_transpose_ijk(const int, float *, float *, float *);
	void     t_transpose_ijk(const int, float *, float *, float *);

	void     c_timereverse(float *, float *, float *);
        int	 timereverse_size();

        /* t array operators */
        float * t_alloc();
        void     t_dealloc(float *);
	void     t_read(const int, float *, const int);
};

#endif

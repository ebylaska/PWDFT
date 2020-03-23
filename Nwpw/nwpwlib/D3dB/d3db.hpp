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

   double *tmpx,*tmpy,*tmpz;

public:
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
        void     rr_copy(const double *,double *);
        void     tt_copy(const double *,double *);
        void     rr_SMul(const double, const double *, double *);
        void	 r_SMul(const double, double *);
        double	 r_dsum(const double *);
        double	 rr_dot(const double *, const double *);
        void  	 r_zero_ends(double *);
        void  	 r_abs(double *);
        void  	 rrr_Sum(const double *, const double *, double *);
        void  	 rrr_Mul(const double *, const double *, double *);
        void  	 rrr_Minus(const double *, const double *, double *);
        void  	 rrr_Divide(const double *, const double *, double *);
        void  	 rr_daxpy(const double, const double *, double *);
        //void  	 r_read(const int, const int, double *);
	void     c_read(const int, double *, const int);
	void     c_write(const int, double *, const int);

	void     cr_fft3d(double *);
	void     rc_fft3d(double *);

	void     c_transpose_jk(double *, double *, double *);
	void     t_transpose_jk(double *, double *, double *);

	void     c_transpose_ijk(const int, double *, double *, double *);
	void     t_transpose_ijk(const int, double *, double *, double *);

	void     c_timereverse(double *, double *, double *);
        int	 timereverse_size();

        /* t array operators */
        double * t_alloc();
        void     t_dealloc(double *);
	void     t_read(const int, double *, const int);
};

#endif

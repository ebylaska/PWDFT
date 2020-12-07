#ifndef _PNEB_H_
#define _PNEB_H_
/* Pneb.h
   Author - Eric Bylaska
*/

#pragma once

#include	"gdevice.hpp"
#include	"Control2.hpp"
#include	"Lattice.hpp"
#include	"Parallel.hpp"
#include	"PGrid.hpp"
#include	"d1db.hpp"
#include	"util.hpp"
#include	"nwpw_timing.hpp"

class Pneb : public PGrid, public d1db  {

   //int ispin,ne[2],neq[2];
   int parallelized;

#ifdef NWPW_SYCL
   double *s22_dev, *s21_dev, *s12_dev, *s11_dev, *sa1_dev, *sa0_dev, *st1_dev; // device-side
   double *s22, *s21, *s12, *s11, *sa1, *sa0, *st1; // host_side

   // index, adiff is required for ggm_lambda_sycl()
   std::int64_t* index = nullptr;
   double* adiff = nullptr;
#else
   double *s22, *s21, *s12, *s11, *sa1, *sa0, *st1;
#endif

   int *ma[2],*ma1[2],*ma2[2],*mc[2],*na[2],*nc[2];
   int mcq[2],ncq[2];
   int  ncqmax;

public:

        /* constructors */
	Pneb(Parallel *, Lattice *, Control2&, int, int *);

        /* destructor */
        ~Pneb()
        {
#ifdef NWPW_SYCL
	    cl::sycl::free(index, *get_syclQue());
	    cl::sycl::free(adiff, *get_syclQue());

            cl::sycl::free(s22_dev, *get_syclQue());
            delete [] s22;
#else
            delete [] s22;
#endif
            if (parallelized)
               for (int ms=0; ms<ispin; ++ms)
               {
                  delete [] ma[ms];
                  delete [] ma1[ms];
                  delete [] ma2[ms];
                  delete [] mc[ms];
                  delete [] na[ms];
                  delete [] nc[ms];
               }
        }

        void g_generate_random(double *);
        void g_read(const int, double *);
        void g_write(const int, double *);



        double *g_allocate(const int nb) {
           double *ptr;
           ptr = new double [2*(neq[0]+neq[1])*npack(nb)];
           return ptr;
        }
        void g_deallocate(double *ptr) { delete [] ptr;}


        double *h_allocate() {
           double *ptr;
           ptr = new double [(neq[0]+neq[1])*n2ft3d];
           return ptr;
        }
        void h_deallocate(double *ptr) { delete [] ptr;}

        int m_size(const int mb) {
           int nsize;
           if (mb==-1) nsize = ne[0]*ne[0] + ne[1]*ne[1];
           else nsize = ne[mb]*ne[mb];
           return nsize;
        }
        double *m_allocate(const int mb, const int nblock) {
           double *ptr;
           int nsize;
           if (mb==-1)
              nsize = ne[0]*ne[0] + ne[1]*ne[1];
           else
              nsize = ne[mb]*ne[mb];

           ptr = new double [nblock*nsize];
           return ptr;
        }
        void m_deallocate(double *ptr) { delete [] ptr;}

        double gg_traceall(double *, double *);
        void gg_copy(double *, double *);
        void g_zero(double *);
        void hr_aSumSqr(const double, double *, double *);
        void ggm_sym_Multiply(double *, double *, double *);
        void ffm_sym_Multiply(const int, double *, double *, double *);

        void ffm3_sym_Multiply(const int, double *, double *, double*, double*, double*);
        void m_scale_s22_s21_s11(const int, const double, double *s22, double *s21, double *s11);
        void fmf_Multiply(const int, double *, double *, double, double *, double);

        void m_scal(const double, double *);
        double m_trace(double *);
        void m_diagonalize(double *, double *);
        void mmm_Multiply(const int, double *, double *, double, double*, double);


        void gh_fftb(double *, double *);
        void ggm_lambda(double, double *, double *, double *);
        void g_ortho(double *);

        void gg_SMul(double, double *, double *);
        void gg_Sum2(double *, double *);
        void ggg_Minus(double *, double *, double *);

#ifdef NWPW_SYCL
        void ggm_lambda_sycl(double, double *, double *, double *);
        void ffm3_sym_Multiply_sycl(const int, const double *, const double *, double*, double*, double* );
        void m_scale_s22_s21_s11_sycl(const int, const double, double *s22, double *s21, double *s11);
        void fmf_Multiply_sycl(const int, double *, double *, double, double *, double);
        void mmm_Multiply_sycl(const int, double *, double *, double, double*, double);
        void ggm_sym_Multiply_sycl(double *psi1, double *psi2, double *hml);
#endif

};

#endif

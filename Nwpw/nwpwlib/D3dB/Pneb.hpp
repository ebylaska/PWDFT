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
   float *s22_dev, *s21_dev, *s12_dev, *s11_dev, *sa1_dev, *sa0_dev, *st1_dev; // device-side
   float *s22, *s21, *s12, *s11, *sa1, *sa0, *st1; // host_side

   // index, adiff is required for ggm_lambda_sycl()
   std::int64_t* index = nullptr;
   float* adiff = nullptr;
#else
   float *s22, *s21, *s12, *s11, *sa1, *sa0, *st1;
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

        void g_generate_random(float *);
        void g_read(const int, float *);
        void g_write(const int, float *);



        float *g_allocate(const int nb) {
           float *ptr;
           ptr = new float [2*(neq[0]+neq[1])*npack(nb)];
           return ptr;
        }
        void g_deallocate(float *ptr) { delete [] ptr;}


        float *h_allocate() {
           float *ptr;
           ptr = new float [(neq[0]+neq[1])*n2ft3d];
           return ptr;
        }
        void h_deallocate(float *ptr) { delete [] ptr;}

        int m_size(const int mb) {
           int nsize;
           if (mb==-1) nsize = ne[0]*ne[0] + ne[1]*ne[1];
           else nsize = ne[mb]*ne[mb];
           return nsize;
        }
        float *m_allocate(const int mb, const int nblock) {
           float *ptr;
           int nsize;
           if (mb==-1)
              nsize = ne[0]*ne[0] + ne[1]*ne[1];
           else
              nsize = ne[mb]*ne[mb];

           ptr = new float [nblock*nsize];
           return ptr;
        }
        void m_deallocate(float *ptr) { delete [] ptr;}

        float gg_traceall(float *, float *);
        void gg_copy(float *, float *);
        void g_zero(float *);
        void hr_aSumSqr(const float, float *, float *);
        void ggm_sym_Multiply(float *, float *, float *);
        void ffm_sym_Multiply(const int, float *, float *, float *);

        void ffm3_sym_Multiply(const int, float *, float *, float*, float*, float*);
        void m_scale_s22_s21_s11(const int, const float, float *s22, float *s21, float *s11);
        void fmf_Multiply(const int, float *, float *, float, float *, float);

        void m_scal(const float, float *);
        float m_trace(float *);
        void m_diagonalize(float *, float *);
        void mmm_Multiply(const int, float *, float *, float, float*, float);


        void gh_fftb(float *, float *);
        void ggm_lambda(float, float *, float *, float *);
        void g_ortho(float *);

        void gg_SMul(float, float *, float *);
        void gg_Sum2(float *, float *);
        void ggg_Minus(float *, float *, float *);

#ifdef NWPW_SYCL
        void ggm_lambda_sycl(float, float *, float *, float *);
        void ffm3_sym_Multiply_sycl(const int, const float *, const float *, float*, float*, float* );
        void m_scale_s22_s21_s11_sycl(const int, const float, float *s22, float *s21, float *s11);
        void fmf_Multiply_sycl(const int, float *, float *, float, float *, float);
        void mmm_Multiply_sycl(const int, float *, float *, float, float*, float);
        void ggm_sym_Multiply_sycl(float *psi1, float *psi2, float *hml);
#endif

};

#endif

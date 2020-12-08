#ifndef _PNEB_H_
#define _PNEB_H_
/* Pneb.h
   Author - Eric Bylaska
*/

#include	"Parallel.hpp"
#include	"PGrid.hpp"
#include	"d1db.hpp"

class Pneb : public PGrid, public d1db  {

   //int ispin,ne[2],neq[2];
   int parallelized;
   float *s22,*s21,*s12,*s11,*sa1,*sa0,*st1;
   int *ma[2],*ma1[2],*ma2[2],*mc[2],*na[2],*nc[2];
   int mcq[2],ncq[2];
   int  ncqmax;

public:

        /* constructor */
	Pneb(Parallel *, int, int *);

        /* destructor */
        ~Pneb() 
        { 
            delete [] s22;

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

        void g_read(const int, float *);



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
        void fmf_Multiply(const int, float *, float *, float, float *, float);

        void m_scal(const float, float *);
        float m_trace(float *);
        void m_diagonalize(float *, float *);
        void mmm_Multiply(const int, float *, float *, float, float*, float);
        float m_dmax(const int, float *);
        void m_scale_s22(const int, const float, float *);
        void m_scale_s21(const int, const float, float *);
        void m_scale_s11(const int, const float, float *);


        void gh_fftb(float *, float *);
        void ggm_lambda(float, float *, float *, float *);

        void gg_SMul(float, float *, float *);
        void gg_Sum2(float *, float *);
        void ggg_Minus(float *, float *, float *);

};

#endif

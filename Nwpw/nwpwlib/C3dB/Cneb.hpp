#ifndef _CNEB_HPP_
#define _CNEB_HPP_

#pragma once

/* Pneb.h
   Author - Eric Bylaska
*/

#pragma once

#include "Control2.hpp"
#include "Lattice.hpp"
#include "Brillouin.hpp"
#include "CGrid.hpp"
#include "Parallel.hpp"
#include "c1db.hpp"
#include "k1db.hpp"
//#include "gdevice.hpp"
#include "nwpw_timing.hpp"
#include "util.hpp"

namespace pwdft {


/**
 * @brief Class representing the PNEB (Parallel NEB) calculation.
 *
 * The `Pneb` class is responsible for managing the Parallel NEB (PNEB) calculation.
 * It inherits properties and methods from several base classes such as `PGrid` and `d1db`.
 * The PNEB calculation involves complex operations related to parallelization,
 * matrix manipulations, and more.
 */

class Cneb : public CGrid, public c1db {

   // int ispin,ne[2],neq[2];
   bool mparallelized=false;
   int parallelized;
   double *s22, *s21, *s12, *s11, *sa1, *sa0, *st1;
   double *mat_tmp, *work1, *work2;
   double *bcolwork, *bwork2, *rwork1, *rwork2;
 
   int *ma[2], *ma1[2], *ma2[2], *mc[2], *na[2], *nc[2];
   int *m2c[2], *n2c[2];
   int mcq[2], ncq[2];
   int m2cq[2],n2cq[2];
   int mcqmax[2],ncqmax[2];
   int mall[3] ,mpack[3], *mindx[3];
   int ncqmax0, npack1_all, nida1_all, n2ft3d_all;
   int g_rnd_algorithm = 1;

   int io_norbs_max = 10;
   bool io_buffer = true;

public:
   /* constructors */
   Cneb(Parallel *, Lattice *, Control2 &, int, int *, Brillouin *);
 
   /* destructor */
   ~Cneb() {
      delete[] s22;
      if (parallelized)
      {
         delete [] mindx[0];
         for (int ms=0; ms<ispin; ++ms) {
            delete[] ma[ms];
            delete[] ma1[ms];
            delete[] ma2[ms];
            delete[] mc[ms];
            delete[] na[ms];
            delete[] nc[ms];
            //delete[] m2c[ms];
            //delete[] n2c[ms];
            delete [] mindx[ms+1];
         }
         delete [] mat_tmp;
         delete [] work1;
         delete [] work2;
         delete [] rwork1;
         delete [] rwork2;
         delete [] bcolwork;
         delete [] bwork2;
      }
   }
 
   void g_generate_random(double *);
   void g_generate1_random(double *);
   void g_generate2_random(double *);
   void g_read(const int, double *);
   void g_read_ne(const int, const int *, const int, double *);
   void g_write(const int, double *);

   void h_read(const int, const int, double *);
   void h_write(const int, const int, const double *);


 
   double *g_allocate(const int nb) {
     double *ptr;
     ptr = new (std::nothrow) double[2*(neq[0]+neq[1]) * CGrid::npack(nb)]();
     return ptr;
   }
   void g_deallocate(double *ptr) { delete[] ptr; }
 
   double *g_nallocate(const int nb, const int nblock) {
     double *ptr;
     ptr = new (std::nothrow) double[nblock * 2*(neq[0]+neq[1])*CGrid::npack(nb)]();
     return ptr;
   }

   double *g_allocate_nbrillq_all() {
     double *ptr;
     ptr = new (std::nothrow) double[nbrillq*2*(neq[0]+neq[1]) * CGrid::npack1_max()]();
     return ptr;
   }
 
   double *h_allocate() 
   {
      double *ptr;
      ptr = new (std::nothrow) double[(neq[0]+neq[1]) * n2ft3d]();
      return ptr;
   }
 
   void h_deallocate(double *ptr) { delete[] ptr; }

   double *h_allocate_nbrillq_all() 
   {
      double *ptr;
      ptr = new (std::nothrow) double[nbrillq*(neq[0]+neq[1]) * n2ft3d]();
      return ptr;
   }

 
   int m_size(const int mb) 
   {
      int nsize;
      if (mb == -1)
        nsize = ne[0]*ne[0] + ne[1]*ne[1];
      else
        nsize = ne[mb]*ne[mb];
      return nsize;
   }
 
   double *m_allocate(const int mb, const int nblock) {
     double *ptr;
     int nsize;
     if (mb == -1)
       nsize = ne[0]*ne[0] + ne[1]*ne[1];
     else
       nsize = ne[mb]*ne[mb];
 
     ptr = new (std::nothrow) double[nblock * nsize]();
     return ptr;
   }
   void m_deallocate(double *ptr) { delete[] ptr; }
 
   double *w_allocate(const int mb, const int nblock) 
   {
      double *ptr;
      int nsize;
      if (mb == -1)
        nsize = 2*(ne[0]*ne[0]+ne[1]*ne[1]);
      else
        nsize = 2*ne[mb]*ne[mb];
     
      ptr = new (std::nothrow) double[nblock * nsize]();
      return ptr;
   }

   void w_deallocate(double *ptr) { delete[] ptr; }


   double *w_allocate_nbrillq_all() 
   {
      double *ptr;
      int nsize;
      nsize = 2*(ne[0]*ne[0]+ne[1]*ne[1]);
      
      ptr = new (std::nothrow) double[nbrillq*nsize]();
      std::memset(ptr,0,nbrillq*nsize*sizeof(double));
      return ptr;
   }

 
 
   double gg_traceall(double *, double *);
   void gg_copy(double *, double *);
   void g_zero(double *);
   void hr_aSumSqr(const double, double *, double *);
   void hhr_aSumMul(const double, const double *, const double *, double *);
 
   void ffw_sym_Multiply(const int, double *, double *, double *);
   void ffw_Multiply(const int, double *, double *, double *);
   void ggw_SVD(double *, double *, double *, double *);

   void ggw_sym_Multiply(double *, double *, double *);
   void ggw_Multiply(double *, double *, double *);

 
   void ffw4_sym_Multiply(const int, const int, double *, double *, double *, double *, double *, double *);
   void ffw3_sym_Multiply(const int, const int, double *, double *, double *, double *, double *);
   void w_scale_s22_s21_s12_s11(const int, const double, double *s22,
                                double *s21, double *s12, double *s11);
   void w_scale_s22_s21_s11(const int, const double, double *s22, double *s21,
                            double *s11);
   void fwf_Multiply(const int, double *, double *, double *, double *, double *);
 
   void fm_QR(const int, double *, double *);
 
   void ww_SCtimesVtrans(const int, const double, double *, double *, double *,
                         double *, double *, double *);
   void ww_SCtimesVtrans2(const int, const double, double *, double *, double *,
                          double *, double *, double *);
   void ww_SCtimesVtrans3(const int, const double, double *, double *, double *,
                          double *, double *, double *);
 
   void m_scal(const double, double *);
   void w_scal(const double, double *);
   double w_trace(double *);
   void w_diagonalize(double *, double *);
   void m_diagonalize(double *, double *);
   void mmm_Multiply(const int, double *, double *, double, double *, double);
   void mmm_Multiply2(const int, double *, double *, double, double *, double);
   void mm_transpose(const int, double *, double *);
   void mm_Kiril_Btransform(const int, double *, double *);
 
   void gh_fftb(double *, double *);
   void gh_fftb0(double *, double *);

   void ggw_lambda(double, double *, double *, double *);
   // void ggm_lambda2(double, double *, double *, double *);
   void ggw_lambda_sic(double, double *, double *, double *);
   void g_ortho(double *);
 
   void gg_SMul(double, double *, double *);
   // void g_SMul1(double, double *);
   void g_Scale(double, double *);
   void gg_Sum2(double *, double *);
   void gg_Minus2(double *, double *);
   void ggg_Minus(double *, double *, double *);
 
   void gg_daxpy(double, double *, double *);
   };
} // namespace pwdft

#endif

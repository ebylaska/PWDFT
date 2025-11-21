#ifndef _PNEB_HPP_
#define _PNEB_HPP_

#pragma once

/* Pneb.h
   Author - Eric Bylaska
*/

#pragma once

#include "Control2.hpp"
#include "Lattice.hpp"
#include "PGrid.hpp"
#include "Parallel.hpp"
#include "d1db.hpp"
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

class Pneb : public PGrid, public d1db {

   // int ispin,ne[2],neq[2];
   bool mparallelized=false;
   int parallelized;
   double *s22, *s21, *s12, *s11, *sa1, *sa0, *st1, *st2;
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
   Pneb(Parallel *, Lattice *, Control2 &, int, int *);
 
   /* destructor */
   ~Pneb() {
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
   void g_generate_excited_random(const int *, double *);
   void g_generate_extra_random(const int, double *);

   void g_read(const int, double *);
   void g_read_ne(const int, const int *, double *);
   void g_read_ne_reverse(const int, const int *, double *);
   void g_read_excited(const int, const int *, double *);
   void g_write(const int, double *);
   void g_write_excited(const int, const int *, double *);
 
   double *g_allocate(const int nb) {
     double *ptr;
     ptr = new (std::nothrow) double[2 * (neq[0] + neq[1]) * PGrid::npack(nb)]();
     return ptr;
   }
   void g_deallocate(double *ptr) { delete[] ptr; }

 
   double *g_nallocate(const int nb, const int nblock) {
     double *ptr;
     ptr = new (std::nothrow) double[nblock * 2 * (neq[0] + neq[1]) *
                                     PGrid::npack(nb)]();
     return ptr;
   }

   
   double *g_allocate_excited(const int ne_excited[], const int nb) 
   {
      double *ptr;
      ptr = new (std::nothrow) double[2 * (ne_excited[0] + ne_excited[1]) * PGrid::npack(nb)]();
      return ptr;
   }

   double *m_nex_allocate(const int mb, const int ne_excited[]) 
   {
      double *ptr;
      int nsize;
      if (mb == -1)
        nsize = ne_excited[0]*ne_excited[0] + ne_excited[1]*ne_excited[1];
      else
        nsize = ne_excited[mb]*ne_excited[mb];
     
      ptr = new (std::nothrow) double[nsize]();
      return ptr;
   }
 
   double *h_allocate() 
   {
      double *ptr;
      ptr = new (std::nothrow) double[(neq[0] + neq[1]) * n2ft3d]();
      return ptr;
   }
 
   void h_deallocate(double *ptr) { delete[] ptr; }
 
   int m_size(const int mb) 
   {
      int nsize;
      if (mb == -1)
        nsize = ne[0] * ne[0] + ne[1] * ne[1];
      else
        nsize = ne[mb] * ne[mb];
      return nsize;
   }
 
   double *m_allocate(const int mb, const int nblock) {
     double *ptr;
     int nsize;
     if (mb == -1)
       nsize = ne[0] * ne[0] + ne[1] * ne[1];
     else
       nsize = ne[mb] * ne[mb];
 
     ptr = new (std::nothrow) double[nblock * nsize]();
     return ptr;
   }
   void m_deallocate(double *ptr) { delete[] ptr; }
 
   double *w_allocate(const int mb, const int nblock) {
     double *ptr;
     int nsize;
     if (mb == -1)
       nsize = 2 * (ne[0] * ne[0] + ne[1] * ne[1]);
     else
       nsize = 2 * ne[mb] * ne[mb];
 
     ptr = new (std::nothrow) double[nblock * nsize]();
     return ptr;
   }
   void w_deallocate(double *ptr) { delete[] ptr; }
 
   double *m4_allocate(const int mb, const int nblock) {
     double *ptr;
     int nsize;
     if (mb == -1)
       nsize = 4 * (ne[0] * ne[0] + ne[1] * ne[1]);
     else
       nsize = 4 * ne[mb] * ne[mb];
 
     ptr = new (std::nothrow) double[nblock * nsize]();
     return ptr;
   }
   void m4_deallocate(double *ptr) { delete[] ptr; }

   // initialize occupations
   /**
       @brief Initializes the occupation array based on spin and electron configuration.
      
       @param nextra Array specifying the number of extra states for each spin channel.
       @param ptr Pointer to a pre-allocated array for storing occupation data.
      
       Occupation values are set to:
       - 0.0 for extra states.
       - 1.0 for regular states.
      
       Example:
       - For ispin = 2, ne = {5, 5}, nextra = {2, 3}, the output array would be:
         ptr = {0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0}.
   */
   void initialize_occupations(const int nextra[], double *ptr)
   {
      for (int ms = 0; ms < ispin; ++ms) 
      {
         int offset = ms * ne[0]; // Precompute offset for indexing
         for (int n = 0; n < ne[ms]; ++n) 
            ptr[offset + n] = (n < nextra[ms]) ? 0.0 : 1.0;
      }
   }
   double* initialize_occupations_with_allocation(const int nextra[]) 
   {
       double* ptr = new double[ne[0] + ne[1]];
       initialize_occupations(nextra, ptr);
       return ptr;
   }


 
   double gg_traceall_excited(const int *, double *, double *);
   double gg_traceall(double *, double *);
   double gg_traceall_occ(double *, double *, double *);
   void gg_copy(double *, double *);
   void g_zero(double *);
   void hr_aSumSqr(const double, double *, double *);
   void hr_aSumSqr_occ(const double, const double *, double *, double *);
   void hhr_aSumMul(const double, const double *, const double *, double *);
 
   void ggm_sym_Multiply(double *, double *, double *);
   void ggm_Multiply(double *, double *, double *);
   void ffm_sym_Multiply(const int, double *, double *, double *);
   void ffm_Multiply(const int, double *, double *, double *);
   void ggm_SVD(double *, double *, double *, double *);
 
   void ffm4_sym_Multiply(const int, double *, double *, double *, double *, double *, double *);
   void ffm3_sym_Multiply(const int, double *, double *, double *, double *, double *);
   void ffm3_Fulls21_sym_Multiply(const int, double *, double *, double *, double *, double *);
   void m_scale_s22_s21_s12_s11(const int, const double, double *s22,
                                double *s21, double *s12, double *s11);
   void m_scale_s22_s21_s11(const int, const double, double *s22, double *s21,
                            double *s11);
   void fmf_Multiply(const int, double *, double *, double, double *, double);
 
   void fm_QR(const int, double *, double *);
 
   void mmm4_AR_to_T4(const int, const double *, const double *, double *);
   void m4_FactorSkew(const int, double *, double *, double *, double *S);
   void m4_RotationSkew(const int, const double, double *, double *, double *,
                        double *, double *, double *);
   void m4_R4_to_MN(const int, const double *, double *, double *);
 
   void mm_SCtimesVtrans(const int, const double, double *, double *, double *,
                         double *, double *, double *);
   void mm_SCtimesVtrans2(const int, const double, double *, double *, double *,
                          double *, double *, double *);
   void mm_SCtimesVtrans3(const int, const double, double *, double *, double *,
                          double *, double *, double *);
 
   void m_scal(const double, double *);
   void m_diag_scal(const double *, double *);
   void m_diag_scal_inv(const double *, double *);
   double m_trace(double *);
   double m_trace_occ(double *, double *);
   void m_diagonalize(double *, double *);
   void m_diagonalize(const int, double *, double *);
   void mmm_Multiply(const int, double *, double *, double, double *, double);
   void mmm_Multiply2(const int, double *, double *, double, double *, double);
   void mmm_Multiply3(const int, double *, double *, double, double *, double);
   void mm_transpose(const int, const double *, double *);
   void mm_Kiril_Btransform(const int, double *, double *);
 
   void gh_fftb(double *, double *);
   void hg_fftf(double *, double *);
   void nbgh_fftb(const int, double *, double *);
   void hgnb_fftf(double *, double *, const int);
   void ggm_lambda(double, double *, double *, double *);
   void ggm_occ_lambda(double, double *, double *, double *, double *);
   // void ggm_lambda2(double, double *, double *, double *);
   void ggm_lambda_sic(double, double *, double *, double *);
   void g_ortho(const int, double *);

   void g_ortho_excited(const int, double *, const int *, double *);
   void g_project_out_filled(double *, const int, double *);
   void g_project_out_filled_extra(const int, const int *, double *);
   void g_project_out_virtual(const int, const int *, const int, double *,  double *);
   void g_project_out_filled_below(double *, const int, const int, double *);
   void g_project_out_filled_above(double *, const int, const int, double *);
   void g_project_out_filled_from_k_up(double *, const int, const int, double *);

   void g_norm(double *);
 
   void gg_SMul(double, double *, double *);
   // void g_SMul1(double, double *);
   void g_Scale(double, double *);
   void gg_Sum2(double *, double *);
   void gg_Minus2(double *, double *);
   void ggg_Minus(double *, double *, double *);
 
   void gg_daxpy(double, double *, double *);

   void gen_Ba_Bs(const int, double *, double *, double *);
   void gen_UD(const int, double *, double *);
   void gen_X(const int, const double *, double *);
   void gen_X(const int, double *, double *, double *, double *, double *, double *, double *, double *, double *, bool *);

   void fnm_to_X(const int, double *, double *, double *, double *, double *);
   void m_eye(const int, double *, double);
   void m_HmldivideDplusD(const int, double *, double *);
   void m_Hmlfweightscale(const int, double *, double *);
   void printNNMatrix(const std::string&, const int, double *);

   void m_0define_occupation(const double, const bool,
                             const int,
                             const double, const double,
                             double *, double *, double *,
                             const int, const double, double *, double *);
   double define_smearfermi(const int, const double *, const double *);
   double add_smearcorrection(const int, const int, const double *, const double *, const double, const double);

   };
} // namespace pwdft

#endif

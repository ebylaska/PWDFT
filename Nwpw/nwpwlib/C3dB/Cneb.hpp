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
   void g_generate_excited_random(const int *, double *);
   void g_generate_extra_random(const int, double *);
   void g_generate_atomic_guess(double *);

   void g_read(const int, double *);
   void g_read_reverse(const int, double *);
   void g_read_excited(const int, const int *, const int,  double *);
   void g_write(const int, double *);
   void g_write_excited(const int, const int *, const int, double *);
   void g_read_ne(const int, const int *, const int, double *);         // probably won't use
   void g_read_ne_reverse(const int, const int *, const int, double *); // probably won't use

   void g_read_occ(const int, double *);
   void g_write_occ(const int, double *);
   //void g_write_occ_old(const int, double *);

   void r_read_occ(const int, double *, const int, const int);
   void r_write_occ(const int, double *, const int, const int);

   void h_read(const int, const int, double *);
   void h_write(const int, const int, const double *);

   double *g_allocate(const int nb) 
   {
      double *ptr;
      ptr = new (std::nothrow) double[2*(neq[0]+neq[1]) * CGrid::npack(nb)]();
      return ptr;
   }
   void g_deallocate(double *ptr) { delete[] ptr; }
 
   double *g_nallocate(const int nb, const int nblock) 
   {
      double *ptr;
      ptr = new (std::nothrow) double[nblock * 2*(neq[0]+neq[1])*CGrid::npack(nb)]();
      return ptr;
   }

   double *g_allocate_nbrillq_all() 
   {
      double *ptr;
      ptr = new (std::nothrow) double[nbrillq*2*(neq[0]+neq[1]) * CGrid::npack1_max()]();
      return ptr;
   }

   double *g_allocate_excited_nbrillq_all(const int nex[]) 
   {
      double *ptr;
      ptr = new (std::nothrow) double[nbrillq*2*(nex[0]+nex[1]) * CGrid::npack1_max()]();
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
      int nsize = 2*(ne[0]*ne[0]+ne[1]*ne[1]);
      
      double *ptr;
      ptr = new (std::nothrow) double[nbrillq*nsize]();
      std::memset(ptr,0,nbrillq*nsize*sizeof(double));
      return ptr;
   }

   double *w_nex_allocate_nbrillq_all(const int nex[0])
   {
      double *ptr;
      int nsize;
      nsize = 2*(nex[0]*nex[0]+nex[1]*nex[1]);
      
      ptr = new (std::nothrow) double[nbrillq*nsize]();
      std::memset(ptr,0,nbrillq*nsize*sizeof(double));
      return ptr;
   }

   // initialize occupations
   /**
    * @brief Initializes occupation numbers for a system with spin channels and extra states.
    *
    * @param nextra Array specifying the number of extra (unoccupied) states for each spin channel.
    * @param ptrb   Pointer to a pre-allocated array of size `nbrillq * (ne[0] + ne[1])`.
    *
    * The occupation array is filled as follows:
    * - For each spin channel (`ispin` = 1 or 2),
    * - Occupation is 0.0 for the first `nextra[ms]` states, and 1.0 for the rest.
    *
    * Example (ispin=2, ne={5,5}, nextra={2,3}):
    *   Output: ptr = {0,0,1,1,1,  0,0,0,1,1}
    */
   void initialize_occupations(const int nextra[], double *ptrb)
   {
      for (int nb=0; nb<nbrillq; ++nb)
      {
         double *ptr = ptrb + nb*(ne[0]+ne[1]);
         for (int ms = 0; ms < ispin; ++ms)
         {
            int offset = ms*ne[0]; // Precompute offset for indexing
            for (int n = 0; n < ne[ms]; ++n)
               ptr[offset + n] = (n < nextra[ms]) ? 0.0 : 1.0;
         }
      }
   }
   /**
    * @brief Allocates and initializes the occupation array.
    *
    * @param nextra Array of extra unoccupied states per spin channel.
    * @return Pointer to a newly allocated array (caller is responsible for deleting it).
    *
    * Internally calls `initialize_occupations` to fill values.
    */
   double* initialize_occupations_with_allocation(const int nextra[])
   {  
       double* ptr = new double[nbrillq*(ne[0] + ne[1])];
       initialize_occupations(nextra, ptr);
       return ptr;
   }     

 
 
   double gg_traceall_excited(const int *, double *, double *);
   double gg_traceall(double *, double *);
   void gg_copy(double *, double *);
   void g_zero(double *);
   void hr_aSumSqr(const double, double *, double *);
   void hr_aSumSqr_occ(const double, double *, double *, double *);
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
   double w_trace_occ(double *,double *);
   void w_diagonalize(double *, double *);
   void m_diagonalize(double *, double *);
   void mmm_Multiply(const int, double *, double *, double, double *, double);
   void mmm_Multiply2(const int, double *, double *, double, double *, double);
   void www_Multiply2(const int, double *, double *, double *, double *, double *);
   void mm_transpose(const int, double *, double *);
   void ww_transpose(const int, double *, double *);
   void ww_hermit_copy(const int, double *, double *);
   void ww_hermit_transpose(const int, double *, double *);
   void mm_Kiril_Btransform(const int, double *, double *);
 
   void gh_fftb(double *, double *);
   void gh_fftb0(double *, double *);

   void ggw_lambda(double, double *, double *, double *);
   // void ggm_lambda2(double, double *, double *, double *);
   void ggw_lambda_sic(double, double *, double *, double *);
   void g_ortho(double *);
   void g_ortho_excited(const int, double *, const int *, double *);
   void g_project_out_filled(const int, double *, const int, double *);
   void g_project_out_virtual(const int, const int, const int *, const int,  double *,  double *);
   void g_project_out_filled_below(const int, double *, const int, const int, double *);
   void g_project_out_filled_above(const int, double *, const int, const int, double *);
   void g_project_out_filled_from_k_up(const int, double *, const int, const int, double *);
   void g_project_out_filled_extra(const int, const int *, double *);

   void g_norm(const int, double *);
 
   void gg_SMul(double, double *, double *);
   // void g_SMul1(double, double *);
   void g_Scale(double, double *);
   void gg_Sum2(double *, double *);
   void gg_Minus2(double *, double *);
   void ggg_Minus(double *, double *, double *);
 
   void gg_daxpy(double, double *, double *);

   void m_0define_occupation(const double, const bool, const int,
                          const double, const double, double *, double *, double *,
                          const int, const double, double *, double *);

   double define_smearfermi(const int, const double *, const double *);
   double add_smearcorrection(const int, const int, const double *, const double *, const double, const double);

   };
} // namespace pwdft

#endif

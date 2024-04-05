#ifndef _GDEVICES_HPP_
#define _GDEVICES_HPP_

#ifdef NWPW_SYCL
#include "gdevices_sycl.hpp"
#elif defined NWPW_CUDA
#include "gdevices_cuda.hpp"
#elif defined NWPW_OPENCL
#include "gdevices_opencl.hpp"
#elif defined NWPW_HIP
#include "gdevices_hip.hpp"
#else
#include "blas.h"
#endif

#include <cstring>   //memset()
#include <stdexcept> // runtime_error()

namespace pwdft {


   /**************************************
    *                                    *
    *          eigsrt_device             *
    *                                    *
    **************************************/
static void eigsrt_device(double *D, double *V, int n) {
   int i, j, k;
   double p;
 
   for (i = 0; i < (n - 1); ++i) {
     k = i;
     p = D[i];
     for (j = i + 1; j < n; ++j)
       if (D[j] >= p) {
         k = j;
         p = D[j];
       }
 
     if (k != i) {
       D[k] = D[i];
       D[i] = p;
       for (j = 0; j < n; ++j) {
         p = V[j + i * n];
         V[j + i * n] = V[j + k * n];
         V[j + k * n] = p;
       }
     }
   }
 }

// just HOST side calls
#if !defined(NWPW_CUDA) && !defined(NWPW_HIP) && !defined(NWPW_SYCL) &&        \
    !defined(NWPW_OPENCL)

#include "fft.h"

class Gdevices {

public:
  int  typegpu = 0;
  bool hasgpu = false;



   /**************************************
    *                                    *
    *              TN4_dgemm             *
    *                                    *
    **************************************/
   void TN4_dgemm(int npack, int ne, double alpha, double *host_a,
                  double *host_b, double beta, double *host_caa,
                  double *host_cab, double *host_cba, double *host_cbb) {
     int one = 1;
     int shift1 = 0;
     int mshift1 = 0;
 
     for (auto k = 1; k <= ne; ++k) {
       DGEMM_PWDFT((char *)"T", (char *)"N", k, one, npack, alpha, host_a, npack,
                   host_a + shift1, npack, beta, host_caa + mshift1, k);
       DGEMM_PWDFT((char *)"T", (char *)"N", k, one, npack, alpha, host_a, npack,
                   host_b + shift1, npack, beta, host_cab + mshift1, k);
       DGEMM_PWDFT((char *)"T", (char *)"N", k, one, npack, alpha, host_b, npack,
                   host_a + shift1, npack, beta, host_cba + mshift1, k);
       DGEMM_PWDFT((char *)"T", (char *)"N", k, one, npack, alpha, host_b, npack,
                   host_b + shift1, npack, beta, host_cbb + mshift1, k);
       shift1 += npack;
       mshift1 += ne;
     }
 
     // DGEMM_PWDFT((char *) "T",(char *)
     // "N",ne,ne,npack,alpha,host_a,npack,host_a,npack,beta,host_caa,ne);
     // DGEMM_PWDFT((char *) "T",(char *)
     // "N",ne,ne,npack,alpha,host_a,npack,host_b,npack,beta,host_cab,ne);
     // DGEMM_PWDFT((char *) "T",(char *)
     // "N",ne,ne,npack,alpha,host_b,npack,host_b,npack,beta,host_cbb,ne);
   }

   /**************************************
    *                                    *
    *              TN3_dgemm             *
    *                                    *
    **************************************/
   void TN3_dgemm(int npack, int ne, double alpha, double *host_a,
                  double *host_b, double beta, double *host_caa,
                  double *host_cab, double *host_cbb) {
     int one = 1;
     int shift1 = 0;
     int mshift1 = 0;
 
     for (auto k = 1; k <= ne; ++k) {
       DGEMM_PWDFT((char *)"T", (char *)"N", k, one, npack, alpha, host_a, npack,
                   host_a + shift1, npack, beta, host_caa + mshift1, k);
       DGEMM_PWDFT((char *)"T", (char *)"N", k, one, npack, alpha, host_a, npack,
                   host_b + shift1, npack, beta, host_cab + mshift1, k);
       DGEMM_PWDFT((char *)"T", (char *)"N", k, one, npack, alpha, host_b, npack,
                   host_b + shift1, npack, beta, host_cbb + mshift1, k);
       shift1 += npack;
       mshift1 += ne;
     }
 
     // DGEMM_PWDFT((char *) "T",(char *)
     // "N",ne,ne,npack,alpha,host_a,npack,host_a,npack,beta,host_caa,ne);
     // DGEMM_PWDFT((char *) "T",(char *)
     // "N",ne,ne,npack,alpha,host_a,npack,host_b,npack,beta,host_cab,ne);
     // DGEMM_PWDFT((char *) "T",(char *)
     // "N",ne,ne,npack,alpha,host_b,npack,host_b,npack,beta,host_cbb,ne);
   }

   /**************************************
    *                                    *
    *              TN1_dgemm             *
    *                                    *
    **************************************/
   void TN1_dgemm(int npack, int ne, double alpha, double *host_a,
                  double *host_b, double beta, double *host_c) {
     DGEMM_PWDFT((char *)"T", (char *)"N", ne, ne, npack, alpha, host_a, npack,
                 host_b, npack, beta, host_c, ne);
   }

   /**************************************
    *                                    *
    *              TN_dgemm              *
    *                                    *
    **************************************/
   void TN_dgemm(int ne, int nprj, int npack, double alpha, double *host_a,
                 double *host_b, double beta, double *host_c) {
     DGEMM_PWDFT((char *)"T", (char *)"N", ne, nprj, npack, alpha, host_a, npack,
                 host_b, npack, beta, host_c, ne);
   }


   /**************************************
    *                                    *
    *              T_free                *
    *                                    *
    **************************************/
   void T_free() {}

   /**************************************
    *                                    *
    *              NN_dgemm              *
    *                                    *
    **************************************/
   void NN_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b,
                 double beta, double *host_c) {
     DGEMM_PWDFT((char *)"N", (char *)"N", npack, ne, ne, alpha, host_a, npack,
                 host_b, ne, beta, host_c, npack);
   }


   /**************************************
    *                                    *
    *              NN_dgemm1             *
    *                                    *
    **************************************/
   void NN_dgemm1(int m, int n, int k, 
                  double alpha, 
                  double *host_a, int lda, 
                  double *host_b, int ldb,
                  double beta, 
                  double *host_c,int ldc) {
      DGEMM_PWDFT((char *)"N", (char *)"N", m, n, k, alpha, host_a, lda, host_b, ldb, beta, host_c, ldc);
   }


   /**************************************
    *                                    *
    *              TN_dgemm2             *
    *                                    *
    **************************************/
   void TN_dgemm2(int m, int n, int k, 
                  double alpha, 
                  double *host_a, int lda, 
                  double *host_b, int ldb,
                  double beta, 
                  double *host_c,int ldc) {
      DGEMM_PWDFT((char *)"T", (char *)"N", m, n, k, alpha, host_a, lda, host_b, ldb, beta, host_c, ldc);
   }

   /**************************************
    *                                    *
    *              TN_dgemm2c            *
    *                                    *
    **************************************/
   void TN_dgemm2c(int n, int m, int npack2, int nida2, 
                  double *host_a, double *host_b, double *host_c) {
      double rtwo  = 2.0;
      double rone  = 1.0;
      double rmone = -1.0;
      double rzero = 0.0;
 
      DGEMM_PWDFT((char *)"T", (char *)"N", n,m,npack2,rtwo,host_a,npack2,host_b,npack2,rzero,host_c,n);
      if (nida2>0)
         DGEMM_PWDFT((char *)"T", (char *)"N", n,m,nida2,rmone,host_a,npack2,host_b,npack2,rone,host_c,n);
   }


   /**************************************
    *                                    *
    *              NT_dgemm3             *
    *                                    *
    **************************************/
   void NT_dgemm3(int m, int n, int k, 
                  double alpha, 
                  double *host_a, int lda, 
                  double *host_b, int ldb,
                  double beta, 
                  double *host_c,int ldc) {
      DGEMM_PWDFT((char *)"N", (char *)"T", m, n, k, alpha, host_a, lda, host_b, ldb, beta, host_c, ldc);
   }


   /**************************************
    *                                    *
    *              NT_dgemm              *
    *                                    *
    **************************************/
   void NT_dgemm(int npack, int ne, int nprj, double alpha, double *host_a,
                 double *host_b, double beta, double *host_c) {
 
      DGEMM_PWDFT((char *)"N", (char *)"T", npack, ne, nprj, alpha, host_a, npack,
                  host_b, ne, beta, host_c, npack);
   }

   /**************************************
    *                                    *
    *              MM6_dgemm             *
    *                                    *
    **************************************/
   void MM6_dgemm(int ne, double *host_s21, double *host_s12, double *host_s11,
                  double *host_sa0, double *host_sa1, double *host_st1) {
      double rzero = 0.0;
      double rone = 1.0;
     
      // mmm_Multiply(ms, s21, sa0, 1.0, sa1, 1.0);
      DGEMM_PWDFT((char *)"N", (char *)"N", ne, ne, ne, rone, host_s21, ne,
                  host_sa0, ne, rone, host_sa1, ne);
     
      // mmm_Multiply(ms, sa0, s12, 1.0, sa1, 1.0);
      DGEMM_PWDFT((char *)"N", (char *)"N", ne, ne, ne, rone, host_sa0, ne,
                  host_s12, ne, rone, host_sa1, ne);
     
      // mmm_Multiply(ms, s11, sa0, 1.0, st1, 0.0);
      DGEMM_PWDFT((char *)"N", (char *)"N", ne, ne, ne, rone, host_s11, ne,
                  host_sa0, ne, rzero, host_st1, ne);
     
      // mmm_Multiply(ms, sa0, st1, 1.0, sa1, 1.0);
      DGEMM_PWDFT((char *)"N", (char *)"N", ne, ne, ne, rone, host_sa0, ne,
                  host_st1, ne, rone, host_sa1, ne);
   }

   /**************************************
    *                                    *
    *              NN_eigensolver        *
    *                                    *
    **************************************/
   void NN_eigensolver(int ispin, int ne[], double *host_hml, double *host_eig) {
      int n, ierr;
      int nn = ne[0] * ne[0] + 14;
      double xmp1[nn];
      // double *xmp1 = new (std::nothrow) double[nn]();
     
      int shift1 = 0;
      int shift2 = 0;
      for (int ms = 0; ms < ispin; ++ms) 
      {
         n = ne[ms];
        
         // eigen_(&n,&n,&hml[shift2],&eig[shift1],xmp1,&ierr);
         //  d3db::parall->Barrier();
         EIGEN_PWDFT(n, host_hml + shift2, host_eig + shift1, xmp1, nn, ierr);
         // if (ierr != 0) throw std::runtime_error(std::string("NWPW Error:
         // EIGEN_PWDFT failed!"));
        
         eigsrt_device(host_eig + shift1, host_hml + shift2, n);
         shift1 += ne[0];
         shift2 += ne[0] * ne[0];
      }
   }



  /// DOUBLE COMPLEX BLAS

   /**************************************
    *                                    *
    *              NN1_zgemm             *
    *                                    *
    **************************************/
   void NN1_zgemm(int npack1, int npack, int ne, double *alpha, double *host_a, double *host_b,
                 double *beta, double *host_c) {
     ZGEMM_PWDFT((char *)"N", (char *)"N", npack, ne, ne, alpha, host_a, npack1,
                 host_b, ne, beta, host_c, npack1);
   }              

   /**************************************
    *                                    *
    *              CN1_zgemm             *
    *                                    *
    **************************************/
   void CN1_zgemm(int npack1, int npack, int ne, double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) {
     ZGEMM_PWDFT((char *)"C", (char *)"N", ne, ne, npack, alpha, host_a, npack1,
                 host_b, npack1, beta, host_c, ne);
   }  

   /**************************************
    *                                    *
    *              CN2_zgemm             *
    *                                    *
    **************************************/
   void CN2_zgemm(int ne, int nprj, int npack, int npack1_max, double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) {
      ZGEMM_PWDFT((char *)"C", (char *)"N", ne, nprj, npack, alpha, host_a, npack1_max,
                  host_b, npack1_max, beta, host_c, ne);
   }

   /**************************************
    *                                    *
    *           CN2_stride_zgemm         *
    *                                    *
    **************************************/
   void CN2_stride_zgemm(int ne, int nprj, int npack, int npack1_max, double *alpha, double *host_a,
                         double *host_b, double *beta, double *host_c) {
      ZGEMM_PWDFT((char *)"C", (char *)"N", ne, nprj, npack, alpha, host_a, npack1_max,
                  host_b, npack1_max, beta, host_c, ne);
   }

   /**************************************
    *                                    *
    *              NC2_zgemm             *
    *                                    *
    **************************************/
   void NC2_zgemm(int npack1, int npack, int ne, int nprj,  double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) {
      ZGEMM_PWDFT((char *)"N", (char *)"C", npack,ne,nprj, alpha, host_a, npack1,
                  host_b, ne, beta, host_c, npack1);
   }

   /**************************************
    *                                    *
    *             NC2_stride_zgemm       *
    *                                    *
    **************************************/
   void NC2_stride_zgemm(int npack1, int npack, int ne, int nprj,  double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) {
      ZGEMM_PWDFT((char *)"N", (char *)"C", npack,ne,nprj, alpha, host_a, npack1,
                  host_b, ne, beta, host_c, npack1);
   }

                 
   /**************************************
    *                                    *
    *              NN_zgemm              *
    *                                    *
    **************************************/
   void NN_zgemm(int m, int n, int k,
                  double *alpha,
                  double *host_a, int lda,
                  double *host_b, int ldb,
                  double *beta,
                  double *host_c,int ldc) {
      ZGEMM_PWDFT((char *)"N", (char *)"N", m, n, k, alpha, host_a, lda, host_b, ldb, beta, host_c, ldc);
   }             

   /**************************************
    *                                    *
    *              CN_zgemm              *
    *                                    *
    **************************************/
   void CN_zgemm(int m, int n, int k,
                  double *alpha,
                  double *host_a, int lda,
                  double *host_b, int ldb,
                  double *beta,
                  double *host_c,int ldc) {
      ZGEMM_PWDFT((char *)"C", (char *)"N", m, n, k, alpha, host_a, lda, host_b, ldb, beta, host_c, ldc);
   } 


   /**************************************
    *                                    *
    *              NC_zgemm              *
    *                                    *
    **************************************/
   void NC_zgemm(int m, int n, int k,
                  double *alpha,
                  double *host_a, int lda,
                  double *host_b, int ldb,
                  double *beta,
                  double *host_c,int ldc) {
      ZGEMM_PWDFT((char *)"N", (char *)"C", m, n, k, alpha, host_a, lda, host_b, ldb, beta, host_c, ldc);
   }

 
   /**************************************
    *                                    *
    *              CN3_zgemm             *
    *                                    *
    **************************************/
   void CN3_zgemm(int npack1, int npack, int ne, double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_caa,
                  double *host_cab, double *host_cbb) {
     int one = 1;
     int shift1 = 0;
     int mshift1 = 0;
       
     for (auto k = 1; k <= ne; ++k) 
     {
        ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack, 
                    alpha, 
                    host_a, npack1,
                    host_a + shift1, npack1, 
                    beta,  
                    host_caa + mshift1, k);
        ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack, 
                    alpha, 
                    host_a, npack1,
                    host_b + shift1, npack1, 
                    beta, 
                    host_cab + mshift1, k);
        ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack, 
                    alpha, 
                    host_b, npack1,
                    host_b + shift1, npack1,
                    beta, 
                    host_cbb + mshift1, k);
        shift1 += 2*npack1; 
        mshift1 += 2*ne;
     }
   }

   /**************************************
    *                                    *
    *              CN4_zgemm             *
    *                                    *
    **************************************/
   void CN4_zgemm(int npack1, int npack, int ne, double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_caa,
                  double *host_cab, double *host_cba, double *host_cbb) {
     int one = 1;
     int shift1 = 0;
     int mshift1 = 0;
         
     for (auto k = 1; k <= ne; ++k) 
     {
        ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack, 
                    alpha, 
                    host_a, npack1,
                    host_a + shift1, npack1, 
                    beta, 
                    host_caa + mshift1, k);
        ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack, 
                    alpha, 
                    host_a, npack1,
                    host_b + shift1, npack1, 
                    beta, 
                    host_cab + mshift1, k);
        ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack, 
                    alpha, 
                    host_b, npack1,
                    host_a + shift1, npack1, 
                    beta, 
                    host_cba + mshift1, k);
        ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack, 
                    alpha, 
                    host_b, npack1,
                    host_b + shift1, npack1, 
                    beta, 
                    host_cbb + mshift1, k);
        shift1  += 2*npack1;
        mshift1 += 2*ne;
     }
 
   }

   /**************************************
    *                                    *
    *              WW6_zgemm             *
    *                                    *
    **************************************/
   void WW6_zgemm(int ne, double *host_s21, double *host_s12, double *host_s11,
                  double *host_sa0, double *host_sa1, double *host_st1) {
      double rone[2]  = {1.0,0.0};
      double rzero[2] = {0.0,0.0};
    
      // www_Multiply1(ms, s21, sa0, 1.0, sa1, 1.0);
      ZGEMM_PWDFT((char *)"N", (char *)"N", ne, ne, ne, rone, host_s21, ne,
                  host_sa0, ne, rone, host_sa1, ne);
                   
      // www_Multiply2(ms, sa0, s12, 1.0, sa1, 1.0);
      ZGEMM_PWDFT((char *)"C", (char *)"N", ne, ne, ne, rone, host_sa0, ne,
                  host_s12, ne, rone, host_sa1, ne);
     
      // www_Multiply3(ms, s11, sa0, 1.0, st1, 0.0);
      ZGEMM_PWDFT((char *)"N", (char *)"C", ne, ne, ne, rone, host_s11, ne,
                  host_sa0, ne, rzero, host_st1, ne);
                   
      // www_Multiply1(ms, sa0, st1, 1.0, sa1, 1.0);
      ZGEMM_PWDFT((char *)"N", (char *)"N", ne, ne, ne, rone, host_sa0, ne,
                  host_st1, ne, rone, host_sa1, ne);
   }  


   /**************************************
    *                                    *
    *              WW_eigensolver        *
    *                                    *
    **************************************/
   void WW_eigensolver(int ispin, int ne[], double *host_hml, double *host_eig) 
   {
      int n, ierr;
      int nn = ne[0] * ne[0] + 14;
      double xmp1[nn];
      double rmp1[nn];
      // double *xmp1 = new (std::nothrow) double[nn]();
     
      int shift1 = 0;
      int shift2 = 0;
      for (int ms=0; ms<ispin; ++ms) 
      {
         n = ne[ms];
        
         // eigen_(&n,&n,&hml[shift2],&eig[shift1],xmp1,&ierr);
         //  d3db::parall->Barrier();
         ZEIGEN_PWDFT(n, host_hml + shift2, host_eig + shift1, xmp1, nn, rmp1, ierr);
         // if (ierr != 0) throw std::runtime_error(std::string("NWPW Error:
         // EIGEN_PWDFT failed!"));
       
         //eigsrt_device(host_eig + shift1, host_hml + shift2, n);
         shift1 += ne[0];
         shift2 += 2*ne[0]*ne[0];
      } 
   }


   /**************************************
    *                                    *
    *          batch_rfftx_tmpx          *
    *                                    *
    **************************************/
   void batch_rfftx_tmpx(bool forward, int nx, int nq, int n2ft3d, double *a, double *tmpx) 
   {
      int nxh2 = nx + 2;
      if (forward) 
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q) 
         {
            drfftf_(&nx, a + indx, tmpx);
            indx += nxh2;
         }
         indx = 1;
         for (auto j = 0; j < (nq); ++j) 
         {
            for (auto i = nx; i >= 2; --i) 
            {
               a[indx + i - 1] = a[indx + i - 2];
            }
            a[indx + 1 - 1] = 0.0;
            a[indx + nx + 1 - 1] = 0.0;
            indx += nxh2;
         }
      } 
      else 
      {
         int indx = 1;
         for (auto j = 0; j < nq; ++j) 
         {
            for (auto i = 2; i <= nx; ++i)
               a[indx + i - 2] = a[indx + i - 1];
            indx += nxh2;
         }
         indx = 0;
         for (auto q = 0; q < nq; ++q) 
         {
            drfftb_(&nx, a + indx, tmpx);
            indx += nxh2;
         }
      }
   }


   /**************************************
    *                                    *
    *          batch_cfftx_tmpx          *
    *                                    *
    **************************************/
   void batch_cfftx_tmpx(bool forward, int nx, int nq, int n2ft3d, double *a, double *tmpx)
   {
      if (forward)
      {
         int indx = 0;
         for (auto q=0; q<nq; ++q)
         {
            dcfftf_(&nx, a + indx, tmpx);
            indx += (2*nx);
         }
      }
      else
      {
         int indx = 0;
         for (auto q=0; q<nq; ++q)
         {
            dcfftb_(&nx, a + indx, tmpx);
            indx += (2*nx);
         }
      }
   }

   /**************************************
    *                                    *
    *          batch_cffty_tmpy          *
    *                                    *
    **************************************/
   void batch_cffty_tmpy(bool forward, int ny, int nq, int n2ft3d, double *a, double *tmpy) 
   {
      if (forward) 
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q) 
         {
            dcfftf_(&ny, a + indx, tmpy);
            indx += (2*ny);
         }
      } 
      else 
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q) 
         {
            dcfftb_(&ny, a + indx, tmpy);
            indx += (2*ny);
         }
      }
   }

   /**************************************
    *                                    *
    *          batch_cffty_tmpy_zero     *
    *                                    *
    **************************************/
   void batch_cffty_tmpy_zero(bool forward, int ny, int nq, int n2ft3d, double *a, double *tmpy, bool *zero) 
   {
      if (forward) 
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q) 
         {
            if (!zero[q])
               dcfftf_(&ny, a + indx, tmpy);
            indx += (2 * ny);
         }
      } 
      else 
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q) 
         {
            if (!zero[q])
               dcfftb_(&ny, a + indx, tmpy);
            indx += (2 * ny);
         }
      }
   }

   /**************************************
    *                                    *
    *          batch_cfftz_tmpz          *
    *                                    *
    **************************************/
   void batch_cfftz_tmpz(bool forward, int nz, int nq, int n2ft3d, double *a, double *tmpz) 
   {
      if (forward) 
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q) 
         {
            dcfftf_(&nz, a + indx, tmpz);
            indx += (2 * nz);
         }
      } 
      else 
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q) 
         {
            dcfftb_(&nz, a + indx, tmpz);
            indx += (2 * nz);
         }
      }
   }

   /**************************************
    *                                    *
    *          batch_cfftz_tmpz_zero     *
    *                                    *
    **************************************/
   void batch_cfftz_tmpz_zero(bool forward, int nz, int nq, int n2ft3d, double *a, double *tmpz, bool *zero) 
   {
      if (forward) 
      {
         int indx = 0;
         for (auto q=0; q<nq; ++q) 
         {
            if (!zero[q])
               dcfftf_(&nz, a + indx, tmpz);
            indx += (2*nz);
         }
      } 
      else 
      {
         int indx = 0;
         for (auto q = 0; q < nq; ++q) 
         {
            if (!zero[q])
               dcfftb_(&nz, a + indx, tmpz);
            indx += (2*nz);
         }
      }
   }

};

#endif // !NWPW_CUDA && !NWPW_HIP && !NWPW_SYCL && !NWPW_OPENCL

} // namespace pwdft
#endif // _GDEVICES_HPP_

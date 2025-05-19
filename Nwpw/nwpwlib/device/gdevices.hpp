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


#include <complex>
#include <cstring>   //memset()
#include <stdexcept> // runtime_error()

#include <iostream> 

namespace pwdft {

   // Define a constants for the radix values
  static constexpr int radix_values[] = {17, 16, 11, 9, 8, 7, 6, 5, 4, 3, 2};
//static constexpr int radix_values[] = {17, 11, 9, 8, 7, 6, 5, 4, 3, 2};
//static constexpr int radix_values[] = {17, 11, 9, 8, 7, 6, 5, 4, 3, 2};

/**************************************
 *                                    *
 *     eigsrt_device_complex          *
 *                                    *
 **************************************/
static void eigsrt_device_complex(double *D, double *V, int n)
{
   int i, j, k;
   double p;

   for (i = 0; i < (n - 1); ++i)
   {
      k = i;
      p = D[i];
      for (j = i + 1; j < n; ++j)
      {
         if (D[j] >= p)
         {
            k = j;
            p = D[j];
         }
      }

      if (k != i)
      {
         // Swap eigenvalues
         std::swap(D[i], D[k]);

         // Swap complex eigenvectors column i and k
         for (j = 0; j < n; ++j)
         {
            int i_idx = 2 * (j + i * n);
            int k_idx = 2 * (j + k * n);
            // Real part
            std::swap(V[i_idx], V[k_idx]);
            // Imaginary part
            std::swap(V[i_idx + 1], V[k_idx + 1]);
         }
      }
   }
}


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
    *           TN3_FullCab_dgemm        *
    *                                    *
    **************************************/
   void TN3_FullCab_dgemm(int npack, int ne, double alpha, double *host_a,
                  double *host_b, double beta, double *host_caa,
                  double *host_cab, double *host_cbb) {
     int one = 1;
     int shift1 = 0;
     int mshift1 = 0;

     for (auto k = 1; k <= ne; ++k) {
       DGEMM_PWDFT((char *)"T", (char *)"N", k, one, npack, alpha, host_a, npack,
                   host_a + shift1, npack, beta, host_caa + mshift1, k);
       //DGEMM_PWDFT((char *)"T", (char *)"N", k, one, npack, alpha, host_a, npack,
       //            host_b + shift1, npack, beta, host_cab + mshift1, k);
       DGEMM_PWDFT((char *)"T", (char *)"N", k, one, npack, alpha, host_b, npack,
                   host_b + shift1, npack, beta, host_cbb + mshift1, k);
       shift1 += npack;
       mshift1 += ne;
     }

     // DGEMM_PWDFT((char *) "T",(char *)
     // "N",ne,ne,npack,alpha,host_a,npack,host_a,npack,beta,host_caa,ne);

      DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_a,npack,host_b,npack,beta,host_cab,ne);

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
   void NN_eigensolver(int ispin, int ne[], double *host_hml, double *host_eig) 
   {
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

   /**************************************
    *                                    *
    *          NN_eigensolver0           *
    *                                    *
    **************************************/
   void NN_eigensolver0(int n, double *host_hml, double *host_eig) 
   {
      int ierr;
      int nn = n*n + 14;
      double xmp1[nn];

      EIGEN_PWDFT(n, host_hml, host_eig, xmp1, nn, ierr);
      eigsrt_device(host_eig, host_hml, n);
   }



  /// DOUBLE COMPLEX BLAS

   /**************************************
    *                                    *
    *              NN1_zgemm             *
    *                                    *
    **************************************/
   void NN1_zgemm(int npack1_max, int npack, int ne, double *alpha, double *host_a, double *host_b,
                 double *beta, double *host_c) {
     ZGEMM_PWDFT((char *)"N", (char *)"N", npack, ne, ne, alpha, host_a, npack1_max,
                 host_b, ne, beta, host_c, npack1_max);
   }              

   /**************************************
    *                                    *
    *              CN1_zgemm             *
    *                                    *
    **************************************/
   void CN1_zgemm(int npack1_max, int npack, int ne, double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) {
     //std::cout << "CN1_zgemm: ne=" << ne << " npack=" << npack << " npack1_max=" << npack1_max << " alpha=" << alpha[0] << " " << alpha[1]  << " beta=" << beta[0] << " " << beta[1] << std::endl;
     ZGEMM_PWDFT((char *)"C", (char *)"N", ne, ne, npack, alpha, host_a, npack1_max,
                 host_b, npack1_max, beta, host_c, ne);
     //std::cout << "hosta=" << host_a[0] << " " << host_a[1] << std::endl;
     //std::cout << "hostb=" << host_b[0] << " " << host_b[1] << std::endl;
     //std::cout << "hostc=" << host_c[0] << " " << host_c[1] << std::endl;
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
   void NC2_zgemm(int npack1_max, int npack, int ne, int nprj,  double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) {
      ZGEMM_PWDFT((char *)"N", (char *)"C", npack,ne,nprj, alpha, host_a, npack1_max,
                  host_b, ne, beta, host_c, npack1_max);
   }

   /**************************************
    *                                    *
    *             NC2_stride_zgemm       *
    *                                    *
    **************************************/
   void NC2_stride_zgemm(int npack1_max, int npack, int ne, int nprj,  double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_c) {
      ZGEMM_PWDFT((char *)"N", (char *)"C", npack,ne,nprj, alpha, host_a, npack1_max,
                  host_b, ne, beta, host_c, npack1_max);
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
   void CN3_zgemm(int npack1_max, int npack, int ne, double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_caa,
                  double *host_cab, double *host_cbb) {
     int one = 1;
     int shift1 = 0;
     int mshift1 = 0;
       
     for (auto k = 1; k <= ne; ++k) 
     {
        ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack, 
                    alpha, 
                    host_a, npack1_max,
                    host_a + shift1, npack1_max, 
                    beta,  
                    host_caa + mshift1, k);
        ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack, 
                    alpha, 
                    host_a, npack1_max,
                    host_b + shift1, npack1_max, 
                    beta, 
                    host_cab + mshift1, k);
        ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack, 
                    alpha, 
                    host_b, npack1_max,
                    host_b + shift1, npack1_max,
                    beta, 
                    host_cbb + mshift1, k);
        shift1 += 2*npack1_max; 
        mshift1 += 2*ne;
     }
   }

   /**************************************
    *                                    *
    *              CN4_zgemm             *
    *                                    *
    **************************************/
   void CN4_zgemm(int npack1_max, int npack, int ne, double *alpha, double *host_a,
                  double *host_b, double *beta, double *host_caa,
                  double *host_cab, double *host_cba, double *host_cbb) {
      int one = 1;
      int shift1 = 0;
      int mshift1 = 0;
          
      for (auto k = 1; k <= ne; ++k) 
      {
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack, 
                     alpha, 
                     host_a, npack1_max,
                     host_a + shift1, npack1_max, 
                     beta, 
                     host_caa + mshift1, k);
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack, 
                     alpha, 
                     host_a, npack1_max,
                     host_b + shift1, npack1_max, 
                     beta, 
                     host_cab + mshift1, k);
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack, 
                     alpha, 
                     host_b, npack1_max,
                     host_a + shift1, npack1_max, 
                     beta, 
                     host_cba + mshift1, k);
         ZGEMM_PWDFT((char *)"C", (char *)"N", k, one, npack, 
                     alpha, 
                     host_b, npack1_max,
                     host_b + shift1, npack1_max,
                     beta, 
                     host_cbb + mshift1, k);
         shift1  += 2*npack1_max;
         mshift1 += 2*ne;
      }
   }

   /**************************************
    *                                    *
    *         computeTrans3_Mult         *
    *                                    *
    **************************************/
   /**
    * @brief Computes the 3D transformation of projection data using transformation matrices.
    *
    * This function computes transformation sums for multiple projection (`prj`) functions and
    * `psi` functions. It uses projection data (`prj`) and input data (`psi`) to calculate
    * intermediate values, and then computes the transformation sums (`sum3`) using the 
    * transformation matrices (`Gx`, `Gy`, `Gz`).
    *
    * @param ne    Number of `psi` functions to process.
    * @param nprj  Number of `prj` functions to process.
    * @param psi   Pointer to the input data array containing the `psi` functions.
    * @param prj   Pointer to the projection data array containing the `prj` functions.
    * @param ng    Grid size for the main computation (number of complex numbers).
    * @param ng0   Reduced grid size for the secondary computation (number of complex numbers).
    * @param Gx    Transformation matrix in the x-direction.
    * @param Gy    Transformation matrix in the y-direction.
    * @param Gz    Transformation matrix in the z-direction.
    * @param xtmp1 Buffer for intermediate computations.
    * @param sum3  Array to store the computed transformation sums.
    * 
    * @note This function can be optimized using SIMD instructions, OpenMP for parallelization,
    *       and GPU acceleration to enhance performance.
    */
/*
   void computeTrans3_Mult(const int ne, const int nprj, 
                           const double *psi, const double *prj,
                           int ng, int ng0,
                           double *Gx, double *Gy, double *Gz, double *xtmp1,
                           double *sum3)
   {
      int one = 1;
      int count3 = 0;
      int nshift = 2*ng;

      for (auto l=0; l<nprj; ++l)
      for (auto n=0; n<ne; ++n)
      {
         // Perform cct_pack_iconjgMul
         const double *a = prj + l*nshift;
         const double *b = psi + n*nshift;
         for (int i=0; i<ng; ++i)
            xtmp1[i] = a[2*i]*b[2*i+1] - a[2*i+1]*b[2*i];
      
         double tsumx = 2.0*DDOT_PWDFT(ng,Gx,one,xtmp1,one);
         double tsumy = 2.0*DDOT_PWDFT(ng,Gy,one,xtmp1,one);
         double tsumz = 2.0*DDOT_PWDFT(ng,Gz,one,xtmp1,one);
         tsumx -= DDOT_PWDFT(ng0,Gx,one,xtmp1,one);
         tsumy -= DDOT_PWDFT(ng0,Gy,one,xtmp1,one);
         tsumz -= DDOT_PWDFT(ng0,Gz,one,xtmp1,one);
      
         sum3[count3]   = tsumx;
         sum3[count3+1] = tsumy;
         sum3[count3+2] = tsumz;
         count3 += 3;
      }
   }
   */

   /**************************************
    *                                    *
    *           computeTrans_Mult        *
    *                                    *
    **************************************/
   /**
    * @brief Computes the transformation of projection data using matrix multiplication.
    *
    * This function computes matrix multiplication between projection (`prj`) and input
    * (`psi`) data using a BLAS DGEMM routine, which is optimized for double-precision
    * general matrix multiplication. It first performs the matrix multiplication with the
    * full grid size (`npack2`), and if a reduced grid size (`npack0`) is specified,
    * subtracts the result using another matrix multiplication.
    *
    * The BLAS DGEMM function can be optimized for various architectures, leveraging
    * specialized hardware and libraries to enhance performance.
    *
    * @param ne    Number of `psi` functions to process.
    * @param nprj  Number of `prj` functions to process.
    * @param psi   Pointer to the input data array containing the `psi` functions.
    * @param prj   Pointer to the projection data array containing the `prj` functions.
    * @param ng    Grid size for the main computation (number of complex numbers).
    * @param ng0   Reduced grid size for the secondary computation (number of complex numbers).
    * @param sum1  Array to store the computed matrix multiplication results.
    */
/*
   void computeTrans_Mult(int ne, int nprj, double alpha, double alpha1, int ng, int ng0,
                          double *psi, double *prj, double beta, double beta1, double *sum1) 
   {
      int npack2 = 2*ng;
      int npack0 = 2*ng0;
      double rtwo  = 2.0;
      double rzero = 0.0;
      double rone  = 1.0;
      double rmone = -1.0;

      DGEMM_PWDFT((char *)"T", (char *)"N",ne,nprj,npack2,alpha,psi,npack2,prj,npack2,beta,sum1,ne);
      if (npack0 > 0) 
      {
         DGEMM_PWDFT((char *)"T", (char *)"N",ne,nprj,npack0,alpha1,psi,npack2,prj,npack2,beta1,sum1,ne);
      }
   }
*/


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
       
         eigsrt_device_complex(host_eig + shift1, host_hml + shift2, n);
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
         for (auto q=0; q<nq; ++q) 
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
   void batch_cffty_tmpy_zero(bool forward, int ny, int nq, int nffts, int n2ft3d, double *a, double *tmpy, bool *zero) 
   {
      if (forward) 
      {
         for (auto i = 0; i<nffts; ++i) 
         {
            int indx = i*n2ft3d;
            for (auto q = 0; q<nq; ++q) 
            {
               if (!zero[q])
                  dcfftf_(&ny, a + indx, tmpy);
               indx += (2 * ny);
            }
         }
      } 
      else 
      {
         for (auto i = 0; i<nffts; ++i) 
         {
            int indx = i*n2ft3d;
            for (auto q = 0; q<nq; ++q) 
            {
               if (!zero[q])
                  dcfftb_(&ny, a + indx, tmpy);
               indx += (2 * ny);
            }
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
   void batch_cfftz_tmpz_zero(bool forward, int nz, int nq, int nffts, int n2ft3d, double *a, double *tmpz, bool *zero) 
   {
      if (forward) 
      {  
         for (auto i = 0; i<nffts; ++i)    
         {
            int indx = i*n2ft3d;
            for (auto q=0; q<nq; ++q) 
            {
               if (!zero[q])
                  dcfftf_(&nz, a + indx, tmpz);
               indx += (2*nz);
            }
         }
      } 
      else 
      {
         for (auto i = 0; i<nffts; ++i)    
         {
            int indx = i*n2ft3d;
            for (auto q = 0; q < nq; ++q) 
            {
              if (!zero[q])
                  dcfftb_(&nz, a + indx, tmpz);
               indx += (2*nz);
            }
         }
      }
   }



   ////////////////////////// special complex-complex fft ////////////////////////////

typedef std::complex<double> complex_t;
         
   /**************************************
    *                                    *
    *             fft_radix              *
    *                                    *
    **************************************/
   inline static void fft_radix(const int n, const int s, bool eo, const int radix, const complex_t* twiddle, complex_t* x, complex_t *y)
   {
      if (n == 1) 
      { 
         // Use std::copy to efficiently copy elements from x to y
         if (eo) std::copy(x, x + s, y); 
         //if (eo) std::memcpy(y, x, s * sizeof(double));
         return;
      }

      const int m = n/radix;
      complex_t Atwiddle[radix*radix];
 
      // Precompute twiddle factors for matrix multiplication
      for (int r2=0; r2<radix; ++r2)
      for (int r1=0; r1<radix; ++r1)
         Atwiddle[r1+r2*radix] = twiddle[(r1*r2)%radix];
 
      for (int p=0; p<m; ++p) 
      {
         for (int q=0; q<s; ++q) 
         {
            complex_t* x1 = x+q+s*p;
            complex_t* y1 = y+q+s*radix*p;
 
            // Initialize the output y1 vector
            for (int r1=0; r1<radix; ++r1) 
               y1[s * r1] = 0.0;
 
            // Matrix-vector multiplication
            for (int r2 = 0; r2 < radix; ++r2)
            for (int r1 = 0; r1 < radix; ++r1)
               y1[s*r1] += x1[s*m*r2] * Atwiddle[r1 + r2*radix];
               //y1[s*r1] += x1[s*m*r2] * twiddle[(r1*r2)%radix];
 
            // Apply phase factor to each result
            for (int r1 = 0; r1 < radix; ++r1)
            {
               //complex_t fac = std::pow(twiddle[radix + r1], p);
               complex_t fac = complex_t(1.0, 0.0);
               for (int ps = 0; ps < p; ++ps)
                  fac *= twiddle[radix + r1];
 
               y1[s * r1] *= fac;
            }
         }
      }
   }

   /**************************************
    *                                    *
    *            fft_twiddle             *
    *                                    *
    **************************************/
   inline static void fft_twiddle(const int n, const complex_t* twiddle, complex_t* x) // Fourier transform
   {
      //complex_t* y = new complex_t[n];
      complex_t y[n];
      int eo = 0;
      int s  = 1;
      int nn = n;
      int nsize = 0;
      while (s<=n) 
      {
         // Identify the largest radix applicable for current nn
         int radix = 2;  // Default to radix-2
         for (int r : radix_values) {
            if (nn % r == 0) {
               radix = r;
               break;
            }
         }
 
         // Perform FFT with the determined radix
         if (eo)
            fft_radix(nn, s, eo, radix, twiddle + nsize, y, x);
         else
            fft_radix(nn, s, eo, radix, twiddle + nsize, x, y);

         nsize += 2*radix;
         nn /= radix;
         s *= radix;
         eo = !eo;  // Toggle the 'even-odd' flag
      }
   }


   /**************************************
    *                                    *
    *         set_sub_fft_twiddle        *
    *                                    *
    **************************************/
   static void set_sub_fft_twiddle(const int isgn, const int n, const int radix, complex_t *twiddle)
   {
      const double theta0 = 2*M_PI/((double) n);
      const double theta_radix = 2*M_PI/((double) radix);
 
      // Calculate radix-specific twiddle factors
      for (int r=0; r<radix; ++r)
         twiddle[r] = complex_t(cos(r*theta_radix), isgn*sin(r*theta_radix));
 
      // Calculate the main twiddle factors for the FFT
      for (int r=0; r<radix; ++r)
         twiddle[radix+r] = complex_t(cos(r*theta0), isgn*sin(r*theta0));
   }



   /**************************************
    *                                    *
    *            set_fft_twiddle         *
    *                                    *
    **************************************/
   void set_fft_twiddle(const int isgn, const int n, double *twiddle) 
   {
      complex_t* complex_twiddle = reinterpret_cast<complex_t*>(twiddle);
      int nsize = 0;
      int s = 1;
      int nn = n;
 
      while (s <= n) 
      {
         bool found = false;
         // Loop through possible radix values to find the largest factor of nn
         for (int radix : radix_values) 
         {
            if (nn % radix == 0) 
            {
               set_sub_fft_twiddle(isgn, nn, radix, complex_twiddle + nsize);
               nsize += 2*radix;
               nn /= radix;
               s *= radix;
               found = true;
               break;
            }
         }
         if (!found) break;
      }
   }

   /**************************************
    *                                    *
    *            size_fft_twiddle        *
    *                                    *
    **************************************/
   int size_fft_twiddle(const int n) 
   {
      int nsize = 0;
      int s = 1;
      int nn = n;
 
      while (s <= n) 
      {
         bool found = false;
         // Loop through possible radix values to find the largest factor of nn
         for (int radix : radix_values) 
         {
            if (nn % radix == 0) 
            {
               nsize += 2 * radix;
               nn /= radix;
               s *= radix;
               found = true;
               break;
            }
         }
         if (!found) break;
      }
 
      return nsize;
   }

   /**************************************
    *                                    *
    *              batch_cfft            *
    *                                    *
    **************************************/
   void batch_cfft(const bool forward, const int nz, const int nq, const int nffts, const int nfft3d, double *a, const double *twiddle, const double *tmpz)
   {
      // Ensure the function processes the right type of data
      // If twiddle is indeed of type complex_t, this cast is necessary
      //complex_t*       complex_a       = reinterpret_cast<complex_t*>(a);
      //const complex_t* complex_twiddle = reinterpret_cast<const complex_t*>(twiddle);

      //int shift = nfft3d;
      int shift = 2*nfft3d;
      int indx = 0;

      // Process each FFT batch
      for (auto q=0; q<nq; ++q)
      {
         (forward ? dcfftf_(&nz, a + indx, tmpz) : dcfftb_(&nz, a + indx, tmpz));
         //fft_twiddle(nz, complex_twiddle, complex_a+indx);

         indx += (shift);
      }
   }

   /**************************************
    *                                    *
    *          batch_cfft_zero           *
    *                                    *
    **************************************/
   void batch_cfft_zero(const bool forward, int nz, int nq, int nffts, int nfft3d, double *a, const double *twiddle, const double *tmpz, const bool *zero) 
   {
      // Ensure the function processes the right type of data
      // If twiddle is indeed of type complex_t, this cast is necessary
      //complex_t*       complex_a       = reinterpret_cast<complex_t*>(a);
      //const complex_t* complex_twiddle = reinterpret_cast<const complex_t*>(twiddle);

      //int shift = nfft3d;
      int shift = 2*nfft3d;
      int indx = 0;

      // Process each FFT batch
      for (auto q=0; q<nq; ++q)
      {
         if (!zero[q])
            (forward ? dcfftf_(&nz, a + indx, tmpz) : dcfftb_(&nz, a + indx, tmpz));
            //fft_twiddle(nz, complex_twiddle, complex_a+indx);

         indx += (shift);
      }
   }


   ////////////////////////// special complex-complex fft ////////////////////////////

};

#endif // !NWPW_CUDA && !NWPW_HIP && !NWPW_SYCL && !NWPW_OPENCL

} // namespace pwdft
#endif // _GDEVICES_HPP_

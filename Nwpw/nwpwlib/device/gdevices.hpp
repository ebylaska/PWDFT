#ifndef _GDEVICES_HPP_
#define _GDEVICES_HPP_

#ifdef NWPW_SYCL
#include        "gdevices_sycl.hpp"
#elif defined NWPW_CUDA
#include        "gdevices_cuda.hpp"
#elif defined NWPW_OPENCL
#include        "gdevices_opencl.hpp"
#elif defined NWPW_HIP
#include        "gdevices_hip.hpp"
#else
#include        "blas.h"
#endif

namespace pwdft {

  // just HOST side calls
#if !defined(NWPW_CUDA) && !defined(NWPW_HIP) && !defined(NWPW_SYCL) && !defined(NWPW_OPENCL)
  
#include "fft.h"

class Gdevices {

public:
     bool hasgpu = false;

     void TN3_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_caa, double *host_cab, double *host_cbb)
     {
        int one = 1;
        int shift1  = 0;
        int mshift1 = 0;

        for (auto k=1; k<=ne; ++k)
        {
           DGEMM_PWDFT((char *) "T",(char *) "N",k,one,npack,alpha,host_a,npack,&host_a[shift1],npack,beta,&host_caa[mshift1],k);
           DGEMM_PWDFT((char *) "T",(char *) "N",k,one,npack,alpha,host_a,npack,&host_b[shift1],npack,beta,&host_cab[mshift1],k);
           DGEMM_PWDFT((char *) "T",(char *) "N",k,one,npack,alpha,host_b,npack,&host_b[shift1],npack,beta,&host_cbb[mshift1],k);
           shift1  += npack;
           mshift1 += ne;
        }

        //DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_a,npack,host_a,npack,beta,host_caa,ne);
        //DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_a,npack,host_b,npack,beta,host_cab,ne);
        //DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_b,npack,host_b,npack,beta,host_cbb,ne);
     }
     void TN1_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
        DGEMM_PWDFT((char *) "T",(char *) "N",ne,ne,npack,alpha,host_a,npack,host_b,npack,beta,host_c,ne);
     }

     void TN_dgemm(int ne, int nprj, int npack, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
        DGEMM_PWDFT((char *) "T",(char *) "N",ne,nprj,npack,alpha,host_a,npack,host_b,npack,beta,host_c,ne);
     }

     void T_free() { }

     void NN_dgemm(int npack, int ne, double alpha, double *host_a, double *host_b, double beta, double *host_c) {
        DGEMM_PWDFT((char *) "N",(char *) "N",npack,ne,ne,alpha,host_a,npack,host_b,ne,beta,host_c,npack);
     }

     void NT_dgemm(int npack, int ne, int nprj, double alpha, double *host_a, double *host_b, double beta, double *host_c) {

        DGEMM_PWDFT((char *) "N",(char *) "T",npack,ne,nprj,alpha,host_a,npack,host_b,ne,beta,host_c,npack);
     }

     void batch_cfftx_tmpx(bool forward, int nx, int nq, int n2ft3d, double *a, double *tmpx) {
        int nxh2 = nx+2;
        if (forward)
        {
           int indx = 0;
           for (auto q=0; q<nq; ++q)
           {
             drfftf_(&nx,a+indx,tmpx);
             indx += nxh2;
           }
           indx = 1;
           for (auto j=0; j<(nq); ++j)
           {
              for (auto i=nx; i>=2; --i)
              {
                 a[indx+i-1] = a[indx+i-2];
              }
              a[indx+1-1]    = 0.0;
              a[indx+nx+1-1] = 0.0;
              indx += nxh2;
           }
        }
        else
        {
           int indx = 1;
           for (auto j=0; j<nq; ++j)
           {
              for (auto i=2; i<=nx; ++i)
                 a[indx+i-2] = a[indx+i-1];
              indx += nxh2;
           }
           indx = 0;
           for (auto q=0; q<nq; ++q)
           {
              drfftb_(&nx,a+indx,tmpx);
              indx += nxh2;
           }
        }
     }
     void batch_cffty_tmpy(bool forward, int ny, int nq, int n2ft3d, double *a, double *tmpy) {
        if (forward)
        {
           int indx = 0;
           for (auto q=0; q<nq; ++q)
           {
              dcfftf_(&ny,a+indx,tmpy);
              indx += (2*ny);
           }
        }
        else
        {
           int indx = 0;
           for (auto q=0; q<nq; ++q)
           {
              dcfftb_(&ny,a+indx,tmpy);
              indx += (2*ny);
           }
        }
     }
     void batch_cffty_tmpy_zero(bool forward, int ny, int nq, int n2ft3d, double *a, double *tmpy, bool *zero) {
        if (forward)
        {
           int indx = 0;
           for (auto q=0; q<nq; ++q)
           {
              if (!zero[q])
                 dcfftf_(&ny,a+indx,tmpy);
              indx += (2*ny);
           }
        }
        else
        {
           int indx = 0;
           for (auto q=0; q<nq; ++q)
           {
              if (!zero[q])
                 dcfftb_(&ny,a+indx,tmpy);
              indx += (2*ny);
           }
        }
     }
     void batch_cfftz_tmpz(bool forward, int nz, int nq, int n2ft3d, double *a, double *tmpz) {
        if (forward)
        {
           int indx = 0;
           for (auto q=0; q<nq; ++q)
           {
              dcfftf_(&nz,a+indx,tmpz);
              indx += (2*nz);
           }
        }
        else
        {
           int indx = 0;
           for (auto q=0; q<nq; ++q)
           {
              dcfftb_(&nz,a+indx,tmpz);
              indx += (2*nz);
           }
        }
     }
     void batch_cfftz_tmpz_zero(bool forward, int nz, int nq, int n2ft3d, double *a, double *tmpz, bool *zero) {
        if (forward)
        {
           int indx = 0;
           for (auto q=0; q<nq; ++q)
           {
              if (!zero[q])
                 dcfftf_(&nz,a+indx,tmpz);
              indx += (2*nz);
           }
        }
        else
        {
           int indx = 0;
           for (auto q=0; q<nq; ++q)
           {
              if (!zero[q])
                 dcfftb_(&nz,a+indx,tmpz);
              indx += (2*nz);
           }
        }
     }


};

#endif // !NWPW_CUDA && !NWPW_HIP && !NWPW_SYCL && !NWPW_OPENCL

}
#endif // _GDEVICES_HPP_

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

};

#endif // !NWPW_CUDA && !NWPW_HIP && !NWPW_SYCL && !NWPW_OPENCL

}
#endif // _GDEVICES_HPP_

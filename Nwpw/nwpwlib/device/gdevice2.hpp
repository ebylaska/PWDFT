#pragma once

#include "gdevices.hpp"

// using namespace pwdft;
namespace pwdft {


class gdevice2 {

Gdevices *mygdevice2;

public:

   gdevice2();

   void TN4_dgemm(int, int, double, double *, double *, double, double *, double *, double *, double *);
   void TN3_dgemm(int, int, double, double *, double *, double, double *, double *, double *);
   void TN1_dgemm(int, int, double, double *, double *, double, double *);
   void TN_dgemm(int, int, int, double, double *, double *, double, double *);
   void T_free();
  
   void NN_dgemm(int, int, double, double *, double *, double, double *);
   void NT_dgemm(int, int, int, double, double *, double *, double, double *);
  
   void MM6_dgemm(int, double *, double *, double *, double *, double *, double *);
  
   void NN_eigensolver(int, int *, double *, double *);
  
   void psi_alloc(int, int, int);
   void psi_dealloc();
  
   void psi_copy_host2gpu(int, int, double *);
   void hpsi_copy_host2gpu(int, int, double *);
  
   void psi_copy_gpu2host(int, int, double *);
   void hpsi_copy_gpu2host(int, int, double *);
  
   void batch_fft_init(int, int, int, int, int, int);
   void batch_fft_end(int, int, int);
   void batch_cfftx(bool, int, int, int, double *);
   void batch_cffty(bool, int, int, int, double *);
   void batch_cfftz(bool, int, int, int, double *);
   void batch_cfftx_tmpx(bool, int, int, int, double *, double *);
   void batch_cffty_tmpy(bool, int, int, int, double *, double *);
   void batch_cfftz_tmpz(bool, int, int, int, double *, double *);
  
   void batch_cffty_tmpy_zero(bool, int, int, int, double *, double *, bool *);
   void batch_cfftz_tmpz_zero(bool, int, int, int, double *, double *, bool *);
};

} // namespace pwdft

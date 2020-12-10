#ifndef _GDEVICE_HPP_
#define _GDEVICE_HPP_

extern void gdevice_TN3_dgemm(int, int, double, double *, double *, double, double *, double *, double *);
extern void gdevice_TN_dgemm(int, int, int, double, double *, double *, double, double *);
extern void gdevice_NN_dgemm(int, int, double, double *, double *, double, double *);
extern void gdevice_NT_dgemm(int, int, int, double, double *, double *, double, double *);

extern void gdevice_psi_alloc(int, int);
extern void gdevice_psi_dealloc();

extern void gdevice_psi_copy_host2gpu(int, int, double *);
extern void gdevice_hpsi_copy_host2gpu(int, int, double *);

extern void gdevice_psi_copy_gpu2host(int, int, double *);
extern void gdevice_hpsi_copy_gpu2host(int, int, double *);

extern void gdevice_batch_fft_init(int, int, int);
extern void gdevice_batch_cfftx(bool,int,int,int,double *);
extern void gdevice_batch_cffty(bool,int,int,int,double *);
extern void gdevice_batch_cfftz(bool,int,int,int,double *);


#ifdef NWPW_SYCL
#include <CL/sycl.hpp>

#define NWPW_SYCL_ERROR( EXPR )                                     \
  try {                                                             \
    EXPR;                                                           \
  }                                                                 \
  catch (cl::sycl::exception const &ex) {                           \
    std::stringstream msg;                                          \
    msg << "SYCL Exception at " << __FILE__ << " : " << __LINE__    \
	<< std::endl;                                               \
    throw(std::runtime_error( ex.what() ));                         \
  }                                                                 \
  catch(std::runtime_error const& ex) {                             \
    std::stringstream msg;                                          \
    msg << "runtime Exception at " << __FILE__ << " : " << __LINE__ \
	<< std::endl;                                               \
    throw(std::runtime_error( ex.what() ));                         \
  }

extern cl::sycl::queue* get_syclQue();


#endif // NWPW_SYCL

#endif

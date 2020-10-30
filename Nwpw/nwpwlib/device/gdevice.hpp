#ifndef _GDEVICE_HPP_
#define _GDEVICE_HPP_

extern void gdevice_TN3_dgemm(int, int, double, double *, double *, double, double *, double *, double *);
extern void gdevice_TN_dgemm(int, int, int, double, double *, double *, double, double *);
extern void gdevice_NN_dgemm(int, int, double, double *, double *, double, double *);
extern void gdevice_NT_dgemm(int, int, int, double, double *, double *, double, double *);

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
extern int get_sycl_mem_index(const size_t memSize);
extern double* get_sycl_mem(const int memIndex);
extern void free_sycl_mem(const int memIndex);
#endif // NWPW_SYCL

#endif

#ifndef _GDEVICE_HPP_
#define _GDEVICE_HPP_

extern void gdevice_TN3_dgemm(int, int, float, float *, float *, float, float *, float *, float *);
extern void gdevice_TN_dgemm(int, int, int, float, float *, float *, float, float *);
extern void gdevice_NN_dgemm(int, int, float, float *, float *, float, float *);
extern void gdevice_NT_dgemm(int, int, int, float, float *, float *, float, float *);

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
extern float* get_sycl_mem(const size_t memInBytes);
extern void free_sycl_mem(float* ptr);
extern float* get_host_mem(const size_t memInBytes);
extern void free_host_mem(float* ptr);

#endif // NWPW_SYCL

#endif

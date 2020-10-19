#ifndef _GDEVICE_HPP_
#define _GDEVICE_HPP_

extern void gdevice_TN3_dgemm(int, int, double, double *, double *, double, double *, double *, double *);
extern void gdevice_TN_dgemm(int, int, int, double, double *, double *, double, double *);
extern void gdevice_NN_dgemm(int, int, double, double *, double *, double, double *);
extern void gdevice_NT_dgemm(int, int, int, double, double *, double *, double, double *);

#ifdef NWPW_SYCL
#include <CL/sycl.hpp>

extern cl::sycl::queue* get_syclQue();
#endif

#endif

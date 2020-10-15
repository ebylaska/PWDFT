#ifndef _GDEVICE_HPP_
#define _GDEVICE_HPP_

#ifdef NWPW_SYCL
#include        <cstdio>
#include        <iostream>
#include        <limits>
#include        <CL/sycl.hpp>
#include        "mkl_blas_sycl.hpp"


class Gdevice {

public:

    auto asyncHandler = [&](cl::sycl::exception_list eL) {
       for (auto& e : eL) {
         try {
           std::rethrow_exception(e);
         } catch (cl::sycl::exception& e) {
           std::cout << e.what() << std::endl;
           std::cout << "fail" << std::endl;
           std::terminate();
         }
       }
    };
    cl::sycl::gpu_selector device_selector;
    cl::sycl::queue device_queue(device_selector,
                                 asyncHandler,
                                 cl::sycl::property_list{cl::sycl::property::queue::in_order{}});
    /* device memory */
    double *dev_a,*dev_b,*dev_c;


    Gdevice(const int npack, const int ne) {
       /* Creating 1D buffers for mkl matrices */
       double* dev_a = cl::sycl::malloc_device<double>(npack*ne, device_queue);
       double* dev_b = cl::sycl::malloc_device<double>(ne*ne,    device_queue);
       double* dev_c = cl::sycl::malloc_device<double>(npack*ne, device_queue);
    }

    ~Gdevice() {
       cl::sycl::free(dev_a, device_queue);
       cl::sycl::free(dev_b, device_queue);
       cl::sycl::free(dev_c, device_queue);
     }
};

#endif
#endif

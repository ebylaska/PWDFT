#ifndef NWPW_GPU_TYPES_H_
#define NWPW_GPU_TYPES_H_

#ifdef NWPW_SYCL
#include <CL/sycl.hpp>
#elif defined(NWPW_CUDA)
#include <cuda/cuda_runtime.h>
#elif defined(NWPW_HIP)
#include <hip/hip_runtime.h>
#endif

namespace Nwpw {

#ifdef NWPW_SYCL
  struct gpuStream_t {
    cl::sycl::queue* queue = nullptr;
    bool operator==(const gpuStream_t& rhs) noexcept { return queue == rhs.queue; }
  };
#elif defined(NWPW_HIP)
  using gpuStream_t = hipStream_t;
#elif defined(NWPW_CUDA)
  using gpuStream_t = cudaStream_t;
#endif

}

#endif

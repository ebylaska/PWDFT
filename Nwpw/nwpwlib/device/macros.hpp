#ifndef NWPW_MACROS_H
#define NWPW_MACROS_H

#include <iostream>
#include <stdexcept>
#include <string>
#include <cstdlib> // abort()

#if defined(NWPW_CUDA)
#include <cuda/cuda_runtime.h>
#elif defined(NWPW_HIP)
#include <hip/hip_runtime.h>
#elif defined(NWPW_SYCL)
#include <CL/sycl.hpp>
#endif

#ifndef NWPW_INLINE
#if NWPW_COMP_ICC
#define NWPW_INLINE __forceinline
#else
#define NWPW_INLINE inline
#endif
#endif

#define NWPW_THROW_X(X) throw X
#define NWPW_THROW throw
#define NWPW_TRY try
#define NWPW_CATCH(X) catch (X)

// convert a token to a string
#define NWPW_MAKESTRING2(a) #a
#define NWPW_MAKESTRING(a) NWPW_MAKESTRING2(a)

// bool copy_bool(bool b) { return b; }
// inline void assert_fail(const char *condition, const char *function, const char *file, int line)
// {
//   std::cerr << "assertion failed: " << condition << " in function " << function
//             << " at " << file << ":"
//             << line << std::endl;
//   abort();
// }

#if defined(NWPW_HIP)
#define NWPW_HIP_OR_CUDA_OR_SYCL(a,b,c) a
#elif defined(NWPW_CUDA)
#define NWPW_HIP_OR_CUDA_OR_SYCL(a,b,c) b
#elif defined(NWPW_SYCL)
#define NWPW_HIP_OR_CUDA_OR_SYCL(a,b,c) c
#else
#define NWPW_HIP_OR_CUDA_OR_SYCL(a,b,c) ((void)0);
#endif

// #define NWPW_ASSERT(x)                                                  \
//   do {                                                                  \
//     if(!copy_bool(x))                                                   \
//       assert_fail(NWPW_MAKESTRING(x), __PRETTY_FUNCTION__, __FILE__, __LINE__); \
//   } while(false)
// #endif

// Define a macro for catching SYCL exceptions
#ifdef NWPW_SYCL
#define NWPW_SYCL_TRY_CATCH(X)                                          \
  do {                                                                  \
    NWPW_TRY {X;}                                                       \
    NWPW_CATCH(const cl::sycl::exception& e) {                          \
      NWPW_THROW_X(std::runtime_error("SYCL exception at " +            \
                                      std::string(__FILE__) + ":" +     \
                                      std::to_string(__LINE__) + "\n" + \
                                      e.what()));                       \
    }                                                                   \
  } while (false)
#endif //NWPW_SYCL

#ifdef NWPW_CUDA
#define NWPW_CUDA_SAFE_CALL(call) {                                     \
    cudaError_t nwpw_i_err = call;                                      \
    if (cudaSuccess != nwpw_i_err) {                                    \
      std::string errStr(std::string("CUDA error ") + std::to_string(nwpw_i_err) \
                         + std::string(" in file ") + __FILE__          \
                         + " line " + std::to_string(__LINE__)          \
                         + ": " + cudaGetErrorString(nwpw_i_err));      \
      NWPW_THROW_X(std::runtime_error(errStr));                         \
    }}
#endif // NWPW_CUDA

#ifdef NWPW_HIP
#define NWPW_HIP_SAFE_CALL(call) {                                      \
    hipError_t nwpw_i_err = call;                                       \
    if (hipSuccess != nwpw_i_err) {                                     \
      std::string errStr(std::string("HIP error in file ") + __FILE__   \
                         + " line " + std::to_string(__LINE__)          \
                         + " " + hipGetErrorString(nwpw_i_err));        \
      NWPW_THROW_X(std::runtime_error(errStr));                         \
    }}
#endif // NWPW_HIP



#endif

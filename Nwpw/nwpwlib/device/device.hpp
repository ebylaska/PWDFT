#ifndef NWPW_GPU_DEVICE_H_
#define NWPW_GPU_DEVICE_H_

#include <cstdlib>
#include <memory>
#include <array>
#include "device_types.hpp"
#include "macros.hpp"

#include "Parallel.hpp"

namespace Nwpw {

#ifdef NWPW_SYCL
    struct gpuDeviceProp_t {
        std::string device_name;
        std::size_t totalGlobalMem;
        std::size_t sharedMemPerBlock;
        int multiProcessorCount;
        int maxThreadsPerMultiProcessor;
        int maxThreadsPerBlock;
        int maxThreadsDim[3];
        int maxGridSize[3]; // NOT yet supported in SYCL/DPC++
        int warpSize;
        cl::sycl::info::local_mem_type sharedMemType;
        long maxMemAllocSize; // oneAPI only
        int managedMemory;
        int concurrentManagedAccess;
        int maxParameterSize;
    };
#elif defined(NWPW_HIP)
  using gpuDeviceProp_t = hipDeviceProp_t;
#elif defined(NWPW_CUDA)
  using gpuDeviceProp_t = cudaDeviceProp;
#endif

  class Device {

  public:
    enum DEVICETYPE {
      HOST,
      DEVICE,
      UNIFIED,
      CPU
    };

    Device(Parallel* myParallel);
    ~Device();

    template <typename T>
    T* allocate(std::size_t N, DEVICETYPE devType);
    void deallocate();

    static gpuStream_t gpuStream() noexcept { return gpu_stream; }
    static gpuStream_t nullStream() noexcept { return gpu_default_stream; }
#ifdef NWPW_SYCL
    static cl::sycl::queue& nullQueue() noexcept { return *(gpu_default_stream.queue); }
    static cl::sycl::queue& syclQueue() noexcept { return *(gpu_stream.queue); }
    static cl::sycl::queue& syclQueue(int i) noexcept { return *(gpu_streams[i].queue); }
    static bool onNullStream() noexcept { return gpu_stream == gpu_default_stream; }
    static bool onNullStream(gpuStream_t stream) noexcept { return stream == gpu_default_stream; }
#endif

    /* ABB: 08/17/2020: STILL EXPERIMENTAL
       static int numGpuStreams() noexcept { return max_gpu_streams; }
       static void setStreamIndex(const int idx) noexcept;
       static void resetStreamIndex() noexcept { setStreamIndex(-1); }

       static gpuStream_t setStream(gpuStream_t s) noexcept;
       static gpuStream_t resetStream() noexcept;
    */
    static int deviceId() noexcept;

    static void deviceSynchronize() noexcept;
    static void streamSynchronize() noexcept;
#ifdef NWPW_SYCL
    static void nonNullSyclQueueSynchronize() noexcept;
#endif

#if defined(NWPW_CUDA) || defined(NWPW_HIP)
    static std::size_t totalGlobalMem() noexcept      { return device_prop.totalGlobalMem; }
    static std::size_t sharedMemPerBlock() noexcept   { return device_prop.sharedMemPerBlock; }
    static int numMultiProcessors() noexcept          { return device_prop.multiProcessorCount; }
    static int maxThreadsPerMultiProcessor() noexcept { return device_prop.maxThreadsPerMultiProcessor; }
    static int maxThreadsPerBlock() noexcept          { return device_prop.maxThreadsPerBlock; }
    static int maxThreadsPerBlock(int dir) noexcept   { return device_prop.maxThreadsDim[dir]; }
    static int maxBlocksPerGrid(int dir) noexcept     { return device_prop.maxGridSize[dir]; }
    static std::string deviceName() noexcept          { return std::string(device_prop.device_name); }

    static int maxBlocksPerLaunch () noexcept { return max_blocks_per_launch; }
#endif

    static constexpr int warp_size = NWPW_HIP_OR_CUDA_OR_SYCL(64,32,16);

#ifdef NWPW_SYCL
    static long maxMemAllocSize() noexcept { return device_prop.maxMemAllocSize; }
    static cl::sycl::context& syclContext() { return *sycl_context; }
    static cl::sycl::device& syclDevice() { return *sycl_device; }
#endif

  private:

    static void initialize_gpu();

#ifdef NWPW_SYCL
    static cl::sycl::vector_class<cl::sycl::device> get_sycl_supported_devices();
#endif

    static int device_id;
    int m_num_devices;

#if defined(NWPW_CUDA) || defined(NWPW_HIP)
    static constexpr int max_gpu_streams = 4;
#elif defined(NWPW_SYCL)
    // Equivalent to "single dependent stream". Fits best with math this is used in ("x/max_streams").
    static constexpr int max_gpu_streams = 1;
#endif

    //static dim3 numThreadsMin;

    static std::array<gpuStream_t, max_gpu_streams> gpu_streams;
    static gpuStream_t gpu_default_stream;
    static gpuStream_t gpu_stream;
    static gpuDeviceProp_t device_prop;
    static int max_blocks_per_launch;

#ifdef NWPW_SYCL
    static std::unique_ptr<cl::sycl::context> sycl_context;
    static std::unique_ptr<cl::sycl::device> sycl_device;
#endif
  }; // class Device

#if defined(NWPW_CUDA) || defined(NWPW_HIP)
  NWPW_INLINE gpuStream_t
  gpuStream() noexcept {
    return Device::gpuStream();
  }

  NWPW_INLINE gpuStream_t
  nullStream() noexcept {
    return Device::nullStream();
  }
#endif

  /*
    NWPW_INLINE int
    numGpuStreams() noexcept {
    return Device::numGpuStreams();
    }
  */

  template <typename T>
  T*
  Device::allocate(std::size_t N, DEVICETYPE devType)
  {
    T* p = nullptr;

#ifdef NWLW_SYCL
    switch(devType) {
        case HOST: {
            p = cl::sycl::malloc_host(N * sizeof(T), syclDevice(), syclContext());
            break;
        }
        case DEVICE: {
            p = cl::sycl::malloc_device(N * sizeof(T), syclDevice(), syclContext());
            break;
        }
        case UNIFIED: {
            p = cl::sycl::malloc_shared(N * sizeof(T), syclDevice(), syclContext());
            break;
        }
        case CPU: {
            p = new T[N];
            break;
        }
        default: {
            throw std::runtime_error(std::string("UNKNOWN deviceType for SYCL Memory allocation! \n"));
            break;
        }
    }
#elif defined(NWPW_HIP)
    switch(devType) {
        case HOST: {
            NWPW_HIP_SAFE_CALL(hipHostMalloc(&p, N * sizeof(T), hipHostMallocMapped));
            break;
        }
        case DEVICE: {
            NWPW_HIP_SAFE_CALL(hipMalloc(&p, N * sizeof(T)));
            break;
        }
        case UNIFIED: {
            NWPW_HIP_SAFE_CALL(hipMallocManaged(&p, N * sizeof(T)));
            break;
        }
        case CPU: {
            p = new T[N];
            break;
        }
        default: {
            throw std::runtime_error(std::string("UNKNOWN deviceType for HIP Memory allocation! \n"));
            break;
        }
    }
#elif defined(NWPW_CUDA)
    switch(devType) {
        case HOST: {
            NWPW_CUDA_SAFE_CALL(cudaHostAlloc(&p, N * sizeof(T), cudaHostAllocMapped));
            break;
        }
        case DEVICE: {
            NWPW_CUDA_SAFE_CALL(cudaMalloc(&p, N * sizeof(T)));
            break;
        }
        case UNIFIED: {
            NWPW_CUDA_SAFE_CALL(cudaMallocManaged(&p, N * sizeof(T)));
            break;
        }
        case CPU: {
            p = new T[N];
            break;
        }
        default: {
            throw std::runtime_error(std::string("UNKNOWN deviceType for CUDA Memory allocation! \n"));
            break;
        }
    }
#endif

    if(p == nullptr)
        throw std::bad_alloc{};
    return p;
  }

  NWPW_INLINE void
  deviceSynchronize() noexcept {
    Device::deviceSynchronize();
  }

  NWPW_INLINE void
  streamSynchronize () noexcept {
    Device::streamSynchronize();
  }

  NWPW_INLINE void
  memcpyHostToDevice(void* p_d, const void* p_h, const std::size_t sz) noexcept {
#ifdef NWPW_SYCL
      auto& q = Device::nullQueue();
      q.submit([&] (sycl::handler& h) { h.memcpy(p_d, p_h, sz); });
      try {
          q.wait_and_throw();
      } catch (sycl::exception const& ex) {
          NWPW_THROW_X(std::runtime_error(std::string("SYCL htod_memcpy: ")+ex.what()+"!!!!!"));
      }
#elif defined(NWPW_CUDA)
      NWPW_CUDA_SAFE_CALL(cudaMemcpy(p_d, p_h, sz, cudaMemcpyHostToDevice));
#elif defined(NWPW_HIP)
      NWPW_HIP_SAFE_CALL(hipMemcpy(p_d, p_h, sz, hipMemcpyHostToDevice));
#endif
  }

  NWPW_INLINE void
  memcpyDeviceToHost(void* p_h, const void* p_d, const std::size_t sz) noexcept {
#ifdef NWPW_SYCL
      auto& q = Device::nullQueue();
      q.submit([&] (sycl::handler& h) { h.memcpy(p_h, p_d, sz); });
      try {
          q.wait_and_throw();
      } catch (sycl::exception const& ex) {
          NWPW_THROW_X(std::runtime_error(std::string("SYCL dtoh_memcpy: ")+ex.what()+"!!!!!"));
      }
#elif defined(NWPW_CUDA)
      NWPW_CUDA_SAFE_CALL(cudaMemcpy(p_h, p_d, sz, cudaMemcpyDeviceToHost));
#elif defined(NWPW_HIP)
      NWPW_HIP_SAFE_CALL(hipMemcpy(p_h, p_d, sz, hipMemcpyDeviceToHost));
#endif
  }

  NWPW_INLINE void
  memcpyDeviceToDevice(void* p_d_dst, const void* p_d_src, const std::size_t sz) noexcept {
#ifdef NWPW_SYCL
      auto& q = Device::nullQueue();
      q.submit([&] (sycl::handler& h) { h.memcpy(p_d_dst, p_d_src, sz); });
      try {
          q.wait_and_throw();
      } catch (sycl::exception const& ex) {
          NWPW_THROW_X(std::runtime_error(std::string("SYCL dtod_memcpy: ")+ex.what()+"!!!!!"));
      }
#elif defined(NWPW_CUDA)
      NWPW_CUDA_SAFE_CALL(cudaMemcpy(p_d_dst, p_d_src, sz, cudaMemcpyDeviceToDevice));
#elif defined(NWPW_HIP)
      NWPW_HIP_SAFE_CALL(hipMemcpy(p_d_dst, p_d_src, sz, hipMemcpyDeviceToDevice));
#endif
  }

  NWPW_INLINE void
  memcpyHostToDeviceAsync(void* p_d, const void* p_h, const std::size_t sz) noexcept {
#ifdef NWPW_SYCL
      auto& q = Device::syclQueue();
      q.submit([&] (sycl::handler& h) { h.memcpy(p_d, p_h, sz); });
#elif defined(NWPW_CUDA)
      NWPW_CUDA_SAFE_CALL(cudaMemcpyAsync(p_d, p_h, sz, cudaMemcpyHostToDevice, gpuStream()));
#elif defined(NWPW_HIP)
      NWPW_HIP_SAFE_CALL(hipMemcpyAsync(p_d, p_h, sz, hipMemcpyHostToDevice, gpuStream()));
#endif
  }

  NWPW_INLINE void
  memcpyDeviceToHostAsync(void* p_h, const void* p_d, const std::size_t sz) noexcept {
#ifdef NWPW_SYCL
      auto& q = Device::syclQueue();
      q.submit([&] (sycl::handler& h) { h.memcpy(p_h, p_d, sz); });
#elif defined(NWPW_CUDA)
      NWPW_CUDA_SAFE_CALL(cudaMemcpyAsync(p_h, p_d, sz, cudaMemcpyDeviceToHost, gpuStream()));
#elif defined(NWPW_HIP)
      NWPW_HIP_SAFE_CALL(hipMemcpyAsync(p_h, p_d, sz, hipMemcpyDeviceToHost, gpuStream()));
#endif
  }

  NWPW_INLINE void
  memcpyDeviceToDeviceAsync(void* p_d_dst, const void* p_d_src, const std::size_t sz) noexcept {
#ifdef NWPW_SYCL
      auto& q = Device::syclQueue();
      q.submit([&] (sycl::handler& h) { h.memcpy(p_d_dst, p_d_src, sz); });
#elif defined(NWPW_CUDA)
      NWPW_CUDA_SAFE_CALL(cudaMemcpyAsync(p_d_dst, p_d_src, sz, cudaMemcpyDeviceToDevice, gpuStream()));
#elif defined(NWPW_HIP)
      NWPW_HIP_SAFE_CALL(hipMemcpyAsync(p_d_dst, p_d_src, sz, cudaMemcpyDeviceToDevice, gpuStream()));
#endif
  }

#ifdef NWPW_SYCL
  NWPW_INLINE bool
  onNullStream() {
      return Device::onNullStream();
  }

  NWPW_INLINE bool
  onNullStream(gpuStream_t stream) {
      return Device::onNullStream(stream);
  }
#endif

} // Nwpw

#endif // NWPW_GPU_DEVICE_H_

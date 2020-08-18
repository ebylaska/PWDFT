#include <iostream>
#include <algorithm>
#include <string>
#include <unordered_set>
#include <exception>
#include <stdexcept>

#include "device.hpp"

#define NWPW_SYCL_USE_DEFAULT_SELECTOR (0)
#define NWPW_AURORA_GPUS_PER_SOCKET (3)

namespace Nwpw {

#ifdef NWPW_SYCL
    int Device::device_id = 0;
    constexpr int Device::max_gpu_streams;
    int Device::max_blocks_per_launch = 640;
    std::unique_ptr<cl::sycl::context> Device::sycl_context;
    std::unique_ptr<cl::sycl::device> Device::sycl_device;

    std::array<gpuStream_t,Device::max_gpu_streams> Device::gpu_streams;
    gpuStream_t                                     Device::gpu_default_stream;
    gpuStream_t                                     Device::gpu_stream;
    gpuDeviceProp_t                                 Device::device_prop;
#endif

#ifdef NWPW_SYCL
  auto nwpw_sycl_error_handler = [] (cl::sycl::exception_list exceptions) {
    for (std::exception_ptr const& e : exceptions) {
      try {
        std::rethrow_exception(e);
      } catch (cl::sycl::exception const& ex) {
        NWPW_THROW_X(std::runtime_error(std::string("Async SYCL exception: ") + ex.what() + "!!!!!"));
      }
    }
  };
#endif // NWPW_SYCL

  //_____________________________________________________________________________
  //

  Device::Device(Parallel* myParallel)
  {
    // Count the number of GPU devices.
    int numDevices = 0;
#ifdef NWPW_SYCL
    {
      auto sycl_devices_list = get_sycl_supported_devices();
      numDevices = sycl_devices_list.size();
      m_num_devices = numDevices;
      if (numDevices <= 0) {
        std::runtime_error(std::string("No DPCPP GPU device found!\n"));
      } else if (numDevices > 1) {
        std::runtime_error(std::string("DPCPP TODO: more than one device not supported yet!\n"));
      }
    }
#elif defined(NWPW_CUDA)
    NWPW_CUDA_SAFE_CALL(cudaGetDeviceCount(&numDevices));
    m_num_devices = numDevices;
    if (numDevices <= 0) {
      std::runtime_error(std::string("No CUDA GPU device found!\n"));
    }
#elif defined(NWPW_HIP)
    NWPW_HIP_SAFE_CALL (hipGetDeviceCount(&numDevices));
    m_num_devices = numDevices;
    if (numDevices <= 0) {
      std::runtime_error(std::string("No HIP GPU device found!\n"));
    }
#endif

    // Now, assign ranks to GPUs. If we only have one GPU,
    // or only one MPI rank, this is easy. Otherwise, we
    // need to do a little more work.

    if (myParallel->np() == 1) {
      device_id = 0;
    }
    else if (m_num_devices == 1) {
      device_id = 0;
    }
    else {
      // ifdef the following against MPI so it compiles, but note
      // that we can only get here if using more than one processor,
      // which requires MPI.

      // Create a communicator out of only the ranks sharing GPUs.
      // The default assumption is that this is all the ranks on the
      // same node, and to get that we'll use the MPI-3.0 split that
      // looks for shared memory communicators (and we'll error out
      // if that standard is unsupported).

      int MPIVersion, MPIsubVersion;
      MPI_Get_version(&MPIVersion, &MPIsubVersion);
      // abb: check here if it is being printed by all processes ??
      if (MPIVersion < 3)
        std::runtime_error(std::string("When using GPUs with MPI, if multiple devices are visible to each rank, MPI-3.0 must be supported!\n"));

      // However, it's possible that the ranks sharing GPUs will be
      // confined to a single socket rather than a full node. Indeed,
      // this is often the optimal configuration; for example, on Summit,
      // a good configuration using jsrun is one resource set per
      // socket (two per node), with three GPUs per resource set.
      // To deal with this where we can, we'll take advantage of OpenMPI's
      // specialized split by socket. However, we only want to do this
      // if in fact our resource set is confined to the socket.
      // To make this determination we need to have system information,
      // which is provided by the build system for the systems
      // we know about. The simple heuristic we'll use to determine
      // this is if the number of visible devices is smaller than
      // the known number of GPUs per socket.

      MPI_Comm local_comm;

      int split_type;
      split_type = MPI_COMM_TYPE_SHARED;

      // We have no preference on how ranks get ordered within this communicator.
      int key = 0;
      MPI_Comm_split_type(myParallel->communicator(), split_type, key, MPI_INFO_NULL, &local_comm);

      // Get rank within the local communicator, and number of ranks.
      int n_procs;
      MPI_Comm_size(local_comm, &n_procs);
      int my_rank;
      MPI_Comm_rank(local_comm, &my_rank);
      // Free the local communicator.
      MPI_Comm_free(&local_comm);

      // For each rank that shares a GPU, use round-robin assignment to assign MPI ranks to GPUs.
      // We will arbitrarily assign ranks to GPUs, assuming that socket awareness has already
      // been handled.
      device_id = my_rank % numDevices;

      // // If we detect more ranks than visible GPUs, warn the user that this will fail in the
      // // case where the devices are set to exclusive process mode and MPS is not enabled.
      if (n_procs > numDevices) {
        std::cout << "NWPW Warning: Mapping more than one rank per GPU. This will fail if the GPUs are in exclusive process mode\n"
                  << "              and MPS is not enabled. \n\n";
      }
    }

    // Set the gpu-device here!
#ifdef NWPW_CUDA

    int can_access = 0;
    for (int i = 0; i < numDevices; i++) {
      NWPW_CUDA_SAFE_CALL( cudaSetDevice(i) );
      NWPW_CUDA_SAFE_CALL( cudaSetDeviceFlags(cudaDeviceMapHost) );

      for (int j = 0; j < numDevices; j++) {
        if (i != j) {
          NWPW_CUDA_SAFE_CALL( cudaDeviceCanAccessPeer(&can_access, i, j) );
          if (can_access) {
            printf("CUDA GPU device #%d can access GPU device #%d\n", i, j);
            NWPW_CUDA_SAFE_CALL( cudaDeviceEnablePeerAccess(j, 0) );
          } else {
            //std::runtim_error(std::string("ERROR\n GPU device #%d cannot access GPU device #%d\n.  NWPW is not yet configured to work with multiple GPUs in different NUMA regions.  For now, use the environment variable CUDA_VISIBLE_DEVICES and don't list GPU device #%d\n." , i, j, j));
            std::runtime_error(std::string("NWPW ERROR** CUDA GPUs in multiple NUMA regions are currently unsupported.", __FILE__, __LINE__));
          }
        }
      }
    }

#elif defined(NWPW_HIP)

    int can_access = 0;
    for (int i = 0; i < numDevices; i++) {
      NWPW_HIP_SAFE_CALL( hipSetDevice(i) );
      NWPW_HIP_SAFE_CALL( hipSetDeviceFlags(hipDeviceMapHost) );

      for (int j = 0; j < numDevices; j++) {
        if (i != j) {
          NWPW_HIP_SAFE_CALL( hipDeviceCanAccessPeer(&can_access, i, j) );
          if (can_access) {
            printf("GOOD\n HIP GPU device #%d can access GPU device #%d\n", i, j);
            NWPW_HIP_SAFE_CALL( hipDeviceEnablePeerAccess(j, 0) );
          } else {
            //std::runtim_error(std::string("ERROR\n GPU device #%d cannot access GPU device #%d\n.  NWPW is not yet configured to work with multiple GPUs in different NUMA regions.  For now, use the environment variable HIP_VISIBLE_DEVICES and don't list GPU device #%d\n." , i, j, j));
            std::runtime_error(std::string("NWPW ERROR** HIP GPUs in multiple NUMA regions are currently unsupported.", __FILE__, __LINE__));
          }
        }
      }
    }

    //NWPW_HIP_SAFE_CALL(hipSetDevice(device_id));
    //NWPW_HIP_SAFE_CALL(hipSetDeviceFlags(hipDeviceMapHost));
#endif

    initialize_gpu();
  }

  //_____________________________________________________________________________
  //

  Device::~Device()
  {
#ifdef NWPW_CUDA
    for (int i = 0; i < max_gpu_streams; i++) {
      NWPW_CUDA_SAFE_CALL(cudaStreamDestroy(gpu_streams[i]));
    }
    NWPW_CUDA_SAFE_CALL(cudaDeviceReset());
#elif defined(NWPW_HIP)
    for (int i = 0; i < max_gpu_streams; i++) {
      NWPW_HIP_SAFE_CALL(hipStreamDestroy(gpu_streams[i]));
    }
    NWPW_HIP_SAFE_CALL(hipDeviceReset());
#endif

#ifdef NWPW_SYCL
    sycl_context.reset();
    sycl_device.reset();
    for (auto& s : gpu_streams) {
      delete s.queue;
      s.queue = nullptr;
    }
    gpu_stream.queue = nullptr;
    delete gpu_default_stream.queue;
    gpu_default_stream.queue = nullptr;
#endif
  }

  //_____________________________________________________________________________
  //

  void
  Device::initialize_gpu()
  {
#ifdef NWPW_SYCL
    { // create device, context and queues
      cl::sycl::gpu_selector device_selector;
      sycl_device.reset(new cl::sycl::device(device_selector));
      sycl_context.reset(new cl::sycl::context(*sycl_device, nwpw_sycl_error_handler));
      gpu_default_stream.queue = new cl::sycl::queue(*sycl_context, device_selector);
      for (int i = 0; i < max_gpu_streams; ++i) {
        gpu_streams[i].queue = new cl::sycl::queue(*sycl_context, device_selector);
      }
    }

    { // device property
      auto const& d = *sycl_device;
      device_prop.device_name         = d.get_info<cl::sycl::info::device::name>();
      device_prop.totalGlobalMem      = d.get_info<cl::sycl::info::device::global_mem_size>();
      device_prop.sharedMemType       = d.get_info<cl::sycl::info::device::local_mem_type>();
      device_prop.sharedMemPerBlock   = d.get_info<cl::sycl::info::device::local_mem_size>();
      device_prop.multiProcessorCount = d.get_info<cl::sycl::info::device::max_compute_units>();
      device_prop.maxThreadsPerMultiProcessor = -1; // SYCL todo: d.get_info<cl::sycl::info::device::max_work_items_per_compute_unit>();
      device_prop.maxThreadsPerBlock = d.get_info<cl::sycl::info::device::max_work_group_size>();
      auto mtd = d.get_info<cl::sycl::info::device::max_work_item_sizes>();
      device_prop.maxThreadsDim[0] = mtd[0];
      device_prop.maxThreadsDim[1] = mtd[1];
      device_prop.maxThreadsDim[2] = mtd[2];

      device_prop.maxGridSize[0] = -1; // xxxxx SYCL todo: unknown
      device_prop.maxGridSize[0] = -1; // unknown
      device_prop.maxGridSize[0] = -1; // unknown
      device_prop.warpSize = warp_size;

      //auto sgss = d.get_info<cl::sycl::info::device::sub_group_sizes>();
      device_prop.maxMemAllocSize         = d.get_info<cl::sycl::info::device::max_mem_alloc_size>();
      device_prop.managedMemory           = d.get_info<cl::sycl::info::device::host_unified_memory>();
      device_prop.concurrentManagedAccess = d.get_info<cl::sycl::info::device::usm_shared_allocations>();
      device_prop.maxParameterSize        = d.get_info<cl::sycl::info::device::max_parameter_size>();
      {
        std::cout << "Device Properties:\n"
                  << "  name               : " << device_prop.device_name << "\n"
                  << "  totalGlobalMem     : " << device_prop.totalGlobalMem << "\n"
                  << "  sharedMemPerBlock  : " << device_prop.sharedMemPerBlock << "\n"
                  << "  multiProcessorCount: " << device_prop.multiProcessorCount << "\n"
                  << "  maxThreadsPerBlock : " << device_prop.maxThreadsPerBlock << "\n"
                  << "  maxThreadsDim      : (" << device_prop.maxThreadsDim[0] << ", " << device_prop.maxThreadsDim[1] << ", " << device_prop.maxThreadsDim[2] << ")\n"
                  << "  warpSize           : " << device_prop.warpSize << "\n"
                  << "  maxMemAllocSize    : " << device_prop.maxMemAllocSize << "\n"
                  << std::endl;
      }
    }

#elif defined(NWPW_HIP)
    NWPW_HIP_SAFE_CALL(hipGetDeviceProperties(&device_prop, device_id));
    gpu_default_stream = 0;
    for (int i = 0; i < max_gpu_streams; i++) {
      NWPW_HIP_SAFE_CALL(hipStreamCreate(&gpu_streams[i]));
    }
#elif defined(NWPW_CUDA)
    NWPW_CUDA_SAFE_CALL(cudaGetDeviceProperties(&device_prop, device_id));
    gpu_default_stream = 0;
    for (int i = 0; i < max_gpu_streams; i++) {
      NWPW_CUDA_SAFE_CALL(cudaStreamCreate(&gpu_streams[i]));
    }
#endif

    gpu_stream = gpu_default_stream;

#ifdef NWPW_SYCL
    max_blocks_per_launch = 1000000; // SYCL todo
#else
    max_blocks_per_launch = numMultiProcessors() * maxThreadsPerMultiProcessor() / NWPW_GPU_MAX_THREADS;
#endif
  }

  //_____________________________________________________________________________
  //

#ifdef NWPW_SYCL
  auto
  Device::get_sycl_supported_devices() -> decltype(cl::sycl::device::get_devices())
  {
#ifdef NWPW_SYCL_USE_DEFAULT_SELECTOR
    return {cl::sycl::device(cl::sycl::default_selector())};
#else
    std::vector<cl::sycl::device> supported_devices;
    auto platform_list = cl::sycl::platform::get_platforms();
    for (const auto &platform : platform_list) {
      auto device_list = platform.get_devices();
      auto platform_name = platform.template get_info<cl::sycl::info::platform::name>();
      std::transform(platform_name.begin(), platform_name.end(), platform_name.begin(), ::tolower);

      for (const auto &device : device_list) {
        bool supported_condition = (device.is_gpu() && platform_name.find("Intel") != std::string::npos)
          if (supported_condition) {
            supported_devices.push_back(device);
          }
      }
    }
    return supported_devices;
#endif
  }
#endif //NWPW_SYCL

  //_____________________________________________________________________________
  //

  /* todo: DPCPP still needs some work
     void
     Device::deallocate(void* p)
     {
     #ifdef NWLW_SYCL
     switch(devType) {
     case HOST:
     case DEVICE:
     case UNIFIED: {
     p = cl::sycl::free(p, syclContext());
     break;
     }
     case CPU: {
     std::free(p);
     break;
     }
     }
     #elif defined(NWPW_HIP)
     switch(devType) {
     case HOST: {
     NWPW_HIP_SAFE_CALL(hipHostFree(p));
     break;
     }
     case DEVICE:
     case UNIFIED: {
     NWPW_HIP_SAFE_CALL(hipFree(p));
     break;
     }
     case CPU: {
     std::free(p);
     break;
     }
     }
     #elif defined(NWPW_CUDA)
     switch(devType) {
     case HOST: {
     NWPW_CUDA_SAFE_CALL(cudaHostFree(p));
     break;
     }
     case DEVICE:
     case UNIFIED: {
     NWPW_CUDA_SAFE_CALL(cudaFree(p));
     break;
     }
     case CPU: {
     std::free(p);
     break;
     }
     }
     #endif
     }
  */

  //_____________________________________________________________________________
  //

  int
  Device::deviceId() noexcept {
    return device_id;
  }

  /* ABB 04/20/2020: STILL EXPERIMENTAL
     void
     Device::setStreamIndex(const int idx) noexcept {
     if (idx < 0) {
     gpu_stream = gpu_default_stream;
     } else {
     gpu_stream = gpu_streams[idx % max_gpu_streams];
     }
     }

     gpuStream_t
     Device::resetStream() noexcept {
     gpuStream_t r = gpu_stream;
     gpu_stream = gpu_default_stream;
     return r;
     }

     gpuStream_t
     Device::setStream(gpuStream_t s) noexcept {
     gpuStream_t r = gpu_stream;
     gpu_stream = s;
     return r;
     }
  */

  //_____________________________________________________________________________
  //

  void
  Device::deviceSynchronize() noexcept {
#ifdef NWPW_SYCL
    nonNullSyclQueueSynchronize();
    try {
      gpu_default_stream.queue->wait_and_throw();
    } catch (sycl::exception const& ex) {
      NWPW_THROW_X(std::runtime_error(std::string("NWPW deviceSynchronize: ") + ex.what() + "!!!!!"));
    }
#elif defined(NWPW_CUDA)
    NWPW_CUDA_SAFE_CALL(cudaDeviceSynchronize());
#elif defined(NWPW_HIP)
    NWPW_HIP_SAFE_CALL(hipDeviceSynchronize());
#endif
  }

  //_____________________________________________________________________________
  //

  void
  Device::streamSynchronize() noexcept {
#ifdef NWPW_SYCL
    auto& q = syclQueue();
    try {
      q.wait_and_throw();
    } catch (sycl::exception const& ex) {
      NWPW_THROW_X(std::runtime_error(std::string("NWPW streamSynchronize: ") + ex.what() + "!!!!!"));
    }
#elif defined(NWPW_CUDA)
    NWPW_CUDA_SAFE_CALL(cudaStreamSynchronize(gpu_stream));
#elif defined(NWPW_HIP)
    NWPW_HIP_SAFE_CALL(hipStreamSynchronize(gpu_stream));
#endif
  }

  //_____________________________________________________________________________
  //

#ifdef NWPW_SYCL
  void
  Device::nonNullSyclQueueSynchronize() noexcept {
    for (auto const& s : gpu_streams) {
      try {
        s.queue->wait_and_throw();
      } catch (sycl::exception const& ex) {
        NWPW_THROW_X(std::runtime_error(std::string("NWPW nonNullSyclQueueSynchronize: ") + ex.what() + "!!!!!"));
      }
    }
  }
#endif

  //_____________________________________________________________________________
  //

  /*
    std::size_t
    Device::freeMemAvailable()
    {
    #ifdef NWPW_CUDA
    std::size_t f, t;
    NWPW_CUDA_SAFE_CALL(cudaMemGetInfo(&f,&t));
    f = device_prop.totalGlobalMem;
    return f;
    #elif defined(NWPW_HIP)
    std::size_t f, t;
    NWPW_HIP_SAFE_CALL(hipMemGetInfo(&f,&t));
    f = device_prop.totalGlobalMem;
    return f;
    #else   // SYCL (not yet supported)
    return 0;
    #endif
    }
  */

  //_____________________________________________________________________________
  //

}// Nwpw

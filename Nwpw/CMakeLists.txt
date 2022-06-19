cmake_minimum_required(VERSION 3.17)

# set the project name
project(PWDFT VERSION 1.0)
option( NWPW_SYCL "Enable SYCL Bindings" OFF )
option( NWPW_CUDA "Enable CUDA Bindings" OFF )
option( NWPW_HIP  "Enable HIP Bindings" OFF )
option( NWPW_OPENCL "Enable OpenCL Bindings" OFF )
option( NWPW_OPENMP "Enable OpenMP Bindings" OFF )
option( NWPW_SYCL_ENABLE_PROFILE "Enable SYCL Queue Profiling Bindings" OFF )

string(TIMESTAMP PWDFT_BUILD_TIMESTAMP "\"%a %b %d %H:%M:%S %Y\"")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable Fortran
enable_language(Fortran)

#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_DEBUG "-O3")
#set(CMAKE_CXX_FLAGS_RELEASE "-Ofast")
#set(CMAKE_Fortran_FLAGS "-Ofast")

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_Fortran_FLAGS "-O3")

#Configure Sycll
if (NWPW_SYCL)
  message("-- Using NWPW_SYCL")

  unset(CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)

  add_definitions(-DNWPW_SYCL=1)
  #SET(SYCL_CXX_COMPILE_FLAGS "-sycl-std=2020 -Wsycl-strict -fsycl-targets=spir64_gen -Xsycl-target-backend \'-device xehp'")
  SET(SYCL_CXX_COMPILE_FLAGS "-fsycl -sycl-std=2020 -Wsycl-strict")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SYCL_CXX_COMPILE_FLAGS}")

  include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindMKL.cmake)

  if( NWPW_SYCL_ENABLE_PROFILE )
    add_definitions(-DNWPW_SYCL_ENABLE_PROFILE=1)
  endif(NWPW_SYCL_ENABLE_PROFILE)
endif (NWPW_SYCL)

if (NWPW_CUDA)
   message("-- Using NWPW_CUDA")

   enable_language(CUDA)
   find_package( CUDAToolkit REQUIRED )
   set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_60,sm_70,sm_80")

   message("-- Using cublas, cufft: " ${CUDAToolkit_INCLUDE_DIRS})
   include_directories(${CUDAToolkit_INCLUDE_DIRS})

   add_definitions(-DNWPW_CUDA=1)
   set(MKL_FOUND True)
endif (NWPW_CUDA)


if (NWPW_HIP)
   message("-- Using NWPW_HIP")

   list(APPEND CMAKE_PREFIX_PATH ${ROCM_PATH}/hip ${ROCM_PATH})

   find_package(hip REQUIRED)
   find_package(rocblas REQUIRED)
   find_package(rocfft REQUIRED)
   set(GPU_TARGETS "gfx906,gfx908,gfx90a" CACHE STRING "GPU targets to compile for")

   set(ROCBLAS_INCLUDE_DIRS ${ROCM_PATH}/rocblas)
   set(ROCFFT_INCLUDE_DIRS ${ROCM_PATH}/rocfft)
   message("-- Using rocblas, rocfft: " ${ROCBLAS_INCLUDE_DIRS}, ${ROCFFT_INCLUDE_DIRS})
   include_directories(${ROCBLAS_INCLUDE_DIRS} ${ROCFFT_INCLUDE_DIRS})

   add_definitions(-DNWPW_HIP=1)
   set(MKL_FOUND True)
endif (NWPW_HIP)


set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")


#Configure OpenCl
if(NWPW_OPENCL)
find_package(OpenCL)
if(OpenCL_FOUND)
   message("-- Using OpenCL libs")
   add_definitions(-DNWPW_OPENCL=1)
   include_directories(${OpenCL_INCLUDE_DIRS})
   link_directories(${OpenCL_LIBRARY})
endif(OpenCL_FOUND)
endif(NWPW_OPENCL)


#Configure OpenMP
if(NWPW_OPENMP)
find_package(OpenMP)
if(OPENMP_FOUND)
   message("-- Using OpenMP")
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif(OPENMP_FOUND)
endif(NWPW_OPENMP)

#Configure MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})


if(MKL_FOUND)
   if(NWPW_SYCL)
      message("-- Using oneMKL libs")
      include_directories(${MKL_INCLUDE_DIRS})
   endif(NWPW_SYCL)
else()
   message("-- Using internal blas and lapack libs")
endif()


# configure a header file to pass the version number only
configure_file(NwpwConfig.h.in NwpwConfig.h)
configure_file("${PROJECT_SOURCE_DIR}/NwpwConfig.h.in" "${PROJECT_SOURCE_DIR}/nwpwlib/psp_library/NwpwLibrarypsConfig.hpp")

# create the nwpwlib and pspw libraries
add_subdirectory(nwpwlib)
add_subdirectory(pspw)


if(MAKE_LIBRARY)
   # create the pwdftlib shared library
   message("-- Creating pwdft shared library!")
   add_library(pwdft SHARED nwpw.cpp nwpwlib pspw)
   #add_library(pwdft STATIC nwpw.cpp)
   target_include_directories(pwdft PUBLIC "${PROJECT_BINARY_DIR}")
else()
   # create the pwdft executable
   message("-- Creating pwdft executable!")
   add_executable(pwdft nwpw.cpp)
   target_include_directories(pwdft PUBLIC "${PROJECT_BINARY_DIR}")
endif()


if(OpenCL_FOUND)
   target_link_libraries(pwdft pspw nwpwlib ${MPI_LIBRARIES} ${OpenCL_LIBRARY})
elseif(NWPW_CUDA)
   target_link_libraries(pwdft pspw nwpwlib ${MPI_LIBRARIES} CUDA::cudart)
elseif(NWPW_HIP)
   target_link_libraries(pwdft pspw nwpwlib ${MPI_LIBRARIES} hip::device)
else()
   target_link_libraries(pwdft pspw nwpwlib ${MPI_LIBRARIES})
endif()


if(MPI_COMPILE_FLAGS)
  set_target_properties(pwdft PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS} ${MKL_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(pwdft PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS} ${MKL_LINK_FLAGS}")
endif()

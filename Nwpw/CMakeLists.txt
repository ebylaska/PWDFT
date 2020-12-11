cmake_minimum_required(VERSION 3.10)

# set the project name
project(PWDFT VERSION 1.0)
option( NWPW_SYCL "Enable SYCL Bindings" OFF )
option( NWPW_OPENCL "Enable OpenCL Bindings" OFF )
option( NWPW_SYCL_ENABLE_PROFILE "Enable SYCL Queue Profiling Bindings" OFF )

string(TIMESTAMP PWDFT_BUILD_TIMESTAMP "\"%a %b %d %H:%M:%S %Y\"")

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindMKL.cmake)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable Fortran
enable_language(Fortran)

#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if( NWPW_SYCL )
  add_definitions(-DNWPW_SYCL=1)
  #SET(SYCL_CXX_COMPILE_FLAGS "-sycl-std=121 -fsycl-targets=spir64_gen-unknown-unknown-sycldevice -Xs \'-device skl\' ")
  SET(SYCL_CXX_COMPILE_FLAGS "-sycl-std=121 -fsycl-targets=spir64_gen-unknown-unknown-sycldevice -Xsycl-target-backend \'-device skl\'")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SYCL_CXX_COMPILE_FLAGS}")
  
  if( NWPW_SYCL_ENABLE_PROFILE )
    add_definitions(-DNWPW_SYCL_ENABLE_PROFILE=1)
  endif(NWPW_SYCL_ENABLE_PROFILE)
endif(NWPW_SYCL)

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")


#Configure OpenCl
if( NWPW_OPENCL )
find_package(OpenCL)
if(OpenCL_FOUND)
   message("-- Using OpenCL libs")
   add_definitions(-DNWPW_OPENCL=1)
   include_directories(${OpenCL_INCLUDE_DIRS})
   link_directories(${OpenCL_LIBRARY})
endif(OPENCL_FOUND)
endif(NWPW_OPENCL)

#Configure MPI
find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})


if(MKL_FOUND)
   include_directories(${MKL_INCLUDE_DIRS})
   message("-- Using MKL libs")
else()
   message("-- Using internal blas and lapack libs")
endif()

# configure a header file to pass the version number only
configure_file(NwpwConfig.h.in NwpwConfig.h)
configure_file("${PROJECT_SOURCE_DIR}/NwpwConfig.h.in" "${PROJECT_SOURCE_DIR}/nwpwlib/psp_library/NwpwLibrarypsConfig.hpp")

# create the nwpwlib and pspw libraries
add_subdirectory(nwpwlib)
add_subdirectory(pspw)

# create the pwdft executable
add_executable(pwdft nwpw.cpp)
target_include_directories(pwdft PUBLIC
                           "${PROJECT_BINARY_DIR}"
                           )
if(OpenCL_FOUND)
   target_link_libraries(pwdft pspw nwpwlib ${MPI_LIBRARIES} ${OpenCL_LIBRARY})
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

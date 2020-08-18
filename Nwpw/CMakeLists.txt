cmake_minimum_required(VERSION 3.10)

# set the project name
project(PWDFT VERSION 1.0)
add_definitions(-DNWPW_SYCL=1)

string(TIMESTAMP PWDFT_BUILD_TIMESTAMP "\"%a %b %d %H:%M:%S %Y\"")

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindMKL.cmake)


# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Enable Fortran
enable_language(Fortran)

#SET(SYCL_CXX_COMPILE_FLAGS "-fsycl")

#if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
#    set(bounds "-fbounds-check")
#endif()
#if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
#    set(dialect "-stand f08 -free -implicitnone")
#    set(bounds "-check bounds")
#endif()
#if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
#    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
#    set(bounds "-C")
#endif()

#SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${SYCL_CXX_COMPILE_FLAGS}")

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")

#Configure OpenMP
find_package(OpenMP)

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

if(OpenMP_CXX_FOUND)
   target_link_libraries(pwdft pspw nwpwlib ${MPI_LIBRARIES} OpenMP::OpenMP_CXX)
else(NOT OpenMP_CXX_FOUND)
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


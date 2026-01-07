
# files in the nwpwlib library
file(GLOB src_blas      blas/*.f)
file(GLOB src_lapack    lapack/*.f)
file(GLOB_RECURSE src_fftpack   fftpack/*.f)
file(GLOB_RECURSE src_fextra    fextra/*.f)
file(GLOB_RECURSE src_nwpwxc    nwpwxc/*.F)
file(GLOB_RECURSE src_Control   Control/*.hpp Control/*.cpp Control/*.h Control/*.c)
file(GLOB_RECURSE src_Parallel  Parallel/*.hpp Parallel/*.cpp)
file(GLOB_RECURSE src_io        io/*.hpp io/*.cpp)
file(GLOB_RECURSE src_ion       ion/*.hpp ion/*.cpp)
file(GLOB_RECURSE src_symmetry  symmetry/*.hpp symmetry/*.cpp)
file(GLOB_RECURSE src_brillouin brillouin/*.hpp brillouin/*.cpp)
file(GLOB_RECURSE src_efield    efield/*.hpp efield/*.cpp)
file(GLOB_RECURSE src_apc       apc/*.hpp apc/*.cpp)
file(GLOB_RECURSE src_dipole    dipole/*.hpp dipole/*.cpp)
#file(GLOB_RECURSE src_dispersion dispersion/*.hpp dispersion/*.cpp)
file(GLOB_RECURSE src_dplot     dplot/*.hpp dplot/*.cpp)
file(GLOB_RECURSE src_aimd_analysis aimd_analysis/*.hpp aimd_analysis/*.cpp)
file(GLOB_RECURSE src_thermostats   thermostats/*.hpp thermostats/*.cpp)
file(GLOB_RECURSE src_utilities     utilities/*.hpp utilities/*.cpp)
file(GLOB_RECURSE src_paw_utilities paw_utilities/*.hpp paw_utilities/*.cpp)
file(GLOB_RECURSE src_xcfunctions xcfunctions/*.hpp xcfunctions/*.cpp)
file(GLOB_RECURSE src_vdw       vdw/*.hpp vdw/*.cpp)
file(GLOB_RECURSE src_lattice   lattice/*.hpp lattice/*.cpp)
file(GLOB_RECURSE src_D3dB      D3dB/*.hpp D3dB/*.cpp)
file(GLOB_RECURSE src_C3dB      C3dB/*.hpp C3dB/*.cpp)
file(GLOB_RECURSE src_pspgen    pseudopotential/*.h  pseudopotential/*.c)
file(GLOB_RECURSE src_psplib    psp_library/*.hpp  psp_library/*.cpp)
file(GLOB_RECURSE src_parse     parse/*.hpp parse/*.cpp)
file(GLOB_RECURSE src_device    device/*.hpp device/*.cpp)

# spglib (space-group detection)
#file(GLOB src_spglib external/spglib/src/*.c)

#add_compile_options(-Wno-maybe-uninitialized)


# create the nwpwlib library
if(MKL_FOUND)
add_library(nwpwlib ${src_D3dB} ${src_C3dB} ${src_ion} ${src_symmetry} ${src_brillouin} ${src_efield} ${src_apc} ${src_dipole} ${src_dplot} ${src_aimd_analysis}  ${src_thermostats} ${src_xcfunctions}  ${src_vdw} ${src_lattice} ${src_Control} ${src_io} ${src_paw_utilities} ${src_utilities} ${src_Parallel} ${src_lapack} ${src_blas} ${src_fftpack} ${src_fextra} ${src_nwpwxc} ${src_pspgen} ${src_psplib} ${src_parse} ${src_device})
else()
   message("-- including interal blas routines contained in project")
add_library(nwpwlib ${src_D3dB} ${src_C3dB} ${src_ion} ${src_symmetry} ${src_brillouin} ${src_efield} ${src_apc} ${src_dipole} ${src_dplot} ${src_aimd_analysis} ${src_thermostats} ${src_xcfunctions} ${src_vdw} ${src_lattice} ${src_Control} ${src_io} ${src_paw_utilities} ${src_utilities} ${src_Parallel} ${src_lapack} ${src_blas} ${src_fftpack} ${src_fextra} ${src_nwpwxc}  ${src_pspgen} ${src_psplib} ${src_parse} ${src_device})
endif()


target_include_directories(nwpwlib
    PUBLIC
        ${PROJECT_BINARY_DIR}    # <-- gives access to NwpwConfig.h
        ${PROJECT_SOURCE_DIR}
)

#target_compile_options(nwpwlib PRIVATE -Wno-maybe-uninitialized)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    add_compile_options(-Wno-maybe-uninitialized)
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    add_compile_options(-Wno-uninitialized)
endif()



# adding target_include_directories to nwpwlib
add_subdirectory(parse)
add_subdirectory(psp_library)
add_subdirectory(pseudopotential)
add_subdirectory(pseudopotential/paw_atom)

#if(NOT MKL_FOUND)
add_subdirectory(blas)
#endif()

add_subdirectory(nwpwxc)
add_subdirectory(fftpack)
add_subdirectory(Control)
add_subdirectory(Parallel)
add_subdirectory(utilities)
add_subdirectory(paw_utilities)
add_subdirectory(io)
add_subdirectory(xcfunctions)
add_subdirectory(vdw)
add_subdirectory(lattice)
add_subdirectory(thermostats)
add_subdirectory(aimd_analysis)
add_subdirectory(dplot)
#add_subdirectory(dispersion)
add_subdirectory(dipole)
add_subdirectory(apc)
add_subdirectory(efield)
add_subdirectory(brillouin)
add_subdirectory(symmetry)
add_subdirectory(ion)
add_subdirectory(C3dB)
add_subdirectory(D3dB)
add_subdirectory(device)

add_subdirectory(external/spglib)
target_link_libraries(nwpwlib PUBLIC spglib)





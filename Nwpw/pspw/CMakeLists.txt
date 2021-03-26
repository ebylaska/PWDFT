
# files in the pspw library
file(GLOB_RECURSE src_minimizer minimizer/*.hpp minimizer/*.cpp)
file(GLOB_RECURSE src_cpsd      cpsd/*.hpp cpsd/*.cpp)
file(GLOB_RECURSE src_psi       lib/psi/*.hpp lib/psi/*.cpp)
file(GLOB_RECURSE src_electron  lib/electron/*.hpp lib/electron/*.cpp)
file(GLOB_RECURSE src_coulomb   lib/coulomb/*.hpp lib/coulomb/*.cpp)
file(GLOB_RECURSE src_kinetic   lib/kinetic/*.hpp lib/kinetic/*.cpp)
file(GLOB_RECURSE src_psp       lib/psp/*.hpp lib/psp/*.cpp)
file(GLOB_RECURSE src_exchange-correlation  lib/exchange-correlation/*.hpp lib/exchange-correlation/*.cpp)


# create to the pspw library
add_library(pspw ${src_minimizer} ${src_cpsd} ${src_psi} ${src_coulomb} ${src_kinetic} ${src_psp} ${src_exchange-correlation} )

# interface nwpwlib library to pspw library
target_link_libraries(pspw nwpwlib)

# add target_include_directories to the pspw library
add_subdirectory(minimizer)
add_subdirectory(cpsd)
add_subdirectory(lib/psi)
add_subdirectory(lib/electron)
add_subdirectory(lib/coulomb)
add_subdirectory(lib/kinetic)
add_subdirectory(lib/psp)
add_subdirectory(lib/exchange-correlation)


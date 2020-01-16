
# add the nwpwlib library
file(GLOB_RECURSE src_cpsd      cpsd/*.hpp cpsd/*.cpp)
file(GLOB_RECURSE src_psi       lib/psi/*.hpp lib/psi/*.cpp)
file(GLOB_RECURSE src_coulomb   lib/coulomb/*.hpp lib/coulomb/*.cpp)
file(GLOB_RECURSE src_kinetic   lib/kinetic/*.hpp lib/kinetic/*.cpp)
file(GLOB_RECURSE src_psp       lib/psp/*.hpp lib/psp/*.cpp)
file(GLOB_RECURSE src_exchange-correlation  lib/exchange-correlation/*.hpp lib/exchange-correlation/*.cpp)


# add to the nwpwlib library
add_library(pspw ${src_cpsd} ${src_psi} ${src_coulomb} ${src_kinetic} ${src_psp} ${src_exchange-correlation} )
target_link_libraries(pspw nwpwlib)

add_subdirectory(cpsd)
add_subdirectory(lib/psi)
add_subdirectory(lib/coulomb)
add_subdirectory(lib/kinetic)
add_subdirectory(lib/psp)
add_subdirectory(lib/exchange-correlation)


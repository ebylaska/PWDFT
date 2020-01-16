
# add the nwpwlib library
file(GLOB_RECURSE src_cpsd      cpsd/*.hpp cpsd/*.cpp)
file(GLOB_RECURSE src_coulomb   lib/coulomb/*.hpp lib/coulomb/*.cpp)


# add to the nwpwlib library
add_library(pspw ${src_cpsd} ${src_coulomb} )

add_subdirectory(cpsd)
add_subdirectory(lib/coulomb)



# files in the band library
file(GLOB_RECURSE src_band_cpsd cpsd/*.hpp cpsd/*.cpp)
file(GLOB_RECURSE src_cKinetic  lib/cKinetic/*.hpp lib/cKinetic/*.cpp)
file(GLOB_RECURSE src_cCoulomb  lib/cCoulomb/*.hpp lib/cCoulomb/*.cpp)
file(GLOB_RECURSE src_cExchange-Correlation  lib/cExchange-Correlation/*.hpp lib/cExchange-Correlation/*.cpp)
file(GLOB_RECURSE src_cpsp      lib/cpsp/*.hpp lib/cpsp/*.cpp)
file(GLOB_RECURSE src_cpsi      lib/cpsi/*.hpp lib/cpsi/*.cpp)

# create to the band library
add_library(band ${src_band_cpsd} ${src_cKinetic} ${src_cCoulomb} ${src_cExchange-Correlation} ${src_cpsp} ${src_cpsi} )

# interface nwpwlib library to pspw library
target_link_libraries(band nwpwlib)

# add target_include_directories to the pspw library
add_subdirectory(cpsd)
add_subdirectory(lib/cKinetic)
add_subdirectory(lib/cCoulomb)
add_subdirectory(lib/cExchange-Correlation)
add_subdirectory(lib/cpsp)
add_subdirectory(lib/cpsi)



# files in the band library
file(GLOB_RECURSE src_band_cpsd cpsd/*.hpp cpsd/*.cpp)
file(GLOB_RECURSE src_cpsp      lib/cpsp/*.hpp lib/cpsp/*.cpp)


# create to the band library
add_library(band ${src_band_cpsd} ${src_cpsp} )

# interface nwpwlib library to pspw library
target_link_libraries(band nwpwlib)

# add target_include_directories to the pspw library
add_subdirectory(cpsd)
add_subdirectory(lib/cpsp)


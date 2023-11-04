
# files in the band library
file(GLOB_RECURSE src_band_cpsd cpsd/*.hpp cpsd/*.cpp)


# create to the pspw library
add_library(band ${src_band_cpsd} )

# interface nwpwlib library to pspw library
target_link_libraries(band nwpwlib)

# add target_include_directories to the pspw library
add_subdirectory(cpsd)


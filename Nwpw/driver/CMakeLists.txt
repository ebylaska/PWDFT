# files in the driver library
file(GLOB_RECURSE src_optimize   optimize/*.hpp optimize/*.cpp)
#file(GLOB_RECURSE src_vibrations vibrations/*.hpp vibrations/*.cpp)


# create to the pspw library
#add_library(driver ${src_optimize} ${src_vibrations} )
add_library(driver ${src_optimize} )

# interface pspw and band libraries to driver library
target_link_libraries(driver pspw band)

# add target_include_directories to the pspw library
add_subdirectory(optimize)
#add_subdirectory(vibrations)

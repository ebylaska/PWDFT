
# files in the file library
file(GLOB_RECURSE src_file_generate   file_generate/*.hpp file_generate/*.cpp)

# create the file library
add_library(file ${src_file_generate} )

# interface nwpwlib library to file library
target_link_libraries(file nwpwlib)

# add target_include_directories to the pspw library
add_subdirectory(file_generate)


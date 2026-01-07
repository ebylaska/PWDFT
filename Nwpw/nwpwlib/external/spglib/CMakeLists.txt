cmake_minimum_required(VERSION 3.15)
project(spglib C)

message(STATUS ">>> ENTERING SPGLIB (vendored) <<<")

file(GLOB SPGLIB_SRC
    src/*.c
)

add_library(spglib STATIC ${SPGLIB_SRC})

target_include_directories(spglib
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)

set_target_properties(spglib PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

# Silence spglib warnings without polluting PWDFT
target_compile_options(spglib PRIVATE
    -Wno-unused-parameter
    -Wno-unused-function
    -Wno-maybe-uninitialized
)


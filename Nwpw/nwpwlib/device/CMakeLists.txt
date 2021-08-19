
target_include_directories(nwpwlib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(NWPW_CUDA)
  target_link_libraries( nwpwlib PUBLIC CUDA::cublas )
  target_link_libraries( nwpwlib PUBLIC CUDA::cufft )
endif()


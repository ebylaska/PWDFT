
target_include_directories(nwpwlib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(NWPW_CUDA)
  target_link_libraries( nwpwlib PUBLIC CUDA::cublas )
  target_link_libraries( nwpwlib PUBLIC CUDA::cufft )
  target_link_libraries( nwpwlib PUBLIC CUDA::cusolver )
endif()

if(NWPW_HIP)
  target_link_libraries( nwpwlib PUBLIC roc::rocblas )
  target_link_libraries( nwpwlib PUBLIC roc::rocfft )
  target_link_libraries( nwpwlib PUBLIC roc::rocsolver )
endif()

target_sources(hemelb-tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/CollisionTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/KernelTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/StreamerTests.cc
  )
add_subdirectory(iolets)

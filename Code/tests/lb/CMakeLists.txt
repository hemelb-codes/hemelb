target_sources(hemelb-tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/CollisionTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/IncompressibilityCheckerTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/KernelTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/LatticeTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/RheologyModelTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/StreamerTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/VirtualSiteIoletStreamerTests.cc
  )
add_subdirectory(iolets)

target_sources(hemelb-tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/GeometryReaderTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/LatticeDataTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/NeedsTests.cc
  )
add_subdirectory(neighbouring)

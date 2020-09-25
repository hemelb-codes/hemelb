target_sources(hemelb-tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/LabelledRequest.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/MpiTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/RecordingNet.cc
)
add_subdirectory(phased)

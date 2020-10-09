target_sources(hemelb-tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/MockIntercommunicand.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/MockIntercommunicator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/MockIntercommunicatorTests.cc
)

if (HEMELB_BUILD_MULTISCALE)
  add_subdirectory(mpwide)
endif()

if (HEMELB_BUILD_MULTISCALE)
  target_sources(hemelb-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/MPWideIntercommunicatorTests.cc
  )
endif()

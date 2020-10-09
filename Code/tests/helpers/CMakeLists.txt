target_sources(hemelb-tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/FolderTestFixture.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/FourCubeBasedTestFixture.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/FourCubeLatticeData.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/HasCommsTestFixture.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/MockCommunicator.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/MockCommsHelper.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/OneInOneOutSimConfig.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/RandomSourceTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/MockCommunicatorTests.cc
  )

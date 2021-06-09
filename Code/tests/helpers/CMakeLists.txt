# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
target_sources(hemelb-tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/FolderTestFixture.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/FourCubeBasedTestFixture.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/FourCubeLatticeData.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/HasCommsTestFixture.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/LatticeDataAccess.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/MockNetHelper.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/OneInOneOutSimConfig.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/RandomSourceTests.cc
  )

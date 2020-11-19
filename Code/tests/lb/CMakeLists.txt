# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
target_sources(hemelb-tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/BroadcastMocks.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/CollisionTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/IncompressibilityCheckerTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/KernelTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/LatticeTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/RheologyModelTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/StreamerTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/VirtualSiteIoletStreamerTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/GuoForcingTests.cc
  )
add_subdirectory(iolets)

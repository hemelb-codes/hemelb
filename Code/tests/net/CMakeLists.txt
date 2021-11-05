# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
target_sources(hemelb-tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/LabelledRequest.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/MpiTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/RecordingNet.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/NeighborCommTests.cc
)
add_subdirectory(phased)

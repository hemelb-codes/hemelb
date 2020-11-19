# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
target_sources(hemelb-tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/GeometryReaderTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/LatticeDataTests.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/NeedsTests.cc
  )
add_subdirectory(neighbouring)

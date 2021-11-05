# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
if (HEMELB_BUILD_MULTISCALE)
  target_sources(hemelb-tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/MPWideIntercommunicatorTests.cc
  )
endif()

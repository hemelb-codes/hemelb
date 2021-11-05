# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_reporting Reporter.cc Timers.cc Dict.cc)
target_link_libraries(hemelb_reporting PRIVATE CTemplate::CTemplate)

configure_file (
  "${PROJECT_SOURCE_DIR}/reporting/BuildInfo.h.in"
  "${PROJECT_BINARY_DIR}/reporting/BuildInfo.h"
  )

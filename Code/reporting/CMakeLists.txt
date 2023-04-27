# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_reporting OBJECT BuildInfo.cc Reporter.cc Timers.cc Dict.cc)
target_link_libraries(hemelb_reporting PRIVATE CTemplate::CTemplate)

# Ah CMake... We can only reference a variable in configure file and
# var names are case sensitive. We also can't use generator
# expressions at configure time because of multiconfig stuff.
string(TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TYPE)

configure_file (
  "${PROJECT_SOURCE_DIR}/reporting/BuildInfo.h.in"
  "${PROJECT_BINARY_DIR}/reporting/BuildInfo.h"
  )

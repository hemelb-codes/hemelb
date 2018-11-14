# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

add_library(hemelb_reporting STATIC
  Reporter.cc
  Timers.cc
  Dict.cc
)

hemelb_add_target_dependency_ctemplate(hemelb_reporting)

configure_file (
  "${PROJECT_SOURCE_DIR}/reporting/BuildInfo.h.in"
  "${PROJECT_BINARY_DIR}/reporting/BuildInfo.h"
  )

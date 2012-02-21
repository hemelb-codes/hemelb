add_library(hemelb_reporting Reporter.cc Timers.cc)
configure_file (
  "${PROJECT_SOURCE_DIR}/reporting/BuildInfo.h.in"
  "${PROJECT_BINARY_DIR}/reporting/BuildInfo.h"
  )
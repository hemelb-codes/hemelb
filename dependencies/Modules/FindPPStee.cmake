#
# Find the PPStee includes and libraries
#
# PPSTEE_INCLUDE_DIRS - where to find ppstee.hpp
# PPSTEE_LIBRARIES    - List of fully qualified libraries to link against.
# PPSTEE_FOUND        - Do not attempt to use if "no" or undefined.

find_path(PPSTEE_INCLUDE_DIR ppstee.hpp)

set(PPSTEE_INCLUDE_DIRS ${PPSTEE_INCLUDE_DIR})

find_library(PPSTEE_LIB_PPSTEE ppstee )

set( PPSTEE_LIBRARIES
  ${PPSTEE_LIB_PPSTEE}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PPStee  DEFAULT_MSG
                                  PPSTEE_INCLUDE_DIR
                                  PPSTEE_INCLUDE_DIRS
                                  PPSTEE_LIB_PPSTEE
                                  PPSTEE_LIBRARIES)

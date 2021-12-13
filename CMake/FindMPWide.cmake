# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

include_guard()

# - Find MPWide
# Find the native MPWide includes and library
#
#   MPWide_FOUND       - True if MPWide found.
#   MPWide_INCLUDE_DIR - where to find MPWide.h, etc.
#   MPWide_LIBRARIES   - List of libraries when using MPWide.
#

IF( MPWide_INCLUDE_DIR )
  # Already in cache, be silent
  SET( MPWide_FIND_QUIETLY TRUE )
ENDIF( MPWide_INCLUDE_DIR )

FIND_PATH( MPWide_INCLUDE_DIR "MPWide.h"
  PATH_SUFFIXES "MPWide" ".mpwide" )

FIND_LIBRARY( MPWide_LIBRARIES
  NAMES "MPW"
  PATH_SUFFIXES "MPWide" ".mpwide" )

# handle the QUIETLY and REQUIRED arguments and set MPWide_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE( "FindPackageHandleStandardArgs" )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( "MPWide" DEFAULT_MSG MPWide_INCLUDE_DIR MPWide_LIBRARIES )

MARK_AS_ADVANCED( MPWide_INCLUDE_DIR MPWide_LIBRARIES )

if (MPWide_FOUND)
  add_library(MPWide::MPWide INTERFACE IMPORTED)
  set_target_properties(MPWide::MPWide PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${MPWide_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${MPWide_LIBRARIES}"
    )
endif()

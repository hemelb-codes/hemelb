# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

include_guard()

# - Find TinyXML
# Find the native TinyXML includes and library
#
#   TINYXML_FOUND       - True if TinyXML found.
#   TINYXML_INCLUDE_DIR - where to find tinyxml.h, etc.
#   TINYXML_LIBRARIES   - List of libraries when using TinyXML.
#

IF( TINYXML_INCLUDE_DIR )
  # Already in cache, be silent
  SET( TinyXML_FIND_QUIETLY TRUE )
ENDIF( TINYXML_INCLUDE_DIR )

FIND_PATH( TINYXML_INCLUDE_DIR "tinyxml.h"
  PATH_SUFFIXES "tinyxml" )

FIND_LIBRARY( TINYXML_LIBRARIES
  NAMES "tinyxml"
  PATH_SUFFIXES "tinyxml" )

# handle the QUIETLY and REQUIRED arguments and set TINYXML_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE( "FindPackageHandleStandardArgs" )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( "TinyXML" DEFAULT_MSG TINYXML_INCLUDE_DIR TINYXML_LIBRARIES )

MARK_AS_ADVANCED( TINYXML_INCLUDE_DIR TINYXML_LIBRARIES )

if(TINYXML_FOUND)
  add_library(TinyXML::TinyXML INTERFACE IMPORTED)
  set_target_properties(TinyXML::TinyXML PROPERTIES
    # STL use is required by HemeLB and the build options
    INTERFACE_COMPILE_DEFINITIONS "TIXML_USE_STL"
    INTERFACE_INCLUDE_DIRECTORIES "${TINYXML_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${TINYXML_LIBRARIES}"
    )
endif()

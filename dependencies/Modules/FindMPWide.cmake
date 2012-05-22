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

# - Find CTEMPLATE
# Find the native CTEMPLATE includes and library
#
#   CTEMPLATE_FOUND       - True if CTEMPLATE found.
#   CTEMPLATE_INCLUDE_DIR - where to find CTEMPLATE.h, etc.
#   CTEMPLATE_LIBRARIES   - List of libraries when using CTEMPLATE.
#

IF( CTEMPLATE_INCLUDE_DIR )
    # Already in cache, be silent
    SET( CTEMPLATE_FIND_QUIETLY TRUE )
ENDIF( CTEMPLATE_INCLUDE_DIR )

FIND_PATH( CTEMPLATE_INCLUDE_DIR "ctemplate/template.h")

FIND_LIBRARY( CTEMPLATE_LIBRARIES
              NAMES "ctemplate"
              PATH_SUFFIXES "ctemplate" )

# handle the QUIETLY and REQUIRED arguments and set CTEMPLATE_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE( "FindPackageHandleStandardArgs" )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( "CTemplate" DEFAULT_MSG CTEMPLATE_INCLUDE_DIR CTEMPLATE_LIBRARIES )

MARK_AS_ADVANCED( CTEMPLATE_INCLUDE_DIR CTEMPLATE_LIBRARIES )
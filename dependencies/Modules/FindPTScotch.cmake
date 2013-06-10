#
# Find the PTScotch includes and libraries
#
# http://hdf.ncsa.uiuc.edu
#
# PTSCOTCH_INCLUDE_DIRS - where to find ptscotch.h
# PTSCOTCH_LIBRARIES    - List of fully qualified libraries to link against.
# PTSCOTCH_FOUND        - Do not attempt to use if "no" or undefined.

find_path(PTSCOTCH_INCLUDE_DIR ptscotch.h)

set(PTSCOTCH_INCLUDE_DIRS ${PTSCOTCH_INCLUDE_DIR})

find_library(PTSCOTCH_LIB_PTSCOTCH ptscotch )
find_library(PTSCOTCH_LIB_PTSCOTCHERR ptscotcherr )
find_library(PTSCOTCH_LIB_PTSCOTCHERREXIT ptscotcherrexit )
find_library(PTSCOTCH_LIB_PTSCOTCHPARMETIS ptscotchparmetis )
#find_library(PTSCOTCH_LIB_SCOTCH scotch )
#find_library(PTSCOTCH_LIB_SCOTCHERR scotcherr )
#find_library(PTSCOTCH_LIB_SCOTCHERREXIT scotcherrexit )
#find_library(PTSCOTCH_LIB_SCOTCHMETIS scotchmetis )

set( PTSCOTCH_LIBRARIES
  #${PTSCOTCH_LIB_SCOTCH}
  #${PTSCOTCH_LIB_SCOTCHERR}
  #${PTSCOTCH_LIB_SCOTCHERREXIT}
  #${PTSCOTCH_LIB_SCOTCHMETIS}
  ${PTSCOTCH_LIB_PTSCOTCH}
  ${PTSCOTCH_LIB_PTSCOTCHERR}
  ${PTSCOTCH_LIB_PTSCOTCHERREXIT}
  #${PTSCOTCH_LIB_PTSCOTCHPARMETIS}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PTScotch  DEFAULT_MSG
                                  PTSCOTCH_INCLUDE_DIR
                                  PTSCOTCH_INCLUDE_DIRS
                                  PTSCOTCH_LIB_PTSCOTCH
                                  PTSCOTCH_LIB_PTSCOTCHERR
                                  PTSCOTCH_LIB_PTSCOTCHERREXIT
                                  PTSCOTCH_LIBRARIES)

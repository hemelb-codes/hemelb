# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

include_guard()

# Find the METIS includes and libraries
#
# METIS is a library that implements a variety of algorithms for 
# partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of 
# sparse matrices. It can be found at:
# 	http://www-users.cs.umn.edu/~karypis/metis/parmetis/index.html
#
# PARMETIS_INCLUDE_DIR - where to find autopack.h
# PARMETIS_LIBRARIES   - List of fully qualified libraries to link against.
# PARMETIS_FOUND       - Do not attempt to use if "no" or undefined.

if(Parmetis_FIND_QUIETLY AND Parmetis_FIND_REQUIRED)
  FIND_PACKAGE(Metis QUIET REQUIRED)
elseif(Parmetis_FIND_QUIETLY)
  FIND_PACKAGE(Metis QUIET)
elseif(Parmetis_FIND_REQUIRED)
  FIND_PACKAGE(Metis REQUIRED)
else()
  FIND_PACKAGE(Metis)
endif()

FIND_PATH(PARMETIS_INCLUDE_DIR parmetis.h HINTS ${METIS_INCLUDE_DIR} ${PARMETIS_ROOT}/include ${PARMETIS_DIR}/include)
get_filename_component(METIS_LIB_DIR ${METIS_LIBRARY} DIRECTORY)
FIND_LIBRARY(PARMETIS_LIBRARY parmetis HINTS ${METIS_LIB_DIR} ${PARMETIS_ROOT}/lib ${PARMETIS_DIR}/lib)

list(APPEND PARMETIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR})
list(APPEND PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR})
list(REMOVE_DUPLICATES PARMETIS_INCLUDE_DIRS)
set(PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIRS} CACHE STRING "All include dirs needed by ParMETIS")

list(APPEND PARMETIS_LIBRARIES ${METIS_LIBRARY})
list(APPEND PARMETIS_LIBRARIES ${PARMETIS_LIBRARY})
list(REMOVE_DUPLICATES PARMETIS_LIBRARIES)
set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARIES} CACHE STRING "All libraries needed by ParMETIS")

include("FindPackageHandleStandardArgs")
FIND_PACKAGE_HANDLE_STANDARD_ARGS("Parmetis" DEFAULT_MSG PARMETIS_INCLUDE_DIR PARMETIS_LIBRARY)
MARK_AS_ADVANCED(PARMETIS_INCLUDE_DIR PARMETIS_LIBRARY PARMETIS_INCLUDE_DIRS PARMETIS_LIBRARIES)

if(PARMETIS_FOUND)
  add_library(Parmetis::Parmetis INTERFACE IMPORTED)
  set_target_properties(Parmetis::Parmetis PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${PARMETIS_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${PARMETIS_LIBRARY}"
  )
  target_link_libraries(Parmetis::Parmetis INTERFACE Metis::Metis)
endif()

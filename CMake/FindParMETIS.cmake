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
cmake_policy(SET CMP0074 NEW)

if(ParMETIS_FIND_QUIETLY AND ParMETIS_FIND_REQUIRED)
  find_package(METIS QUIET REQUIRED)
elseif(ParMETIS_FIND_QUIETLY)
  find_package(METIS QUIET)
elseif(ParMETIS_FIND_REQUIRED)
  find_package(METIS REQUIRED)
else()
  find_package(METIS)
endif()

find_package(PkgConfig)
pkg_check_modules(PC_PARMETIS QUIET parmetis)

find_path(ParMETIS_INCLUDE_DIR parmetis.h HINTS ${METIS_INCLUDE_DIR} ${ParMETIS_ROOT}/include ${ParMETIS_DIR}/include ${PC_PARMETIS_INCLUDE_DIRS})
get_filename_component(METIS_LIB_DIR ${METIS_LIBRARY} DIRECTORY)
find_library(ParMETIS_LIBRARY NAMES parmetis ${PC_PARMETIS_LIBRARIES} HINTS ${METIS_LIB_DIR} ${ParMETIS_ROOT}/lib ${ParMETIS_DIR}/lib ${PC_PARMETIS_LIBRARY_DIRS})

include("FindPackageHandleStandardArgs")
find_package_handle_standard_args("ParMETIS" DEFAULT_MSG ParMETIS_INCLUDE_DIR ParMETIS_LIBRARY)
mark_as_advanced(ParMETIS_INCLUDE_DIR ParMETIS_LIBRARY)

if(ParMETIS_FOUND)
  add_library(ParMETIS::ParMETIS INTERFACE IMPORTED)
  set_target_properties(ParMETIS::ParMETIS PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${ParMETIS_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${ParMETIS_LIBRARY}"
  )
  target_link_libraries(ParMETIS::ParMETIS INTERFACE METIS::METIS)
endif()

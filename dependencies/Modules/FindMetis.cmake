# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

include_guard()

# Find the PARMETIS includes and libraries
#
# ParMETIS is an MPI-based parallel library that implements a variety of algorithms for 
# partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of 
# sparse matrices. It can be found at:
# 	http://www-users.cs.umn.edu/~karypis/metis/parmetis/index.html
#
# You can help this script by setting METIS_ROOT, if you know it
#
# METIS_INCLUDE_DIR - where to find autopack.h
# METIS_LIBRARY     - List of fully qualified libraries to link against.
# METIS_FOUND       - Do not attempt to use if "no" or undefined.

FIND_PATH(METIS_INCLUDE_DIR metis.h HINTS ${METIS_ROOT}/include ${METIS_DIR}/include)
FIND_LIBRARY(METIS_LIBRARY metis HINTS ${METIS_ROOT}/lib ${METIS_DIR}/include)

include("FindPackageHandleStandardArgs")
FIND_PACKAGE_HANDLE_STANDARD_ARGS("Metis" DEFAULT_MSG METIS_INCLUDE_DIR METIS_LIBRARY)
MARK_AS_ADVANCED(METIS_INCLUDE_DIR METIS_LIBRARY)

if (METIS_FOUND)
  add_library(Metis::Metis INTERFACE IMPORTED)
  set_target_properties(Metis::Metis PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${METIS_LIBRARY}"
  )
endif()

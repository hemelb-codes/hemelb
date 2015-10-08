#
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

if(METIS_FOUND)
  FIND_PATH(PARMETIS_INCLUDE_DIR parmetis.h)
  FIND_LIBRARY(PARMETIS_LIBRARY parmetis)
endif()

if(PARMETIS_INCLUDE_DIR)
  set(PARMETIS_INCLUDE_DIRS
     ${PARMETIS_INCLUDE_DIR} ${METIS_INCLUDE_DIRS})
endif()
if(PARMETIS_LIBRARY)
  set(PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARIES})
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PARMETIS
   REQUIRED_VARS PARMETIS_INCLUDE_DIR PARMETIS_LIBRARY)

#
# Find the PARMETIS includes and libraries
#
# ParMETIS is an MPI-based parallel library that implements a variety of algorithms for 
# partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of 
# sparse matrices. It can be found at:
# 	http://www-users.cs.umn.edu/~karypis/metis/parmetis/index.html
#
# METIS_INCLUDE_DIR - where to find autopack.h
# METIS_LIBRARY     - List of fully qualified libraries to link against.
# METIS_FOUND       - Do not attempt to use if "no" or undefined.

FIND_PATH(METIS_INCLUDE_DIR metis.h)
FIND_LIBRARY(METIS_LIBRARY metis)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS
   REQUIRED_VARS METIS_INCLUDE_DIR METIS_LIBRARY)

set(METIS_INCLUDE_DIRS ${METIS_INCLUDE_DIR})
set(METIS_LIBRARIES ${METIS_LIBRARY})

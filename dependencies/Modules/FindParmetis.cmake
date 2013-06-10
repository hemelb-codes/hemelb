#
# Find the ParMETIS includes and libraries
#
# ParMETIS is an MPI-based parallel library that implements a variety of algorithms for 
# partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of 
# sparse matrices. It can be found at:
# 	http://www-users.cs.umn.edu/~karypis/metis/parmetis/index.html
#
# PARMETIS_INCLUDE_DIR  - where to find parmetis.h
# METIS_INCLUDE_DIR     - where to find metis.h
# PARMETIS_INCLUDE_DIRS - List of include dir.
# PARMETIS_LIBRARIES    - List of fully qualified libraries to link against.
# PARMETIS_FOUND        - Do not attempt to use if "no" or undefined.

find_path(PARMETIS_INCLUDE_DIR parmetis.h)
find_path(METIS_INCLUDE_DIR metis.h)

set(PARMETIS_INCLUDE_DIRS ${PARMETIS_INCLUDE_DIR} ${METIS_INCLUDE_DIR})

find_library(PARMETIS_LIB_PARMETIS parmetis)
find_library(PARMETIS_LIB_METIS metis)

set(PARMETIS_LIBRARIES ${PARMETIS_LIB_PARMETIS} ${PARMETIS_LIB_METIS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Parmetis  DEFAULT_MSG
                                  PARMETIS_INCLUDE_DIR METIS_INCLUDE_DIR
                                  PARMETIS_INCLUDE_DIRS
                                  PARMETIS_LIB_PARMETIS PARMETIS_LIB_METIS
                                  PARMETIS_LIBRARIES)

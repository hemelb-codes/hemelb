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

find_package(PkgConfig)
pkg_check_modules(PC_METIS QUIET metis)

find_path(METIS_INCLUDE_DIR metis.h HINTS ${METIS_ROOT}/include ${METIS_DIR}/include ${PC_METIS_INCLUDE_DIRS})
find_library(METIS_LIBRARY NAMES metis ${PC_METIS_LIBRARIES} HINTS ${METIS_ROOT}/lib ${METIS_ROOT}/lib64 ${PC_METIS_LIBRARY_DIRS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args("METIS" DEFAULT_MSG METIS_INCLUDE_DIR METIS_LIBRARY)
mark_as_advanced(METIS_INCLUDE_DIR METIS_LIBRARY)

# Figure out how METIS was compiled (i.e. what types does it use for
# indices and reals).
#
# category = IDX or REAL
function(_get_metis_type category widths types output)
  list(LENGTH widths nw)
  list(LENGTH types nt)
  if(NOT (${nw} EQUAL ${nt}))
    message(FATAL_ERROR "Widths and types have different lengths")
  endif()

  list(APPEND CMAKE_REQUIRED_INCLUDES ${METIS_INCLUDE_DIR})
  math(EXPR loop_max "${nw} - 1")
  foreach(i RANGE ${loop_max})
    list(GET widths ${i} WIDTH)
    CHECK_CXX_SOURCE_COMPILES("#include <metis.h>
int main() {
  static_assert(${category}TYPEWIDTH==${WIDTH}, \"Not ${WIDTH} bit\");
}
" HAVE_METIS_${category}_${WIDTH})
    if(${HAVE_METIS_${category}_${WIDTH}})
      list(GET types ${i} TYPE)
      set(${output} ${TYPE} PARENT_SCOPE)
      return()
    endif()
  endforeach()
  message(FATAL_ERROR "Could not figure out metis ${category} width")
endfunction()

if (METIS_FOUND)
  include(CheckCXXSourceCompiles)
  _get_metis_type(IDX "32;64" "int32_t;int64_t" METIS_IDX_T)
  _get_metis_type(REAL "32;64" "float;double" METIS_REAL_T)
  set(METIS_IDX_T ${METIS_IDX_T} PARENT_SCOPE)
  set(METIS_REAL_T ${METIS_REAL_T} PARENT_SCOPE)

  add_library(METIS::METIS INTERFACE IMPORTED)
  set_target_properties(METIS::METIS PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${METIS_LIBRARY}"
  )
endif()

# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

# ------MPI------------------
# Require MPI for this project:
find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${CMAKE_CXX_LINK_FLAGS}")
include_directories(${MPI_INCLUDE_PATH})

# Figure out if this MPI implementation supports MPI 3
set(CMAKE_REQUIRED_DEFINITIONS ${MPI_COMPILE_FLAGS})
set(CMAKE_REQUIRED_INCLUDES ${MPI_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_LIBRARIES})

CHECK_CXX_SOURCE_COMPILES("#include <mpi.h>
int main(int argc, char* argv[]) {
  static_assert(MPI_VERSION >= 3, \"\");
}" HAVE_MPI_STANDARD_VERSION_3)

if (NOT HAVE_MPI_STANDARD_VERSION_3)
  message(FATAL_ERROR "MPI standard version >= 3.0 required")
endif()

# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

# ------MPI------------------
# Require MPI for this project:
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${CMAKE_CXX_LINK_FLAGS}")
include_directories(${MPI_INCLUDE_PATH})

# Figure out if this MPI implementation has a const-correct API (supports MPI 3)
function(get_constcorrectmpi)
  set(CMAKE_REQUIRED_FLAGS -Werror)
  set(CMAKE_REQUIRED_DEFINITIONS ${MPI_COMPILE_FLAGS})
  set(CMAKE_REQUIRED_INCLUDES ${MPI_INCLUDE_PATH})
  set(CMAKE_REQUIRED_LIBRARIES ${MPI_LIBRARIES})
  CHECK_CXX_SOURCE_COMPILES("#include <mpi.h>
int main(int argc, char* argv[]) {
  const int send = 0;
  int recv;
  MPI_Request req;
  MPI_Init(&argc, &argv);
  MPI_Irecv(&recv, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &req);
  MPI_Send(&send, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
  MPI_Wait(&req, MPI_STATUS_IGNORE);
  MPI_Finalize();
}" HAVE_CONSTCORRECTMPI)
endfunction()

get_constcorrectmpi()

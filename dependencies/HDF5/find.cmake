# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()
cmake_policy(SET CMP0074 NEW)

set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 1.8 ${DEPS_FIND_MODE_HDF5} COMPONENTS C HL)

if(HDF5_FOUND AND NOT HDF5_IS_PARALLEL)
  string(FIND "${DEPS_FIND_MODE_HDF5}" REQUIRED req_ind)
  if(req_ind GREATER -1)
    set(mtype FATAL_ERROR)
  else()
    set(mtype WARNING)
  endif()
  message(${mtype} "HDF5 found but not parallel")
  set(HDF5_FOUND 0)
endif()

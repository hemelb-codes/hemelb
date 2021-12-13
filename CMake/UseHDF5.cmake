# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

include(HemeLbDependency)

find_hemelb_dependency(HDF5 REQUIRED)

if(NOT TARGET hdf5::hdf5)
  message(STATUS "Hacking in a proper target for HDF5")

  add_library(hdf5::hdf5 INTERFACE IMPORTED)
  set_target_properties(hdf5::hdf5 PROPERTIES
    INTERFACE_COMPILE_DEFINITIONS "${HDF5_C_DEFINITIONS}"
    INTERFACE_INCLUDE_DIRECTORIES "${HDF5_C_INCLUDE_DIRS}"
    INTERFACE_LINK_LIBRARIES "${HDF5_C_LIBRARIES}"
    )
endif()

if(NOT TARGET hdf5::hdf5_hl)
  message(STATUS "Hacking in a proper target for HDF5 high-level")

  add_library(hdf5::hdf5_hl INTERFACE IMPORTED)
  set_target_properties(hdf5::hdf5 PROPERTIES
    INTERFACE_LINK_LIBRARIES "hdf5::hdf5"
    )
endif()

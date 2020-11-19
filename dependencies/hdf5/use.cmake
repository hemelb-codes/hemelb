# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
find_package(HDF5 1.8 REQUIRED COMPONENTS C HL)

macro(hemelb_add_target_dependency_hdf5 tgt)
  target_link_libraries(${tgt} PRIVATE hdf5::hdf5_hl)
endmacro()

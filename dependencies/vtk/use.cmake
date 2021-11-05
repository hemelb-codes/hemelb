# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

find_package(VTK 8.2 REQUIRED NO_MODULE)

macro(hemelb_add_target_dependency_vtk tgt)
  target_link_libraries(${tgt} PRIVATE VTK::CommonCore VTK::CommonDataModel VTK::IOXML VTK::FiltersCore)
endmacro()

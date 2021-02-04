# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

include(HemeLbDependency)

hemelb_dependency(tinyxml use)
hemelb_dependency(parmetis use)
hemelb_dependency(boost use)
hemelb_dependency(ctemplate use)
hemelb_dependency(zlib use)

if(HEMELB_BUILD_RBC)
  hemelb_dependency(hdf5 use)
  hemelb_dependency(vtk use)
endif()

if(HEMELB_BUILD_MULTISCALE)
  hemelb_dependency(MPWide use)
endif()

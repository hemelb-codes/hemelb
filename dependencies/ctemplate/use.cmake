# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
hemelb_dependency(ctemplate find)

macro(hemelb_add_target_dependency_ctemplate tgt)
  target_include_directories(${tgt} PRIVATE ${CTEMPLATE_INCLUDE_DIR})
  target_link_libraries(${tgt} PRIVATE ${CTEMPLATE_LIBRARIES})
endmacro()

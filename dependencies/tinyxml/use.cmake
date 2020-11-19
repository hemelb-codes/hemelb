# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
hemelb_dependency(tinyxml find)

macro(hemelb_add_target_dependency_tinyxml tgt)
  # STL use is required by HemeLB and the build options 
  target_compile_definitions(${tgt}  PRIVATE -DTIXML_USE_STL)
  target_include_directories(${tgt} PRIVATE ${TINYXML_INCLUDE_DIR})
  target_link_libraries(${tgt} PRIVATE ${TINYXML_LIBRARIES})
endmacro()


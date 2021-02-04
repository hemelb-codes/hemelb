# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

hemelb_dependency(cppunit find)

macro(hemelb_add_target_dependency_cppunit tgt)
  target_include_directories(${tgt} PRIVATE ${CPPUNIT_INCLUDE_DIR})
  target_link_libraries(${tgt} PRIVATE ${CPPUNIT_LIBRARY}
    ${CMAKE_DL_LIBRARY}) #Because on some systems CPPUNIT needs to be linked to libdl
endmacro()

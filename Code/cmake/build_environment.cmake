# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
include_guard()

#------Capture build environment -------------
find_package(Git REQUIRED)
execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
  RESULT_VARIABLE rev_ok
  OUTPUT_VARIABLE HEMELB_REVISION_NUMBER
  OUTPUT_STRIP_TRAILING_WHITESPACE
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

if (NOT rev_ok EQUAL 0)
  message("Could not get revision number from Git, looking for revision_info.txt")
  execute_process(COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../revision_info.txt
    RESULT_VARIABLE rev_ok
    OUTPUT_VARIABLE HEMELB_REVISION_NUMBER
    OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
message("Revision number: ${HEMELB_REVISION_NUMBER}")

execute_process(COMMAND date
  OUTPUT_VARIABLE HEMELB_BUILD_TIME
  OUTPUT_STRIP_TRAILING_WHITESPACE)


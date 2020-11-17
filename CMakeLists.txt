# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

cmake_minimum_required (VERSION 3.10)

Project(HemeLBMaster)

include(ExternalProject)

set(HEMELB_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# Include the options for the main code
include("Code/cmake/options.cmake")
# and create a list of them to pass to the main configure step
hemelb_pass_cmake_defines(HEMELB_CODE_CMAKE_DEFINES)

#------- Dependencies -----------

add_subdirectory(dependencies)

# ------ HemeLB ------------------
message("${PROJECT_SOURCE_DIR}/Code")

ExternalProject_Add(
  hemelb
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/Code
  CONFIGURE_COMMAND "${CMAKE_COMMAND};${PROJECT_SOURCE_DIR}/Code;${HEMELB_CODE_CMAKE_DEFINES}"
  BUILD_COMMAND make -j${HEMELB_SUBPROJECT_MAKE_JOBS}
  )

# These are handled in dependencies/CMakeLists.txt and needed here
# to ensure they are build in correct order.

add_dependencies(hemelb parmetis tinyxml boost catch2 ctemplate zlib)
if(HEMELB_BUILD_MULTISCALE)
  add_dependencies(hemelb MPWide)
endif()

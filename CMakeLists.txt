# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

cmake_minimum_required (VERSION 3.13)

project(HemeLBSuperBuild)

set(HEMELB_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})
list(INSERT CMAKE_MODULE_PATH 0 "${HEMELB_ROOT_DIR}/CMake")

include(ExternalProject)
# Include all the options for the main code
include(GlobalOptions)
include(HemeLbOptions)
include(DependencyOptions)

#------- Dependencies -----------

add_subdirectory(dependencies)

# ------ HemeLB ------------------
# Create a list of them to pass to the main configure step
set(HEMELB_CODE_CMAKE_DEFINES "")
pass_cmake_defines(GLOBAL HEMELB_CODE_CMAKE_DEFINES)
pass_cmake_defines(HEMELB HEMELB_CODE_CMAKE_DEFINES)
message("${CMAKE_COMMAND};${PROJECT_SOURCE_DIR}/Code;${HEMELB_CODE_CMAKE_DEFINES}")
ExternalProject_Add(
  hemelb
  SOURCE_DIR ${PROJECT_SOURCE_DIR}/Code
  CONFIGURE_COMMAND "${CMAKE_COMMAND};${PROJECT_SOURCE_DIR}/Code;${HEMELB_CODE_CMAKE_DEFINES}"
  BUILD_COMMAND make -j${HEMELB_SUBPROJECT_MAKE_JOBS}
  )

# These are handled in dependencies/CMakeLists.txt and needed here
# to ensure they are build in correct order.
list(LENGTH HEMELB_PROJECT_DEPENDENCIES ndeps)
if (ndeps GREATER 0)
  add_dependencies(hemelb ${HEMELB_PROJECT_DEPENDENCIES})
endif()

# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

include(PassOptions)

# This file defines the CMake options and cache variables across the
# Superbuild, Code and dependencies project. Basically those HemeLB
# options that affect which dependencies are needed.

pass_group_declare(GLOBAL)

pass_option(GLOBAL HEMELB_BUILD_MULTISCALE "Build HemeLB Multiscale functionality" OFF)
pass_option(GLOBAL HEMELB_BUILD_RBC "Build the resolved red blood cells functionality" OFF)
pass_option(GLOBAL HEMELB_BUILD_TESTS "Build the tests" ON)

pass_cachevar(GLOBAL HEMELB_DEPENDENCIES_PATH "${HEMELB_ROOT_DIR}/dependencies"
  FILEPATH "Path to find dependency find modules")
pass_cachevar(GLOBAL HEMELB_DEPENDENCIES_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}
  FILEPATH "Path to find dependency includes and libraries")
pass_cachevar(GLOBAL HEMELB_SUBPROJECT_MAKE_JOBS 1
  STRING "Number of jobs to use for subproject build steps")

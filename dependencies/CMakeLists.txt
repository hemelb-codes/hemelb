# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

cmake_minimum_required (VERSION 3.10)

Project(HemeLBDependencies)

include(ExternalProject)

if(NOT DEFINED HEMELB_ROOT_DIR)
  set(HEMELB_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)
endif()
if(NOT DEFINED HEMELB_DEPENDENCIES_PATH)
  set(HEMELB_DEPENDENCIES_PATH ${CMAKE_CURRENT_SOURCE_DIR})
endif()
if(NOT DEFINED HEMELB_DEPENDENCIES_INSTALL_PATH)
  set(HEMELB_DEPENDENCIES_INSTALL_PATH ${HEMELB_DEPENDENCIES_PATH} CACHE STRING "Path to install Dependencies")
endif()

list(APPEND CMAKE_PREFIX_PATH "${HEMELB_DEPENDENCIES_INSTALL_PATH}")
set(CMAKE_MODULE_PATH "${HEMELB_DEPENDENCIES_PATH}/Modules/" ${CMAKE_MODULE_PATH})

#option(HEMELB_BUILD_MULTISCALE "Build HemeLB Multiscale functionality" OFF)

find_package(MPI REQUIRED) # Do this *after* we include the patched find MPI module

include(HemeLbDependency)

hemelb_dependency(tinyxml find_or_build)
hemelb_dependency(parmetis find_or_build)
hemelb_dependency(catch2 find_or_build)
hemelb_dependency(boost find_or_build)
hemelb_dependency(ctemplate find_or_build)
hemelb_dependency(zlib find_or_build)

if (HEMELB_BUILD_RBC)
  hemelb_dependency(hdf5 find_or_build)
  hemelb_dependency(vtk find_or_build)
endif()

if (HEMELB_BUILD_MULTISCALE)
  hemelb_dependency(MPWide find_or_build)
endif()

# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
cmake_minimum_required (VERSION 3.10)

project(HemeLB)

set(HEMELB_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#---- OPTION switches accepted by the build -------
include("cmake/options.cmake")

# Add warnings flags to development build types
if (HEMELB_USE_ALL_WARNINGS_GNU)
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
  set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
endif()
if(HEMELB_OPTIMISATION)
  set( CMAKE_CXX_FLAGS_RELEASE ${HEMELB_OPTIMISATION})
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HEMELB_OPTIMISATION}")
endif()


add_definitions(-DHEMELB_CODE)
add_definitions(-DHEMELB_READING_GROUP_SIZE=${HEMELB_READING_GROUP_SIZE})
add_definitions(-DHEMELB_COMPUTE_ARCHITECTURE=${HEMELB_COMPUTE_ARCHITECTURE})

if(HEMELB_VALIDATE_GEOMETRY)
  add_definitions(-DHEMELB_VALIDATE_GEOMETRY)
endif()

if (NOT HEMELB_USE_STREAKLINES)
  add_definitions(-DNO_STREAKLINES)
endif()

if (HEMELB_WAIT_ON_CONNECT)
  add_definitions(-DHEMELB_WAIT_ON_CONNECT)
endif()

if (HEMELB_IMAGES_TO_NULL)
  add_definitions(-DHEMELB_IMAGES_TO_NULL)
endif()

if (HEMELB_USE_SSE3)
  add_definitions(-DHEMELB_USE_SSE3)
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
  set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -msse3")
  set( CMAKE_CXX_FLAGS_RELEASE "${HEMELB_OPTIMISATION} -msse3")
endif()

if (HEMELB_BUILD_RBC)
  add_definitions(-DHEMELB_BUILD_RBC)
endif()

if (HEMELB_CALLGRIND)
  add_definitions(-DHEMELB_CALLGRIND)
endif()

if (HEMELB_USE_VELOCITY_WEIGHTS_FILE)
  add_definitions(-DHEMELB_USE_VELOCITY_WEIGHTS_FILE)
endif()

if (HEMELB_USE_KRUEGER_ORDERING)
  add_definitions(-DHEMELB_USE_KRUEGER_ORDERING)
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" "${HEMELB_DEPENDENCIES_PATH}/Modules/")
list(APPEND CMAKE_PREFIX_PATH ${HEMELB_DEPENDENCIES_INSTALL_PATH})
list(APPEND CMAKE_INCLUDE_PATH ${HEMELB_DEPENDENCIES_INSTALL_PATH}/include)
list(APPEND CMAKE_LIBRARY_PATH ${HEMELB_DEPENDENCIES_INSTALL_PATH}/lib)

#--- Set up runtime search path for DLLs -----
if(HEMELB_DEPENDENCIES_SET_RPATH)
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  SET(CMAKE_INSTALL_RPATH "${HEMELB_DEPENDENCIES_INSTALL_PATH}/lib")
endif()

include(build_environment)
include(platform_checks)

if(NOT UBUNTU_BUG_WORKAROUND)
  if(HAVE_ISNAN)
    add_definitions(-DHAVE_ISNAN)
  endif()
endif()

if(HAVE_STD_ISNAN)
  add_definitions(-DHAVE_STD_ISNAN)
endif()

if(HAVE_RUSAGE)
  add_definitions(-DHAVE_RUSAGE)
endif()

if(LINUX_SCANDIR)
  add_definitions(-DLINUX_SCANDIR)
endif()

include(mpi)
include(dependencies)

#-------------Resources -----------------------

set(BUILD_RESOURCE_PATH ${PROJECT_BINARY_DIR}/resources)
file(MAKE_DIRECTORY ${BUILD_RESOURCE_PATH})
set(INSTALL_RESOURCE_PATH ${CMAKE_INSTALL_PREFIX}/share/hemelb/resources)

function(add_to_resources)
  unset(donotexist)
  set(hasmissing FALSE)
  foreach(resource ${ARGN})
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${resource}")
      set(donotexist "${donotexist}\n - ${CMAKE_CURRENT_SOURCE_DIR}/${resource}")
      set(hasmissing TRUE)
    else()
      configure_file(${CMAKE_CURRENT_SOURCE_DIR}/${resource} ${BUILD_RESOURCE_PATH} COPYONLY)
      get_filename_component(resource_leaf ${resource} NAME)
      install(FILES ${BUILD_RESOURCE_PATH}/${resource_leaf} DESTINATION share/hemelb/resources)
    endif()
  endforeach()
  if(hasmissing)
    message(WARNING "The following resources are missing:${donotexist}")
  endif()
endfunction()

configure_file (
  "${PROJECT_SOURCE_DIR}/resources/path_parameters.h.in"
  "${PROJECT_BINARY_DIR}/resources/path_parameters.h"
  )
configure_file("Traits.in.h" "${PROJECT_BINARY_DIR}/Traits.h")

# add the binary tree to the search path for include files
# so that we will find TutorialConfig.h
include_directories("${PROJECT_BINARY_DIR}")

# ----------- HemeLB ------------------
if (APPLE)
  add_definitions(-DHEMELB_CFG_ON_BSD -DHEMELB_CFG_ON_OSX)
endif()
if (HEMELB_USE_DEBUGGER)
  # Force a debug build, because the debugger makes no sense without debug symbols
  set(CMAKE_BUILD_TYPE DEBUG)
endif()

add_executable(${HEMELB_EXECUTABLE} main.cc)

include_directories(${PROJECT_SOURCE_DIR})
set(package_subdirs
  configuration
  extraction
  reporting
  steering
  vis
  lb
  geometry
  net
  debug
  util
  io
  log
  colloids
  redblood
  )

foreach(subdir ${package_subdirs})
  set(lib "hemelb_${subdir}")
  list(APPEND heme_libraries ${lib})
  add_subdirectory(${subdir})
endforeach()
add_subdirectory(resources)
target_link_libraries(${HEMELB_EXECUTABLE} 
  ${heme_libraries} ${heme_libraries}
  ${MPI_LIBRARIES}
  )

# Config XML checker (confcheck)
add_executable(${HEMELB_EXECUTABLE}-confcheck confcheck.cc)
target_link_libraries(${HEMELB_EXECUTABLE}-confcheck
  ${heme_libraries} ${heme_libraries}
  ${MPI_LIBRARIES}
  )

INSTALL(TARGETS ${HEMELB_EXECUTABLE} ${HEMELB_EXECUTABLE}-confcheck RUNTIME DESTINATION bin)
add_to_resources(resources/report.txt.ctp resources/report.xml.ctp)


# ----------- HemeLB Multiscale ------------------
if (HEMELB_BUILD_MULTISCALE)
  if (APPLE)
    add_definitions(-DHEMELB_CFG_ON_BSD -DHEMELB_CFG_ON_OSX)
  endif()
  if (HEMELB_USE_DEBUGGER)
    # Force a debug build, because the debugger makes no sense without debug symbols
    set(CMAKE_BUILD_TYPE DEBUG)
  endif()
  
  add_executable(multiscale_hemelb mainMultiscale.cc)
  
  include_directories(${PROJECT_SOURCE_DIR})

  add_subdirectory(multiscale)
  target_link_libraries(multiscale_hemelb
    PRIVATE ${heme_libraries} hemelb_multiscale ${MPI_LIBRARIES}
    )
  INSTALL(TARGETS multiscale_hemelb RUNTIME DESTINATION bin)
  add_to_resources(resources/report.txt.ctp resources/report.xml.ctp)
endif()

# ----------- HEMELB unittests ---------------
if(HEMELB_BUILD_TESTS_ALL OR HEMELB_BUILD_TESTS_UNIT)
  enable_testing()
  add_subdirectory(tests)
  add_subdirectory(unittests)
endif()

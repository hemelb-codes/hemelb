# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
cmake_minimum_required (VERSION 3.10)

project(HemeLB)

set(HEMELB_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)

set(CMAKE_CXX_STANDARD 11)
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
add_definitions(-DHEMELB_LATTICE=${HEMELB_LATTICE})
add_definitions(-DHEMELB_KERNEL=${HEMELB_KERNEL})
add_definitions(-DHEMELB_WALL_BOUNDARY=${HEMELB_WALL_BOUNDARY})
add_definitions(-DHEMELB_INLET_BOUNDARY=${HEMELB_INLET_BOUNDARY})
add_definitions(-DHEMELB_OUTLET_BOUNDARY=${HEMELB_OUTLET_BOUNDARY})
add_definitions(-DHEMELB_WALL_INLET_BOUNDARY=${HEMELB_WALL_INLET_BOUNDARY})
add_definitions(-DHEMELB_WALL_OUTLET_BOUNDARY=${HEMELB_WALL_OUTLET_BOUNDARY})
add_definitions(-DHEMELB_COMPUTE_ARCHITECTURE=${HEMELB_COMPUTE_ARCHITECTURE})
add_definitions(-DHEMELB_LOG_LEVEL=${HEMELB_LOG_LEVEL})

if(HEMELB_VALIDATE_GEOMETRY)
  add_definitions(-DHEMELB_VALIDATE_GEOMETRY)
endif()

if (NOT HEMELB_USE_STREAKLINES)
  add_definitions(-DNO_STREAKLINES)
endif()

if (HEMELB_WAIT_ON_CONNECT)
  add_definitions(-DHEMELB_WAIT_ON_CONNECT)
endif()

if (NOT HEMELB_STATIC_ASSERT)
  add_definitions(-DHEMELB_NO_STATIC_ASSERT)
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

if (HEMELB_USE_VELOCITY_WEIGHTS_FILE)
  add_definitions(-DHEMELB_USE_VELOCITY_WEIGHTS_FILE)
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" "${HEMELB_DEPENDENCIES_PATH}/Modules/")
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

configure_file (
  "${PROJECT_SOURCE_DIR}/resources/path_parameters.h.in"
  "${PROJECT_BINARY_DIR}/resources/path_parameters.h"
  )

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

set(root_sources SimulationMaster.cc)
add_executable(${HEMELB_EXECUTABLE} main.cc ${root_sources})

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
  )
foreach(subdir ${package_subdirs})
  set(lib "hemelb_${subdir}")
  list(APPEND heme_libraries ${lib})
  add_subdirectory(${subdir})
endforeach()
add_subdirectory(resources)
target_link_libraries(${HEMELB_EXECUTABLE} 
  ${heme_libraries}
  ${MPI_LIBRARIES}
  ${Boost_LIBRARIES}
  )
INSTALL(TARGETS ${HEMELB_EXECUTABLE} RUNTIME DESTINATION bin)
list(APPEND RESOURCES resources/report.txt.ctp resources/report.xml.ctp)

# ----------- HemeLB Multiscale ------------------
if (HEMELB_BUILD_MULTISCALE)
  if (APPLE)
    add_definitions(-DHEMELB_CFG_ON_BSD -DHEMELB_CFG_ON_OSX)
  endif()
  if (HEMELB_USE_DEBUGGER)
    # Force a debug build, because the debugger makes no sense without debug symbols
    set(CMAKE_BUILD_TYPE DEBUG)
  endif()
  
  set(root_sources SimulationMaster.cc multiscale/MultiscaleSimulationMaster.h)
  add_executable(multiscale_hemelb mainMultiscale.cc ${root_sources})
  hemelb_add_target_dependency_mpwide(multiscale_hemelb)
  
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
    multiscale
    )
  foreach(subdir ${package_subdirs})
    set(lib "hemelb_${subdir}")
    list(APPEND heme_libraries ${lib})
  endforeach()
  add_subdirectory(multiscale)
  target_link_libraries(multiscale_hemelb
    ${heme_libraries}
    ${MPI_LIBRARIES}
    ${Boost_LIBRARIES}
    )
  INSTALL(TARGETS multiscale_hemelb RUNTIME DESTINATION bin)
  list(APPEND RESOURCES resources/report.txt.ctp resources/report.xml.ctp)
endif()

# ----------- HEMELB unittests ---------------
if(HEMELB_BUILD_TESTS_ALL OR HEMELB_BUILD_TESTS_UNIT)
  #------CPPUnit ---------------
  hemelb_dependency(cppunit use)
  add_executable(unittests_hemelb ${root_sources})
  add_subdirectory(unittests)
  target_link_libraries(unittests_hemelb 
    hemelb_unittests 
    ${heme_libraries}
    ${MPI_LIBRARIES}
    ${Boost_LIBRARIES})
  INSTALL(TARGETS unittests_hemelb RUNTIME DESTINATION bin)
  list(APPEND RESOURCES unittests/resources/four_cube.gmy unittests/resources/four_cube.xml unittests/resources/four_cube_multiscale.xml
    unittests/resources/config.xml unittests/resources/config0_2_0.xml
    unittests/resources/config_file_inlet.xml unittests/resources/iolet.txt 
    unittests/resources/config-velocity-iolet.xml unittests/resources/config_new_velocity_inlets.xml
    unittests/resources/velocity_inlet.txt.weights.txt
    unittests/resources/xmltest.xml unittests/resources/config_file_velocity_inlet.xml unittests/resources/velocity_inlet.txt 
    )
endif()

# ----------- HEMELB functionaltests ---------------
if(HEMELB_BUILD_TESTS_ALL OR HEMELB_BUILD_TESTS_FUNCTIONAL)
  #------CPPUnit ---------------
  hemelb_dependency(cppunit use)
  add_executable(functionaltests_hemelb ${root_sources})
  add_subdirectory(functionaltests/cpptests)
  target_link_libraries(functionaltests_hemelb 
    hemelb_functionaltests 
    ${heme_libraries}
    ${MPI_LIBRARIES}
    ${Boost_LIBRARIES})
  INSTALL(TARGETS functionaltests_hemelb RUNTIME DESTINATION bin)
endif()

#-------- Copy and install resources --------------

foreach(resource ${RESOURCES})
  configure_file(${PROJECT_SOURCE_DIR}/${resource} ${BUILD_RESOURCE_PATH} COPYONLY)
  get_filename_component(resource_leaf ${resource} NAME)
  INSTALL(FILES ${BUILD_RESOURCE_PATH}/${resource_leaf} DESTINATION share/hemelb/resources)
endforeach()

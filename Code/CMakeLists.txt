# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.
cmake_minimum_required (VERSION 3.10)

project(HemeLB
  VERSION 0.0.0
  DESCRIPTION "A high performance parallel lattice-Boltzmann code for large scale fluid flow in complex geometries"
  LANGUAGES C CXX CUDA
)

include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

set(HEMELB_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/..)

#---- OPTION switches accepted by the build -------
include("cmake/options.cmake")

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set( CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++11")

if (HEMELB_USE_ALL_WARNINGS_GNU)
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

if(HEMELB_OPTIMISATION)
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HEMELB_OPTIMISATION}")
endif()

if(HEMELB_PROFILING)
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HEMELB_PROFILING}")
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
endif()

if (HEMELB_USE_VELOCITY_WEIGHTS_FILE)
  add_definitions(-DHEMELB_USE_VELOCITY_WEIGHTS_FILE)
endif()

if (HEMELB_LATTICE_INCOMPRESSIBLE)
  add_definitions(-DHEMELB_LATTICE_INCOMPRESSIBLE)
endif()

if (HEMELB_CUDA_AWARE_MPI)
  add_definitions(-DHEMELB_CUDA_AWARE_MPI)
endif()

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" "${HEMELB_DEPENDENCIES_PATH}/Modules/")
list(APPEND CMAKE_INCLUDE_PATH ${HEMELB_DEPENDENCIES_INSTALL_PATH}/include)
list(APPEND CMAKE_LIBRARY_PATH ${HEMELB_DEPENDENCIES_INSTALL_PATH}/lib)

#--- Set up runtime search path for DLLs -----
if(HEMELB_DEPENDENCIES_SET_RPATH)
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  set(CMAKE_INSTALL_RPATH "${HEMELB_DEPENDENCIES_INSTALL_PATH}/lib")
endif()

include(build_environment)
include(platform_checks)

if (APPLE)
  add_definitions(-DHEMELB_CFG_ON_BSD -DHEMELB_CFG_ON_OSX)
endif()

if(HAVE_ISNAN)
  add_definitions(-DHAVE_ISNAN)
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

include_directories(${PROJECT_SOURCE_DIR})
set(package_subdirs
  configuration
  extraction
  reporting
  steering
  vis
  geometry
  lb
  net
  debug
  util
  io
  logging
  colloids
)
foreach(subdir ${package_subdirs})
  set(lib "hemelb_${subdir}")
  list(APPEND heme_libraries ${lib})
  add_subdirectory(${subdir})
endforeach()
add_subdirectory(resources)

list(APPEND RESOURCES
  resources/report.txt.ctp
  resources/report.xml.ctp
)

# ----------- HemeLB ------------------
if (HEMELB_BUILD_MAIN)
  set(root_sources SimulationMaster.cc)
  add_executable(${HEMELB_EXECUTABLE} main.cc ${root_sources})

  target_link_libraries(${HEMELB_EXECUTABLE}
    ${heme_libraries}
    ${MPI_LIBRARIES}
    ${Boost_LIBRARIES}
    )
  install(TARGETS ${HEMELB_EXECUTABLE} RUNTIME DESTINATION bin)
endif()

# ----------- HemeLB Multiscale ------------------
if (HEMELB_BUILD_MULTISCALE)
  set(root_sources SimulationMaster.cc multiscale/MultiscaleSimulationMaster.h)
  add_executable(multiscale_hemelb mainMultiscale.cc ${root_sources})
  hemelb_add_target_dependency_mpwide(multiscale_hemelb)

  add_subdirectory(multiscale)

  target_link_libraries(multiscale_hemelb
    ${heme_libraries}
    ${MPI_LIBRARIES}
    ${Boost_LIBRARIES}
    )
  install(TARGETS multiscale_hemelb RUNTIME DESTINATION bin)
endif()

# ----------- HEMELB unittests ---------------
if(HEMELB_BUILD_TESTS_ALL OR HEMELB_BUILD_TESTS_UNIT)
  #------CPPUnit ---------------
  hemelb_dependency(cppunit use)
  add_executable(hemelb-unittests ${root_sources})
  add_subdirectory(unittests)

  target_link_libraries(hemelb-unittests
    hemelb_unittests
    ${heme_libraries}
    ${MPI_LIBRARIES}
    ${Boost_LIBRARIES})
  install(TARGETS hemelb-unittests RUNTIME DESTINATION bin)

  list(APPEND RESOURCES
    unittests/resources/four_cube.gmy
    unittests/resources/four_cube.xml
    unittests/resources/four_cube_multiscale.xml
    unittests/resources/config.xml
    unittests/resources/config0_2_0.xml
    unittests/resources/config_file_inlet.xml
    unittests/resources/iolet.txt
    unittests/resources/config-velocity-iolet.xml
    unittests/resources/config_new_velocity_inlets.xml
    unittests/resources/velocity_inlet.txt.weights.txt
    unittests/resources/xmltest.xml
    unittests/resources/config_file_velocity_inlet.xml
    unittests/resources/velocity_inlet.txt
  )
endif()

# ----------- HEMELB functionaltests ---------------
if(HEMELB_BUILD_TESTS_ALL OR HEMELB_BUILD_TESTS_FUNCTIONAL)
  #------CPPUnit ---------------
  hemelb_dependency(cppunit use)
  add_executable(hemelb-functionaltests ${root_sources})
  add_subdirectory(functionaltests/cpptests)

  target_link_libraries(hemelb-functionaltests
    hemelb_functionaltests
    ${heme_libraries}
    ${MPI_LIBRARIES}
    ${Boost_LIBRARIES})
  install(TARGETS hemelb-functionaltests RUNTIME DESTINATION bin)
endif()

#-------- Copy and install resources --------------
foreach(resource ${RESOURCES})
  configure_file(${PROJECT_SOURCE_DIR}/${resource} ${BUILD_RESOURCE_PATH} COPYONLY)
  get_filename_component(resource_leaf ${resource} NAME)
  install(FILES ${BUILD_RESOURCE_PATH}/${resource_leaf} DESTINATION share/hemelb/resources)
endforeach()

cmake_minimum_required (VERSION 2.8)

project(HemeLB)

#---- OPTION switches accepted by the build -------

option(HEMELB_USE_DEBUGGER "Use built in hemelb debugger" OFF)
option(HEMELB_USE_MULTIMACHINE "Use multi-level parallelism support" OFF)
option(HEMELB_BUILD_UNITTESTS "Build the unit-tests" ON)
option(HEMELB_USE_ALL_WARNINGS_GNU "Show all compiler warnings on development builds (gnu-style-compilers)" ON)
option(HEMELB_USE_STREAKLINES "Calculate streakline images" ON)
option(HEMELB_USE_BOOST "Use Boost" OFF)
option(CTEMPLATE_USE_STATIC "Prefer Static CTemplate library" OFF)
set(HEMELB_LOG_LEVEL info
	CACHE STRING "Log level, choose 'debug', 'warning', or 'info'" )
set(HEMELB_STEERING_LIB basic
	CACHE STRING "Steering library, choose 'basic' or 'none'" )
set(HEMELB_DEPENDENCIES_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies"
	CACHE FILEPATH "Path to find dependency find modules")
set(HEMELB_DEPENDENCIES_INSTALL_PATH ${HEMELB_DEPENDENCIES_PATH}
	CACHE FILEPATH "Path to find dependency includes and libraries")
set(HEMELB_OPTIMISATION "-O4" CACHE STRING "Optimisation level (can be blank or -O1 to -O4)")
	
# Add warnings flags to development build types
if (HEMELB_USE_ALL_WARNINGS_GNU)
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
	set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
endif()
if(HEMELB_OPTIMISATION)
	set( CMAKE_CXX_FLAGS_RELEASE ${HEMELB_OPTIMISATION})
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HEMELB_OPTIMISATION}")
endif()

if (NOT HEMELB_USE_STREAKLINES)
	add_definitions(-DNO_STREAKLINES)
endif()



#------Dependencies --------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${HEMELB_DEPENDENCIES_PATH}/Modules/")
list(APPEND CMAKE_INCLUDE_PATH ${HEMELB_DEPENDENCIES_INSTALL_PATH}/include)
list(APPEND CMAKE_LIBRARY_PATH ${HEMELB_DEPENDENCIES_INSTALL_PATH}/lib)

#--- Set up runtime search path for DLLs -----

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
SET(CMAKE_INSTALL_RPATH "${HEMELB_DEPENDENCIES_PATH}/lib")

#-------Check for platform dependent features ----

include(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES("#include <cmath>\n int main(int c,char** v){ return isnan(1.0); }" HAVE_ISNAN)
CHECK_CXX_SOURCE_COMPILES("#include <cmath>\n int main(int c,char** v){ return std::isnan(1.0); }" HAVE_STD_ISNAN)

if(HAVE_ISNAN)
  add_definitions(-DHAVE_ISNAN)
endif()

if( HAVE_STD_ISNAN)
  add_definitions(-D HAVE_STD_ISNAN)
endif()

# ------MPI------------------
# Require MPI for this project:
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${CMAKE_CXX_LINK_FLAGS}")
include_directories(${MPI_INCLUDE_PATH})

#------Parmetis  ------------
find_package(Parmetis REQUIRED)
include_directories(${PARMETIS_INCLUDE_DIR})
MESSAGE(STATUS "Found Parmetis: ${PARMETIS_LIBRARY}")

#------TinyXML ----------------
find_package(TinyXML REQUIRED)
OPTION(TIXML_USE_STL "Use STL with TIXML" ON)
if(TIXML_USE_STL)
	add_definitions(-DTIXML_USE_STL)
endif()
include_directories(${TINYXML_INCLUDE_DIR})
if(HEMELB_USE_BOOST)
	#------BOOST ------------------
	SET(Boost_ADDITIONAL_VERSIONS "1.48" "1.48.0")
	SET(BOOST_LIBS mpi regex)
	find_package(Boost 1.48 REQUIRED)
	include_directories(${Boost_INCLUDE_DIRS})
endif()
#------CTemplate ----------------
find_package(CTemplate REQUIRED)
include_directories(${CTEMPLATE_INCLUDE_DIR})

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

# ----------- HemeLB ------------------s
if (APPLE)
	add_definitions(-DHEMELB_CFG_ON_BSD -DHEMELB_CFG_ON_OSX)
endif()
if (HEMELB_USE_DEBUGGER)
	# Force a debug build, because the debugger makes no sense without debug symbols
	set(CMAKE_BUILD_TYPE DEBUG)
endif()

set(root_sources SimulationMaster.cc mpiInclude.cc D3Q15.cc)
add_executable(hemelb main.cc ${root_sources})
include_directories(${PROJECT_SOURCE_DIR})
set(package_subdirs
	configuration
	reporting
	steering
	vis
	lb
	net
	debug
	topology
	util
	geometry
	io
	log
	)
foreach(subdir ${package_subdirs})
	set(lib "hemelb_${subdir}")
	list(APPEND heme_libraries ${lib})
	add_subdirectory(${subdir})
endforeach()
add_subdirectory(resources)
target_link_libraries(hemelb 
	${heme_libraries}
	${MPI_LIBRARIES}
	${PARMETIS_LIBRARIES}
	${TINYXML_LIBRARIES}
	${Boost_LIBRARIES}
	${CTEMPLATE_LIBRARIES}
	)
INSTALL(TARGETS hemelb RUNTIME DESTINATION bin)
list(APPEND RESOURCES resources/report.txt.ctp resources/report.xml.ctp)
# ----------- HEMELB unittests ---------------
if(HEMELB_BUILD_UNITTESTS)
	#------CPPUnit ---------------
		find_package(CPPUnit REQUIRED)
		include_directories(${CPPUNIT_INCLUDE_DIR})
	add_executable(unittests_hemelb ${root_sources})
	add_subdirectory(unittests)
	target_link_libraries(unittests_hemelb 
		hemelb_unittests 
		${heme_libraries}
		${MPI_LIBRARIES}
		${PARMETIS_LIBRARIES}
		${TINYXML_LIBRARIES}
		${CPPUNIT_LIBRARY}
		${Boost_LIBRARIES}
		${CTEMPLATE_LIBRARIES})
	INSTALL(TARGETS unittests_hemelb RUNTIME DESTINATION bin)
	list(APPEND RESOURCES unittests/resources/four_cube.dat unittests/resources/four_cube.xml unittests/resources/config.xml)
endif()

#-------- Copy and install resources --------------

foreach(resource ${RESOURCES})
	configure_file(${PROJECT_SOURCE_DIR}/${resource} ${BUILD_RESOURCE_PATH} COPYONLY)
endforeach()

file(GLOB ${BUILD_RESOURCE_PATH}/* resource_files)
INSTALL(TARGETS ${resource_files} RESOURCE DESTINATION ${BUILD_RESOURCE_PATH})

cmake_minimum_required (VERSION 2.8)

project(HemeLB)

#---- OPTION switches accepted by the build -------

option(HEMELB_USE_DEBUGGER "Use built in hemelb debugger" OFF)
option(HEMELB_VALIDATE_GEOMETRY "Validate geometry" OFF)
option(HEMELB_USE_MULTIMACHINE "Use multi-level parallelism support" OFF)
option(HEMELB_BUILD_TESTS_ALL "Build all the tests" ON)
option(HEMELB_BUILD_TESTS_UNIT "Build the unit-tests (HEMELB_BUILD_TESTS_ALL takes precedence)" ON)
option(HEMELB_BUILD_TESTS_FUNCTIONAL "Build the functional tests (HEMELB_BUILD_TESTS_ALL takes precedence)" ON)
option(HEMELB_USE_ALL_WARNINGS_GNU "Show all compiler warnings on development builds (gnu-style-compilers)" ON)
option(HEMELB_USE_STREAKLINES "Calculate streakline images" ON)
option(HEMELB_USE_BOOST "Use Boost" ON)
option(HEMELB_DEPENDENCIES_SET_RPATH "Set runtime RPATH" ON)
option(HEMELB_STATIC_ASSERT "Use simple compile-time assertions" ON)
option(HEMELB_WAIT_ON_CONNECT "Wait for steering client" OFF)
option(HEMELB_BUILD_MULTISCALE "Build HemeLB Multiscale functionality" OFF)
option(HEMELB_IMAGES_TO_NULL "Write images to null" OFF)

set(HEMELB_READING_GROUP_SIZE 5
  CACHE INTEGER "Number of cores to use to read geometry file.")
set(HEMELB_LOG_LEVEL info
	CACHE STRING "Log level, choose 'debug', 'warning', or 'info'" )
set(HEMELB_STEERING_LIB basic
	CACHE STRING "Steering library, choose 'basic' or 'none'" )
set(HEMELB_DEPENDENCIES_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../dependencies"
	CACHE FILEPATH "Path to find dependency find modules")
set(HEMELB_DEPENDENCIES_INSTALL_PATH ${HEMELB_DEPENDENCIES_PATH}
	CACHE FILEPATH "Path to find dependency includes and libraries")
set(HEMELB_OPTIMISATION "-O4" CACHE STRING "Optimisation level (can be blank or -O1 to -O4)")
set(HEMELB_LATTICE "D3Q15"
  CACHE STRING "Select the lattice type to use (D3Q15,D3Q19,D3Q27)")
set(HEMELB_KERNEL "LBGK"
  CACHE STRING "Select the kernel to use (LBGK,EntropicAnsumali,EntropicChik,MRT,NNCY,NNC,NNTPL)")
set(HEMELB_WALL_BOUNDARY "SIMPLEBOUNCEBACK"
  CACHE STRING "Select the boundary conditions to be used at the walls (FINTERPOLATION,GZS,REGULARISED,SIMPLEBOUNCEBACK,JUNKYANG)")
set(HEMELB_POINTPOINT_IMPLEMENTATION Coalesce
	CACHE STRING "Point to point comms implementation, choose 'Coalesce', 'Separated', or 'Immediate'" )
set(HEMELB_GATHERS_IMPLEMENTATION Separated
	CACHE STRING "Gather comms implementation, choose 'Separated', or 'ViaPointPoint'" )
set(HEMELB_ALLTOALL_IMPLEMENTATION Separated
	CACHE STRING "Alltoall comms implementation, choose 'Separated', or 'ViaPointPoint'" )
option(HEMELB_SEPARATE_CONCERNS "Communicate for each concern separately" OFF)
	
# Add warnings flags to development build types
if (HEMELB_USE_ALL_WARNINGS_GNU)
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
	set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
endif()
if(HEMELB_OPTIMISATION)
	set( CMAKE_CXX_FLAGS_RELEASE ${HEMELB_OPTIMISATION})
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${HEMELB_OPTIMISATION}")
endif()

add_definitions(-DHEMELB_READING_GROUP_SIZE=${HEMELB_READING_GROUP_SIZE})
add_definitions(-DHEMELB_LATTICE=${HEMELB_LATTICE})
add_definitions(-DHEMELB_KERNEL=${HEMELB_KERNEL})
add_definitions(-DHEMELB_WALL_BOUNDARY=${HEMELB_WALL_BOUNDARY})

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

#------Capture build environment -------------

execute_process(COMMAND hg id -i RESULT_VARIABLE rev_ok OUTPUT_VARIABLE HEMELB_REVISION_NUMBER OUTPUT_STRIP_TRAILING_WHITESPACE)
if (NOT rev_ok EQUAL 0)
  message("Could not get revision number from mercurial, looking for revision_info.txt")
  execute_process(COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../revision_info.txt RESULT_VARIABLE rev_ok OUTPUT_VARIABLE HEMELB_REVISION_NUMBER OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()
message("Revision number: ${HEMELB_REVISION_NUMBER}")

execute_process(COMMAND date OUTPUT_VARIABLE HEMELB_BUILD_TIME OUTPUT_STRIP_TRAILING_WHITESPACE)
#------Dependencies --------------------------
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${HEMELB_DEPENDENCIES_PATH}/Modules/")
list(APPEND CMAKE_INCLUDE_PATH ${HEMELB_DEPENDENCIES_INSTALL_PATH}/include)
list(APPEND CMAKE_LIBRARY_PATH ${HEMELB_DEPENDENCIES_INSTALL_PATH}/lib)

#--- Set up runtime search path for DLLs -----
if(HEMELB_DEPENDENCIES_SET_RPATH)
	SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
	SET(CMAKE_INSTALL_RPATH "${HEMELB_DEPENDENCIES_INSTALL_PATH}/lib")
endif()

#-------Check for platform dependent features ----

include(CheckCXXSourceCompiles)
CHECK_CXX_SOURCE_COMPILES("#include <cmath>\n int main(int c,char** v){ return isnan(1.0); }" HAVE_ISNAN)
CHECK_CXX_SOURCE_COMPILES("#include <cmath>\n int main(int c,char** v){ return std::isnan(1.0); }" HAVE_STD_ISNAN)
CHECK_CXX_SOURCE_COMPILES("#include <sys/time.h>\n#include <sys/resource.h>\nint main(int c,char** v){ rusage usage;\ngetrusage(RUSAGE_SELF, &usage);\nreturn usage.ru_maxrss; }" HAVE_RUSAGE)

if(HAVE_ISNAN)
  add_definitions(-DHAVE_ISNAN)
endif()

if(HAVE_STD_ISNAN)
  add_definitions(-DHAVE_STD_ISNAN)
endif()

if(HAVE_RUSAGE)
  add_definitions(-DHAVE_RUSAGE)
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
if(HEMELB_BUILD_MULTISCALE)
  #------MPWide ----------------
  find_package(MPWide REQUIRED)
  include_directories(${MPWide_INCLUDE_DIR})
  add_definitions(-DHEMELB_BUILD_MULTISCALE)
endif()

#------zlib ----------------
find_package(ZLIB REQUIRED)
include_directories(${ZLIB_INCLUDE_DIR})

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

set(root_sources SimulationMaster.cc mpiInclude.cc)
add_executable(hemelb main.cc ${root_sources})

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
	topology
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
target_link_libraries(hemelb 
	${heme_libraries}
	${MPI_LIBRARIES}
	${PARMETIS_LIBRARIES}
	${TINYXML_LIBRARIES}
	${Boost_LIBRARIES}
	${CTEMPLATE_LIBRARIES}
	${ZLIB_LIBRARIES}
    ${MPWide_LIBRARIES}
	)
INSTALL(TARGETS hemelb RUNTIME DESTINATION bin)
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
	
	set(root_sources SimulationMaster.cc mpiInclude.cc multiscale/MultiscaleSimulationMaster.h)
	add_executable(hemelb_multiscale mainMultiscale.cc ${root_sources})
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
		topology
		util
		io
		log
	    colloids
		)
	foreach(subdir ${package_subdirs})
		set(lib "hemelb_${subdir}")
		list(APPEND heme_libraries ${lib})
#		add_subdirectory(${subdir})
	endforeach()
#	add_subdirectory(resources)
	target_link_libraries(hemelb_multiscale 
		${heme_libraries}
		${MPI_LIBRARIES}
		${PARMETIS_LIBRARIES}
		${TINYXML_LIBRARIES}
		${Boost_LIBRARIES}
		${CTEMPLATE_LIBRARIES}
		${ZLIB_LIBRARIES}
	    ${MPWide_LIBRARIES}
		)
	INSTALL(TARGETS hemelb_multiscale RUNTIME DESTINATION bin)
	list(APPEND RESOURCES resources/report.txt.ctp resources/report.xml.ctp)
endif()

# ----------- HEMELB unittests ---------------
if(HEMELB_BUILD_TESTS_ALL OR HEMELB_BUILD_TESTS_UNIT)
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
		${CTEMPLATE_LIBRARIES}
		${ZLIB_LIBRARIES}
		${MPWide_LIBRARIES}
		${CMAKE_DL_LIBS}) #Because on some systems CPPUNIT needs to be linked to libdl
	INSTALL(TARGETS unittests_hemelb RUNTIME DESTINATION bin)
	list(APPEND RESOURCES unittests/resources/four_cube.gmy unittests/resources/four_cube.xml unittests/resources/four_cube_multiscale.xml
		unittests/resources/config.xml unittests/resources/config0_2_0.xml
		unittests/resources/config_file_inlet.xml unittests/resources/iolet.txt
		)
endif()

# ----------- HEMELB functionaltests ---------------
if(HEMELB_BUILD_TESTS_ALL OR HEMELB_BUILD_TESTS_FUNCTIONAL)
	#------CPPUnit ---------------
	find_package(CPPUnit REQUIRED)
	include_directories(${CPPUNIT_INCLUDE_DIR})
	add_executable(functionaltests_hemelb ${root_sources})
	add_subdirectory(functionaltests/cpptests)
	target_link_libraries(functionaltests_hemelb 
		hemelb_functionaltests 
		${heme_libraries}
		${MPI_LIBRARIES}
		${PARMETIS_LIBRARIES}
		${TINYXML_LIBRARIES}
		${CPPUNIT_LIBRARY}
		${Boost_LIBRARIES}
		${CTEMPLATE_LIBRARIES}
		${MPWide_LIBRARIES}
		${ZLIB_LIBRARIES}
		${CMAKE_DL_LIBS}) #Because on some systems CPPUNIT needs to be linked to libdl
	INSTALL(TARGETS functionaltests_hemelb RUNTIME DESTINATION bin)
endif()

#-------- Copy and install resources --------------

foreach(resource ${RESOURCES})
	configure_file(${PROJECT_SOURCE_DIR}/${resource} ${BUILD_RESOURCE_PATH} COPYONLY)
	get_filename_component(resource_leaf ${resource} NAME)
	INSTALL(FILES ${BUILD_RESOURCE_PATH}/${resource_leaf} DESTINATION share/hemelb/resources)
endforeach()

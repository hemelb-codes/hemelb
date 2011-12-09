cmake_minimum_required (VERSION 2.8)
Project(HemeLBMaster)
include(ExternalProject)

#------- Conditional compile options ----------

OPTION(TIXML_USE_STL "Use STL with TIXML" ON)
if(TIXML_USE_STL)
	add_definitions(-DTIXML_USE_STL)
endif(TIXML_USE_STL)

option(USE_HEMELB_DEBUGGER "Use built in hemelb debugger" OFF)
set(HEMELB_LOG_LEVEL info
	CACHE STRING "Log level, choose 'debug', 'warning', or 'info'" )
set(HEMELB_STEERING_LIB basic
	CACHE STRING "Steering library, choose 'basic' or 'none'" )
option(USE_MULTIMACHINE "Use multi-level parallelism support" OFF)
option(HEMELB_BUILD_UNITTESTS "Build the unit-tests" ON)
option(USE_STREAKLINES "Calculate streakline images" ON)
set(OPTIMISATION "-O4" CACHE STRING "Optimisation level (can be blank or -O1 to -O4)")

#------- Dependencies -----------

add_subdirectory(dependencies)
	
# ------ HemeLB ------------------

ExternalProject_Add(
   hemelb
	URL ${PROJECT_SOURCE_DIR}/Code	
	CONFIGURE_COMMAND cmake <SOURCE_DIR> 
	-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> 
	-DDEPENDENCIES_PATH=${CMAKE_CURRENT_SOURCE_DIR}/dependencies
	-DTIXML_USE_STL=${TIXML_USE_STL}
	-DUSE_HEMELB_DEBUGGER=${USE_HEMELB_DEBUGGER}
	-DHEMELB_LOG_LEVEL=${HEMELB_LOG_LEVEL}
	-DHEMELB_STEERING_LIB=${HEMELB_STEERING_LIB}
	-DUSE_MULTIMACHINE=${USE_MULTIMACHINE}
	-DHEMELB_BUILD_UNITTESTS=${HEMELB_BUILD_UNITTESTS}
	-DUSE_STREAKLINES=${USE_STREAKLINES}
	-DOPTIMISATION=${OPTIMISATION}
)
ExternalProject_Get_Property(hemelb install_dir)
FILE(GLOB results ${install_dir}/bin/*)
INSTALL(FILES ${results} DESTINATION bin)

add_dependencies(hemelb cppunit parmetis tinyxml)
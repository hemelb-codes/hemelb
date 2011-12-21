cmake_minimum_required (VERSION 2.8)
Project(HemeLBMaster)
include(ExternalProject)

#------- Conditional compile options ----------

OPTION(TIXML_USE_STL "Use STL with TIXML" ON)
if(TIXML_USE_STL)
	add_definitions(-DTIXML_USE_STL)
endif(TIXML_USE_STL)

option(HEMELB_USE_DEBUGGER "Use built in hemelb debugger" OFF)
set(HEMELB_LOG_LEVEL info
	CACHE STRING "Log level, choose 'debug', 'warning', or 'info'" )
set(HEMELB_STEERING_LIB basic
	CACHE STRING "Steering library, choose 'basic' or 'none'" )
option(HEMELB_USE_MULTIMACHINE "Use multi-level parallelism support" OFF)
option(HEMELB_BUILD_UNITTESTS "Build the unit-tests" ON)
option(HEMELB_USE_STREAKLINES "Calculate streakline images" ON)
option(HEMELB_USE_ALL_WARNINGS_GNU "Show all compiler warnings on development builds (gnu-style-compilers)" ON)
option(HEMELB_USE_BOOST "Use Boost" OFF)
set(HEMELB_OPTIMISATION "-O4" CACHE STRING "Optimisation level (can be blank or -O1 to -O4)")
set(HEMELB_STEERING_HOST "CCS" CACHE STRING "Use a default host suffix for steering? (CCS, NGS2Leeds, NGS2Manchester, LONI, NCSA or blank)")

#------- Dependencies -----------

add_subdirectory(dependencies)
	
# ------ HemeLB ------------------
message("${PROJECT_SOURCE_DIR}/Code")
ExternalProject_Add(
   hemelb
	URL ${PROJECT_SOURCE_DIR}/Code	
	CONFIGURE_COMMAND cmake <SOURCE_DIR> 
	-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> 
	-DHEMELB_DEPENDENCIES_PATH=${CMAKE_CURRENT_SOURCE_DIR}/dependencies
	-DHEMELB_DEPENDENCIES_INSTALL_PATH=${HEMELB_DEPENDENCIES_INSTALL_PATH}
	-DTIXML_USE_STL=${TIXML_USE_STL}
	-DHEMELB_USE_DEBUGGER=${HEMELB_USE_DEBUGGER}
	-DHEMELB_LOG_LEVEL=${HEMELB_LOG_LEVEL}
	-DHEMELB_STEERING_LIB=${HEMELB_STEERING_LIB}
	-DHEMELB_USE_MULTIMACHINE=${HEMELB_USE_MULTIMACHINE}
	-DHEMELB_BUILD_UNITTESTS=${HEMELB_BUILD_UNITTESTS}
	-DHEMELB_USE_STREAKLINES=${HEMELB_USE_STREAKLINES}
	-DHEMELB_OPTIMISATION=${HEMELB_OPTIMISATION}
	-DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
	-DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
	-DMPI_C_COMPILER=${MPI_C_COMPILER}
	-DMPI_CXX_COMPILER=${MPI_CXX_COMPILER}
	-DMPI_C_NO_INTERROGATE=${MPI_C_NO_INTERROGATE}
	-DMPI_CXX_NO_INTERROGATE=${MPI_CXX_NO_INTERROGATE}
	-DHEMELB_USE_ALL_WARNINGS_GNU=${HEMELB_USE_ALL_WARNINGS_GNU}
	-DHEMELB_USE_BOOST=${HEMELB_USE_BOOST}
	-DHEMELB_STEERING_HOST=${HEMELB_STEERING_HOST}
)
ExternalProject_Get_Property(hemelb install_dir)

INSTALL(FILES ${install_dir}/bin/hemelb DESTINATION bin)
if (HEMELB_BUILD_UNITTESTS)
	INSTALL(FILES ${install_dir}/bin/unittests_hemelb DESTINATION bin)
endif()

add_dependencies(hemelb cppunit parmetis tinyxml boost ctemplate)

cmake_minimum_required (VERSION 2.8)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
project(HemeLB_Code)
# ------MPI------------------
# Require MPI for this project:
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${CMAKE_CXX_LINK_FLAGS}")
include_directories(MPI_INCLUDE_PATH)

#------Parmetis  ------------
find_package(Parmetis REQUIRED)
include_directories(${PARMETIS_INCLUDE_DIR})
MESSAGE(STATUS "Found Parmetis: ${PARMETIS_LIBRARY}")
#------CPPUnit ---------------
find_package(CPPUnit REQUIRED)
include_directories(${CPPUNIT_INCLUDE_DIR})

#------TinyXML ----------------
find_package(TinyXML REQUIRED)
add_definitions(-DTIXML_USE_STL)
include_directories(${TINYXML_INCLUDE_DIR})

# ----------- HemeLB ------------------s

add_definitions(-DHEMELB_CFG_ON_BSD -DHEMELB_CFG_ON_OSX)
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
	log)
foreach(subdir ${package_subdirs})
	set(lib "hemelb_${subdir}")
	list(APPEND heme_libraries ${lib})
	add_subdirectory(${subdir})
endforeach(${subdir})

target_link_libraries(hemelb 
	${heme_libraries}
	${MPI_LIBRARIES}
	${PARMETIS_LIBRARIES}
	${TINYXML_LIBRARIES}
	${CPPUNIT_LIBRARY}
	)

# ----------- HEMELB unittests ---------------
add_executable(unittests_hemelb ${root_sources})
add_subdirectory(unittests)
target_link_libraries(unittests_hemelb 
	hemelb_unittests 
	${heme_libraries}
	${MPI_LIBRARIES}
	${PARMETIS_LIBRARIES}
	${TINYXML_LIBRARIES}
	${CPPUNIT_LIBRARY})

INSTALL(TARGETS hemelb unittests_hemelb RUNTIME DESTINATION .)

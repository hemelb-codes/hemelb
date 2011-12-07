cmake_minimum_required (VERSION 2.8)
project(HemeLB_Code)

# ------MPI------------------
# Require MPI for this project:
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${CMAKE_CXX_LINK_FLAGS}")
include_directories(MPI_INCLUDE_PATH)
target_link_libraries(hemelb ${MPI_LIBRARIES})
target_link_libraries(unittests_hemelb )

#------Parmetis  ------------
find_package(Parmetis REQUIRED)
include_directories(${PARMETIS_INCLUDE_DIRS})

#------CPPUnit ---------------
find_package(CPPUnit REQUIRED)
include_directories(${CPPUNIT_INCLUDE_DIRS})

#------TinyXML ----------------
find_package(TinyXML REQUIRED)
include_directories(${TINYXML_INCLUDE_DIRS})

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
	${CPPUNIT_LIBRARIES}
	)

# ----------- HEMELB unittests ---------------
add_executable(unittests_hemelb ${root_sources})
add_subdirectory(unittests)
add_dependencies(hemelb_unittests ${externals})
target_link_libraries(unittests_hemelb 
	hemelb_unittests 
	${heme_libraries}
	${MPI_LIBRARIES}
	${PARMETIS_LIBRARIES}
	${TINYXML_LIBRARIES}
	${CPPUNIT_LIBRARIES})



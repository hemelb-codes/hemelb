cmake_minimum_required (VERSION 2.6)
project(HemeLB)
add_definitions(-DTIXML_USE_STL -DHEMELB_CFG_ON_BSD -DHEMELB_CFG_ON_OSX)
add_executable(hemelb SimulationMaster.cc main.cc mpiInclude.cc D3Q15.cc)
add_executable(heme_unittests SimulationMaster.cc mpiInclude.cc D3Q15.cc)
# Require MPI for this project:
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${CMAKE_CXX_LINK_FLAGS}")
include_directories(MPI_INCLUDE_PATH)
include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/parmetis/include 
	${PROJECT_SOURCE_DIR}/parmetis/metis/include)

target_link_libraries(hemelb ${MPI_LIBRARIES})
target_link_libraries(heme_unittests ${MPI_LIBRARIES})

target_link_libraries(heme_unittests cppunit)
#include_directories(unittests configuration reporting steering vis lb net debug topology xml util geometry io log parmetis)
add_subdirectory(unittests)
add_subdirectory(configuration)
add_subdirectory(reporting)
add_subdirectory(steering)
add_subdirectory(vis)
add_subdirectory(lb)
add_subdirectory(net)
add_subdirectory(debug)
add_subdirectory(topology)
add_subdirectory(xml)
add_subdirectory(util)
add_subdirectory(geometry)
add_subdirectory(io)
add_subdirectory(log)
add_library(parmetis STATIC IMPORTED)
add_library(metis STATIC IMPORTED)
set_property(TARGET parmetis PROPERTY IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/parmetis/build/Darwin-i386/libparmetis/libparmetis.a)
set_property(TARGET metis PROPERTY IMPORTED_LOCATION ${PROJECT_SOURCE_DIR}/parmetis/build/Darwin-i386/libmetis/libmetis.a)
#add_subdirectory(parmetis)
target_link_libraries(hemelb hemelb_configuration hemelb_reporting hemelb_steering 
	hemelb_vis hemelb_lb hemelb_net hemelb_debug hemelb_topology hemelb_xml hemelb_util hemelb_geometry
	hemelb_io hemelb_log parmetis metis)
target_link_libraries(heme_unittests hemelb_unittests hemelb_configuration hemelb_reporting hemelb_steering 
		hemelb_vis hemelb_lb hemelb_net hemelb_debug hemelb_topology hemelb_xml hemelb_util hemelb_geometry
	hemelb_io hemelb_log parmetis metis)


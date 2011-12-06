cmake_minimum_required (VERSION 2.8)
project(HemeLB)
add_definitions(-DTIXML_USE_STL -DHEMELB_CFG_ON_BSD -DHEMELB_CFG_ON_OSX)
add_executable(hemelb SimulationMaster.cc main.cc mpiInclude.cc D3Q15.cc)
add_executable(heme_unittests SimulationMaster.cc mpiInclude.cc D3Q15.cc)
# ------MPI------------------
# Require MPI for this project:
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${CMAKE_CXX_LINK_FLAGS}")
include_directories(MPI_INCLUDE_PATH)
target_link_libraries(hemelb ${MPI_LIBRARIES})
target_link_libraries(heme_unittests ${MPI_LIBRARIES})

# -------CPPUNIT -------------
target_link_libraries(heme_unittests cppunit)

# ------- parmetis --------------------
include(ExternalProject)
ExternalProject_Add(
   parmetis_project
	URL http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.2.tar.gz
   CONFIGURE_COMMAND make config
	BUILD_IN_SOURCE 1
	BUILD_COMMAND make
	INSTALL_COMMAND ""
)
ExternalProject_Get_Property(parmetis_project binary_dir)
ExternalProject_Get_Property(parmetis_project source_dir) 
add_library(parmetis STATIC IMPORTED)
add_library(metis STATIC IMPORTED)
add_dependencies(parmetis parmetis_project)
add_dependencies(metis parmetis_project)
include_directories(${source_dir}/include 
	${source_dir}/metis/include)
set_property(
	TARGET parmetis 
	PROPERTY IMPORTED_LOCATION ${binary_dir}/build/Darwin-i386/libparmetis/libparmetis.a)
set_property(
	TARGET metis 
	PROPERTY IMPORTED_LOCATION ${binary_dir}/build/Darwin-i386/libmetis/libmetis.a)
	
# ----------- HemeLB ------------------	
include_directories(${PROJECT_SOURCE_DIR})
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
target_link_libraries(hemelb hemelb_configuration hemelb_reporting hemelb_steering 
	hemelb_vis hemelb_lb hemelb_net hemelb_debug hemelb_topology hemelb_xml hemelb_util hemelb_geometry
	hemelb_io hemelb_log parmetis metis)
target_link_libraries(heme_unittests hemelb_unittests hemelb_configuration hemelb_reporting hemelb_steering 
		hemelb_vis hemelb_lb hemelb_net hemelb_debug hemelb_topology hemelb_xml hemelb_util hemelb_geometry
	hemelb_io hemelb_log parmetis metis)


cmake_minimum_required (VERSION 2.8)
Project(HemeLB)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/Code/cmake/Modules/")
include(ExternalProject)

# --------- tinyxml -------------------
find_package(TinyXML QUIET)
if (${TINYXML_FOUND})
	message("TinyXML already installed, no need to download")
else(${TINYXML_FOUND})
	message("TinyXML not installed, will download and build")
	ExternalProject_Add(
	   tinyxml
		URL http://sourceforge.net/projects/tinyxml/files/tinyxml/2.6.2/tinyxml_2_6_2.tar.gz
		CONFIGURE_COMMAND cmake <SOURCE_DIR> -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
		PATCH_COMMAND cp ${CMAKE_SOURCE_DIR}/Code/cmake/patches/tinyxml.cmake CMakeLists.txt
	)
	ExternalProject_Get_Property(tinyxml install_dir)
	list(APPEND libraries ${install_dir}/lib)
	list(APPEND includes ${install_dir}/include)
endif(${TINYXML_FOUND})

# ------- parmetis --------------------
find_package(Parmetis QUIET)
if (${PARMETIS_FOUND})
	message("Parmetis already installed, no need to download")
else(${PARMETIS_FOUND})
	message("Parmetis not installed, will download and build")
	ExternalProject_Add(
	   parmetis
		URL http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.2.tar.gz
	   CONFIGURE_COMMAND make config prefix=<INSTALL_DIR> && cd metis && make config prefix=<INSTALL_DIR>
		BUILD_IN_SOURCE 1
		BUILD_COMMAND make MAKEFLAGS= #Workaround for bug in parmetis makefile treating switches as targets
		INSTALL_COMMAND make install MAKEFLAGS= && cd metis && make install MAKEFLAGS=
	)
	ExternalProject_Get_Property(parmetis install_dir)
	list(APPEND libraries ${install_dir}/lib)
	list(APPEND includes ${install_dir}/include)
endif(${PARMETIS_FOUND})

# -------CPPUNIT -------------
find_package(CPPUnit QUIET)
if (${CPPUNIT_FOUND})
	message("CPPUnit already installed, no need to download")
else(${CPPUNIT_FOUND})
	message("CPPUnit not installed, will download and build")
	ExternalProject_Add(
	   cppunit
		URL http://downloads.sourceforge.net/cppunit/cppunit-1.12.1.tar.gz
		PATCH_COMMAND autoreconf -fvi #autoreconf required on osx - based on contents of portfile
	   CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --mandir=<INSTALL_DIR>/share/man --disable-doxygen --disable-dot
		BUILD_IN_SOURCE 1
	)
	ExternalProject_Get_Property(cppunit install_dir)
	list(APPEND libraries ${install_dir}/lib)
	list(APPEND includes ${install_dir}/include)
endif(${CPPUNIT_FOUND})

# ---------BOOST -----------------------
#ExternalProject_Add(
#   boost_project
#	URL http://gitorious.org/boost/cmake/archive-tarball/cmake-1.48.0.tar.gz
#)
	
# ------ HemeLB ------------------

string(REPLACE ";" "::" incpaths "${includes}")
string(REPLACE ";" "::" libpaths "${libraries}")

ExternalProject_Add(
   hemelb
	URL ${PROJECT_SOURCE_DIR}/Code
	LIST_SEPARATOR "::"
	CONFIGURE_COMMAND cmake <SOURCE_DIR> -DCMAKE_INCLUDE_PATH=${incpaths} 
		-DCMAKE_LIBRARY_PATH=${libpaths} 
		-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
)

add_dependencies(hemelb cppunit parmetis tinyxml)
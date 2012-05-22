cmake_minimum_required (VERSION 2.8)
Project(HemeLBDependencies)
include(ExternalProject)
set(HEMELB_DEPENDENCIES_PATH ${CMAKE_CURRENT_SOURCE_DIR})
set(HEMELB_DEPENDENCIES_INSTALL_PATH ${HEMELB_DEPENDENCIES_PATH} CACHE STRING "Path to install Dependencies")
set(HEMELB_SUBPROJECT_MAKE_JOBS 1 CACHE INTEGER "Number of jobs to use for subproject build steps")
set(CMAKE_MODULE_PATH "${HEMELB_DEPENDENCIES_PATH}/Modules/" ${CMAKE_MODULE_PATH})
set(CPPUNIT_CONFIGURE_OPTIONS "" CACHE STRING "Extra configuration options for CPPUNIT")
option(CPPUNIT_PATCH_LDL "Add -ldl option to cppunit" OFF)
set(CTEMPLATE_CONFIGURE_OPTIONS "" CACHE STRING "Extra configurations options for CTEMPLATE")
option(HEMELB_BUILD_MULTISCALE "Build HemeLB Multiscale functionality" OFF)
find_package(MPI REQUIRED) # Do this *after* we include the patched find MPI module

# --------- tinyxml -------------------
find_package(TinyXML QUIET)
if (TINYXML_FOUND)
	message("TinyXML already installed, no need to download")
else()
	message("TinyXML not installed, wil build from source")
	OPTION(TIXML_USE_STL "Use STL with TIXML" ON)
	if(TIXML_USE_STL)
		add_definitions(-DTIXML_USE_STL)
	endif()
	find_file(TINYXML_TARBALL tinyxml_2_6_2.tar.gz 
		DOC "Path to download TinyXML (can be url http://)"
		PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
	)
	if(NOT TINYXML_TARBALL)
		message("No tinyxml source found, will download.")
		set(TINYXML_TARBALL http://sourceforge.net/projects/tinyxml/files/tinyxml/2.6.2/tinyxml_2_6_2.tar.gz
			CACHE STRING "Path to download TinyXML (can be local file://)" FORCE)
	endif()
	ExternalProject_Add(
	   tinyxml
		INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PATH}
		URL ${TINYXML_TARBALL}
		CONFIGURE_COMMAND cmake <SOURCE_DIR>
		-DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
		-DTIXML_USE_STL=${TIXML_USE_STL}
		-DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
		-DCMAKE_CXX_FLAGS_RELEASE=${CMAKE_CXX_FLAGS_RELEASE}
		PATCH_COMMAND cp ${HEMELB_DEPENDENCIES_PATH}/patches/tinyxml.cmake CMakeLists.txt
		BUILD_COMMAND make -j${HEMELB_SUBPROJECT_MAKE_JOBS}
	)
endif()

# --------- mpwide -------------------
if (HEMELB_BUILD_MULTISCALE)
    find_package(MPWide QUIET)
    if (MPWide_FOUND)
	    message("MPWide already installed, no need to download")
    else()
	    message("MPWide not installed, will build from source")
	    find_file(MPWIDE_TARBALL MPWide-1.2b.tar.gz 
		    DOC "Path to download MPWide (can be url http://)"
		    PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
	    )
	    if(NOT MPWIDE_TARBALL)
		    message("No MPWide source found, will download.")
		    set(MPWIDE_TARBALL http://castle.strw.leidenuniv.nl/~derek/MPWide-1.2.1.tar.gz
			    CACHE STRING "Path to download MPWide (can be local file://)" FORCE)
	    endif()
	    ExternalProject_Add(
	       MPWide
		    INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PATH}
		    URL ${MPWIDE_TARBALL}
		    BUILD_IN_SOURCE 1
            CONFIGURE_COMMAND echo "No configure step for MPWIDE"
		    BUILD_COMMAND make -j${HEMELB_SUBPROJECT_MAKE_JOBS} 
            INSTALL_COMMAND make install INSTALL_PREFIX=${HEMELB_DEPENDENCIES_INSTALL_PATH} ${MWIDE_BUILD_OPTIONS}
	    )
    endif()
endif()

# ------- parmetis --------------------
find_package(Parmetis QUIET)
if (PARMETIS_FOUND)
	message("Parmetis already installed, no need to download")
else()
	message("Parmetis not installed, will build from source")
	find_file(PARMETIS_TARBALL parmetis-4.0.2.tar.gz 
		DOC "Path to download Parmetis (can be url http://)"
		PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
	)
	if(NOT PARMETIS_TARBALL)
		message("No parmetis source found, will download.")
		set(PARMETIS_TARBALL http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/parmetis-4.0.2.tar.gz 
		CACHE STRING "Path to download Parmetis (can be local file://)" FORCE)
	endif()
	set(PARMETIS_CC ${MPI_C_COMPILER} CACHE STRING "MPI Compiler to use for Parmetis, leave blank to let parmetis guess")
	set(PARMETIS_CXX ${MPI_CXX_COMPILER} CACHE STRING "MPI Compiler to use for Parmetis, leave blank to let parmetis guess")
	if(PARMETIS_CC)
		set(PARMETIS_CC_OPTION cc=${PARMETIS_CC})
	endif()
	if(PARMETIS_CXX)
		set(PARMETIS_CXX_OPTION cxx=${PARMETIS_CXX})
	endif()
	ExternalProject_Add(
	   parmetis
		INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PATH}
		URL ${PARMETIS_TARBALL}
	   CONFIGURE_COMMAND make config prefix=<INSTALL_DIR> ${PARMETIS_CC_OPTION} ${PARMETIS_CXX_OPTION} && 
				cd metis && 
				make config prefix=<INSTALL_DIR> ${PARMETIS_CC_OPTION} ${PARMETIS_CXX_OPTION}
		BUILD_IN_SOURCE 1
		BUILD_COMMAND make -j${HEMELB_SUBPROJECT_MAKE_JOBS} MAKEFLAGS= #Workaround for bug in parmetis makefile treating switches as targets
		INSTALL_COMMAND make install MAKEFLAGS= && cd metis && make install MAKEFLAGS=
	)
endif()

# -------CPPUNIT -------------
find_package(CPPUnit QUIET)
if (CPPUNIT_FOUND)
	message("CPPUnit already installed, no need to download")
else()
	message("CPPUnit not installed, will build from source")
	find_file(CPPUNIT_TARBALL cppunit-1.12.1.tar.gz 
		DOC "Path to download CPPUNIT (can be url http://)"
		PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
	)
	if(NOT CPPUNIT_TARBALL)
		message("No cppunit source found, will download.")
		set(CPPUNIT_TARBALL http://downloads.sourceforge.net/cppunit/cppunit-1.12.1.tar.gz
			CACHE STRING "Path to download CPPUNIT (can be local file://)" FORCE)
	endif()
	if(APPLE)
		set(RECONFIGURE_CPPUNIT_DEFAULT ON)
	else()
		set(RECONFIGURE_CPPUNIT_DEFAULT OFF)
	endif()
	option(RECONFIGURE_CPPUNIT "Include Reconfigure step for CPPUNIT" ${RECONFIGURE_CPPUNIT_DEFAULT})
	if (CPPUNIT_PATCH_LDL)
    	set(PATCH_COMMAND_LDL patch -p1 < ${HEMELB_DEPENDENCIES_PATH}/patches/cppunit_ldl.diff)
    else()
        set(PATCH_COMMAND_LDL echo noldl)
    endif()
    if (RECONFIGURE_CPPUNIT)
    	set(PATCH_COMMAND_RECONFIGURE autoreconf -fvi) #autoreconf required on osx - based on contents of portfile)
    else()
        set(PATCH_COMMAND_RECONFIGURE echo noreconf)
    endif()
	ExternalProject_Add(
	   cppunit
		INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PATH}
		URL ${CPPUNIT_TARBALL}
	   CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR> --mandir=<INSTALL_DIR>/share/man --disable-doxygen --disable-dot "${CPPUNIT_CONFIGURE_OPTIONS}"
		BUILD_IN_SOURCE 1
		BUILD_COMMAND make -j${HEMELB_SUBPROJECT_MAKE_JOBS}
		PATCH_COMMAND ${PATCH_COMMAND_LDL} && ${PATCH_COMMAND_RECONFIGURE}
	)
endif()

if(HEMELB_USE_BOOST)
	# ---------BOOST -----------------------
	set(Boost_ADDITIONAL_VERSIONS "1.48" "1.48.0")
	set(BOOST_LIBS mpi regex)
	string(REPLACE ";" "," BOOST_LIBS_COMMA "${BOOST_LIBS}")
	find_package(Boost 1.48 COMPONENTS ${BOOST_LIBS})
	if (BOOST_FOUND)
		message("Boost 1.48 already installed, no need to download")
	else()
		message("Boost not installed, will build from source")
		find_file(BOOST_TARBALL boost_1_48_0.tar.gz
			DOC "Path to download BOOST (can be url http://)"
			PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
		)
		if(NOT BOOST_TARBALL)
			message("No boost source found, will download.")
			set(BOOST_TARBALL http://sourceforge.net/projects/boost/files/boost/1.48.0/boost_1_48_0.tar.gz
				CACHE STRING "Path to download BOOST (can be local file://)" FORCE)
		endif()
		ExternalProject_Add(
	   	boost
			INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PATH}
			URL ${BOOST_TARBALL}
			BUILD_COMMAND <SOURCE_DIR>/b2
			INSTALL_COMMAND <SOURCE_DIR>/b2 install
			CONFIGURE_COMMAND <SOURCE_DIR>/bootstrap.sh --prefix=<INSTALL_DIR> "--with-libraries=${BOOST_LIBS_COMMA}"
			BUILD_IN_SOURCE 1
		)
		ExternalProject_Add_Step(
			boost mpi
			COMMAND cat ${HEMELB_DEPENDENCIES_PATH}/patches/boost_options.jam >> project-config.jam
			DEPENDEES configure
			DEPENDERS build
			WORKING_DIRECTORY <SOURCE_DIR>
		)
	endif()
endif()

#----------------Google C++ template ----------------

option(CTEMPLATE_PATCH_VACOPY "Define va_copy macro through patch" OFF)
option(CTEMPLATE_PATCH_ALIGN "Define GTL align macros as gnu" OFF)
find_package(CTemplate)
if (CTEMPLATE_FOUND)
	message("Ctemplate already installed, no need to download")
else()
	message("Ctemplate not installed, will build from source")
	find_file(CTEMPLATE_TARBALL ctemplate-2.2.tar.gz
		DOC "Path to download CTemplate (can be url http://)"
		PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
	)
	if(NOT CTEMPLATE_TARBALL)
		message("No ctemplate source found, will download.")
		set(CTEMPLATE_TARBALL http://ctemplate.googlecode.com/files/ctemplate-2.2.tar.gz
			CACHE STRING "Path to download CTemplate (can be local file://)" FORCE)
	endif()
	if (CTEMPLATE_PATCH_VACOPY)
    	set(PATCH_COMMAND_VACOPY patch -p1 < ${HEMELB_DEPENDENCIES_PATH}/patches/ctemplate_vacopy.diff)
    else()
        set(PATCH_COMMAND_VACOPY echo novacopy)
    endif()
    if (PATCH_ALIGN)
    	set(PATCH_COMMAND_ALIGN patch -p1 < ${HEMELB_DEPENDENCIES_PATH}/patches/ctemplate_align.diff)
    else()
        set(PATCH_COMMAND_ALIGN echo noreconf)
    endif()
	ExternalProject_Add(
   	ctemplate
		INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PATH}
		URL ${CTEMPLATE_TARBALL}
	   CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>  ${CTEMPLATE_CONFIGURE_OPTIONS}
		BUILD_COMMAND make -j${HEMELB_SUBPROJECT_MAKE_JOBS}
		BUILD_IN_SOURCE 1
		PATCH_COMMAND ${PATCH_COMMAND_ALIGN} && ${PATCH_COMMAND_VACOPY}
	)
endif()

#-----------------ZLIB --------------------------------
#http://sourceforge.net/projects/libpng/files/zlib/1.2.6/zlib-1.2.6.tar.gz
find_package(ZLIB)
if (ZLIB_FOUND)
	message("ZLIB already installed, no need to download")
else()
	message("ZLIB not installed, will build from source")
	find_file(ZLIB_TARBALL zlib-1.2.6.tar.gz
		DOC "Path to download ZLIB (can be url http://)"
		PATHS ${HEMELB_DEPENDENCIES_PATH}/distributions
	)
	if(NOT ZLIB_TARBALL)
		message("No ctemplate source found, will download.")
		set(ZLIB_TARBALL http://sourceforge.net/projects/libpng/files/zlib/1.2.6/zlib-1.2.6.tar.gz
			CACHE STRING "Path to download ZLIB (can be local file://)" FORCE)
	endif()
	ExternalProject_Add(
   	zlib
		INSTALL_DIR ${HEMELB_DEPENDENCIES_INSTALL_PATH}
		URL ${ZLIB_TARBALL}
	   CONFIGURE_COMMAND <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
		BUILD_COMMAND make -j${HEMELB_SUBPROJECT_MAKE_JOBS}
		BUILD_IN_SOURCE 1
	)
endif()

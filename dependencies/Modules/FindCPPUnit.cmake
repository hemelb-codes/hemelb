#
# http://root.cern.ch/viewvc/trunk/cint/reflex/cmake/modules/FindCppUnit.cmake
#
# - Find CppUnit
# This module finds an installed CppUnit package.
#
# It sets the following variables:
#  CPPUNIT_FOUND       - Set to false, or undefined, if CppUnit isn't found.
#  CPPUNIT_INCLUDE_DIR - The CppUnit include directory.
#  CPPUNIT_LIBRARY     - The CppUnit library to link against.

FIND_PATH(CPPUNIT_INCLUDE_DIR cppunit/Test.h)
option(CPPUNIT_USE_STATIC "Prefer Static CPPUNIT library" OFF)
if(CPPUNIT_USE_STATIC)
set(__old_cmake_find_lib_suffixes ${CMAKE_FIND_LIBRARY_SUFFIXES})
	set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()
FIND_LIBRARY(CPPUNIT_LIBRARY NAMES cppunit)
if(CPPUNIT_USE_STATIC)
	set(CMAKE_FIND_LIBRARY_SUFFIXES ${__old_cmake_find_lib_suffixes})
endif()

IF (CPPUNIT_INCLUDE_DIR AND CPPUNIT_LIBRARY)
   SET(CPPUNIT_FOUND TRUE)
ENDIF (CPPUNIT_INCLUDE_DIR AND CPPUNIT_LIBRARY)

IF (CPPUNIT_FOUND)

   # show which CppUnit was found only if not quiet
   IF (NOT CppUnit_FIND_QUIETLY)
      MESSAGE(STATUS "Found CppUnit: ${CPPUNIT_LIBRARY}")
   ENDIF (NOT CppUnit_FIND_QUIETLY)

ELSE (CPPUNIT_FOUND)

   # fatal error if CppUnit is required but not found
   IF (CppUnit_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find CppUnit")
   ENDIF (CppUnit_FIND_REQUIRED)

ENDIF (CPPUNIT_FOUND)
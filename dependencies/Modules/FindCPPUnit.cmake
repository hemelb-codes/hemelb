# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

include_guard()

find_package(PkgConfig QUIET)
if (PkgConfig_FOUND)
  pkg_check_modules(_CPPUnit cppunit)
  if (_CPPUnit_FOUND)
    # Forward consistently
    set(CPPUnit_FOUND "${_CPPUnit_FOUND}")
    set(CPPUnit_INCLUDE_DIR "${_CPPUnit_INCLUDE_DIRS}")
    set(CPPUnit_LIBRARY "${_CPPUnit_LINK_LIBRARIES}")
    add_library(CPPUnit::CPPUnit INTERFACE IMPORTED)
    set_target_properties(CPPUnit::CPPUnit PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${_CPPUnit_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES "${_CPPUnit_LINK_LIBRARIES}"
      )
  endif()
endif()

if (NOT PkgConfig_FOUND OR NOT _CPPUnit_FOUND)
  find_path(CPPUnit_INCLUDE_DIR cppunit/Test.h)

  option(CPPUnit_USE_STATIC "Prefer Static CPPUnit library" OFF)
  if(CPPUnit_USE_STATIC)
    set(__old_cmake_find_lib_suffixes ${CMAKE_FIND_LIBRARY_SUFFIXES})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()

  find_library(CPPUnit_LIBRARY NAMES cppunit)
  if(CPPUnit_USE_STATIC)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${__old_cmake_find_lib_suffixes})
  endif()

  add_library(CPPUnit::CPPUnit INTERFACE IMPORTED)
  set_target_properties(CPPUnit::CPPUnit PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${CPPUnit_INCLUDE_DIR}"
    INTERFACE_LINK_LIBRARIES "${CPPUnit_LIBRARY}"
    )

endif()

# handle the QUIETLY and REQUIRED arguments
include( "FindPackageHandleStandardArgs" )
find_package_handle_standard_args(
  "CPPUnit"
  DEFAULT_MSG
  CPPUnit_INCLUDE_DIR CPPUnit_LIBRARY
  )


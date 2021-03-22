# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

cmake_minimum_required(VERSION 3.2)

project(HemeLbSetupTool)

macro(dvar name)
  message("${name} = ${${name}}")
endmacro()

get_filename_component(HEMELB_CODE_DIR ${PROJECT_SOURCE_DIR}/../../Code ABSOLUTE)

find_package(Boost REQUIRED COMPONENTS system filesystem)
find_package(VTK REQUIRED COMPONENTS vtkCommonCore)
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})
find_package(HDF5 REQUIRED)

# This macro does list append but ensures the paths are canonicalised
# so that the remove duplicates works (badly done
# FindAnonymousPackage.cmake)
macro(append listname)
  foreach(dirname ${ARGN})
    get_filename_component(canonical ${dirname} ABSOLUTE)
    list(APPEND ${listname} ${canonical})
  endforeach()
endmacro()

append(SETUPTOOL_INCLUDE_DIRS ${HEMELB_CODE_DIR})
# HemeLB objects compiled 

append(SETUPTOOL_INCLUDE_DIRS ${Boost_INCLUDE_DIRS})
append(SETUPTOOL_LIBRARY_DIRS ${Boost_LIBRARY_DIRS})
append(SETUPTOOL_LIBRARIES ${Boost_LIBRARIES})

append(SETUPTOOL_INCLUDE_DIRS ${VTK_INCLUDE_DIRS})
# VTK libraries covered by this .so being loaded by the interpreter
# which has already loaded that module

append(SETUPTOOL_INCLUDE_DIRS ${HDF5_INCLUDE_DIRS})
append(SETUPTOOL_LIBRARY_DIRS ${HDF5_LIBRARY_DIRS})
append(SETUPTOOL_LIBRARIES ${HDF5_LIBRARIES})

append(SETUPTOOL_INCLUDE_DIRS ${CGAL_INCLUDE_DIRS})
append(SETUPTOOL_LIBRARY_DIRS ${CGAL_LIBRARIES_DIRS})
append(SETUPTOOL_LIBRARIES ${CGAL_LIBRARIES})
foreach(_cgal_lib ${CGAL_ESSENTIAL_3RD_PARTY_LIBRARIES})
  append(SETUPTOOL_LIBRARIES ${${_cgal_lib}_LIBRARIES})
endforeach()


list(REMOVE_DUPLICATES SETUPTOOL_INCLUDE_DIRS)
list(REMOVE_DUPLICATES SETUPTOOL_LIBRARY_DIRS)
list(REMOVE_DUPLICATES SETUPTOOL_LIBRARIES)

# Python distutils expect colon separated path lists for include_dirs
# and library_dirs
string(REPLACE ";" ":" SETUPTOOL_INCLUDE_DIRS "${SETUPTOOL_INCLUDE_DIRS}")
string(REPLACE ";" ":" SETUPTOOL_LIBRARY_DIRS "${SETUPTOOL_LIBRARY_DIRS}")

# Because CMake loves full paths and distutils expects a list of
# search dirs and library names suitable for the usual compiler
# arguments -L COLON:SEPARATED:PATHS and -lfoo -lbar, we just give the
# libraries to the link_objects argument instead. Just to be difficult
# distutils expects comma separated list here.
string(REPLACE ";" "," SETUPTOOL_LIBRARIES "${SETUPTOOL_LIBRARIES}")

configure_file(setup.cfg.in setup.cfg)

# find Zoltan includes and libraries
#
# http://www.cs.sandia.gov/Zoltan/Zoltan.html
#
# ZOLTAN_INCLUDE_DIR_C - where to find zoltan.h
# ZOLTAN_INCLUDE_DIRS - List of include dir.
# ZOLTAN_LIBRARIES    - List of fully qualified libraries to link against.
# ZOLTAN_FOUND        - Do not attempt to use if "no" or undefined.

# Does NOT work
# Try to find Zoltan using Trilinos recommendations
#if( DEFINED TRILINOS_HOME AND NOT DEFINED ZOLTAN_HOME )
#    set( ZOLTAN_HOME ${TRILINOS_HOME} )
#endif()
#
#if( DEFINED ZOLTAN_HOME )
#    find_package(Zoltan PATHS ${ZOLTAN_HOME}/lib/cmake/Zoltan ${ZOLTAN_HOME}/include )
#endif()
#
#if( DEFINED DEPS_ROOT )
#    find_package(Zoltan PATHS ${DEPS_ROOT}/lib/cmake/Zoltan ${DEPS_ROOT}/include )
#endif()
#
#if(Zoltan_FOUND)
#
#  set( ZOLTAN_INCLUDE_DIRS "" )
#
#  list( APPEND ZOLTAN_INCLUDE_DIRS ${Zoltan_INCLUDE_DIRS})
#  list( APPEND ZOLTAN_INCLUDE_DIRS ${Zoltan_TPL_INCLUDE_DIRS})
#
#  set( ZOLTAN_LIBRARIES "" )
#
#  foreach( test_lib ${Zoltan_LIBRARIES} )
#    find_library( ${test_lib}_lib ${test_lib} PATHS ${Zoltan_LIBRARY_DIRS} NO_DEFAULT_PATH)
#    find_library( ${test_lib}_lib ${test_lib})
#    mark_as_advanced( ${test_lib}_lib )
#    list( APPEND ZOLTAN_LIBRARIES ${${test_lib}_lib} )
#  endforeach()
#
#  list( APPEND ZOLTAN_LIBRARIES ${Zoltan_TPL_LIBRARIES} )
#
#endif()

find_path(ZOLTAN_INCLUDE_DIR_CXX zoltan_cpp.h)
find_path(ZOLTAN_INCLUDE_DIR_C zoltan.h)

set(ZOLTAN_INCLUDE_DIRS ${ZOLTAN_INCLUDE_DIR_CXX} ${ZOLTAN_INCLUDE_DIR_C})

find_library(ZOLTAN_LIB zoltan)

set(ZOLTAN_LIBRARIES ${ZOLTAN_LIB})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Zoltan DEFAULT_MSG
                                  ZOLTAN_INCLUDE_DIR_C
                                  ZOLTAN_INCLUDE_DIRS
                                  ZOLTAN_LIB
                                  ZOLTAN_LIBRARIES)

# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

include_guard()

# Get all properties that cmake supports
if(NOT CMAKE_PROPERTY_LIST)
  execute_process(COMMAND cmake --help-property-list OUTPUT_VARIABLE CMAKE_PROPERTY_LIST)
  # Convert command output into a CMake list
  string(REGEX REPLACE ";" "\\\\;" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
  string(REGEX REPLACE "\n" ";" CMAKE_PROPERTY_LIST "${CMAKE_PROPERTY_LIST}")
endif()

# Debug print a variable nicely
macro(dvar name)
  message("${name} = '${${name}}'")
endmacro()

function(dvar_re regex)
  get_cmake_property(all_var_names VARIABLES)
  list (SORT all_var_names)
  foreach (varname ${all_var_names})
    string(REGEX MATCH "${regex}" matched "${varname}")
    if (matched)
      dvar(${varname})
    endif()
  endforeach()
endfunction()

function(dprop tgt)
  foreach(pname ${ARGN})
    get_target_property(pval ${tgt} ${pname})
    message("${tgt}.${pname} = '${pval}'")
  endforeach()
endfunction()

function(dprop_all target)
  if(NOT TARGET ${target})
    message(STATUS "There is no target named '${target}'")
    return()
  endif()

  foreach(property ${CMAKE_PROPERTY_LIST})
    string(REPLACE "<CONFIG>" "${CMAKE_BUILD_TYPE}" property ${property})

    # Fix https://stackoverflow.com/questions/32197663/how-can-i-remove-the-the-location-property-may-not-be-read-from-target-error-i
    if(property STREQUAL "LOCATION" OR property MATCHES "^LOCATION_" OR property MATCHES "_LOCATION$")
      continue()
    endif()

    get_property(was_set TARGET ${target} PROPERTY ${property} SET)
    if(was_set)
      get_target_property(value ${target} ${property})
      message("${target}.${property} = '${value}'")
    endif()
  endforeach()
endfunction()

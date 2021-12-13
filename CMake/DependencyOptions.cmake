# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

include(PassOptions)

pass_group_declare(DEPS)

file(GLOB _dep_findcmake_list RELATIVE "${HEMELB_ROOT_DIR}/dependencies" */find.cmake)
foreach(_dep_findcmake ${_dep_findcmake_list})
  string(REPLACE "/find.cmake" "" _dep_name "${_dep_findcmake}")
  string(TOUPPER ${_dep_name} _up_name)
  pass_cachevar_choice(DEPS "DEPS_${_up_name}" Auto
  STRING "How to get ${_dep_name}. System: use an existing one. Build: build ourselves. Auto: try a system one and fall back to building." System;Auto;Build)
endforeach()
  

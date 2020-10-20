# This file is part of HemeLB and is Copyright (C)
# the HemeLB team and/or their institutions, as detailed in the
# file AUTHORS. This software is provided under the terms of the
# license in the file LICENSE.

if(HEMELB_USE_DEBUGGER)
  # We have been asked to compile a debugger helper

  # Different platforms need different things
  if(APPLE)
    set(debugger_implementation OSX/OsxDebugger.cc common/ActiveDebugger.cc)
    set(HEMELB_DEBUGGER_IMPLEMENTATION OSX)

  elseif(UNIX)
    set(debugger_implementation linux/LinuxDebugger.cc common/ActiveDebugger.cc)
    set(HEMELB_DEBUGGER_IMPLEMENTATION linux)

  else()
    message("No hemelb built-in debugger available for system ${CMAKE_SYSTEM_NAME}")
    set(debugger_implementation none/NullDebugger.cc)
    set(HEMELB_DEBUGGER_IMPLEMENTATION none)
  endif()

else()
  # We have not been asked for a debugger; compile the no-op one.
  set(debugger_implementation none/NullDebugger.cc)
  set(HEMELB_DEBUGGER_IMPLEMENTATION none)

endif()

# Set a var for use in PlatformDebugger.h
set(HEMELB_DEBUGGER_IMPLEMENTATION_${HEMELB_DEBUGGER_IMPLEMENTATION} ON)

# This file gets configured by CMake to include the platform-specific
# debugger header and alias it to a common name.
configure_file (
  "${PROJECT_SOURCE_DIR}/debug/PlatformDebugger.h.in"
  "${PROJECT_BINARY_DIR}/debug/PlatformDebugger.h"
  )

add_library(hemelb_debug Debugger.cc none/NullDebugger.cc ${debugger_implementation})

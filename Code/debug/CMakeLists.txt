
if(HEMELB_BUILD_DEBUGGER)
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
	set(HEMELB_DEBUGGER_IMPLEMENTATION none)	
endif()
set(HEMELB_DEBUGGER_IMPLEMENTATION_${HEMELB_DEBUGGER_IMPLEMENTATION} ON)

configure_file (
  "${PROJECT_SOURCE_DIR}/debug/PlatformDebugger.h.in"
  "${PROJECT_BINARY_DIR}/debug/PlatformDebugger.h"
  )

add_library(hemelb_debug Debugger.cc none/NullDebugger.cc ${debugger_implementation})
#No debugging for now - pending platform dependence

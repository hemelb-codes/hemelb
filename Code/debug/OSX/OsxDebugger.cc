#include <string>

#include "debug/OSX/OsxDebugger.h"

namespace hemelb
{
  namespace debug
  {
    OsxDebugger::OsxDebugger(char* executable) : ActiveDebugger(executable) {}
    
    const std::string OsxDebugger::GetPlatformInterpreter(void) const {
      return std::string("osascript");
    }
    
    const std::string OsxDebugger::GetPlatformScript(void) const {
      std::string include (__FILE__);
      std::string debugOsxDir = include.substr(0, include.rfind('/'));
      
      return debugOsxDir + "/MPIdebug.scpt";
    }
    
    const std::string OsxDebugger::GetPlatformGdbScript(void) const {
      std::string include (__FILE__);
      std::string debugOsxDir = include.substr(0, include.rfind('/'));
      
      return debugOsxDir + "/resume.gdb";
    }

    Debugger* PlatformDebuggerFactory(char *executable) {
      return new OsxDebugger(executable);
    }
 
  } // namespace debug
} // namespace hemelb

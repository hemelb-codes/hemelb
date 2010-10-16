#include <string>

#include "dbg/OSX/OsxDebugger.h"

namespace hemelb
{
  namespace dbg
  {
    OsxDebugger::OsxDebugger(char* executable) : ActiveDebugger(executable) {}
    
    const std::string OsxDebugger::GetPlatformInterpreter(void) const {
      return std::string("osascript");
    }
    
    const std::string OsxDebugger::GetPlatformScript(void) const {
      std::string include (__FILE__);
      std::string dbgOsxDir = include.substr(0, include.rfind('/'));
      
      return dbgOsxDir + "/MPIdebug.scpt";
    }
    
    Debugger* PlatformDebuggerFactory(char *executable) {
      return new OsxDebugger(executable);
    }
 
  }
}

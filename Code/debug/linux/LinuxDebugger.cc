#include <string>

#include "debug/linux/LinuxDebugger.h"

namespace hemelb
{
  namespace debug
  {
    LinuxDebugger::LinuxDebugger(char* executable) :
      ActiveDebugger(executable) {}
    
    const std::string LinuxDebugger::GetPlatformInterpreter(void) const {
      return std::string("bash");
    }
    
    const std::string LinuxDebugger::GetPlatformScript(void) const {
      std::string include (__FILE__);
      std::string debugLinuxDir = include.substr(0, include.rfind('/'));
      
      return debugLinuxDir + "/launchGdbs.sh";
    }
    
    Debugger* PlatformDebuggerFactory(char *executable) {
      return new LinuxDebugger(executable);
    }
 
  } // namespace debug
} // namespace hemelb

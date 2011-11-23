#include <string>

#include "debug/linux/LinuxDebugger.h"

namespace hemelb
{
  namespace debug
  {
    LinuxDebugger::LinuxDebugger(const char* const executable) :
      ActiveDebugger(executable) {}
    
    const std::string LinuxDebugger::GetPlatformInterpreter(void) const {
      return std::string("bash");
    }
    
    const std::string LinuxDebugger::GetPlatformScript(void) const {
      std::string include (__FILE__);
      std::string debugLinuxDir = include.substr(0, include.rfind('/'));
      
      return debugLinuxDir + "/launchGdbs.sh";
    }
    
    const std::string LinuxDebugger::GetPlatformGdbScript(void) const {
      std::string include (__FILE__);
      std::string debugLinuxDir = include.substr(0, include.rfind('/'));
      
      return debugLinuxDir + "/resume.gdb";
    }
    
    Debugger* PlatformDebuggerFactory(const char * const executable) {
      return new LinuxDebugger(executable);
    }
 
  } // namespace debug
} // namespace hemelb

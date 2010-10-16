#ifndef HEMELB_DEBUG_OSX_OSXDEBUGGER_H
#define HEMELB_DEBUG_OSX_OSXDEBUGGER_H

#include <debug/common/ActiveDebugger.h>

namespace hemelb
{
  namespace debug
  {
    
    class OsxDebugger : public ActiveDebugger {
    protected:
      // Platform specific getters
      const std::string GetPlatformInterpreter(void) const;
      const std::string GetPlatformScript(void) const;
      
      // C'tor...
      OsxDebugger(char* executable);
      // ... which the factory function needs to be able to get at.
      friend Debugger* PlatformDebuggerFactory(char* executable);  
      
    };
    
    // Factory. Don't be calling this.
    Debugger* PlatformDebuggerFactory(char* executable);
  }
}

#endif // HEMELB_DEBUG_OSX_OSXDEBUGGER_H

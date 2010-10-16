#include <iostream>
#include "dbg/none/NullDebugger.h"

namespace hemelb
{
  namespace dbg
  {
    // Do nothing!
    void NullDebugger::Attach(void) {}
    void NullDebugger::BreakHere(void){}
    NullDebugger::NullDebugger(char* executable) {
      mExecutable = std::string(executable);
    }
    
    Debugger* PlatformDebuggerFactory(char *executable) {
      std::cerr << "Construct OsxDebugger" << std::endl;
      return new NullDebugger(executable);
    }
  }
}

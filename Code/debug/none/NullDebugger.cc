#include <iostream>
#include "debug/none/NullDebugger.h"

namespace hemelb
{
  namespace debug
  {
    // Do nothing!
    void NullDebugger::Attach(void)
    {
    }
    void NullDebugger::BreakHere(void)
    {
    }
    void NullDebugger::Print(const char* iFormat, ...)
    {
    }

    NullDebugger::NullDebugger(const char* const executable) :
      Debugger(executable)
    {
    }

    Debugger* PlatformDebuggerFactory(const char * const executable)
    {
      return new NullDebugger(executable);
    }
  }
}

#ifndef HEMELB_DEBUG_NONE_DEBUGGER_H
#define HEMELB_DEBUG_NONE_DEBUGGER_H

#include "debug/Debugger.h"

namespace hemelb
{
  namespace debug
  {
    class NullDebugger : public Debugger
    {
      public:
        void BreakHere(void);
        void Print(const char* iFormat, ...);

      protected:
        void Attach(void);
        NullDebugger(const char* const executable);
        friend Debugger* PlatformDebuggerFactory(const char* const executable);
    };

    Debugger* PlatformDebuggerFactory(const char* const executable);
  }
}

#endif // HEMELB_DEBUG_NONE_DEBUGGER_H

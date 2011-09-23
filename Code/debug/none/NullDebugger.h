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
        NullDebugger(char* executable);
        friend Debugger* PlatformDebuggerFactory(char* executable);
      };

    Debugger* PlatformDebuggerFactory(char* executable);
  }
}

#endif // HEMELB_DEBUG_NONE_DEBUGGER_H

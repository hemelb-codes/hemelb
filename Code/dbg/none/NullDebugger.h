#ifndef HEMELB_DBG_NONE_DEBUGGER_H
#define HEMELB_DBG_NONE_DEBUGGER_H

#include "dbg/Debugger.h"

namespace hemelb
{
  namespace dbg
  {
    class NullDebugger : public Debugger {
    public:
       void BreakHere(void);
      
    protected:
      void Attach(void);
      NullDebugger(char* executable);
      friend Debugger* PlatformDebuggerFactory(char* executable);  
    };
    
    Debugger* PlatformDebuggerFactory(char* executable);  
  }
}

#endif // HEMELB_DBG_NONE_DEBUGGER_H

#ifndef HEMELB_DEBUG_DEBUGGER_H
#define HEMELB_DEBUG_DEBUGGER_H

#include <string>

namespace hemelb
{
  namespace debug
  {

    class Debugger {
      /* Interface for debugger control.
	 
       */
    public:
      // the singleton pattern
      static Debugger* Init(char *);
      static Debugger* Get(void);
      
      virtual void BreakHere(void) = 0;

    protected:
      Debugger(char* executable);
      virtual ~Debugger();
      
      virtual void Attach() = 0;
      
      std::string mExecutable;
      
      // Singleton pattern
      static bool isSingletonCreated;
      static Debugger* singleton;
      
    };
    
    Debugger* PlatformDebuggerFactory(char* executable);
    
    
  }
}

#endif // HEMELB_DEBUG_DEBUGGER_H

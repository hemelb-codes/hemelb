#include "dbg/Debugger.h"
#include "dbg/platform.h"

namespace hemelb
{
  namespace dbg
  {
    
    Debugger* Debugger::Init(char *executable) {
      /* Static member function that implements the singleton pattern.
       * Use the namespace function PlatformDebuggerFactory to
       * actually construct the instance. It should be defined in the
       * appropriate platform subdirectory.
       */
      if (!Debugger::isSingletonCreated) {
	Debugger::singleton =  PlatformDebuggerFactory(executable);
      }
      Debugger::singleton->Attach();
    }
    
    Debugger* Debugger::Get(void) {
      // Get the single instance.
      return Debugger::singleton;
    }
    
    // Init static members
    bool Debugger::isSingletonCreated = false;
    Debugger* Debugger::singleton = 0;
    
    Debugger::Debugger(char* executable) {
      // Ctor
      mExecutable = std::string(executable);
    }
    
    // Dtor
    Debugger::~Debugger() {}
    
  } // namespace dbg
} // namespace hemelb

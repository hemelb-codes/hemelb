#include "debug/Debugger.h"
#include "debug/platform.h"

namespace hemelb
{
  namespace debug
  {

    Debugger* Debugger::Init(char *executable)
    {
      /* Static member function that implements the singleton pattern.
       * Use the namespace function PlatformDebuggerFactory to
       * actually construct the instance. It should be defined in the
       * appropriate platform subdirectory.
       */
      if (Debugger::singleton == NULL)
      {
        Debugger::singleton = PlatformDebuggerFactory(executable);
      }
      Debugger::singleton->Attach();
      return Debugger::singleton;
    }

    Debugger* Debugger::Get(void)
    {
      // Get the single instance.
      return Debugger::singleton;
    }

    // Init static members
    Debugger* Debugger::singleton = NULL;

    Debugger::Debugger(char* executable)
    {
      // Ctor
      mExecutable = std::string(executable);
    }

    // Dtor
    Debugger::~Debugger()
    {
    }

  } // namespace debug
} // namespace hemelb

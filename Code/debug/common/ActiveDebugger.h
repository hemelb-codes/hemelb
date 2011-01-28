#ifndef HEMELB_DEBUG_ACTIVEDEBUGGER_H
#define HEMELB_DEBUG_ACTIVEDEBUGGER_H

#include <vector>

#include <debug/Debugger.h>

namespace hemelb
{
  namespace debug
  {
    
    class ActiveDebugger : public Debugger {
      /* This class does the bulk of the pause/attach process. It
       * still has some pure virtual functions that must be overriden
       * to give platform-specific data.
       *
       */
      
    public:
      void BreakHere(void);
      void Print(const char* iFormat, ...);
      
    protected:
      typedef std::vector<int> VoI; // Vector of Ints
      typedef std::vector<std::string> VoS; // Vector of Strings
      
      // C'tor
      ActiveDebugger(char* executable);
      
      bool mAmAttached; // Indicate attachment state
      VoI* mPIds; // vector of process IDs
      
      void Attach(void);
      void GatherProcessIds(void);
      void SpawnDebuggers(void);
      
      // Platform specific stuff
      virtual const std::string GetPlatformInterpreter(void) const = 0;
      virtual const std::string GetPlatformScript(void) const = 0;
      virtual const std::string GetPlatformGdbScript(void) const = 0;
      
      // Helper function
      static std::string ConvertIntToString(int i);

    };
    
  } // namespace debug
} // namespace hemelb

#endif // HEMELB_DEBUG_ACTIVE_DEBUGGER_H

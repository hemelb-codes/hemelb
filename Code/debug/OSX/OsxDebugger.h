// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_DEBUG_OSX_OSXDEBUGGER_H
#define HEMELB_DEBUG_OSX_OSXDEBUGGER_H

#include <debug/common/ActiveDebugger.h>

namespace hemelb
{
  namespace debug
  {

    class OsxDebugger : public ActiveDebugger
    {
      protected:
        // Platform specific getters
        const std::string GetBinaryPath(void) const;
        const std::string GetPlatformInterpreter(void) const;
        const std::string GetPlatformScript(void) const;

        // C'tor...
        OsxDebugger(const char* const executable, const net::MpiCommunicator& comm);
        // ... which the factory function needs to be able to get at.
        friend class Debugger;

    };

  }
}

#endif // HEMELB_DEBUG_OSX_OSXDEBUGGER_H

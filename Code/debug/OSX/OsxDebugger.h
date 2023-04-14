// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_DEBUG_OSX_OSXDEBUGGER_H
#define HEMELB_DEBUG_OSX_OSXDEBUGGER_H

#include <debug/common/ActiveDebugger.h>

namespace hemelb::debug
{

    class OsxDebugger : public ActiveDebugger
    {
      protected:
        // Platform specific getters
        [[nodiscard]] std::string GetBinaryPath() const override;
        [[nodiscard]] std::string GetPlatformInterpreter() const override;
        [[nodiscard]] std::string GetPlatformScript() const override;

        // C'tor...
        OsxDebugger(const char* executable, const net::MpiCommunicator& comm);
        // ... which the factory function needs to be able to get at.
        friend class Debugger;

    };

}

#endif // HEMELB_DEBUG_OSX_OSXDEBUGGER_H

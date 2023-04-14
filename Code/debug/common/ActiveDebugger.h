// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_DEBUG_COMMON_ACTIVEDEBUGGER_H
#define HEMELB_DEBUG_COMMON_ACTIVEDEBUGGER_H

#include <vector>

#include <debug/Debugger.h>

namespace hemelb::debug
{

    class ActiveDebugger : public Debugger
    {
        /* This class does the bulk of the pause/attach process. It
         * still has some pure virtual functions that must be overridden
         * to give platform-specific data.
         */

      public:
        void BreakHere() override;
        void Print(const char* iFormat, ...) override;

      protected:
        using VoI = std::vector<int>; // Vector of Ints
        using VoS = std::vector<std::string>; // Vector of Strings

        // C'tor
        ActiveDebugger(const char* executable, const net::MpiCommunicator& comm);

        bool mAmAttached; // Indicate attachment state
        VoI mPIds; // vector of process IDs

        void Attach() override;
        void GatherProcessIds();
        void SpawnDebuggers();

        // Platform specific stuff
        [[nodiscard]] virtual std::string GetBinaryPath() const = 0;
        [[nodiscard]] virtual std::string GetPlatformInterpreter() const = 0;
        [[nodiscard]] virtual std::string GetPlatformScript() const = 0;
    };
}
#endif

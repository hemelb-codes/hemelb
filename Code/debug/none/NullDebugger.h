// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_DEBUG_NONE_NULLDEBUGGER_H
#define HEMELB_DEBUG_NONE_NULLDEBUGGER_H

#include "debug/Debugger.h"

namespace hemelb::debug
{
    class NullDebugger : public Debugger
    {
      public:
        void BreakHere() override;
        void Print(const char* iFormat, ...) override;

      protected:
        void Attach() override;
        NullDebugger(const char* executable, const net::MpiCommunicator& comm);
        friend class Debugger;
    };

}

#endif

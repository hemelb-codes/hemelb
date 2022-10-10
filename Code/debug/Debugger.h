// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_DEBUG_DEBUGGER_H
#define HEMELB_DEBUG_DEBUGGER_H

#include <string>
#include "debug.h"
#include "net/MpiCommunicator.h"

namespace hemelb::debug
{

    // Interface for debugger control.
    // This is a singleton and that makes sense.
    class Debugger
    {
      public:
        // the singleton pattern
        static Debugger* Init(bool active, const char *, const net::MpiCommunicator& comm);
        static Debugger* Get();

        virtual void BreakHere() = 0;
        virtual void Print(const char* iFormat, ...) = 0;

      protected:
        Debugger(const char* executable, net::MpiCommunicator comm);
        virtual ~Debugger() = default;
        virtual void Attach() = 0;

        std::string mExecutable;
        const net::MpiCommunicator mCommunicator;
        // Singleton pattern
        static Debugger* singleton;
    };
}

#endif // HEMELB_DEBUG_DEBUGGER_H

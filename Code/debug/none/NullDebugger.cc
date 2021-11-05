// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iostream>
#include "debug/none/NullDebugger.h"

namespace hemelb
{
  namespace debug
  {
    // Do nothing!
    void NullDebugger::Attach(void)
    {
    }
    void NullDebugger::BreakHere(void)
    {
    }
    void NullDebugger::Print(const char* iFormat, ...)
    {
    }

    NullDebugger::NullDebugger(const char* const executable, const net::MpiCommunicator& comm) :
        Debugger(executable, comm)
    {
    }

  }
}

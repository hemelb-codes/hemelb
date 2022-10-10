// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "debug/Debugger.h"
#include "debug/PlatformDebugger.h"
#include "debug/none/NullDebugger.h"

namespace hemelb::debug
{
    void Break() {
        Debugger::Get()->BreakHere();
    }

    void Init(bool active, const char * const executable,
              const net::MpiCommunicator& comm)
    {
        Debugger::Init(active, executable, comm);
    }

    Debugger* Debugger::Init(bool active, const char * const executable,
                             const net::MpiCommunicator& comm)
    {
      /* Static member function that implements the singleton pattern.
       * Use the namespace function PlatformDebuggerFactory to
       * actually construct the instance. It should be defined in the
       * appropriate platform subdirectory.
       */
      if (Debugger::singleton == nullptr)
      {
        if (active)
          Debugger::singleton = new PlatformDebugger(executable, comm);
        else
          Debugger::singleton = new NullDebugger(executable, comm);
      }
      Debugger::singleton->Attach();
      return Debugger::singleton;
    }

    Debugger* Debugger::Get()
    {
      // Get the single instance.
      return Debugger::singleton;
    }

    // Init static members
    Debugger* Debugger::singleton = nullptr;

    Debugger::Debugger(const char* const executable, net::MpiCommunicator comm) :
        mExecutable(executable), mCommunicator(std::move(comm))
    {
    }

} // namespace


// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "debug/Debugger.h"
#include "debug/PlatformDebugger.h"
#include "debug/none/NullDebugger.h"

namespace hemelb
{
  namespace debug
  {

    Debugger* Debugger::Init(bool active, const char * const executable, comm::Communicator::ConstPtr comm)
    {
      /* Static member function that implements the singleton pattern.
       *
       * PlatformDebugger.h is configured to have an appropriate
       * typedef for the current platform.
       */
      if (Debugger::singleton == NULL)
      {
        if (active)
          Debugger::singleton = new PlatformDebugger(executable, comm);
        else
          Debugger::singleton = new NullDebugger(executable, comm);
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

    Debugger::Debugger(const char* const executable, comm::Communicator::ConstPtr comm) :
      mExecutable(executable), mCommunicator(comm)
    {
    }

    // Dtor
    Debugger::~Debugger()
    {
    }

  } // namespace debug
} // namespace hemelb

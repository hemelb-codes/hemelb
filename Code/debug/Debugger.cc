// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "debug/Debugger.h"
#include "debug/PlatformDebugger.h"
#include "debug/none/NullDebugger.h"

namespace hemelb
{
  namespace debug
  {

    Debugger* Debugger::Init(bool active, const char * const executable, const net::MpiCommunicator& comm)
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

    Debugger* Debugger::Get(void)
    {
      // Get the single instance.
      return Debugger::singleton;
    }

    // Init static members
    Debugger* Debugger::singleton = nullptr;

    Debugger::Debugger(const char* const executable, const net::MpiCommunicator& comm) :
      mExecutable(executable), mCommunicator(comm)
    {
    }

    // Dtor
    Debugger::~Debugger()
    {
    }

  } // namespace debug
} // namespace hemelb

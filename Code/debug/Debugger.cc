// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include "debug/Debugger.h"
#include "debug/platform.h"

namespace hemelb
{
  namespace debug
  {

    Debugger* Debugger::Init(const char * const executable)
    {
      /* Static member function that implements the singleton pattern.
       * Use the namespace function PlatformDebuggerFactory to
       * actually construct the instance. It should be defined in the
       * appropriate platform subdirectory.
       */
      if (Debugger::singleton == NULL)
      {
        Debugger::singleton = PlatformDebuggerFactory(executable);
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

    Debugger::Debugger(const char* const executable)
    {
      // Ctor
      mExecutable = std::string(executable);
    }

    // Dtor
    Debugger::~Debugger()
    {
    }

  } // namespace debug
} // namespace hemelb

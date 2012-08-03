// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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

    NullDebugger::NullDebugger(const char* const executable) :
      Debugger(executable)
    {
    }

    Debugger* PlatformDebuggerFactory(const char * const executable)
    {
      return new NullDebugger(executable);
    }
  }
}

// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_DEBUG_NONE_NULLDEBUGGER_H
#define HEMELB_DEBUG_NONE_NULLDEBUGGER_H

#include "debug/Debugger.h"

namespace hemelb
{
  namespace debug
  {
    class NullDebugger : public Debugger
    {
      public:
        void BreakHere(void);
        void Print(const char* iFormat, ...);

      protected:
        void Attach(void);
        NullDebugger(const char* const executable, const net::MpiCommunicator& comm);
        friend Debugger* PlatformDebuggerFactory(const char* const executable, const net::MpiCommunicator& comm);
    };

    Debugger* PlatformDebuggerFactory(const char* const executable, const net::MpiCommunicator& comm);
  }
}

#endif // HEMELB_DEBUG_NONE_DEBUGGER_H

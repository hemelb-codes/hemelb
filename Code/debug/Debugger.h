// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_DEBUG_DEBUGGER_H
#define HEMELB_DEBUG_DEBUGGER_H

#include <string>
#include "net/MpiCommunicator.h"

namespace hemelb
{
  namespace debug
  {

    class Debugger
    {
        /* Interface for debugger control.

         */
      public:
        // the singleton pattern
        static Debugger* Init(bool active, const char * const, const net::MpiCommunicator& comm);
        static Debugger* Get(void);

        virtual void BreakHere(void) = 0;
        virtual void Print(const char* iFormat, ...) = 0;

      protected:
        Debugger(const char* const executable, const net::MpiCommunicator& comm);
        virtual ~Debugger();

        virtual void Attach() = 0;

        std::string mExecutable;
        const net::MpiCommunicator mCommunicator;
        // Singleton pattern
        static Debugger* singleton;

    };

  }
}

#endif // HEMELB_DEBUG_DEBUGGER_H

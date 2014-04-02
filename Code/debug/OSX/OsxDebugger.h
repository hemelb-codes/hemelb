// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_DEBUG_OSX_OSXDEBUGGER_H
#define HEMELB_DEBUG_OSX_OSXDEBUGGER_H

#include <debug/common/ActiveDebugger.h>

namespace hemelb
{
  namespace debug
  {
    
    class OsxDebugger : public ActiveDebugger {
    protected:
      // Platform specific getters
      const std::string GetBinaryPath(void) const;
      const std::string GetPlatformInterpreter(void) const;
      const std::string GetPlatformScript(void) const;

      // C'tor...
      OsxDebugger(const char* const executable, const net::MpiCommunicator& comm);
      // ... which the factory function needs to be able to get at.
      friend Debugger* PlatformDebuggerFactory(const char* const executable, const net::MpiCommunicator& comm);
      
    };
    
    // Factory. Don't be calling this.
    Debugger* PlatformDebuggerFactory(const char* const executable, const net::MpiCommunicator& comm);
  }
}

#endif // HEMELB_DEBUG_OSX_OSXDEBUGGER_H

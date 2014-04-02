// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include <string>

#include "debug/OSX/OsxDebugger.h"

namespace hemelb
{
  namespace debug
  {
    OsxDebugger::OsxDebugger(const char* const executable, const net::MpiCommunicator& comm) :
      ActiveDebugger(executable, comm)
    {
    }

    const std::string OsxDebugger::GetPlatformInterpreter(void) const
    {
      return std::string("osascript");
    }

    const std::string OsxDebugger::GetPlatformScript(void) const
    {
      std::string include(__FILE__);
      std::string debugOsxDir = include.substr(0, include.rfind('/'));

      return debugOsxDir + "/MPIdebug.applescript";
    }

    Debugger* PlatformDebuggerFactory(const char * const executable, const net::MpiCommunicator& comm)
    {
      return new OsxDebugger(executable, comm);
    }

  } // namespace debug
} // namespace hemelb

// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include <string>

#include "debug/linux/LinuxDebugger.h"

namespace hemelb
{
  namespace debug
  {
    LinuxDebugger::LinuxDebugger(const char* const executable, const net::MpiCommunicator& comm) :
      ActiveDebugger(executable, comm) {}

    const std::string LinuxDebugger::GetPlatformInterpreter(void) const {
      return std::string("bash");
    }

    const std::string LinuxDebugger::GetPlatformScript(void) const {
      std::string include (__FILE__);
      std::string debugLinuxDir = include.substr(0, include.rfind('/'));

      return debugLinuxDir + "/launchGdbs.sh";
    }

    Debugger* PlatformDebuggerFactory(const char * const executable, const net::MpiCommunicator& comm) {
      return new LinuxDebugger(executable, comm);
    }

  } // namespace debug
} // namespace hemelb

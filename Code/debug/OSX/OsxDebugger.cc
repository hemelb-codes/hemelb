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
#include <mach-o/dyld.h>
namespace hemelb
{
  namespace debug
  {
    OsxDebugger::OsxDebugger(const char* const executable, const net::MpiCommunicator& comm) :
      ActiveDebugger(executable, comm)
    {
    }
    const std::string OsxDebugger::GetBinaryPath(void) const
    {
      char* path = nullptr;
      uint32_t size = 0;
      _NSGetExecutablePath(path, &size);
      path = new char[size];
      int ret = _NSGetExecutablePath(path, &size);
      if (ret != 0)
        // error
        throw Exception() << "Error getting executable path.";

      std::string ans(path, size);
      delete[] path;
      return ans;
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

  } // namespace debug
} // namespace hemelb

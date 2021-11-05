// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

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

// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <string>

#include "debug/linux/LinuxDebugger.h"
#include <unistd.h>
#include <cerrno>

namespace hemelb
{
  namespace debug
  {
    LinuxDebugger::LinuxDebugger(const char* const executable, const net::MpiCommunicator& comm) :
        ActiveDebugger(executable, comm)
    {
    }

    const std::string LinuxDebugger::GetBinaryPath(void) const
    {
      char buf[1024];
      ssize_t len;
      len = readlink("/proc/self/exe", buf, sizeof (buf) - 1);

      if (len == -1)
        // error
        throw Exception() << "Error getting executable path. Error code: " << errno;
      buf[len] = '\0';
      return std::string(buf, len);
    }

    const std::string LinuxDebugger::GetPlatformInterpreter(void) const
    {
      return std::string("bash");
    }

    const std::string LinuxDebugger::GetPlatformScript(void) const
    {
      std::string include(__FILE__);
      std::string debugLinuxDir = include.substr(0, include.rfind('/'));

      return debugLinuxDir + "/launchGdbs.sh";
    }

  } // namespace debug
} // namespace hemelb

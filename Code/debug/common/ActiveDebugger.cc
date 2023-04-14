// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include <iostream>
#include <optional>
#include <string>
#include <sstream>
#include <vector>

#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <cstdio>
#include <cstdarg>

#include <unistd.h>
#include <sys/wait.h>

#include "net/mpi.h"

#include "debug/common/ActiveDebugger.h"

namespace hemelb::debug
{

    ActiveDebugger::ActiveDebugger(const char* const executable, const net::MpiCommunicator& comm) :
        Debugger(executable, comm), mAmAttached(false), mPIds()
    {
    }

    namespace {
        std::string ConvertIntToString(int i) {
            // Convert an int to a string.
            std::stringstream ss;
            ss << i;
            return ss.str();
        }

        // Get an environment variable if it exists
        std::optional<std::string> GetEnv(char const* vname) {
            char const* val = std::getenv(vname);
            if (val == nullptr)
                return std::nullopt;
            return std::make_optional<std::string>(val);
        }
    }

    void ActiveDebugger::BreakHere()
    {
      return;
    }

    void ActiveDebugger::Print(const char* iFormat, ...)
    {
      std::va_list args;
      va_start(args, iFormat);
      std::vprintf(iFormat, args);
      va_end(args);
    }

    void ActiveDebugger::Attach()
    {
      // Start up the debugger(s), tell them to attach to the
      // processes and wait for them to attach. This function forks
      // another process on the rank 0 task and waits for it.
      if (mAmAttached)
        return;

      volatile int amWaiting = 1;

      int childPid = -1;
      // To rank 0
      GatherProcessIds();

      if (mCommunicator.Rank() == 0)
      {
        childPid = fork();
        // Fork gives the PID of the child to the parent and zero to the child
        if (childPid == 0)
        {
          // This function won't return.
          SpawnDebuggers();
        }
      } // if rank 0

      // Sit here waiting for GDB to attach
      while (amWaiting)
        sleep(5);

      if (mCommunicator.Rank() == 0)
      {
        // Reap the spawner
        int deadPid = waitpid(childPid, nullptr, 0);
        if (deadPid != childPid)
          std::cerr << "Error in waitpid, code: " << errno << std::endl;
      }

      mAmAttached = true;
    }

    void ActiveDebugger::GatherProcessIds()
    {
      int pId = getpid();

      mPIds = mCommunicator.Gather(pId, 0);
    }

    void ActiveDebugger::SpawnDebuggers()
    {
      // Run the script to Tell the OS to start appropriate
      // terminal(s), launch the debuggers and attach them to our
      // processes.  We're the child, so we are DOOMED.  (Either to
      // exec() or if that fails exit().)
      std::string srcFile(__FILE__);
      std::string debugCommonDir = srcFile.substr(0, srcFile.rfind('/'));

      const std::string binaryPath = GetBinaryPath();

      VoS args;

      args.push_back(GetPlatformInterpreter());

      args.push_back(GetPlatformScript());

      args.push_back(binaryPath);

      auto debugger_commands_file = GetEnv("HEMELB_DEBUG_COMMANDS");
      if (debugger_commands_file.has_value()) {
          args.push_back("-file");
          args.push_back(debugger_commands_file.value());
      }

      for (VoI::iterator i = mPIds.begin(); i < mPIds.end(); ++i)
      {
        args.push_back(ConvertIntToString(*i));
      }

      // +1 to include required nullptr pointer for execvp()
      char **argv = new char *[args.size() + 1];

      // convert to C array of char arrays.
      for (unsigned int i = 0; i < args.size(); ++i)
      {
        argv[i] = new char[args[i].length() + 1]; // for terminating null
        std::strcpy(argv[i], args[i].c_str());
      }

      // Terminating nullptr
      argv[args.size()] = nullptr;

      // Exec to replace hemelb with osascript
      int code = execvp(argv[0], argv);

      // OK- that didn't work if we get here, better die (since we're
      // the extra process). Print the error code too.
      std::cerr << "Couldn't exec() script to launch debuggers. Return value: " << code
          << "; error code: " << errno << std::endl;
      // Now print the command we wanted to exec()
      for (VoS::iterator it = args.begin(); it < args.end(); it++)
      {
        std::cerr << *it << " ";
      }
      std::cerr.flush();
      // Die
      exit(1);
    }

}

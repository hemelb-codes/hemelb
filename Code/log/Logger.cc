// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cstdarg>
#include <sys/time.h>
#include <sys/resource.h>

#include "util/utilityFunctions.h"
#include "net/mpi.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace log
  {
    const LogLevel Logger::currentLogLevel = HEMELB_LOG_LEVEL;
    // Use negative value to indicate uninitialised.
    int Logger::thisRank = -1;
    double Logger::startTime = -1.0;

    void Logger::Init()
    {
      // If Logger uninitialised
      if (thisRank < 0)
      {
        // Check that MPI is ready
        if (net::MpiEnvironment::Initialized())
        {
          thisRank = net::MpiCommunicator::World().Rank();
        }
        startTime = util::myClock();
      }
    }

    template<>
    void Logger::LogInternal<OnePerCore>(std::string format, std::va_list args)
    {
      std::stringstream output;

      // Set the fill digit to be 0, so the integer 1 renders as 0000001
      output.fill('0');

      output << "[Rank " << std::setw(7) << thisRank << ", " << std::setiosflags(std::ios::fixed)
          << std::setprecision(1) << (util::myClock() - startTime) << "s";

#ifdef HAVE_RUSAGE
      rusage usage;
      getrusage(RUSAGE_SELF, &usage);

      output << ", mem: " << std::setw(7) << usage.ru_maxrss;
#endif

      output << "]: " << format << '\n';

      std::string overFormat(output.str());

      std::vprintf(overFormat.c_str(), args);
    }

    template<>
    void Logger::LogInternal<Singleton>(std::string format, std::va_list args)
    {
      if (thisRank == 0)
      {
        char lead[20];
        std::sprintf(lead, "![%.1fs]", util::myClock() - startTime);

        std::string newFormat = std::string(lead);
        std::vprintf(newFormat.append(format).append("\n").c_str(), args);
      }
    }

  }
}

#include <iostream>
#include <cstdio>
#include <cstdarg>
#include <sys/time.h>
#include <sys/resource.h>

#include "mpiInclude.h"
#include "util/utilityFunctions.h"
#include "mpiInclude.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace log
  {
    int Logger::thisRank = -1;
    double Logger::startTime = -1.0;

    template<>
    void Logger::LogInternal<OnePerCore>(std::string format, std::va_list args)
    {
      if (thisRank < 0)
      {
        // need to be able to log, even if MPI not initialised (for testability)
        int initialized;
        MPI_Initialized(&initialized);
        if (initialized) {
          MPI_Comm_rank(MPI_COMM_WORLD, &thisRank);
        }
        startTime = util::myClock();
      }

      rusage usage;
      getrusage(RUSAGE_SELF, &usage);

      char lead[60];
      std::sprintf(lead, "[Rank %.6i, %.1fs, mem: %li]: ", thisRank, util::myClock() - startTime, usage.ru_maxrss);

      std::string overFormat(lead);
      overFormat.append(format).append("\n");

      std::vprintf(overFormat.c_str(), args);
    }

    template<>
    void Logger::LogInternal<Singleton>(std::string format, std::va_list args)
    {
      if (thisRank < 0)
      {
        int initialized;
        MPI_Initialized(&initialized);
        if (initialized) {
          MPI_Comm_rank(MPI_COMM_WORLD, &thisRank);
        }
        startTime = util::myClock();
      }

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

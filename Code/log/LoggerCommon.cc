#include <iostream>
#include <cstdio>
#include <cstdarg>

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
        MPI_Comm_rank(MPI_COMM_WORLD, &thisRank);
        startTime = util::myClock();
      }

      char lead[40];
      std::sprintf(lead, "[Rank %.6i, %.1fs]: ", thisRank, util::myClock() - startTime);

      std::string overFormat(lead);
      overFormat.append(format).append("\n");

      std::vprintf(overFormat.c_str(), args);
    }

    template<>
    void Logger::LogInternal<Singleton>(std::string format, std::va_list args)
    {
      if (thisRank < 0)
      {
        MPI_Comm_rank(MPI_COMM_WORLD, &thisRank);
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

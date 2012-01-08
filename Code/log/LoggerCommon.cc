#include <cstdlib>
#include <iostream>
#include <stdio.h>

#include "mpiInclude.h"
#include "util/utilityFunctions.h"
#include "log/Logger.h"

namespace hemelb
{
  namespace log
  {
    int Logger::thisRank = -1;
    double Logger::startTime = -1.0;

    template<>
    void Logger::LogInternal<OnePerCore>(std::string format, va_list args)
    {
      if (thisRank < 0)
      {
        MPI_Comm_rank(MPI_COMM_WORLD, &thisRank);
        startTime = util::myClock();
      }

      char lead[40];
      sprintf(lead, "[Rank %.6i, %.1fs]: ", thisRank, util::myClock() - startTime);

      std::string overFormat(lead);
      overFormat.append(format).append("\n");

      vprintf(overFormat.c_str(), args);
    }

    template<>
    void Logger::LogInternal<Singleton>(std::string format, va_list args)
    {
      if (thisRank < 0)
      {
        MPI_Comm_rank(MPI_COMM_WORLD, &thisRank);
        startTime = util::myClock();
      }

      if (thisRank == 0)
      {
        char lead[20];
        sprintf(lead, "![%.1fs]", util::myClock() - startTime);

        std::string newFormat = std::string(lead);
        vprintf(newFormat.append(format).append("\n").c_str(), args);
      }
    }
  }
}

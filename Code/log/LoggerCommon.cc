#include <cstdlib>
#include <iostream>
#include <stdio.h>

#include "log/Logger.h"

namespace hemelb
{
  namespace log
  {
    int Logger::thisRank = -1;

    template<>
    void Logger::LogInternal<OnePerCore>(std::string format, va_list args)
    {
      if (thisRank < 0)
      {
        MPI_Comm_rank(MPI_COMM_WORLD, &thisRank);
      }

      char lead[20];
      sprintf(lead, "[Rank %.6i]: ", thisRank);

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
      }

      if (thisRank == 0)
      {
        std::string newFormat = "!";
        vprintf(newFormat.append(format).append("\n").c_str(), args);
      }
    }
  }
}

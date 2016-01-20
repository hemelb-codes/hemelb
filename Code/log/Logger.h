
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LOG_LOGGER_H
#define HEMELB_LOG_LOGGER_H

#include <cstdarg>
#include <string>

namespace hemelb
{
  namespace log
  {
    enum LogLevel
    {
      //! Unrecoverable error
      Critical = 0,
      //! Recoverable error or unexpected state
      Error = 1,
      //! Undesirable but not necessarily wrong states
      Warning = 2,
      //! A conservative set of interesting run time events
      Info = 3,
      //! Detailed information on runtime flow
      //! This level should also be used for any non-essential computation.
      Debug = 4,
      //! Low-level trace information, including actual data.
      Trace = 5
    };

    enum LogType
    {
      Singleton,
      OnePerCore
    };

    class Logger
    {
      public:
        template<LogLevel queryLogLevel>
        static bool ShouldDisplay()
        {
          return queryLogLevel <= currentLogLevel;
        }

        static void Init();

        template<LogLevel queryLogLevel, LogType logType>
        static void Log(std::string format, ...)
        {
          if (queryLogLevel <= currentLogLevel)
          {
            va_list args;
            va_start(args, format);
            LogInternal<logType> (format, args);
            va_end(args);
          }
        }

      private:
        template<LogType>
        static void LogInternal(std::string format, va_list args);

        static const LogLevel currentLogLevel;
        static int thisRank;
        static double startTime;
    };

  }
}

#endif /* HEMELB_LOG_LOGGER_H */

// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_LOG_LOGGER_H
#define HEMELB_LOG_LOGGER_H

#include <cstdarg>
#include <sstream>
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
            LogInternal<logType>(format, args);
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
#define HEMELB_CAPTURE(VARIABLE)                                        \
  hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>(\
    #VARIABLE " := %s",                                                 \
    (std::ostringstream() << VARIABLE).str().c_str()                    \
  )
#define HEMELB_CAPTURE2(V0, V1)                                         \
  hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>(\
    #V0 " := %s " #V1 " := %s",                                         \
    (std::ostringstream() << V0).str().c_str(),                         \
    (std::ostringstream() << V1).str().c_str()                          \
  )
#define HEMELB_CAPTURE3(V0, V1, V2)                                     \
  hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>(\
    #V0 " := %s " #V1 " := %s " #V2 " := %s",                           \
    (std::ostringstream() << V0).str().c_str(),                         \
    (std::ostringstream() << V1).str().c_str(),                         \
    (std::ostringstream() << V2).str().c_str()                          \
  )
#define HEMELB_CAPTURE4(V0, V1, V2, V3)                                 \
  hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>(\
    #V0 " := %s " #V1 " := %s " #V2 " := %s",                           \
    (std::ostringstream() << V0).str().c_str(),                         \
    (std::ostringstream() << V1).str().c_str(),                         \
    (std::ostringstream() << V2).str().c_str(),                         \
    (std::ostringstream() << V3).str().c_str()                          \
  )

#endif /* HEMELB_LOG_LOGGER_H */

// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_LOG_LOGGER_H
#define HEMELB_LOG_LOGGER_H

#include <vector>
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

    template<class T>
       std::string __convert_to_string(T && value)
       {
         std::ostringstream sstr;
         sstr << value;
         return sstr.str();
       }

    template<class ... ARGS>
      void capture(std::string message, ARGS &&... args)
      {
         hemelb::log::Logger::Log<hemelb::log::Debug, hemelb::log::OnePerCore>(
             message,
             __convert_to_string(std::forward<ARGS>(args)).c_str()...
         );
      }
  }
}

//! Simple debug macro to print stuff to the log
#define HEMELB_CAPTURE(VARIABLE) hemelb::log::capture(#VARIABLE " := %s", VARIABLE)
//! Simple debug macro to print stuff to the log
#define HEMELB_CAPTURE2(V0, V1) hemelb::log::capture(#V0 " := %s " #V1 " := %s", V0, V1)
//! Simple debug macro to print stuff to the log
#define HEMELB_CAPTURE3(V0, V1, V2)     \
  hemelb::log::capture(#V0 " := %s " #V1 " := %s " #V2 " := %s", V0, V1, V2)
//! Simple debug macro to print stuff to the log
#define HEMELB_CAPTURE4(V0, V1, V2, V3) \
  hemelb::log::capture(#V0 " := %s " #V1 " := %s " #V2 " := %s " #V3 ": = %s", V0, V1, V2, V3)
//! Simple debug macro to print stuff to the log
#define HEMELB_CAPTURE5(V0, V1, V2, V3, V4) \
  hemelb::log::capture(                     \
    #V0 " := %s " #V1 " := %s " #V2 " := %s " #V3 ": = %s " #V4 ": = %s", V0, V1, V2, V3, V4)

#endif /* HEMELB_LOG_LOGGER_H */

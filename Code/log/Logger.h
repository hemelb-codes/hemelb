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
      Warning = 0,
      Info = 1,
      Debug = 2
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

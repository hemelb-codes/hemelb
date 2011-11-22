#ifndef HEMELB_REPORTER_H
#define HEMELB_REPORTER_H

#include <string>
#include <stdarg.h>
#include "configuration/SimConfig.h"
#include "configuration/CommandLine.h"
#include "log/Logger.h"
#include "util/fileutils.h"
#include "reporting/Timers.h"
namespace hemelb{
  namespace reporting {
    // abstracting this policy so we can mock it to test the reporter
    class FileWriterPolicy {
      public:
        FileWriterPolicy(const std::string &path){
          file = fopen(path.c_str(), "w");
        }
        ~FileWriterPolicy()
        {
          fclose(file);
        }
      protected:
        void Print(const char * format, ...){
          std::va_list arg;
          va_start(arg, format);
          vfprintf(file,format, arg);
          va_end(arg);
        }
      private:
        FILE* file;
    };


   template<class TimersPolicy, class WriterPolicy> class ReporterBase : public WriterPolicy {
      public:
        ReporterBase(const std::string &path, const std::string &inputFile, const long int asite_count, TimersPolicy&timers);
        void Cycle();
        void TimeStep(){timestep_count++;}
        void Image();
        void Snapshot();
        void Write();
        void Stability(bool astability){stability=astability;}
      private:
        bool doIo;
        unsigned int cycle_count;
        unsigned int snapshot_count;
        unsigned int image_count;
        unsigned long int timestep_count;
        long int site_count;
        bool stability;
        TimersPolicy &timings;
    };

   typedef ReporterBase<Timers,FileWriterPolicy> Reporter;
  }
}

#endif

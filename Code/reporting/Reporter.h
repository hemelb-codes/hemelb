#ifndef HEMELB_REPORTER_H
#define HEMELB_REPORTER_H

#include <string>
#include "configuration/SimConfig.h"
#include "configuration/CommandLine.h"
#include "log/Logger.h"
#include "util/fileutils.h"
#include "reporting/Timers.h"
namespace hemelb{
  namespace reporting {
   template<class TimersPolicy> class ReporterBase {
      public:
        ReporterBase(const std::string &path, const std::string &inputFile, const long int asite_count, TimersPolicy&timers);
        ~ReporterBase();
        void Cycle();
        void TimeStep(){timestep_count++;}
        void Image();
        void Snapshot();
        void Write();
        void Stability(bool astability){stability=astability;}
      private:
        FILE *ReportFile(){
          return mTimingsFile;
        }
        bool doIo;
        FILE *mTimingsFile;
        unsigned int cycle_count;
        unsigned int snapshot_count;
        unsigned int image_count;
        unsigned long int timestep_count;
        long int site_count;
        bool stability;
        TimersPolicy &timings;
    };

   typedef ReporterBase<Timers> Reporter;
  }
}

#endif

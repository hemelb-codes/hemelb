#ifndef HEMELB_REPORTER_H
#define HEMELB_REPORTER_H

#include <string>
#include "configuration/SimConfig.h"
#include "configuration/CommandLine.h"
#include "log/Logger.h"
#include "util/fileutils.h"
namespace hemelb{
  namespace reporting {
    class Reporter {
      public:
        Reporter(bool doio, const std::string &name, std::string &inputFile);
        ~Reporter();
        void Cycle(long int cycle_id);
        void Phase1(long int site_count, int total_time_steps, long int cycle_id, double simulation_time, bool unstable,
                    unsigned long time_steps_per_cycle,
                    double decomposition_time, double initialise_time, double read_time, double creation_time);
        void ProcessorTimings(std::string *const names,double *const mins,double *const means,double *const maxes);
      private:
        FILE *ReportFile(){
            return mTimingsFile;
        }
        bool doIo;
        FILE *mTimingsFile;
    };
  }
}

#endif

#ifndef HEMELB_REPORTING_REPORTER_H
#define HEMELB_REPORTING_REPORTER_H

#include <string>
#include <stdarg.h>
#include "configuration/SimConfig.h"
#include "configuration/CommandLine.h"
#include "log/Logger.h"
#include "util/fileutils.h"
#include "Timers.h"
#include "Policies.h"
namespace hemelb
{
  namespace reporting
  {
    template<class TimersPolicy, class WriterPolicy, class CommsPolicy> class ReporterBase : public WriterPolicy,
                                                                                             public CommsPolicy
    {
      public:
        ReporterBase(const std::string &path,
                     const std::string &inputFile,
                     const long int aSiteCount,
                     TimersPolicy& timers);
        void Cycle();
        void TimeStep()
        {
          timestepCount++;
        }
        void Image();
        void Snapshot();
        void Write();
        void Stability(bool astability)
        {
          stability = astability;
        }
      private:
        bool doIo;
        unsigned int cycleCount;
        unsigned int snapshotCount;
        unsigned int imageCount;
        unsigned long int timestepCount;
        long int siteCount;
        bool stability;
        TimersPolicy &timings;
    };

    typedef ReporterBase<Timers, FileWriterPolicy, MPICommsPolicy> Reporter;
  }
}

#endif // HEMELB_REPORTING_REPORTER_H

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
        ReporterBase(const std::string &path, const std::string &inputFile, const long int asite_count, TimersPolicy &timers):
          cycle_count(0), timestep_count(0), site_count(asite_count),stability(false),  timings(timers)
       {
         mTimingsFile = fopen(path.c_str(), "w");
         fprintf(mTimingsFile, "***********************************************************\n");
         fprintf(mTimingsFile, "Opening config file:\n %s\n", inputFile.c_str());
       }
        ~ReporterBase(){
          fclose(mTimingsFile);
        }
        void Cycle(){
          cycle_count++;
          fprintf(ReportFile(), "cycle id: %u\n", cycle_count);
        }
        void TimeStep(){timestep_count++;}
        void Image(){
          image_count++;
          fprintf(ReportFile(), "Image written: %u\n", image_count);
        }
        void Snapshot(){
          snapshot_count++;
          fprintf(ReportFile(), "Snapshot written: %u\n", snapshot_count);
        }
        void Write(){

          fprintf(ReportFile(), "\n");
          fprintf(ReportFile(),
                  "threads: %i, machines checked: %i\n\n",
                  hemelb::topology::NetworkTopology::Instance()->GetProcessorCount(),
                  hemelb::topology::NetworkTopology::Instance()->GetMachineCount());
          fprintf(ReportFile(),
                  "topology depths checked: %i\n\n",
                  hemelb::topology::NetworkTopology::Instance()->GetDepths());
          fprintf(ReportFile(), "fluid sites: %li\n\n", site_count);
          fprintf(ReportFile(),
                  "cycles and total time steps: %u, %lu \n\n",
                  cycle_count,
                  timestep_count);
          fprintf(ReportFile(),
                  "time steps per second: %.3f\n\n",
                  timestep_count / timings[TimersPolicy::simulation].Get());

          if (!stability)
          {
            fprintf(ReportFile(),
                    "Attention: simulation unstable with %lu timesteps/cycle\n",
                    timestep_count / cycle_count);
            fprintf(ReportFile(), "Simulation terminated\n");
          }

          fprintf(ReportFile(), "time steps per cycle: %lu\n", timestep_count / cycle_count);
          fprintf(ReportFile(), "\n");

          fprintf(ReportFile(), "\n");

          fprintf(ReportFile(), "\n");
          fprintf(ReportFile(),
                  "total time (s):                            %.3f\n\n",
                  (timings[TimersPolicy::total].Get()));

          fprintf(ReportFile(), "Sub-domains info:\n\n");

          for (hemelb::proc_t n = 0;
              n < hemelb::topology::NetworkTopology::Instance()->GetProcessorCount(); n++)
          {
            fprintf(ReportFile(),
                    "rank: %lu, fluid sites: %lu\n",
                    (unsigned long) n,
                    (unsigned long) hemelb::topology::NetworkTopology::Instance()->FluidSitesOnEachProcessor[n]);
          }

          // Note that CycleId is 1-indexed and will have just been incremented when we finish.
          double cycles = cycle_count;

          double normalisations[TimersPolicy::numberOfTimers] = { 1.0,
                                                            1.0,
                                                            1.0,
                                                            1.0,
                                                            cycles,
                                                            image_count,
                                                            cycles,
                                                            cycles,
                                                            snapshot_count,
                                                            1.0 };

          fprintf(ReportFile(),
                  "\n\nPer-proc timing data (secs per [simulation,simulation,simulation,simulation,cycle,image,cycle,cycle,snapshot,simulation]): \n\n");
          fprintf(ReportFile(), "\t\tLocal \tMin \tMean \tMax\n");
          for (unsigned int ii = 0; ii < TimersPolicy::numberOfTimers; ii++)
          {
            fprintf(ReportFile(),
                    "%s\t\t%.3g\t%.3g\t%.3g\t%.3g\n",
                    timerNames[ii].c_str(),
                    timings[ii].Get(),
                    timings.Mins()[ii] / normalisations[ii],
                    timings.Means()[ii] / normalisations[ii],
                    timings.Maxes()[ii] / normalisations[ii]);
          }
        }
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

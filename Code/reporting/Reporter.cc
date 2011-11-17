#include "Reporter.h"
namespace hemelb{
  namespace reporting {
    Reporter::Reporter(bool doio, const std::string &name, std::string &inputFile):
                          doIo(doio)
    {
      mTimingsFile = fopen(name.c_str(), "w");
      fprintf(mTimingsFile, "***********************************************************\n");
      fprintf(mTimingsFile, "Opening config file:\n %s\n", inputFile.c_str());
    }

    Reporter::~Reporter(){
      fclose(mTimingsFile);
    }

    void Reporter::Cycle(long int cycle_id)
    {
      fprintf(ReportFile(), "cycle id: %li\n", cycle_id);
    }

    void Reporter::ProcessorTimings(std::string *const names,double *const mins,double *const means,double *const maxes){



    }


    void Reporter::Phase1(long int site_count, int total_time_steps, long int cycle_id,
                          bool unstable, unsigned long time_steps_per_cycle, unsigned int image_count, unsigned int snapshot_count,
                          Timers &timings){

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
              "cycles and total time steps: %li, %i \n\n",
              (cycle_id - 1), // Note that the cycle-id is 1-indexed.
              total_time_steps);
      fprintf(ReportFile(), "time steps per second: %.3f\n\n", total_time_steps / timings[Timers::simulation].Get());

      if (unstable)
      {
        fprintf(ReportFile(),
                "Attention: simulation unstable with %li timesteps/cycle\n", time_steps_per_cycle
        );
        fprintf(ReportFile(), "Simulation terminated\n");
      }

      fprintf(ReportFile(),
              "time steps per cycle: %li\n",
              time_steps_per_cycle);
      fprintf(ReportFile(), "\n");

      fprintf(ReportFile(), "\n");

      fprintf(ReportFile(), "\n");
      fprintf(ReportFile(),
              "total time (s):                            %.3f\n\n",
              (timings[Timers::total].Get()));

      fprintf(ReportFile(), "Sub-domains info:\n\n");

      for (hemelb::proc_t n = 0; n
      < hemelb::topology::NetworkTopology::Instance()->GetProcessorCount(); n++)
      {
        fprintf(ReportFile(),
                "rank: %lu, fluid sites: %lu\n",
                (unsigned long) n,
                (unsigned long) hemelb::topology::NetworkTopology::Instance()->FluidSitesOnEachProcessor[n]);
      }


      // Note that CycleId is 1-indexed and will have just been incremented when we finish.
      double cycles = hemelb::util::NumericalFunctions::max(1.0, (double) cycle_id- 1);


      //std::string lNames[5] = { "LBM", "MPISend", "MPIWait", "Images", "Snaps" };

      double normalisations[Timers::numberOfTimers]={1.0,1.0,1.0,1.0,cycles,image_count,cycles,cycles,snapshot_count,1.0};

      fprintf(ReportFile(),
              "\n\nPer-proc timing data (secs per [simulation,simulation,simulation,simulation,cycle,image,cycle,cycle,snapshot,simulation]): \n\n");
      fprintf(ReportFile(), "\t\tLocal \tMin \tMean \tMax\n");
      for (unsigned int ii = 0; ii < Timers::numberOfTimers; ii++)
      {
        fprintf(ReportFile(),
                "%s\t\t%.3g\t%.3g\t%.3g\t%.3g\n",
                timerNames[ii].c_str(),
                timings[ii].Get(),
                timings.Mins()[ii]/ normalisations[ii],
                timings.Means()[ii]/ normalisations[ii],
                timings.Maxes()[ii]/ normalisations[ii]);
      }
    }
  }
}

#ifndef HEMELB_REPORTING_REPORTER_HPP
#define HEMELB_REPORTING_REPORTER_HPP
#include "Reporter.h"
namespace hemelb
{
  namespace reporting
  {
    template<class TimersPolicy, class WriterPolicy, class CommsPolicy> ReporterBase<TimersPolicy, WriterPolicy, CommsPolicy>::ReporterBase(const std::string &name,
                       const std::string &inputFile,
                       const long int asite_count,
                       TimersPolicy& timers) :
        WriterPolicy(name),cycle_count(0), snapshot_count(0), image_count(0), timestep_count(0), site_count(asite_count), stability(true), timings(timers)
    {
      WriterPolicy::Print( "***********************************************************\n");
      WriterPolicy::Print("Opening config file:\n %s\n", inputFile.c_str());
    }

    template<class TimersPolicy, class WriterPolicy, class CommsPolicy> void ReporterBase<TimersPolicy, WriterPolicy, CommsPolicy>::Cycle()
    {
      cycle_count++;
      WriterPolicy::Print("cycle id: %u\n", cycle_count);
    }

    template<class TimersPolicy, class WriterPolicy, class CommsPolicy> void ReporterBase<TimersPolicy, WriterPolicy, CommsPolicy>::Image()
    {
      image_count++;
      WriterPolicy::Print( "Image written: %u\n", image_count);
    }

    template<class TimersPolicy, class WriterPolicy, class CommsPolicy> void ReporterBase<TimersPolicy, WriterPolicy, CommsPolicy>::Snapshot()
    {
      snapshot_count++;
      WriterPolicy::Print( "Snapshot written: %u\n", snapshot_count);
    }

    template<class TimersPolicy, class WriterPolicy, class CommsPolicy> void ReporterBase<TimersPolicy, WriterPolicy, CommsPolicy>::Write()
    {

      WriterPolicy::Print( "\n");
      WriterPolicy::Print(
              "threads: %i, machines checked: %i\n\n",
             CommsPolicy::GetProcessorCount(),
              CommsPolicy::GetMachineCount());
      WriterPolicy::Print(
              "topology depths checked: %i\n\n",
              CommsPolicy::GetDepths());
      WriterPolicy::Print( "fluid sites: %li\n\n", site_count);
      WriterPolicy::Print(
              "cycles and total time steps: %u, %lu \n\n",
              cycle_count,
              timestep_count);
      WriterPolicy::Print(
              "time steps per second: %.3f\n\n",
              timestep_count / timings[TimersPolicy::simulation].Get());

      if (!stability)
      {
        WriterPolicy::Print(
                "Attention: simulation unstable with %lu timesteps/cycle\n",
                timestep_count / cycle_count);
        WriterPolicy::Print( "Simulation terminated\n");
      }

      WriterPolicy::Print( "time steps per cycle: %lu\n", timestep_count / cycle_count);
      WriterPolicy::Print( "\n");

      WriterPolicy::Print( "\n");

      WriterPolicy::Print( "\n");
      WriterPolicy::Print(
              "total time (s):                            %.3f\n\n",
              (timings[TimersPolicy::total].Get()));

      WriterPolicy::Print( "Sub-domains info:\n\n");

      for (hemelb::proc_t n = 0;
          n < CommsPolicy::GetProcessorCount(); n++)
      {
        WriterPolicy::Print(
                "rank: %lu, fluid sites: %lu\n",
                (unsigned long) n,
                (unsigned long) CommsPolicy::FluidSitesOnProcessor(n));
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

      WriterPolicy::Print(
              "\n\nPer-proc timing data (secs per [simulation,simulation,simulation,simulation,cycle,image,cycle,cycle,snapshot,simulation]): \n\n");
      WriterPolicy::Print( "\t\tLocal \tMin \tMean \tMax\n");
      for (unsigned int ii = 0; ii < TimersPolicy::numberOfTimers; ii++)
      {
        WriterPolicy::Print(
                "%s\t\t%.3g\t%.3g\t%.3g\t%.3g\n",
                timerNames[ii].c_str(),
                timings[ii].Get(),
                timings.Mins()[ii] / normalisations[ii],
                timings.Means()[ii] / normalisations[ii],
                timings.Maxes()[ii] / normalisations[ii]);
      }
    }
  }
}
#endif //HEMELB_REPORTING_REPORTER_HPP

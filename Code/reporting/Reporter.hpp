#ifndef HEMELB_REPORTING_REPORTER_HPP
#define HEMELB_REPORTING_REPORTER_HPP
#include "Reporter.h"
#include <iomanip>
namespace hemelb
{
  namespace reporting
  {
    template<class ClockPolicy, class WriterPolicy, class CommsPolicy, class BroadcastPolicy> ReporterBase<
        ClockPolicy, WriterPolicy, CommsPolicy, BroadcastPolicy>::ReporterBase(const std::string &name,
                                                                               const std::string &inputFile,
                                                                               const long int aSiteCount,
                                                                               const TimersBase<
                                                                                   ClockPolicy,
                                                                                   CommsPolicy>& timers,
                                                                               const lb::SimulationState &aState,
                                                                               const lb::IncompressibilityChecker<
                                                                                   BroadcastPolicy> &aChecker) :
        WriterPolicy(name), snapshotCount(0), imageCount(0), siteCount(aSiteCount), stability(true), timings(timers), state(aState), incompressibilityChecker(aChecker)
    {
      WriterPolicy::Print("***********************************************************\n");
      WriterPolicy::Print("Opening config file:\n %s\n", inputFile.c_str());
    }

    template<class ClockPolicy, class WriterPolicy, class CommsPolicy, class BroadcastPolicy> void ReporterBase<
        ClockPolicy, WriterPolicy, CommsPolicy, BroadcastPolicy>::Image()
    {
      imageCount++;
      WriterPolicy::Print("Image written: %u\n", imageCount);
    }

    template<class ClockPolicy, class WriterPolicy, class CommsPolicy, class BroadcastPolicy> void ReporterBase<
        ClockPolicy, WriterPolicy, CommsPolicy, BroadcastPolicy>::Snapshot()
    {
      snapshotCount++;
      WriterPolicy::Print("Snapshot written: %u\n", snapshotCount);
    }

    template<class ClockPolicy, class WriterPolicy, class CommsPolicy, class BroadcastPolicy> void ReporterBase<
        ClockPolicy, WriterPolicy, CommsPolicy, BroadcastPolicy>::Write()
    {

      // Note that CycleId is 1-indexed and will have just been incremented when we finish.
      unsigned long cycles = state.GetCycleId() - 1;
      WriterPolicy::Stream() << std::setprecision(3);
      WriterPolicy::Print("\n");
      WriterPolicy::Print("threads: %i, machines checked: %i\n\n",
                          CommsPolicy::GetProcessorCount(),
                          CommsPolicy::GetMachineCount());
      WriterPolicy::Print("topology depths checked: %i\n\n", CommsPolicy::GetDepths());
      WriterPolicy::Print("fluid sites: %li\n\n", siteCount);
      WriterPolicy::Print("cycles and total time steps: %lu, %lu \n\n",
                          cycles,
                          state.GetTimeStepsPassed() - 1);
      WriterPolicy::Print("time steps per second: %.3f\n\n",
                          (state.GetTimeStepsPassed() - 1) / timings[Timers::simulation].Get());

      if (incompressibilityChecker.AreDensitiesAvailable()
          && !incompressibilityChecker.IsDensityDiffWithinRange())
      {
        WriterPolicy::Print("Maximum relative density difference allowed (%.1f%%) was violated: %.1f%% \n\n",
                            incompressibilityChecker.GetMaxRelativeDensityDifferenceAllowed() * 100,
                            incompressibilityChecker.GetMaxRelativeDensityDifference() * 100);
      }

      if (!stability)
      {
        WriterPolicy::Print("Attention: simulation unstable with %lu timesteps/cycle\n",
                            state.GetTimeStepsPerCycle());
        WriterPolicy::Print("Simulation terminated\n");
      }

      WriterPolicy::Print("time steps per cycle: %lu\n", state.GetTimeStepsPerCycle());
      WriterPolicy::Print("\n");

      WriterPolicy::Print("\n");

      WriterPolicy::Print("\n");
      WriterPolicy::Print("total time (s):                            %.3f\n\n",
                          (timings[Timers::total].Get()));

      WriterPolicy::Print("Sub-domains info:\n\n");

      for (hemelb::proc_t n = 0; n < CommsPolicy::GetProcessorCount(); n++)
      {
        WriterPolicy::Print("rank: %lu, fluid sites: %lu\n",
                            (unsigned long) n,
                            (unsigned long) CommsPolicy::FluidSitesOnProcessor(n));
      }

      std::vector<double> normalisations;
      normalisations.push_back(1.0);
      normalisations.push_back(1.0);
      normalisations.push_back(1.0);
      normalisations.push_back(1.0);
      normalisations.push_back(cycles);
      normalisations.push_back(imageCount);
      normalisations.push_back(cycles);
      normalisations.push_back(cycles);
      normalisations.push_back(cycles);
      normalisations.push_back(snapshotCount);
      normalisations.push_back(1.0);
      // If this assertion trips, it means you have added a new timer and not defined its normalisation factor above.
      assert(normalisations.size() == Timers::numberOfTimers);

      WriterPolicy::Print("\n\nPer-proc timing data (secs per [simulation,simulation,simulation,simulation,cycle,image,cycle,cycle,cycle,snapshot,simulation]): \n\n");
      WriterPolicy::Print("\t\tLocal \tMin \tMean \tMax\n");
      for (unsigned int ii = 0; ii < Timers::numberOfTimers; ii++)
      {
        WriterPolicy::Stream() << timerNames[ii] << "\t\t" << timings[ii].Get() << "\t"
            << timings.Mins()[ii] / normalisations[ii] << "\t"
            << timings.Means()[ii] / normalisations[ii] << "\t"
            << timings.Maxes()[ii] / normalisations[ii] << std::endl;
      }
    }
  }
}
#endif //HEMELB_REPORTING_REPORTER_HPP

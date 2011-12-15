#ifndef HEMELB_REPORTING_REPORTER_HPP
#define HEMELB_REPORTING_REPORTER_HPP
#include "Reporter.h"
#include <iomanip>
namespace hemelb
{
  namespace reporting
  {
    template<class ClockPolicy, class CommsPolicy, class BroadcastPolicy> ReporterBase<
        ClockPolicy,  CommsPolicy, BroadcastPolicy>::ReporterBase(const std::string &apath,
                                                                               const std::string &inputFile,
                                                                               const site_t* fluidSitesOnEachProcessor,
                                                                               const long int aSiteCount,
                                                                               const TimersBase<
                                                                                   ClockPolicy,
                                                                                   CommsPolicy>& timers,
                                                                               const lb::SimulationState &aState,
                                                                               const lb::IncompressibilityChecker<
                                                                                   BroadcastPolicy> &aChecker) :
      path(apath),snapshotCount(0), imageCount(0),
          fluidSitesOnEachProcessor(fluidSitesOnEachProcessor), siteCount(aSiteCount),
          stability(true), timings(timers), state(aState), incompressibilityChecker(aChecker), dictionary("Reporting dictionary")
    {
        dictionary.SetValue("CONFIG",inputFile);
    }

    template<class ClockPolicy,  class CommsPolicy, class BroadcastPolicy> void ReporterBase<
        ClockPolicy,  CommsPolicy, BroadcastPolicy>::Image()
    {
      imageCount++;
      dictionary.AddSectionDictionary("IMAGE")->SetIntValue("COUNT",imageCount);
      ////WriterPolicy::Print("Image written: %u\n", imageCount);
    }

    template<class ClockPolicy,  class CommsPolicy, class BroadcastPolicy> void ReporterBase<
        ClockPolicy,  CommsPolicy, BroadcastPolicy>::Snapshot()
    {
      snapshotCount++;
      dictionary.AddSectionDictionary("SNAPSHOT")->SetIntValue("COUNT",snapshotCount);
    }

    template<class ClockPolicy,  class CommsPolicy, class BroadcastPolicy> void ReporterBase<
    ClockPolicy,  CommsPolicy, BroadcastPolicy>::Write(const std::string &ctemplate,const std::string &as)
    {
        std::string output;
        ctemplate::ExpandTemplate(ctemplate, ctemplate::DO_NOT_STRIP, &dictionary, &output);
        std::string to=path+"/"+as;
        std::fstream file(to.c_str(),std::ios_base::out);
        file<<output<<std::flush;
        file.close();
    };

    template<class ClockPolicy,  class CommsPolicy, class BroadcastPolicy> void ReporterBase<
        ClockPolicy,  CommsPolicy, BroadcastPolicy>::FillDictionary()
    {

      // Note that CycleId is 1-indexed and will have just been incremented when we finish.
      unsigned long cycles = state.GetCycleId() - 1;
      dictionary.SetIntValue("THREADS",CommsPolicy::GetProcessorCount());
      dictionary.SetIntValue("MACHINES",CommsPolicy::GetMachineCount());
      dictionary.SetIntValue("DEPTHS",CommsPolicy::GetDepths());
      dictionary.SetIntValue("SITES",siteCount);
      dictionary.SetIntValue("CYCLES",cycles);
      dictionary.SetIntValue("STEPS",state.GetTimeStepsPassed() - 1);
      dictionary.SetFormattedValue("STEPS_PER_SECOND","%.3f", (state.GetTimeStepsPassed() - 1) / timings[Timers::simulation].Get());
      dictionary.SetIntValue("STEPS_PER_CYCLE",state.GetTimeStepsPerCycle());
      if(!stability){
        dictionary.AddSectionDictionary("UNSTABLE");
      }
      dictionary.SetIntValue("STEPS_PER_CYCLE",state.GetTimeStepsPerCycle());

      if (incompressibilityChecker.AreDensitiesAvailable()
          && !incompressibilityChecker.IsDensityDiffWithinRange())
      {
        ctemplate::TemplateDictionary *incomp=dictionary.AddSectionDictionary("DENSITIES");
        incomp->
            SetFormattedValue("ALLOWED","%.1f%%",incompressibilityChecker.GetMaxRelativeDensityDifferenceAllowed() * 100);
        incomp->
                    SetFormattedValue("ACTUAL","%.1f%%",incompressibilityChecker.GetMaxRelativeDensityDifference() * 100);
      }

      for (hemelb::proc_t n = 0; n < CommsPolicy::GetProcessorCount(); n++)
      {
        ctemplate::TemplateDictionary *proc=dictionary.AddSectionDictionary("PROCESSOR");
        proc->SetIntValue("RANK",n);
        proc->SetIntValue("SITES",fluidSitesOnEachProcessor[n]);
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
      for (unsigned int ii = 0; ii < Timers::numberOfTimers; ii++)
      {
        ctemplate::TemplateDictionary *timer=dictionary.AddSectionDictionary("TIMER");
        timer->SetValue("NAME",timerNames[ii]);
        timer->SetFormattedValue("LOCAL","%.3g",timings[ii].Get()/ normalisations[ii] );
        timer->SetFormattedValue("MIN","%.3g",timings.Mins()[ii] / normalisations[ii] );
        timer->SetFormattedValue("MEAN","%.3g",timings.Means()[ii]/ normalisations[ii] );
        timer->SetFormattedValue("MAX","%.3g",timings.Maxes()[ii] / normalisations[ii] );
      }
    }
  }
}
#endif //HEMELB_REPORTING_REPORTER_HPP

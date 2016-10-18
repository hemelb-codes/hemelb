
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REPORTING_TIMERS_HPP
#define HEMELB_REPORTING_TIMERS_HPP

#include "reporting/Timers.h"

namespace hemelb
{
  namespace reporting
  {
    template<class ClockPolicy>
    void TimersBase<ClockPolicy>::Reduce()
    {
      std::vector<double> timings(numberOfTimers);
      for (unsigned int ii = 0; ii < numberOfTimers; ii++)
      {
        timings[ii] = timers[ii].Get();
      }

      maxes = comms->Reduce(timings, MPI_MAX, 0);
      means = comms->Reduce(timings, MPI_SUM, 0);
      mins = comms->Reduce(timings, MPI_MIN, 0);
      const int np = comms->Size();
      std::transform(means.begin(), means.end(), means.begin(),
		     [&](double x) { return x / np; });
    }

    template<class ClockPolicy>
    void TimersBase<ClockPolicy>::Report(ctemplate::TemplateDictionary& dictionary)
    {
      dictionary.SetIntValue("THREADS", comms->Size());

      for (unsigned int ii = 0; ii < numberOfTimers; ii++)
      {
        ctemplate::TemplateDictionary *timer = dictionary.AddSectionDictionary("TIMER");
        timer->SetValue("NAME", timerNames[ii]);
        timer->SetFormattedValue("LOCAL", "%.3g", timers[ii].Get());
        timer->SetFormattedValue("MIN", "%.3g", Mins()[ii]);
        timer->SetFormattedValue("MEAN", "%.3g", Means()[ii]);
        timer->SetFormattedValue("MAX", "%.3g", Maxes()[ii]);
      }
    }

  }
}

#endif // HEMELB_REPORTING_TIMERS_HPP

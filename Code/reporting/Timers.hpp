// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#ifndef HEMELB_REPORTING_TIMERS_HPP
#define HEMELB_REPORTING_TIMERS_HPP

#include "reporting/Timers.h"

namespace hemelb::reporting
{
    template<class ClockPolicy>
    template <typename Communicator>
    void TimersBase<ClockPolicy>::Reduce(Communicator &&comm)
    {
        n_processes = comm.Size();
        std::vector<double> timings(numberOfTimers);
        for (unsigned int i = 0; i < numberOfTimers; i++)
            timings[i] = timers[i].Get();

        maxes = comm.Reduce(timings, MPI_MAX, 0);
        means = comm.Reduce(timings, MPI_SUM, 0);
        mins =  comm.Reduce(timings, MPI_MIN, 0);

        for (auto& m: means)
            m /= n_processes;
    }

    template<class ClockPolicy>
    void TimersBase<ClockPolicy>::Report(Dict& dictionary)
    {
      dictionary.SetIntValue("THREADS", n_processes);

      for (unsigned int ii = 0; ii < numberOfTimers; ii++)
      {
        Dict timer = dictionary.AddSectionDictionary("TIMER");
        timer.SetValue("NAME", timerNames[ii]);
        timer.SetFormattedValue("LOCAL", "%.3g", timers[ii].Get());
        timer.SetFormattedValue("MIN", "%.3g", Mins()[ii]);
        timer.SetFormattedValue("MEAN", "%.3g", Means()[ii]);
        timer.SetFormattedValue("MAX", "%.3g", Maxes()[ii]);
      }
    }

}
#endif // HEMELB_REPORTING_TIMERS_HPP

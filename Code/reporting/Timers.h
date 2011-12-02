#ifndef HEMELB_REPORTING_TIMERS_H
#define HEMELB_REPORTING_TIMERS_H

#include <vector>
#include "util/utilityFunctions.h"
#include "Policies.h"
namespace hemelb
{
  namespace reporting
  {
    template<class ClockPolicy> class TimerBase : public ClockPolicy
    {
      public:
        TimerBase() :
            start(0), time(0)
        {
        }
        double Get() const
        {
          return time;
        }
        void Set(double t)
        {
          time = t;
        }
        void Start()
        {
          start = ClockPolicy::CurrentTime();
        }
        void Stop()
        {
          time += ClockPolicy::CurrentTime() - start;
        }
      private:
        double start;
        double time;
    };

    /***
     * Manages a set of timings associated with the run
     */

    template<class ClockPolicy, class CommsPolicy> class TimersBase : public CommsPolicy
    {
      public:
        typedef TimerBase<ClockPolicy> Timer;
        enum TimerName
        {
          total,
          domainDecomposition,
          fileRead,
          netInitialise,
          lb,
          visualisation,
          mpiSend,
          mpiWait,
          snapshot,
          simulation,
          last
        };
        static const unsigned int numberOfTimers = last;
        TimersBase() :
            timers(numberOfTimers), maxes(numberOfTimers), mins(numberOfTimers), means(numberOfTimers)
        {
        }
        const std::vector<double> &Maxes() const
        {
          return maxes;
        }
        const std::vector<double> &Mins() const
        {
          return mins;
        }
        const std::vector<double> &Means() const
        {
          return means;
        }
        Timer & operator[](TimerName t)
        {
          return timers[t];
        }
        const Timer & operator[](TimerName t) const
        {
          return timers[t];
        }
        Timer & operator[](unsigned int t)
        {
          return timers[t];
        }
        const Timer & operator[](unsigned int t) const
        {
          return timers[t];
        }
        void Reduce();
      private:
        std::vector<Timer> timers;
        std::vector<double> maxes;
        std::vector<double> mins;
        std::vector<double> means;
    };
    typedef TimerBase<HemeLBClockPolicy> Timer;
    typedef TimersBase<HemeLBClockPolicy, MPICommsPolicy> Timers;
  }

  static const std::string timerNames[hemelb::reporting::Timers::numberOfTimers] =
      { "Total",
        "Domain Decomposition",
        "File Read",
        "Net initialisation",
        "Lattice Boltzmann",
        "Visualisation",
        "MPI Send",
        "MPI Wait",
        "Snapshots",
        "Simulation total" };
}

#endif //HEMELB_REPORTING_TIMERS_H

#ifndef HEMELB_TIMERS_H
#define HEMELB_TIMERS_H

#include <vector>
#include "util/utilityFunctions.h"
namespace hemelb
{
  namespace reporting
  {

    class Timer
    {
      public:
        Timer() :
            start(0), time(0)
        {
        }
        double Get()
        {
          return time;
        }
        void Set(double t)
        {
          time = t;
        }
        void Start()
        {
          start = CurrentTime();
        }
        void Stop()
        {
          time += CurrentTime() - start;
        }
      private:
        double start;
        double time;
        static double CurrentTime()
        {
          return hemelb::util::myClock();
        }
    };

    /***
     * Manages a set of timings associated with the run
     */
    class Timers
    {
      public:
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
        Timers() :
            timers(numberOfTimers), maxes(numberOfTimers), mins(numberOfTimers), means(numberOfTimers)
        {
        }
        std::vector<double> &Maxes()
        {
          return maxes;
        }
        std::vector<double> &Mins()
        {
          return mins;
        }
        std::vector<double> &Means()
        {
          return means;
        }
        Timer & operator[](TimerName t)
        {
          return timers[t];
        }
        Timer & operator[](unsigned int t)
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

#endif

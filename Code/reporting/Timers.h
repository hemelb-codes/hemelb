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
            timers(numberOfTimers)
        {
        }
        Timer & operator[](TimerName t)
        {
          return timers[t];
        }
      private:
        std::vector<Timer> timers;
    };
  }
}

#endif

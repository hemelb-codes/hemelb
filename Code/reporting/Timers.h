#ifndef HEMELB_TIMERS_H
#define HEMELB_TIMERS_H

#include <vector>
#include "util/utilityFunctions.h"
#include "topology/NetworkTopology.h"
namespace hemelb
{
  namespace reporting
  {
    class HemeLBClockPolicy
    {
      protected:
        static double CurrentTime()
        {
          return hemelb::util::myClock();
        }
    };

    class MPICommsPolicy
    {
      protected:
        int Reduce(void *sendbuf,
                   void *recvbuf,
                   int count,
                   MPI_Datatype datatype,
                   MPI_Op op,
                   int root,
                   MPI_Comm comm)
        {
          return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm);
        }

    };

    template<class ClockPolicy> class TimerBase : public ClockPolicy
    {
      public:
        TimerBase() :
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
        void Reduce(){
          double timings[numberOfTimers];
          for (unsigned int ii = 0; ii < numberOfTimers; ii++) {
            timings[ii] = timers[ii].Get();
          }

          CommsPolicy::Reduce(timings, &maxes[0], numberOfTimers, hemelb::MpiDataType<double>(), MPI_MAX, 0, MPI_COMM_WORLD);
          CommsPolicy::Reduce(timings, &means[0], numberOfTimers, hemelb::MpiDataType<double>(), MPI_SUM, 0, MPI_COMM_WORLD);
          CommsPolicy::Reduce(timings, &mins[0],  numberOfTimers, hemelb::MpiDataType<double>(), MPI_MIN, 0, MPI_COMM_WORLD);
          for (unsigned int ii = 0; ii < numberOfTimers; ii++) {
            means[ii] /= (double) (hemelb::topology::NetworkTopology::Instance()->GetProcessorCount());
          }
        };
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

#endif

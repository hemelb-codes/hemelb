#ifndef HEMELB_REPORTING_TIMERS_H
#define HEMELB_REPORTING_TIMERS_H

#include <vector>
#include "util/utilityFunctions.h"
#include "Policies.h"
namespace hemelb
{
  namespace reporting
  {
    /**
     * Timer which manages performance measurement for a single aspect of the code
     * @tparam ClockPolicy Policy defining how to get the current time
     */
    template<class ClockPolicy> class TimerBase : public ClockPolicy
    {
      public:
        /**
         * Starts with the timer stopped.
         */
        TimerBase() :
            start(0), time(0)
        {
        }
        /**
         * Get the current total time spent on this timer
         * @return current total time spent on this timer
         */
        double Get() const
        {
          return time;
        }
        /**
         * Force the current total time to an arbitrary value.
         * Use to reset if the simulation restarts
         * @param t The time to which to set the timer
         */
        void Set(double t)
        {
          time = t;
        }
        /**
         * Start the timer.
         */
        void Start()
        {
          start = ClockPolicy::CurrentTime();
        }
        /**
         * Stop the timer.
         */
        void Stop()
        {
          time += ClockPolicy::CurrentTime() - start;
        }
      private:
        double start; //! Time when the timer was last started.
        double time; //! Current running total time.
    };

    /**
     * Class which manages a set of timers timing aspects of a HemeLB run
     * @tparam ClockPolicy How to get the current time
     * @tparam CommsPolicy How to share information between processes
     */
    template<class ClockPolicy, class CommsPolicy> class TimersBase : public CommsPolicy
    {
      public:
        typedef TimerBase<ClockPolicy> Timer;
        /**
         * The set of possible timers
         */
        enum TimerName
        {
          total, //!< Total time
          domainDecomposition, //!< Time spent in domain decomposition
          fileRead, //!< Time spent in reading the geometry description file
          netInitialise, //!< Time spent initialising the network topology
          lb, //!< Time spent doing the core lattice boltzman simulation
          visualisation, //!< Time spent on visualisation
          mpiSend, //!< Time spent sending MPI data
          mpiWait, //!< Time spent waiting for MPI
          snapshot, //!< Time spent producing snapshots
          simulation, //!< Total time for running the simulation
          last //!< last
        };
        static const unsigned int numberOfTimers = last;
        TimersBase() :
            timers(numberOfTimers), maxes(numberOfTimers), mins(numberOfTimers), means(numberOfTimers)
        {
        }
        /**
         * Max across all processes.
         * Following the sharing of timing data between processes, the max time across all processes for each timer.
         * @return the max time across all processes for each timer.
         */
        const std::vector<double> &Maxes() const
        {
          return maxes;
        }
        /**
         * Min across all processes.
         * Following the sharing of timing data between processes, the minimum time across all processes for each timer.
         * @return the minimum time across all processes for each timer.
         */
        const std::vector<double> &Mins() const
        {
          return mins;
        }
        /**
         * Averages across all processes.
         * Following the sharing of timing data between processes, the average time across all processes for each timer.
         * @return the average time across all processes for each timer.
         */
        const std::vector<double> &Means() const
        {
          return means;
        }
        /**
         * The timer for the given timer name
         * @param t the timer name
         * @return Reference to the given timer
         */
        Timer & operator[](TimerName t)
        {
          return timers[t];
        }
        /**
         * The timer for the given timer name
         * @param t the timer name
         * @return Reference to the given timer
         */
        const Timer & operator[](TimerName t) const
        {
          return timers[t];
        }
        /**
         * The timer for the given timer name
         * @param t the timer name
         * @return Reference to the given timer
         */
        Timer & operator[](unsigned int t)
        {
          return timers[t];
        }
        /**
         * The timer for the given timer name
         * @param t the timer name
         * @return Reference to the given timer
         */
        const Timer & operator[](unsigned int t) const
        {
          return timers[t];
        }
        /**
         * Share timing information across timers
         */
        void Reduce();
      private:
        std::vector<Timer> timers; //! The set of timers
        std::vector<double> maxes; //! Max across processes
        std::vector<double> mins; //! Min across processes
        std::vector<double> means; //! Average across processes
    };
    typedef TimerBase<HemeLBClockPolicy> Timer;
    typedef TimersBase<HemeLBClockPolicy, MPICommsPolicy> Timers;
  }

  /**
   * String message label for each timer for reporting
   */
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

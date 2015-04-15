// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

#ifndef HEMELB_REPORTING_TIMERS_H
#define HEMELB_REPORTING_TIMERS_H

#include <vector>
#include "reporting/Reportable.h"
#include "util/utilityFunctions.h"
#include "reporting/Policies.h"
namespace hemelb
{
  namespace reporting
  {
    /**
     * Timer which manages performance measurement for a single aspect of the code
     * @tparam ClockPolicy Policy defining how to get the current time
     */
    template<class ClockPolicy>
    class TimerBase : public ClockPolicy
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
    template<class ClockPolicy, class CommsPolicy>
    class TimersBase : public CommsPolicy,
                       public Reportable
    {
      public:
        typedef TimerBase<ClockPolicy> Timer;
        /**
         * The set of possible timers
         */
        enum TimerName
        {
          total = 0, //!< Total time
          initialDecomposition, //!< Initial seed decomposition
          domainDecomposition, //!< Time spent in parmetis domain decomposition
          fileRead, //!< Time spent in reading the geometry description file
          reRead, //!< Time spend in re-reading the geometry after second decomposition
          unzip, //!< Time spend in un-zipping
          moves, //!< Time spent moving things around post-parmetis
          parmetis, //!< Time spent in Parmetis
          latDatInitialise, //!< Time spent initialising the lattice data
          lb, //!< Time spent doing the core lattice boltzman simulation
          lb_calc, //!< Time spent doing calculations in the core lattice boltzmann simulation
          visualisation, //!< Time spent on visualisation
          monitoring, //!< Time spent monitoring for stability, compressibility, etc.
          mpiSend, //!< Time spent sending MPI data
          mpiWait, //!< Time spent waiting for MPI
          simulation, //!< Total time for running the simulation,
          readNet,
          readParse,
          readBlock,
          readBlocksPrelim,
          readBlocksAll,
          steeringWait, //!< Time spent waiting for a steering client in wait on connect mode.
          moveForcingNumbers,
          moveForcingData,
          blockRequirements,
          moveCountsSending,
          moveDataSending,
          PopulateOptimisationMovesList,
          InitialGeometryRead,
          colloidInitialisation,
          colloidCommunicatePositions,
          colloidCommunicateVelocities,
          colloidCalculateForces,
          colloidUpdateCalculations,
          colloidOutput,
          extractionWriting,
          last
        //!< last, this has to be the last element of the enumeration so it can be used to track cardinality
        };
        static const unsigned int numberOfTimers = last;

        /**
         * String message label for each timer for reporting
         */
        static const std::string timerNames[TimersBase::numberOfTimers];

        TimersBase(const net::IOCommunicator& comms) :
            CommsPolicy(comms), timers(numberOfTimers), maxes(numberOfTimers), mins(numberOfTimers),
                means(numberOfTimers)
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

        void Report(ctemplate::TemplateDictionary& dictionary);

      private:
        std::vector<Timer> timers; //! The set of timers
        std::vector<double> maxes; //! Max across processes
        std::vector<double> mins; //! Min across processes
        std::vector<double> means; //! Average across processes
    };
    typedef TimerBase<HemeLBClockPolicy> Timer;
    typedef TimersBase<HemeLBClockPolicy, MPICommsPolicy> Timers;

    template<class ClockPolicy, class CommsPolicy>
    const std::string TimersBase<ClockPolicy, CommsPolicy>::timerNames[TimersBase<ClockPolicy,
        CommsPolicy>::numberOfTimers] =

    { "Total",
      "Seed Decomposition",
      "Domain Decomposition",
      "File Read",
      "Re Read",
      "Unzip",
      "Moves",
      "Parmetis",
      "Lattice Data initialisation",
      "Lattice Boltzmann",
      "LB calc only",
      "Visualisation",
      "Monitoring",
      "MPI Send",
      "MPI Wait",
      "Simulation total",
      "Reading communications",
      "Parsing",
      "Read IO",
      "Read Blocks prelim",
      "Read blocks all",
      "Steering Client Wait",
      "Move Forcing Counts",
      "Move Forcing Data",
      "Block Requirements",
      "Move Counts Sending",
      "Move Data Sending",
      "Populating moves list for decomposition optimisation",
      "Initial geometry reading",
      "Colloid initialisation",
      "Colloid position communication",
      "Colloid velocity communication",
      "Colloid force calculations",
      "Colloid calculations for updating",
      "Colloid outputting",
      "Extraction writing" };
  }

}

#endif //HEMELB_REPORTING_TIMERS_H

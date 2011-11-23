#ifndef HEMELB_TIMERS_HPP
#define HEMELB_TIMERS_HPP
#include "Timers.h"
namespace hemelb
{
  namespace reporting
  {
  	template<class ClockPolicy, class CommsPolicy> void TimersBase<ClockPolicy,CommsPolicy>::Reduce(){
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
  	}
  }
}

#endif

#include "Timers.h"
#include "topology/NetworkTopology.h"

namespace hemelb{
  namespace reporting{

    template <> void Timers::Reduce(){

      timers[total].Stop();
      double timings[numberOfTimers];
      for (unsigned int ii = 0; ii < numberOfTimers; ii++) {
        timings[ii] = timers[ii].Get();
      }

      MPI_Reduce(timings, &maxes[0], numberOfTimers, hemelb::MpiDataType<double>(), MPI_MAX, 0, MPI_COMM_WORLD);
      MPI_Reduce(timings, &means[0], numberOfTimers, hemelb::MpiDataType<double>(), MPI_SUM, 0, MPI_COMM_WORLD);
      MPI_Reduce(timings, &mins[0],  numberOfTimers, hemelb::MpiDataType<double>(), MPI_MIN, 0, MPI_COMM_WORLD);
      for (unsigned int ii = 0; ii < numberOfTimers; ii++) {
        means[ii] /= (double) (hemelb::topology::NetworkTopology::Instance()->GetProcessorCount());
      }
    }
  }
}

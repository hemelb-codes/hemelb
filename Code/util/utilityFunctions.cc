#include "mpiInclude.h"
#include "util/utilityFunctions.h"

namespace hemelb
{
  namespace util
  {
    // Returns the number of seconds to 6dp elapsed since the Epoch
    double myClock()
    {
      return MPI_Wtime();
    }

  }

}

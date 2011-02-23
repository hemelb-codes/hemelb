#include "steering/ImageSendComponent.h"
#include "mpiInclude.h"

namespace hemelb
{
  namespace steering
  {

    /**
     * Return seconds since epoch to microsec precision.
     */
    double ImageSendComponent::frameTiming()
    {
      return MPI_Wtime();
    }

  }
}

#include "vis/streaklineDrawer/StreakPixel.h"

namespace hemelb
{
  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::streaklinedrawer::StreakPixel>::RegisterMpiDataType()
  {
    MPI_Datatype ret = vis::streaklinedrawer::StreakPixel::GetMPIType();
    MPI_Type_commit(&ret);
    return ret;
  }
}

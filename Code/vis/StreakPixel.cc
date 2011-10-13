#include "vis/StreakPixel.h"

namespace hemelb
{
  template<>
  MPI_Datatype MpiDataTypeTraits<hemelb::vis::StreakPixel>::RegisterMpiDataType()
  {
    MPI_Datatype ret = vis::StreakPixel::GetMPIType();
    MPI_Type_commit(&ret);
    return ret;
  }
}

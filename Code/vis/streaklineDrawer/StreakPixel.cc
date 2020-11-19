// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "vis/streaklineDrawer/StreakPixel.h"

namespace hemelb
{
  namespace net
  {
    template<>
    MPI_Datatype MpiDataTypeTraits<hemelb::vis::streaklinedrawer::StreakPixel>::RegisterMpiDataType()
    {
      MPI_Datatype ret = vis::streaklinedrawer::StreakPixel::GetMPIType();
      MPI_Type_commit(&ret);
      return ret;
    }
  }
}

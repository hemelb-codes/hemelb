// 
// Copyright (C) University College London, 2007-2012, all rights reserved.
// 
// This file is part of HemeLB and is CONFIDENTIAL. You may not work 
// with, install, use, duplicate, modify, redistribute or share this
// file, or any part thereof, other than as allowed by any agreement
// specifically made by you with University College London.
// 

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

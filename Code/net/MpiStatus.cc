// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/MpiStatus.h"

namespace hemelb
{
  namespace net
  {
    MpiStatus::MpiStatus() :
      statPtr()
    {
    }
    MpiStatus::MpiStatus(MPI_Status stat) : statPtr()
    {
      statPtr.reset(new MPI_Status(stat));
    }
  }
}

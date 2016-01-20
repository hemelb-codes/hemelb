
// This file is part of HemeLB and is Copyright (C)
// the HemeLB team and/or their institutions, as detailed in the
// file AUTHORS. This software is provided under the terms of the
// license in the file LICENSE.

#include "net/mpi.h"
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
